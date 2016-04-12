#ifndef LCRQUEUE_H
#define LCRQUEUE_H
#include "align.h"
#include "hzdptr.h"
#include "primitives.h"

#ifndef LCRQ_RING_SIZE
#define LCRQ_RING_SIZE (1ull << 12)
#endif


class LCRQueue {
private:
#define EMPTY ((void *) -1)
    
    const size_t RING_SIZE = LCRQ_RING_SIZE;
    
    typedef struct RingNode {
	volatile uint64_t val;
	volatile uint64_t idx;
	uint64_t pad[14];
    } RingNode DOUBLE_CACHE_ALIGNED;
    
    typedef struct RingQueue {
	volatile int64_t head DOUBLE_CACHE_ALIGNED;
	volatile int64_t tail DOUBLE_CACHE_ALIGNED;
	struct RingQueue *next DOUBLE_CACHE_ALIGNED;
	RingNode array[LCRQ_RING_SIZE];
    } RingQueue DOUBLE_CACHE_ALIGNED;
    
    typedef struct {
	RingQueue * volatile head DOUBLE_CACHE_ALIGNED;
	RingQueue * volatile tail DOUBLE_CACHE_ALIGNED;
	int nprocs;
    } queue_t;
    
    typedef struct {
	RingQueue * next;
	hzdptr_t hzdptr;
    } handle_t;
    
    queue_t * q = NULL;
    handle_t ** hds;
    
private:

    inline void init_ring(RingQueue *rq) {
	for (size_t i = 0; i < RING_SIZE; i++) {
	    rq->array[i].val = -1;
	    rq->array[i].idx = i;
	}
	
	rq->head = rq->tail = 0;
	rq->next = NULL;
    }
    
    inline int is_empty(uint64_t v)  {
	return (v == (uint64_t)-1);
    }
    
    
    inline uint64_t node_index(uint64_t i) {
	return (i & ~(1ull << 63));
    }
    
    
    inline uint64_t set_unsafe(uint64_t i) {
	return (i | (1ull << 63));
    }
    
    
    inline uint64_t node_unsafe(uint64_t i) {
	return (i & (1ull << 63));
    }
    
    
    inline uint64_t tail_index(uint64_t t) {
	return (t & ~(1ull << 63));
    }
    
    
    inline int crq_is_closed(uint64_t t) {
	return (t & (1ull << 63)) != 0;
    }
    

    inline void fixState(RingQueue *rq) {
	while (1) {
	    uint64_t t = rq->tail;
	    uint64_t h = rq->head;
	    
	    if ((uint64_t)(rq->tail) != t)
		continue;
	    
	    if (h > t) {
		if (CAS(&rq->tail, &t, h)) break;
		continue;
	    }
	    break;
	}
    }

    inline int close_crq(RingQueue *rq, const uint64_t t, const int tries) {
	uint64_t tt = t + 1;
	
	if (tries < 10)
	    return CAS(&rq->tail, &tt, tt|(1ull<<63));
	else
	    return BTAS(&rq->tail, 63);
    }
    
    void lcrq_put(handle_t * handle, uint64_t arg) {
	int try_close = 0;
	
	while (1) {
	    RingQueue *rq = (RingQueue*)hzdptr_setv(&q->tail, &handle->hzdptr, 0);
	    RingQueue *next = rq->next;
	    
	    if (next != NULL) {
		CAS(&q->tail, &rq, next);
		continue;
	    }
	    
	    uint64_t t = FAA(&rq->tail, 1);
	    
	    if (crq_is_closed(t)) {
		RingQueue * nrq;
	    alloc:
		nrq = handle->next;
		
		if (nrq == NULL) {
		    nrq = (RingQueue*)align_malloc(PAGE_SIZE, sizeof(RingQueue));
		    init_ring(nrq);
		}
		
		// Solo enqueue
		nrq->tail = 1;
		nrq->array[0].val = (uint64_t) arg;
		nrq->array[0].idx = 0;
		
		if (CAS(&rq->next, &next, nrq)) {
		    CAS(&q->tail, &rq, nrq);
		    handle->next = NULL;
		    return;
		}
		continue;
	    }
	    
	    RingNode* cell = &rq->array[t & (RING_SIZE-1)];
	    
	    uint64_t idx = cell->idx;
	    uint64_t val = cell->val;
	    
	    if (is_empty(val)) {
		if (node_index(idx) <= t) {
		    if ((!node_unsafe(idx) || (uint64_t)(rq->head) < t) &&
			CAS2(cell, &val, &idx, arg, t)) {
			return;
		    }
		}
	    }
	    
	    uint64_t h = rq->head;
	    
	    if ((int64_t)(t - h) >= (int64_t)RING_SIZE &&
		close_crq(rq, t, ++try_close)) {
		goto alloc;
	    }
	}
	
	hzdptr_clear(&handle->hzdptr, 0);
    }
    
    uint64_t lcrq_get(handle_t * handle) {
	while (1) {
	    RingQueue *rq = (RingQueue *)hzdptr_setv(&q->head, &handle->hzdptr, 0);
	    RingQueue *next;
	    
	    uint64_t h = FAA(&rq->head, 1);
	    
	    RingNode* cell = &rq->array[h & (RING_SIZE-1)];
	    
	    uint64_t tt = 0;
	    int r = 0;
	    
	    while (1) {
		
		uint64_t cell_idx = cell->idx;
		uint64_t unsafe = node_unsafe(cell_idx);
		uint64_t idx = node_index(cell_idx);
		uint64_t val = cell->val;
		
		if (idx > h) break;
		
		if (!is_empty(val)) {
		    if (idx == h) {
			if (CAS2(cell, &val, &cell_idx, -1, (unsafe | h) + RING_SIZE))
			    return val;
		    } else {
			if (CAS2(cell, &val, &cell_idx, val, set_unsafe(idx))) {
			    break;
			}
		    }
		} else {
		    if ((r & ((1ull << 10) - 1)) == 0)
			tt = rq->tail;
		    
		    // Optimization: try to bail quickly if queue is closed.
		    int crq_closed = crq_is_closed(tt);
		    uint64_t t = tail_index(tt);
		    
		    if (unsafe) { // Nothing to do, move along
			if (CAS2(cell, &val, &cell_idx, val, (unsafe | h) + RING_SIZE))
			    break;
		    } else if (t < h + 1 || r > 200000 || crq_closed) {
			if (CAS2(cell, &val, &idx, val, h + RING_SIZE)) {
			    if (r > 200000 && tt > RING_SIZE)
				BTAS(&rq->tail, 63);
			    break;
			}
		    } else {
			++r;
		    }
		}
	    }
	    
	    if (tail_index(rq->tail) <= h + 1) {
		fixState(rq);
		// try to return empty
		next = rq->next;
		if (next == NULL)
		    return -1;  // EMPTY
		if (tail_index(rq->tail) <= h + 1) {
		    if (CAS(&q->head, &rq, next)) {
			hzdptr_retire(&handle->hzdptr, rq);
		    }
		}
	    }
	}
	
	hzdptr_clear(&handle->hzdptr, 0);
    }
    
    void enqueue(handle_t * th, void * val) {
	lcrq_put(th, (uint64_t) val);
    }
    
    void * dequeue(handle_t * th)  {
	return (void *) lcrq_get(th);
    }

    void queue_init(int nprocs)  {
	RingQueue *rq = (RingQueue*)align_malloc(PAGE_SIZE, sizeof(RingQueue));
	init_ring(rq);

	q->head = rq;
	q->tail = rq;
	q->nprocs = nprocs;
    }
    
    void queue_register(handle_t * th, int id)  {
	hzdptr_init(&th->hzdptr, q->nprocs, 1);
    }

    void queue_deregister(handle_t *th, int id) {}
      
public:

    // must be called once before starting using the queue
    bool init(int nprocs) {
	q = (queue_t*)align_malloc(PAGE_SIZE, sizeof(queue_t));
	if (!q) return false;
	queue_init(nprocs);
	hds = (handle_t**)align_malloc(PAGE_SIZE, sizeof(handle_t * [nprocs]));
	if (!hds) {
	    return false;
	}
	return true;
    }

    // called once by a thread before calling push/pop
    void registerq(int pid) {
	hds[pid] = (handle_t*)align_malloc(PAGE_SIZE, sizeof(handle_t));
	queue_register(hds[pid], pid);
    }

    void deregisterq(int pid) {
	queue_deregister(hds[pid], pid);
    }

    bool push(void *data, int pid) {
	enqueue(hds[pid], data);
	//P[pid] += 1;
	return true;
    }

    bool pop(void **val, int pid) {
	*val = dequeue(hds[pid]);
	//if (*val != EMPTY) C[pid] += 1;
	return (*val == EMPTY) ? false: true;
    }

    void printStats(int nthreads) {
	size_t sum = 0;
	for(int i=0;i<nthreads;++i) {
	    printf("TH%d n push = %ld, n pop = %ld\n", i, P[i], C[i]);
	    sum += C[i];
	}
	printf("Total elements read: %ld\n", sum);
    }

private:
    // used only for statistics
    long P[MAX_PROCS] CACHE_ALIGNED;
    long C[MAX_PROCS] CACHE_ALIGNED;
};



#endif /* end of include guard: LCRQUEUE_H */
