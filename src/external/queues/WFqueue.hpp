#ifndef WFQUEUE_H
#define WFQUEUE_H
#include "align.h"
#include "primitives.h"

#ifndef WFQUEUE_NODE_SIZE
#define WFQUEUE_NODE_SIZE ((1 << 10) - 2)
#endif

#ifndef MAX_SPIN
#define MAX_SPIN 100
#endif

#ifndef MAX_PATIENCE
#define MAX_PATIENCE 10
#endif


class WFQueue {
private:
#define EMPTY ((void *) 0)
#define BOT ((void *) 0)
#define TOP ((void *)-1)    
#define MAX_GARBAGE(n) (2 * n)

    const size_t N = WFQUEUE_NODE_SIZE;

    struct _enq_t {
	long volatile id;
	void * volatile val;
    } CACHE_ALIGNED;
    
    struct _deq_t {
	long volatile id;
	long volatile idx;
    } CACHE_ALIGNED;
    
    struct _cell_t {
	void * volatile val;
	struct _enq_t * volatile enq;
	struct _deq_t * volatile deq;
	void * pad[5];
    };
    
    struct _node_t {
	struct _node_t * volatile next CACHE_ALIGNED;
	long id CACHE_ALIGNED;
	struct _cell_t cells[WFQUEUE_NODE_SIZE] CACHE_ALIGNED;
    };
    
    typedef struct DOUBLE_CACHE_ALIGNED {
	/**
	 * Index of the next position for enqueue.
	 */
	volatile long Ei DOUBLE_CACHE_ALIGNED;
	
	/**
	 * Index of the next position for dequeue.
	 */
	volatile long Di DOUBLE_CACHE_ALIGNED;
	
	/**
	 * Index of the head of the queue.
	 */
	volatile long Hi DOUBLE_CACHE_ALIGNED;
	
	/**
	 * Pointer to the head node of the queue.
	 */
	struct _node_t * volatile Hp;
	
	/**
	 * Number of processors.
	 */
	long nprocs;
#ifdef RECORD
	long slowenq;
	long slowdeq;
	long fastenq;
	long fastdeq;
	long empty;
#endif
    } queue_t;
    
    typedef struct _handle_t {
	/**
	 * Pointer to the next handle.
	 */
	struct _handle_t * next;
	
	/**
	 * Hazard pointer.
	 */
	struct _node_t * volatile Hp;
	
	/**
	 * Pointer to the node for enqueue.
	 */
	struct _node_t * volatile Ep;
	
	/**
	 * Pointer to the node for dequeue.
	 */
	struct _node_t * volatile Dp;
	
	/**
	 * Enqueue request.
	 */
	struct _enq_t Er CACHE_ALIGNED;
	
	/**
	 * Dequeue request.
	 */
	struct _deq_t Dr CACHE_ALIGNED;
	
	/**
	 * Handle of the next enqueuer to help.
	 */
	struct _handle_t * Eh CACHE_ALIGNED;
	
	long Ei;
	
	/**
	 * Handle of the next dequeuer to help.
	 */
	struct _handle_t * Dh;
	
	/**
	 * Pointer to a spare node to use, to speedup adding a new node.
	 */
	struct _node_t * spare CACHE_ALIGNED;
	
	/**
	 * Count the delay rounds of helping another dequeuer.
	 */
	int delay;
	
#ifdef RECORD
	long slowenq;
	long slowdeq;
	long fastenq;
	long fastdeq;
	long empty;
#endif
    } handle_t;
    
    typedef struct _enq_t enq_t;
    typedef struct _deq_t deq_t;
    typedef struct _cell_t cell_t;
    typedef struct _node_t node_t;
       
    queue_t * q = NULL;
    handle_t ** hds;
    pthread_barrier_t barrier;

private:

    inline void * spin(void * volatile * p) {
	int patience = MAX_SPIN;
	void * v = *p;
	
	while (!v && patience-- > 0) {
	    v = *p;
	    PAUSEQ();
	}
	
	return v;
    }
    
    inline node_t * new_node() {
	node_t * n = (node_t*)align_malloc(PAGE_SIZE, sizeof(node_t));
	memset(n, 0, sizeof(node_t));
	return n;
    }
    
    node_t * update(node_t * volatile * pPn, node_t * cur,
		    node_t * volatile * pHp) {
	node_t * ptr = *pPn;
	
	if (ptr->id < cur->id) {
	    if (!CAScs(pPn, &ptr, cur)) {
		if (ptr->id < cur->id) cur = ptr;
	    }
	    
	    node_t * Hp = *pHp;
	    if (Hp && Hp->id < cur->id) cur = Hp;
	}
	
	return cur;
    }
    
    void cleanup(handle_t * th) {
	long oid = q->Hi;
	node_t * _new = th->Dp;
	
	if (oid == -1) return;
	if (_new->id - oid < MAX_GARBAGE(q->nprocs)) return;
	if (!CASa(&q->Hi, &oid, -1)) return;
	
	node_t * old = q->Hp;
	handle_t * ph = th;
	handle_t * phs[q->nprocs];
	int i = 0;
	
	do {
	    node_t * Hp = ACQUIRE(&ph->Hp);
	    if (Hp && Hp->id < _new->id) _new = Hp;
	    
	    _new = update(&ph->Ep, _new, &ph->Hp);
	    _new = update(&ph->Dp, _new, &ph->Hp);
	    
	    phs[i++] = ph;
	    ph = ph->next;
	} while (_new->id > oid && ph != th);
	
	while (_new->id > oid && --i >= 0) {
	    node_t * Hp = ACQUIRE(&phs[i]->Hp);
	    if (Hp && Hp->id < _new->id) _new = Hp;
	}
	
	long nid = _new->id;
	
	if (nid <= oid) {
	    RELEASE(&q->Hi, oid);
	} else {
	    q->Hp = _new;
	    RELEASE(&q->Hi, nid);
	    
	    while (old != _new) {
		node_t * tmp = old->next;
		free(old);
		old = tmp;
	    }
	}
    }
    
    cell_t * find_cell(node_t * volatile * p, long i, handle_t * th) {
	node_t * c = *p;
	
	for (size_t j = c->id; j < i / N; ++j) {
	    node_t * n = c->next;
	    
	    if (n == NULL) {
		node_t * t = th->spare;
		
		if (t == NULL) {
		    t = new_node();
		    th->spare = t;
		}
		
		t->id = j + 1;
		
		if (CASra(&c->next, &n, t)) {
		    n = t;
		    th->spare = NULL;
		}
	    }
	    
	    c = n;
	}
	
	*p = c;
	return &c->cells[i % N];
    }
    
    int enq_fast(handle_t * th, void * v, long * id) {
	long i = FAAcs(&q->Ei, 1);
	cell_t * c = find_cell(&th->Ep, i, th);
	void * cv = BOT;
	
	if (CAS(&c->val, &cv, v)) {
#ifdef RECORD
	    th->fastenq++;
#endif
	    return 1;
	} else {
	    *id = i;
	    return 0;
	}
    }
    
    void enq_slow(handle_t * th, void * v, long id) {
	enq_t * enq = &th->Er;
	enq->val = v;
	RELEASE(&enq->id, id);
	
	node_t * tail = th->Ep;
	long i; cell_t * c;
	
	do {
	    i = FAA(&q->Ei, 1);
	    c = find_cell(&tail, i, th);
	    enq_t * ce = (enq_t*)BOT;
	    
	    if (CAScs(&c->enq, &ce, enq) && c->val != TOP) {
		if (CAS(&enq->id, &id, -i)) id = -i;
		break;
	    }
	} while (enq->id > 0);
	
	id = -enq->id;
	c = find_cell(&th->Ep, id, th);
	if (id > i) {
	    long Ei = q->Ei;
	    while (Ei <= id && !CAS(&q->Ei, &Ei, id + 1));
	}
	c->val = v;
	
#ifdef RECORD
	th->slowenq++;
#endif
    }
    
    void * help_enq(handle_t * th, cell_t * c, long i) {
	void * v = spin(&c->val);
	
	if ((v != TOP && v != BOT) ||
	    (v == BOT && !CAScs(&c->val, &v, TOP) && v != TOP)) {
	    return v;
	}
	
	enq_t * e = c->enq;
	
	if (e == BOT) {
	    handle_t * ph; enq_t * pe; long id;
	    ph = th->Eh, pe = &ph->Er, id = pe->id;
	    
	    if (th->Ei != 0 && th->Ei != id) {
		th->Ei = 0;
		th->Eh = ph->next;
		ph = th->Eh, pe = &ph->Er, id = pe->id;
	    }
	    
	    if (id > 0 && id <= i && !CAS(&c->enq, &e, pe))
		th->Ei = id;
	    else
		th->Eh = ph->next;
	    
	    if (e == BOT && CAS(&c->enq, &e, TOP)) e = (enq_t*)TOP;
	}
	
	if (e == TOP) return (q->Ei <= i ? BOT : TOP);
	
	long ei = ACQUIRE(&e->id);
	void * ev = ACQUIRE(&e->val);
	
	if (ei > i) {
	    if (c->val == TOP && q->Ei <= i) return BOT;
	} else {
	    if ((ei > 0 && CAS(&e->id, &ei, -i)) ||
		(ei == -i && c->val == TOP)) {
		long Ei = q->Ei;
		while (Ei <= i && !CAS(&q->Ei, &Ei, i + 1));
		c->val = ev;
	    }
	}
	
	return c->val;
    }

    void help_deq(handle_t * th, handle_t * ph) {
	deq_t * deq = &ph->Dr;
	long idx = ACQUIRE(&deq->idx);
	long id = deq->id;
	
	if (idx < id) return;
	
	node_t * Dp = ph->Dp;
	th->Hp = Dp;
	FENCE();
	idx = deq->idx;
	
	long i = id + 1, old = id, _new = 0;
	while (1) {
	    node_t * h = Dp;
	    for (; idx == old && _new == 0; ++i) {
		cell_t * c = find_cell(&h, i, th);
		
		long Di = q->Di;
		while (Di <= i && !CAS(&q->Di, &Di, i + 1));
		
		void * v = help_enq(th, c, i);
		if (v == BOT || (v != TOP && c->deq == BOT)) _new = i;
		else idx = ACQUIRE(&deq->idx);
	    }
	    
	    if (_new != 0) {
		if (CASra(&deq->idx, &idx, _new)) idx = _new;
		if (idx >= _new) _new = 0;
	    }
	    
	    if (idx < 0 || deq->id != id) break;
	    
	    cell_t * c = find_cell(&Dp, idx, th);
	    deq_t * cd = (deq_t*)BOT;
	    if (c->val == TOP || CAS(&c->deq, &cd, deq) || cd == deq) {
		CAS(&deq->idx, &idx, -idx);
		break;
	    }
	    
	    old = idx;
	    if (idx >= i) i = idx + 1;
	}
    }

    void * deq_fast(handle_t * th, long * id) {
	long i = FAAcs(&q->Di, 1);
	cell_t * c = find_cell(&th->Dp, i, th);
	void * v = help_enq(th, c, i);
	deq_t * cd = (deq_t*)BOT;
	
	if (v == BOT) return BOT;
	if (v != TOP && CAS(&c->deq, &cd, TOP)) return v;
	
	*id = i;
	return TOP;
    }

    void * deq_slow(handle_t * th, long id) {
	deq_t * deq = &th->Dr;
	RELEASE(&deq->id, id);
	RELEASE(&deq->idx, id);
	
	help_deq(th, th);
	long i = -deq->idx;
	cell_t * c = find_cell(&th->Dp, i, th);
	void * val = c->val;
	
#ifdef RECORD
	th->slowdeq++;
#endif
	return val == TOP ? BOT : val;
    }


    
    void enqueue(handle_t * th, void * v) {
	th->Hp = th->Ep;
	
	long id;
	int p = MAX_PATIENCE;
	while (!enq_fast(th, v, &id) && p-- > 0);
	if (p < 0) enq_slow(th, v, id);
	
	RELEASE(&th->Hp, NULL);
    }
    
    void * dequeue(handle_t * th)  {
	th->Hp = th->Dp;
	
	void * v;
	long id;
	int p = MAX_PATIENCE;
	
	do v = deq_fast(th, &id);
	while (v == TOP && p-- > 0);
	if (v == TOP) v = deq_slow(th, id);
	else {
#ifdef RECORD
	    th->fastdeq++;
#endif
	}
	
	if (v != EMPTY) {
	    help_deq(th, th->Dh);
	    th->Dh = th->Dh->next;
	}
	
	RELEASE(&th->Hp, NULL);
	
	if (th->spare == NULL) {
	    cleanup(th);
	    th->spare = new_node();
	}
	
#ifdef RECORD
	if (v == EMPTY) th->empty++;
#endif
	return v;
    }

    void queue_init(int nprocs)  {
	q->Hi = 0;
	q->Hp = new_node();
	
	q->Ei = 1;
	q->Di = 1;
	
	q->nprocs = nprocs;
	
#ifdef RECORD
	q->fastenq = 0;
	q->slowenq = 0;
	q->fastdeq = 0;
	q->slowdeq = 0;
	q->empty = 0;
#endif
	pthread_barrier_init(&barrier, NULL, nprocs);
    }
    
    void queue_register(handle_t * th, int id)  {
	th->next = NULL;
	th->Hp = NULL;
	th->Ep = q->Hp;
	th->Dp = q->Hp;
	
	th->Er.id = 0;
	th->Er.val = BOT;
	th->Dr.id = 0;
	th->Dr.idx = -1;
	
	th->Ei = 0;
	th->spare = new_node();
#ifdef RECORD
	th->slowenq = 0;
	th->slowdeq = 0;
	th->fastenq = 0;
	th->fastdeq = 0;
	th->empty = 0;
#endif
	
	static handle_t * volatile _tail;
	handle_t * tail = _tail;
	
	if (tail == NULL) {
	    th->next = th;
	    if (CASra(&_tail, &tail, th)) {
		th->Eh = th->next;
		th->Dh = th->next;
		return;
	    }
	}
	
	handle_t * next = tail->next;
	do th->next = next;
	while (!CASra(&tail->next, &next, th));
	
	th->Eh = th->next;
	th->Dh = th->next;
	
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



#endif /* end of include guard: WFQUEUE_H */
