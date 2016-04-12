#ifndef CCQUEUE_H
#define CCQUEUE_H
#include "align.h"
#include "primitives.h"
#include "ccsynch.h"


typedef struct _node_t {
    struct _node_t * next CACHE_ALIGNED;
    void * volatile data;
} node_t;


static inline void serialEnqueue(void * state, void * data) {
    node_t * volatile * tail = (node_t **) state;
    node_t * node = (node_t *) data;
 
    (*tail)->next = node;
    *tail = node;
}
    
static inline void serialDequeue(void * state, void * data) {
    node_t * volatile * head = (node_t **) state;
    node_t ** ptr = (node_t **) data;
    
    node_t * node = *head;
    node_t * next = node->next;
    
    if (next) {
	node->data = next->data;
	*head = next;
    } else {
	node = (node_t *) -1;
    }
    
    *ptr = node;
}

class CCQueue {
private:
#define EMPTY (void *) -1
    
    typedef struct _queue_t {
	ccsynch_t enq DOUBLE_CACHE_ALIGNED;
	ccsynch_t deq DOUBLE_CACHE_ALIGNED;
	node_t * head DOUBLE_CACHE_ALIGNED;
	node_t * tail DOUBLE_CACHE_ALIGNED;
    } queue_t DOUBLE_CACHE_ALIGNED;
    
    typedef struct _handle_t {
	ccsynch_handle_t enq;
	ccsynch_handle_t deq;
	node_t * next;
    } handle_t DOUBLE_CACHE_ALIGNED;
    
    
    queue_t * q = NULL;
    handle_t ** hds;
    
private:
    
    void enqueue(handle_t * handle, void * data) {
	node_t * node = handle->next;
	
	if (node) handle->next = NULL;
	else node = (node_t*)align_malloc(CACHE_LINE_SIZE, sizeof(node_t));
	
	node->data = data;
	node->next = NULL;
	
	ccsynch_apply(&q->enq, &handle->enq, &serialEnqueue, &q->tail, node);	
    }
    
    void * dequeue(handle_t * handle)  {
	node_t * node;
	ccsynch_apply(&q->deq, &handle->deq, &serialDequeue, &q->head, &node);
	
	void * data;
	
	if (node == (void *) -1) {
	    data = (void *) -1;
	} else {
	    data = node->data;
	    if (handle->next) free(node);
	    else handle->next = node;
	}
	
	return data;	
    }

    void queue_init(int nprocs)  {
	ccsynch_init(&q->enq);
	ccsynch_init(&q->deq);
	
	node_t * dummy = (node_t*)align_malloc(CACHE_LINE_SIZE, sizeof(node_t));
	dummy->data = 0;
	dummy->next = NULL;
	
	q->head = dummy;
	q->tail = dummy;
    }
    
    void queue_register(handle_t * handle, int id)  {
	ccsynch_handle_init(&handle->enq);
	ccsynch_handle_init(&handle->deq);
	
	handle->next = (node_t*)align_malloc(CACHE_LINE_SIZE, sizeof(node_t));
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

#endif /* end of include guard: CCQUEUE_H */
