#ifndef _HEAP_H
#define _HEAP_H

#include "gibbs_util.h"

typedef struct {
	double score;
	RunNode *rnode;
} node;

typedef struct {
	int size, capacity;
	node **words; 
	node **arr;
} Heap;

Heap* constructHeap(int capacity);
void destructHeap(Heap *h);
void zeroHeap(Heap *h);
void processNode(Heap *h, RunNode *rnode);
RunNode* deleteRootNode(Heap *h);

#endif

