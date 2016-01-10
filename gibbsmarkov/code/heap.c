#include <stdlib.h>
#include "heap.h"


//h->arr is the heap structure
//h->words is a continuous block for faster spatial locality
Heap* constructHeap(int capacity) {
	Heap *h = (Heap*) malloc(sizeof(Heap));
	h->size = 0;
	h->capacity = capacity;
		
	int i;
	h->words = (node**) malloc(capacity * sizeof(node*));
	for ( i =0; i < capacity; i++) {
		h->words[i] = (node*)malloc(sizeof(node));
	}
	h->arr = (node**)malloc(capacity * sizeof(node*));
	for(i = 0; i < capacity; i++) {
		//reset h->arr as a continuous block
		h->arr[i] = h->words[i];
	}
	return h;
}

void destructHeap(Heap *h) {
	int i;
	for ( i =0; i < h->capacity; i++) {
		free(h->words[i]);
	}
	free(h->words);
	free(h->arr);
	free(h);
}


static
void copyNode(node *source, node *dest) {
	dest->rnode = source->rnode;
	dest->score = source->score;
}


static
void upHeap(Heap *h, int index){
	while (index != 0) {
		if (h->arr[index]->score > h->arr[(index-1)/2]->score) {
			node *temp;
			temp = h->arr[index];
			h->arr[index] = h->arr[(index-1)/2];
			h->arr[(index-1)/2] = temp;
			index = (index-1)/2;
		}
		else {
			break;
		}
	}
}

static
void downHeap(Heap *h, int index){
	int left = 2*index+1;
	int right = 2*index+2;
	
	if(left >= h->size) {
		return;
	}
	int maxpos = index;
	if(right < h->size && h->arr[right]->score > h->arr[maxpos]->score) {
		maxpos = right;
	}
	if(h->arr[left]->score > h->arr[maxpos]->score) {
		maxpos = left; 
	}
	if(maxpos != index) {
		//swap node at minpos and index
		node *temp;
		temp = h->arr[index];
		h->arr[index] = h->arr[maxpos];
		h->arr[maxpos] = temp;

		downHeap(h, maxpos);
	}
}

static
void insertNode(Heap *h, node *nptr) {
	if (h->size==0){ /* check if it's empty */
		copyNode(nptr, h->arr[0]);
	}
	else {
		copyNode(nptr, h->arr[h->size]);
		upHeap(h, h->size);
	}
	h->size++;
}



//The heap stores the p elements with lowest score
//arr[0] has the highest score in the heap
void processNode(Heap *h, RunNode *rnode) {
	node nd;
	nd.rnode = rnode;
	nd.score = rnode->score;

	if (h->size != h->capacity){ // check if it's full 
		insertNode(h, &nd);
	}
	else if (nd.score < h->arr[0]->score) {
		copyNode(&nd, h->arr[0]);
		downHeap(h, 0);
	}
}

RunNode* deleteRootNode(Heap *h) {
	if(h->size == 0) {
		return NULL;
	}

	RunNode *rnode = h->arr[0]->rnode;

	node *temp;
	temp = h->arr[0];
	h->arr[0] = h->arr[h->size-1];
	h->arr[h->size-1] = temp;

	h->size--;

	downHeap(h, 0);

	return rnode;
}

void zeroHeap(Heap *h){
	int i;
	for(i = 0; i < h->size; i++) {
		//reset h->arr as a continuous block
		h->arr[i] = h->words[i];
	}
	h->size = 0;
}
