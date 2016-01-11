#ifndef _MULTI_MOTIFS_H
#define _MULTI_MOTIFS_H

#include "stdinc.h"
#include "profile.h"

typedef struct MotifNode_el {
	struct MotifNode_el *next;
	Profile *countmat; 
} MotifNode;


typedef struct {
	int capacity;
	int size;
	MotifNode *head;

	double *log_fact; //log of factorials
	double *log_int;
	int numseqs;
	double ntFreq[NUMALPHAS];
	double minOverlapRatio;
} MultiMotifs;

extern MultiMotifs* initMultiMotifs(int capacity, int numseqs, double ntFreq[]);
extern void nilMultiMotifs(MultiMotifs *multiMotifs);

extern double computeMotifDistance(Profile *countmat1, Profile *countmat2, MultiMotifs *multiMotifs);
extern void prependToMultiMotifsList(Profile *countmat, MultiMotifs *multiMotifs);

#endif

