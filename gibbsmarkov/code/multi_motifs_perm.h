#ifndef _MULTI_MOTIFS_PERM_H
#define _MULTI_MOTIFS_PERM_H

#include "stdinc.h"
#include "profile.h"
#include "gibbs_util.h"

typedef struct MotifNode_el {
	struct MotifNode_el *next;
	Profile *countmat; 
} MotifNode;


typedef struct {
	int capacity;
	int size;
	MotifNode *head;

	double minOverlapRatio;
	int numseqs;
	int maxspan;
	double distCutoff;

	double *log_fact; //log of factorials
	//double *log_int;
	//double ntFreq[NUMALPHAS];

	Profile *revcomplMat; //temporary matrix
} MultiMotifs;

extern MultiMotifs* initMultiMotifs(int capacity, int numseqs, int maxspan);
extern void nilMultiMotifs(MultiMotifs *multiMotifs);
void setDistCutoffForMultiMotifs(MultiMotifs *multiMotifs, RunSet *runset, int numPermPairs, double pvalCutoff);

extern double computeMotifDist(Profile *countmat1, Profile *countmat2, MultiMotifs *multiMotifs);
extern void prependToMultiMotifsList(Profile *countmat, MultiMotifs *multiMotifs);

#endif

