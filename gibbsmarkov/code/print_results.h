#ifndef _PRINT_RESULTS_H
#define _PRINT_RESULTS_H

#include "gibbs_util.h"

extern void printRunNode(FILE *fptr, RunNode *node, enum ScoreMetric metric, Markov *markov, Dataset *data, Zoops *zoops);
//extern void printTopRankRunNodes(FILE *fptr, int numTop, RunSet*, enum ScoreMetric metric, Markov*, Dataset*, Zoops*, double pvalCutoff);


#endif
