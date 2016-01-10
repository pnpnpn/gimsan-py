#ifndef _PRINT_JSON_H
#define _PRINT_JSON_H

#include "gibbs_util.h"

extern void printJsonRunNode(
		FILE *fptr, 
		RunNode *node, 
		enum ScoreMetric metric, 
		Markov *markov, 
		Dataset *data, 
		Zoops *zoops);


#endif
