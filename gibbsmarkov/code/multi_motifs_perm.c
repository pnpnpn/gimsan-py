#include "multi_motifs_perm.h"

static
MotifNode* initMotifNode() {
	MotifNode *node = (MotifNode*) malloc(sizeof(MotifNode));
	node->countmat = NULL;
	node->next = NULL;

	return node;
}

static
void nilMotifNode(MotifNode *node) {
	nilProfile(node->countmat);
	free(node);
}

//Make a copy of the countmat
void prependToMultiMotifsList(Profile *countmat, MultiMotifs *multiMotifs) {
	Profile *copy = initProfile(countmat->maxspan, countmat->maxspan);
	copyProfile(copy, countmat);

	MotifNode *newnode = initMotifNode();
	newnode->countmat = copy;

	newnode->next = multiMotifs->head;
	multiMotifs->head = newnode;
	multiMotifs->size++;
}

MultiMotifs* initMultiMotifs(int capacity, int numseqs, int maxspan) {
	MultiMotifs *mm = (MultiMotifs*) malloc(sizeof(MultiMotifs));
	mm->numseqs = numseqs;
	mm->capacity = capacity;
	mm->minOverlapRatio = 0.50; //minimum overlap ratio to be considered as alignment
	mm->size = 0;
	mm->head = NULL;
	//for(int i = 0; i < NUMALPHAS; i++) {
	//	mm->ntFreq[i] = ntFreq[i];
	//}
	mm->distCutoff = NAN;
	mm->maxspan = maxspan;

	mm->log_fact = (double*) malloc(sizeof(double) * (numseqs+1));
	mm->log_fact[0] = log((double)1.0); // 0.0
	for(int i = 1; i < numseqs + 1; i++) {
		mm->log_fact[i] = mm->log_fact[i-1] + log((double) i );
	}

	//mm->log_int = (double*) malloc(sizeof(double) * (numseqs+1));
	//mm->log_int[0] = LOGZERO; //we want to have 0*log(0) = 0 [-INFINITY * 0 gives NAN]
	//for(int i = 1; i < numseqs + 1; i++) {
	//	mm->log_int[i] = log((double) i );
	//}

	mm->revcomplMat = initProfile(maxspan, maxspan);
	
	if(DEBUG1) {
		for(int i = 1; i < numseqs + 1; i++) {
			fprintf(stderr, "log_fact[%d]: %.5lg, %.5lg\n", i, mm->log_fact[i], exp(mm->log_fact[i]));
		}
	}

	return mm;
}

void nilMultiMotifs(MultiMotifs *multiMotifs) {
	free(multiMotifs->log_fact);
	//free(multiMotifs->log_int);

	MotifNode *garbage;
	while(multiMotifs->head != NULL) {
		garbage = multiMotifs->head;
		multiMotifs->head = multiMotifs->head->next;
		nilMotifNode(garbage);
	}
	nilProfile(multiMotifs->revcomplMat);
	free(multiMotifs);
}

//-----------------------------------------------------------------------------
// Distance
//-----------------------------------------------------------------------------

//compute log(a+b) given log(a) and log(b)
static inline
double logsum(double log_a, double log_b) {
	return (log_a < log_b ? 
		log_b + log((double)(1.0 + exp(log_a - log_b))) : log_a + log((double)(1.0 + exp(log_b - log_a))) 
		);
}

static
void calculateNtLProb(int count[NUMALPHAS], double lp[NUMALPHAS]) {
	const int NUM_COMPS = 5; //5-component mixture
	double hyperParams[NUM_COMPS][NUMALPHAS];
	for(int k = 0; k < 4; k++) {
		for(int i = 0; i < NUMALPHAS; i++) {
			hyperParams[k][i] = ( k == i ? 5 : 1); //(5,1,1,1), (1,5,1,1), etc
		}
	}
	for(int i = 0; i < NUMALPHAS; i++) {
		hyperParams[4][i] = 2; //(2,2,2,2)
	}

	if(DEBUG0) {
		for(int k = 0; k < NUM_COMPS; k++) {
			double sum = 0.0; 
			for(int i = 0; i < NUMALPHAS; i++) {
				sum += hyperParams[k][i];
			}
			if(fabs(sum - 8.0) > 0.00001) {
				fprintf(stderr, "The sum of the hyperparameters should be the same\n");
				exit(1);
			}
		}
	}
	
	//see "A Novel Bayesian DNA Motif..." by Habib, etal Nir Friedman. (p.13)
	double log_pr_n_given_k[NUM_COMPS];
	double sum_log_pr_n_given_k = NINF;
	for(int k = 0; k < NUM_COMPS; k++) {
		double log_prod_term = 0.0;
		for(int i = 0; i < NUMALPHAS; i++) {
			log_prod_term += lgamma(count[i] + hyperParams[k][i]) - lgamma(hyperParams[k][i]);
		}
		log_pr_n_given_k[k] = log_prod_term; //the first term in the paper is cancelled out in Pr(k;n)
		sum_log_pr_n_given_k = logsum(sum_log_pr_n_given_k, log_pr_n_given_k[k]);
	}
	double pr_k_given_n[NUM_COMPS];
	for(int k = 0; k < NUM_COMPS; k++) {
		//Pr(k)/Pr(j) is cancelled out
		pr_k_given_n[k] = exp(log_pr_n_given_k[k] - sum_log_pr_n_given_k);
	}
	if(DEBUG0) {
		double sum = 0.0;
		for(int k = 0; k < NUM_COMPS; k++) {
			sum += pr_k_given_n[k];
		}
		if(fabs(sum - 1.0) > 0.00000001) {
			fprintf(stderr, "Error: Pr(k;n) does not sum up to 1.0 over k\n");
			exit(1);
		}
	}
	//if(DEBUG1) {
	//	for(int k = 0; k < NUM_COMPS; k++) {
	//		fprintf(stderr, "Pr(k = %d ; n) = %.5lg\n", k, pr_k_given_n[k]);
	//	}
	//}

	//finally use the pseudocount alpha
	for(int i = 0; i < NUMALPHAS; i++) {
		double prob = 0.0;
		for(int k = 0; k < NUM_COMPS; k++) {
			double count_sum = 0.0;
			for(int j = 0; j < NUMALPHAS; j++) {
				count_sum += count[j] + hyperParams[k][j];
			}
			prob += (count[i] + hyperParams[k][i]) / count_sum * pr_k_given_n[k];
		}
		lp[i] = log((double)prob);
	}

	if(DEBUG0) {
		double sum = 0.0;
		for(int i = 0; i < NUMALPHAS; i++) {
			sum += exp(lp[i]);
		}
		if(fabs(sum - 1.0) > 0.000000000001) {
			fprintf(stderr, "Error: sum of nucleotide probabilites does not equal 1.0\n");
			for(int i = 0; i < NUMALPHAS; i++) {
				fprintf(stderr, "%c: %.5lg\n", numToChar(i), exp(lp[i]));
			}
			fprintf(stderr, "sum=%.8lg\n", sum);
			exit(1);
		}
	}

}
static
double calculateBinomLProb(int count[NUMALPHAS], double lp[NUMALPHAS], double *log_fact) {
	int numsites = 0;
	for(int i = 0; i < NUMALPHAS; i++) {
		numsites += count[i];
	}

	//double sum = log_fact[numsites];
	//for(int i = 0; i < NUMALPHAS; i++) {
	//	if(count[i] > 0) {
	//		sum += count[i] * lp[i] - log_fact[count[i]];
	//	}
	//}
	double sum = 0.0;
	for(int i = 0; i < NUMALPHAS; i++) {
		if(count[i] > 0) {
			sum += count[i] * lp[i];
		}
	}
	return sum;
}

static
double computeColumnDist(int count1[NUMALPHAS], int count2[NUMALPHAS], double lprob1, double lprob2, MultiMotifs *multiMotifs) {
	int countJ[NUMALPHAS]; //joint count
	for(int a = 0; a < NUMALPHAS; a++) {
		countJ[a] = count1[a] + count2[a];
	}

	//if(DEBUG1) {
	//	//print count1, count2, countJ
	//	fprintf(stderr, "Nucleotides count:\n");
	//	for(int a = 0; a < NUMALPHAS; a++) {
	//		fprintf(stderr, "%c: %4d %4d %4d\n", numToChar(a), count1[a], count2[a], countJ[a]);
	//	}
	//}

	double ntLProbJ[NUMALPHAS];//joint prob
	calculateNtLProb(countJ, ntLProbJ);
	double lprobJ = calculateBinomLProb(countJ, ntLProbJ, multiMotifs->log_fact);

	//if(DEBUG1) {
	//	fprintf(stderr, "Nucleotide prob after Dirichlet prior:\n");
	//	for(int a = 0; a < NUMALPHAS; a++) {
	//		fprintf(stderr, "%c: %.5lf %.5lf %.5lf\n", numToChar(a), exp(ntLProb1[a]), exp(ntLProb2[a]), exp(ntLProbJ[a]));
	//	}
	//	fprintf(stderr, "lprob1: %.5lg\n", lprob1);
	//	fprintf(stderr, "lprob2: %.5lg\n", lprob2);
	//	fprintf(stderr, "lprobJ: %.5lg\n", lprobJ);
	//	fprintf(stderr, "\n");
	//}

	double dist = lprob1 + lprob2 - lprobJ;
	return dist;
}

static
double computeAlignedDist(Profile *cmat1, Profile *cmat2, double lprob1[SPAN_CAPACITY], double lprob2[SPAN_CAPACITY], 
						  int pos1, int pos2, int alignlen, MultiMotifs *multiMotifs) 
{
	double distVect[SPAN_CAPACITY];
	for(int i = 0; i < alignlen; i++) {
		int count1[NUMALPHAS];
		int count2[NUMALPHAS];
			
		for(int a = 0; a < NUMALPHAS; a++) {
			count1[a] = (int) cmat1->mat[pos1 + i][a];
			count2[a] = (int) cmat2->mat[pos2 + i][a];
		}
		distVect[i] = computeColumnDist(count1, count2, lprob1[pos1+i], lprob2[pos2+i], multiMotifs);
	}
	if(DEBUG1) {
		for(int k = 0; k < alignlen; k++) {
			fprintf(stderr, "[%3d] ", k);
		}
		fprintf(stderr, "\n");
		for(int a = 0; a < NUMALPHAS; a++) {
			for(int k = 0; k < alignlen; k++) {
				fprintf(stderr, "%5.0lf ", cmat1->mat[k+pos1][a]);
			}
			fprintf(stderr, "\n");
		}
		fprintf(stderr, "\n");
		for(int a = 0; a < NUMALPHAS; a++) {
			for(int k = 0; k < alignlen; k++) {
				fprintf(stderr, "%5.0lf ", cmat2->mat[k+pos2][a]);
			}
			fprintf(stderr, "\n");
		}

		for(int i = 0; i < alignlen; i++) {
			fprintf(stderr, "%3d: %6.2lf\n", i, distVect[i]);
		}
	}

	//take geometric mean of p-value (arithmetic mean of log(pval))
	double sum = 0.0;
	for(int i = 0; i < alignlen; i++) {
		sum += distVect[i];
	}
	double mean = sum / alignlen;

	if(DEBUG1) {
		fprintf(stderr, "Distance mean: %.8lg\n", mean);
	}

	return mean;
}


static
double computeMotifDistForward(Profile *cmat1, Profile *cmat2, MultiMotifs *multiMotifs) {
	double leastDist = DBL_MAX;
	
	double lprob1[SPAN_CAPACITY];
	double lprob2[SPAN_CAPACITY];

	//precalculate the log-prob 
	for(int i = 0; i < cmat1->span; i++) {
		int count[NUMALPHAS];
		for(int a = 0; a < NUMALPHAS; a++) {
			count[a] = (int) cmat1->mat[i][a];
		}
		double ntLProb[NUMALPHAS];
		calculateNtLProb(count, ntLProb);
		lprob1[i] = calculateBinomLProb(count, ntLProb, multiMotifs->log_fact);//binomial probability
	}
	for(int i = 0; i < cmat2->span; i++) {
		int count[NUMALPHAS];
		for(int a = 0; a < NUMALPHAS; a++) {
			count[a] = (int) cmat2->mat[i][a];
		}
		double ntLProb[NUMALPHAS];
		calculateNtLProb(count, ntLProb);
		lprob2[i] = calculateBinomLProb(count, ntLProb, multiMotifs->log_fact);//binomial probability
	}

	int smallerSpan = (cmat1->span < cmat2->span ? cmat1->span : cmat2->span);
	for(int i = 0; i < cmat1->span; i++) {
		for(int j = 0; j < cmat2->span; j++) {
			//one of the aligned seq must start with 0, or the alignment is extendable.
			if(i == 0 || j == 0) {
				//compute align length
				int alignlen = -1;
				if(cmat1->span - i < cmat2->span - j) {
					alignlen = cmat1->span - i;
				}
				else {
					alignlen = cmat2->span - j;
				}

				if( (alignlen/(double)smallerSpan) >= multiMotifs->minOverlapRatio) {
					double dist = computeAlignedDist(cmat1, cmat2, lprob1, lprob2, i, j, alignlen, multiMotifs);
					if(dist < leastDist) {
						leastDist = dist;
					}
				}
			}
		}
	}
	return leastDist;

}


double computeMotifDist(Profile *countmat1, Profile *countmat2, MultiMotifs *multiMotifs) {
	if(DEBUG1) {
		printCountmat(stderr, countmat1);
		printCountmat(stderr, countmat2);
	}

	//create reverse-complement version of countmat1
	copyProfile(multiMotifs->revcomplMat, countmat1);
	revcomplProfile(multiMotifs->revcomplMat);

	double forward = computeMotifDistForward(countmat1, countmat2, multiMotifs);
	double revcompl = computeMotifDistForward(multiMotifs->revcomplMat, countmat2, multiMotifs);

	double dist = (forward < revcompl ? forward : revcompl); //min()

	if(DEBUG1) {
		printCountmat(stderr, countmat1);
		printCountmat(stderr, countmat2);
		fprintf(stderr, "Distance: %.15lf\n", dist);
		fprintf(stderr, "==========================\n");
	}

	return dist;
}

//-----------------------------------------------------------------------------
// compute distance cutoff
//-----------------------------------------------------------------------------

//change the matrix itself
static
void randPermuteProfile(Profile *profile) {
	if(profile->cols != profile->span) {
		fprintf(stderr, "Error: randPermuteProfile() is only implemented for ungapped profile.\n");
		exit(1);
	}

	for(int i = 0; i < profile->span; i++) {
		int range = profile->span - i;
		int rand = ((int)(Random() * range)) + i;
		
		//swap columns btw i and rand
		double *temp = profile->mat[i];
		profile->mat[i] = profile->mat[rand];
		profile->mat[rand] = temp;
	}
}

static
void swapArrayElements(double *array, int index1, int index2) {
	double temp = array[index1];
	array[index1] = array[index2];
	array[index2] = temp;
}

static
int partitionByPivot(double *array, int pivotIndex, int left, int right) {
	swapArrayElements(array, pivotIndex, left); 

	//"right-1" now has the pivot value, pivotIndex is useless
	double pivot = array[right-1];
	int numLtPivot = 0; //number of elements less than pivot within left to right
	for(int i = left; i < right - 1; i++) { //right-1 is excluded
		if(array[i] < pivot) {
			//invariant: array[left+numLtPivot] >= pivot
			swapArrayElements(array, left+numLtPivot, i);
			numLtPivot++;
		}
	}
	swapArrayElements(array, right-1, left + numLtPivot);
	int newPivotIndex = left + numLtPivot;

	if(DEBUG0) {
		assert(fabs(array[newPivotIndex] - pivot) < 0.000001);
		for(int i = left; i < newPivotIndex; i++) {
			assert(array[i] < pivot);
		}
		for(int i = newPivotIndex + 1; i < right; i++) {
			assert(array[i] >= pivot);
		}
	}
	return newPivotIndex; //new pivot index
}

//right is excluded, thus length is right - left
static
double selectK(double *array, int k, int left, int right) {
	while(true) {
		int pivotIndex = (int)(Random() * (right - left)) + left;
		int newPivotIndex = partitionByPivot(array, pivotIndex, left, right);
		if(newPivotIndex == k) {
			return array[newPivotIndex];
		}
		else if(newPivotIndex < k) {
			left = newPivotIndex + 1;
		}
		else {
			right = newPivotIndex; //because right is excluded from the selection
		}
	}
}

void setDistCutoffForMultiMotifs(MultiMotifs *multiMotifs, RunSet *runset, int numPermPairs, double pvalCutoff) {
	RunNode **nodeArray = (RunNode**) malloc(sizeof(RunNode*) * runset->len);
	int numNodes = 0; 
	for(RunNode *node = runset->head; node != NULL; node = node->next) {
		nodeArray[numNodes++] = node;
	}
	if(numNodes != runset->len) {
		fprintf(stderr, "Error: Inconsistent runset->len!\n");
		exit(1);

	}

	Profile *countmat1 = initProfile(multiMotifs->maxspan, multiMotifs->maxspan);
	Profile *countmat2 = initProfile(multiMotifs->maxspan, multiMotifs->maxspan);

	double *distArray = (double*) malloc(sizeof(double) * numPermPairs);

	for(int i = 0; i < numPermPairs; i++) {
		RunNode* node1 = nodeArray[(int)(Random() * numNodes)];
		RunNode* node2 = nodeArray[(int)(Random() * numNodes)];

		copyProfile(countmat1, node1->countmat);
		copyProfile(countmat2, node2->countmat);

		randPermuteProfile(countmat1);
		randPermuteProfile(countmat2);

		distArray[i] = computeMotifDist(countmat1, countmat2, multiMotifs);
	}

	int indexOfCutoff = (int)(numPermPairs * (1.0 - pvalCutoff)); //floor()
	multiMotifs->distCutoff = selectK(distArray, indexOfCutoff, 0, numPermPairs);

	free(distArray);
	free(nodeArray); //DO NOT deallocate the elements inside the array (they should be retained for runset)

	nilProfile(countmat1);
	nilProfile(countmat2);

	if(DEBUG1) {
		fprintf(stderr, "Distance cutoff value: %.5lg\n", multiMotifs->distCutoff);
		fprintf(stderr, "Index of cutoff= %d\n", indexOfCutoff);
		fprintf(stderr, "Number of permutation pairs= %d\n", numPermPairs);
		fprintf(stderr, "P-value cutoff= %.5lg\n", pvalCutoff);
	}
}

