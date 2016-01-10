#include "multi_motifs.h"

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

MultiMotifs* initMultiMotifs(int capacity, int numseqs, double ntFreq[]) {
	MultiMotifs *mm = (MultiMotifs*) malloc(sizeof(MultiMotifs));
	mm->numseqs = numseqs;
	mm->capacity = capacity;
	mm->minOverlapRatio = 0.50; //minimum overlap ratio to be considered as alignment
	mm->size = 0;
	mm->head = NULL;
	for(int i = 0; i < NUMALPHAS; i++) {
		mm->ntFreq[i] = ntFreq[i];
	}

	mm->log_fact = (double*) malloc(sizeof(double) * (numseqs+1));
	mm->log_fact[0] = log(1.0); // 0.0
	for(int i = 1; i < numseqs + 1; i++) {
		mm->log_fact[i] = mm->log_fact[i-1] + log((double) i );
	}

	mm->log_int = (double*) malloc(sizeof(double) * (numseqs+1));
	mm->log_int[0] = LOGZERO; //we want to have 0*log(0) = 0 [-INFINITY * 0 gives NAN]
	for(int i = 1; i < numseqs + 1; i++) {
		mm->log_int[i] = log((double) i );
	}

	
	if(DEBUG1) {
		for(int i = 1; i < numseqs + 1; i++) {
			fprintf(stderr, "log_fact[%d]: %.5lg, %.5lg\n", i, mm->log_fact[i], exp(mm->log_fact[i]));
		}
	}

	return mm;
}

void nilMultiMotifs(MultiMotifs *multiMotifs) {
	free(multiMotifs->log_fact);
	free(multiMotifs->log_int);

	MotifNode *garbage;
	while(multiMotifs->head != NULL) {
		garbage = multiMotifs->head;
		multiMotifs->head = multiMotifs->head->next;
		nilMotifNode(garbage);
	}
	free(multiMotifs);
}

//-----------------------------------------------------------------------------
// P-values
//-----------------------------------------------------------------------------

//compute log(a+b) given log(a) and log(b)
static inline
double logsum(double log_a, double log_b) {
	return (log_a < log_b ? 
		log_b + log(1.0 + exp(log_a - log_b)) : log_a + log(1.0 + exp(log_b - log_a)) 
		);
}

//The score is a the entropy + N*log(N) [Missing the N*log(N) term]
static 
double computeColumnScore(int fgcol[NUMALPHAS], double lbgfreq[NUMALPHAS], MultiMotifs *multiMotifs ) {
	double score = 0.0;
	for(int a = 0; a < NUMALPHAS; a++) {
		if(fgcol[a] >= 1) {
			score += fgcol[a] * (multiMotifs->log_int[fgcol[a]] - lbgfreq[a]);
		}
	}
	return score;
}

static 
double computeColumnLogProb(int fgcol[NUMALPHAS], double lbgfreq[NUMALPHAS], MultiMotifs *multiMotifs) {
	int numsites = 0;
	for(int i = 0; i < NUMALPHAS; i++) {
		numsites += fgcol[i];
	}

	double lprob = multiMotifs->log_fact[numsites]; //factorial for multinomial term
	for(int i = 0; i < NUMALPHAS; i++) {
		if(fgcol[i] >= 1) {
			lprob += fgcol[i] * lbgfreq[i];
		}
		lprob -= multiMotifs->log_fact[fgcol[i]]; //factorial for multinomial term
	}
	return lprob;
}


static
void computeMultiColumnSimilarityLPval(Profile *fgmat, int bgcount[], bool needDist[], double distVect[], MultiMotifs *multiMotifs) {
	double lbgfreq[NUMALPHAS];
	int bgNumsites = 0; 
	for(int i = 0; i < NUMALPHAS; i++) {
		bgNumsites += bgcount[i];
	}
	double sum = 0.0;
	for(int i = 0; i < NUMALPHAS; i++) {
		double pseudocountTotal = 0.10 * bgNumsites; //pseudocount is needed for background
		lbgfreq[i] = bgcount[i] + pseudocountTotal * multiMotifs->ntFreq[i];
		sum += lbgfreq[i];
	}
	for(int i = 0; i < NUMALPHAS; i++) {
		lbgfreq[i] = log((double)lbgfreq[i]) - log((double)sum);
	}

	double lntFreq[NUMALPHAS];
	for(int a = 0; a< NUMALPHAS; a++) {
		lntFreq[a] = log(multiMotifs->ntFreq[a]);
	}

	int numsites = 0;
	for(int a = 0; a < NUMALPHAS; a++) {
		numsites += (int)fgmat->mat[0][a];
	}
	if(DEBUG0) {
		for(int m = 0; m < fgmat->span; m++) {
			int ns = 0;
			for(int a = 0; a < NUMALPHAS; a++) {
				ns += (int)fgmat->mat[m][a];
			}
			assert(ns == numsites);
		}
	}

	int fgcountmat[SPAN_CAPACITY][NUMALPHAS];
	for(int i = 0; i < fgmat->span; i++) {
		for(int j = 0; j < NUMALPHAS; j++) {
			fgcountmat[i][j] = (int)fgmat->mat[i][j];
		}
	}

	double fgscore[SPAN_CAPACITY];
	double minFgscore = DBL_MAX;
	for(int i = 0; i < fgmat->span; i++) {
		if(needDist[i]) {
			fgscore[i] = computeColumnScore(fgcountmat[i], lbgfreq, multiMotifs);
			if(fgscore[i] < minFgscore) {
				minFgscore = fgscore[i];
			}
		}
		else {
			fgscore[i] = NAN;
		}
	}

	if(DEBUG1) {
		//fprintf(stderr, "fgscore: %.5lg\n", fgscore);
		//fprintf(stderr, "bgfreq: ");
		//for(int i = 0; i < NUMALPHAS; i++) {
		//	fprintf(stderr, "%.4lf ", exp(lbgfreq[i]));
		//}
		//fprintf(stderr, "\n");
	}

	//initialize distMat for correct entries
	int indicesOfNeedDist[SPAN_CAPACITY];
	int numOfNeedDist = 0;
	for(int i = 0; i < fgmat->span; i++) {
		if(needDist[i]) {
			indicesOfNeedDist[numOfNeedDist++] = i;
			distVect[i] = NINF;
		}
		else {
			distVect[i] = NAN;
		}
	}

	double slack = 0.0001;
	double lprob_dp[NUMALPHAS];
	double score_dp[NUMALPHAS];
	lprob_dp[0] = multiMotifs->log_fact[numsites];
	score_dp[0] = 0.0;
	double lprobAll = NINF;

	//Note the less than and equal signs for these cases.
	for(int a = 0; a <= numsites; a++) {
		//lprob_dp[1] = lprob_dp[0] + a * lbgfreq[0] - multiMotifs->log_fact[a] ;
		lprob_dp[1] = lprob_dp[0] + a * lntFreq[0] - multiMotifs->log_fact[a] ;

		//log_int[0] = LOGZERO -- thus 0 * log(0) is 0. Hence, if-else case is not needed.
		score_dp[1] = score_dp[0] + a * (multiMotifs->log_int[a] - lbgfreq[0]);

		for(int b = 0; b <= numsites - a; b++) {
			//lprob_dp[2] = lprob_dp[1] + b * lbgfreq[1] - multiMotifs->log_fact[b] ;
			lprob_dp[2] = lprob_dp[1] + b * lntFreq[1] - multiMotifs->log_fact[b] ;
			score_dp[2] = score_dp[1] + b * (multiMotifs->log_int[b] - lbgfreq[1]);

			for(int c = 0; c <= numsites - a - b; c++) {
				int d = numsites - a - b - c;

				double score = score_dp[2]
					+ c * (multiMotifs->log_int[c] - lbgfreq[2])
					+ d * (multiMotifs->log_int[d] - lbgfreq[3]);

				if(score > minFgscore - slack) {
					//double lprob = lprob_dp[2]
					//	+ c * lbgfreq[2] - multiMotifs->log_fact[c] 
					//	+ d * lbgfreq[3] - multiMotifs->log_fact[d] ;
					double lprob = lprob_dp[2]
						+ c * lntFreq[2] - multiMotifs->log_fact[c] 
						+ d * lntFreq[3] - multiMotifs->log_fact[d] ;
					//There is an error term to avoid having cumlprob to be -INFINITY (includes fgcol[][])
					for(int i = 0; i < numOfNeedDist; i++) {
						if(score > fgscore[indicesOfNeedDist[i]] - slack) { 
							distVect[indicesOfNeedDist[i]] = logsum(distVect[indicesOfNeedDist[i]], lprob);
						}
					}
				}

				if(DEBUG0) {
					double lprob = lprob_dp[2]
						+ c * lntFreq[2] - multiMotifs->log_fact[c] 
						+ d * lntFreq[3] - multiMotifs->log_fact[d] ;

					lprobAll = logsum(lprobAll, lprob);

					int col[] = {a, b, c, d};
					double s1 = computeColumnScore(col, lbgfreq, multiMotifs);
					if(fabs(s1 - score) > 0.0000001) {
						fprintf(stderr, "Incorrect column-score computed using DP.\n");
						fprintf(stderr, "DP: %.2lf, direct: %.2lf\n", score, s1);
						exit(1);
					}
				}

				//if(DEBUG1) {
				//	//fprintf(stderr, "score=%.5lg, lprob=%.5lg, cumlprob=%.5lg\n", score, lprob, cumlprob);
				//}
			}
		}
	}
	if(DEBUG0) {
		if(fabs(exp(lprobAll) - 1.0) > 0.00000001) {
			fprintf(stderr, "Inconsistent pmf with log of sum of all prob %.6lf\n", lprobAll);
			abort();
		}
	}
	if(DEBUG1) {
		for(int i = 0; i < fgmat->span; i++) {
			if(needDist[i]) {
				fprintf(stderr, "fgcol: ");
				for(int a = 0; a < NUMALPHAS; a++) {
					fprintf(stderr, " %5d ", fgcountmat[i][a]);
				}
				fprintf(stderr, "\n");
				fprintf(stderr, "bgcol: ");
				for(int a = 0; a < NUMALPHAS; a++) {
					fprintf(stderr, " %5d ", bgcount[a]);
				}
				fprintf(stderr, "\n");
				fprintf(stderr, "fgscore[%d]: %.8lg\n", i, fgscore[i]);
			}
		}
	}

}

//static
//void computeMultiColumnSimilarityLPval(Profile *fgmat, int bgcount[], bool needDist[], double distVect[], MultiMotifs *multiMotifs) {
//	for(int i = 0; i < fgmat->span; i++) {
//		if(needDist[i]) {
//			int fgcount[NUMALPHAS];
//			for(int a = 0; a < NUMALPHAS; a++) {
//				fgcount[a] = (int)(fgmat->mat[i][a]);
//			}
//			distVect[i] = computeColumnSimilarityLPval(fgcount, bgcount, multiMotifs);
//		}
//	}
//}

static
double computeAlignedSimilarityLPval(double distMat[SPAN_CAPACITY][SPAN_CAPACITY], int fgind, int bgind, int alignlen, MultiMotifs *multiMotifs) {
	double lpvalArr[SPAN_CAPACITY];
	for(int i = 0; i < alignlen; i++) {
		//distMat[bg][fg]
		lpvalArr[i] = distMat[bgind+i][fgind+i];

		if(DEBUG0) {
			if(isnan(distMat[bgind+i][fgind+i])) {
				fprintf(stderr, "Invalid fgind and bgind for distMat[][]\n");
				exit(1);
			}
		}
	}
	if(DEBUG1) {
		for(int i = 0; i < alignlen; i++) {
			fprintf(stderr, "%3d: %6.2lf\n", i, distMat[bgind+i][fgind+i]);
		}
	}

	//take geometric mean of p-value (arithmetic mean of log(pval))
	double sum = 0.0;
	for(int i = 0; i < alignlen; i++) {
		sum += lpvalArr[i];
	}
	double lpval_mean = sum / alignlen;

	if(DEBUG1) {
		fprintf(stderr, "lpval_mean: %.8lg\n", lpval_mean);
	}

	return lpval_mean;
}

static
double computeMotifSimilarityLPval_forward(Profile *fgmat, Profile *bgmat, MultiMotifs *multiMotifs) {
	int minspan = (fgmat->span < bgmat->span ? fgmat->span : bgmat->span); //take min()

	bool needDist[SPAN_CAPACITY][SPAN_CAPACITY]; //[bgmat->span][fgmat->span]
	double distMat[SPAN_CAPACITY][SPAN_CAPACITY]; //[bgmat->span][fgmat->span]

	for(int i = 0; i < SPAN_CAPACITY; i++) {
		for(int j = 0; j < SPAN_CAPACITY; j++) {
			needDist[i][j] = false;
			distMat[i][j] = NAN;
		}
	}

	//decide which distance we need to compute
	for(int bpos = 0; bpos < bgmat->span; bpos++) {
		for(int fpos = 0; fpos < fgmat->span; fpos++) {
			//one of the aligned seq must start with 0, or the alignment is extendable.
			if(bpos == 0 || fpos == 0) {
				//compute align length
				int alignlen = -1;
				if(bgmat->span - bpos < fgmat->span - fpos) {
					alignlen = bgmat->span - bpos;
				}
				else {
					alignlen = fgmat->span - fpos;
				}

				if( (alignlen/(double)minspan) >= multiMotifs->minOverlapRatio) {
					for(int k = 0; k < alignlen; k++) {
						needDist[bpos+k][fpos+k] = true;
					}
				}
			}
		}
	}

	for(int i = 0; i < bgmat->span; i++) {
		int bgcount[NUMALPHAS];
		for(int a = 0; a < NUMALPHAS; a++) {
			bgcount[a] = (int)(bgmat->mat[i][a]);
		}
		computeMultiColumnSimilarityLPval(fgmat, bgcount, needDist[i], distMat[i], multiMotifs);
	}

	//find most similar for all alignments
	double lpvalMostSimilar = -DBL_MAX;
	for(int bpos = 0; bpos < bgmat->span; bpos++) {
		for(int fpos = 0; fpos < fgmat->span; fpos++) {
			//one of the aligned seq must start with 0, or the alignment is extendable.
			if(bpos == 0 || fpos == 0) {
				//compute align length
				int alignlen = -1;
				if(bgmat->span - bpos < fgmat->span - fpos) {
					alignlen = bgmat->span - bpos;
				}
				else {
					alignlen = fgmat->span - fpos;
				}



				if( (alignlen/(double)minspan) >= multiMotifs->minOverlapRatio) {
					if(DEBUG1) {
						for(int k = 0; k < alignlen; k++) {
							fprintf(stderr, "[%3d] ", k);
						}
						fprintf(stderr, "\n");
						for(int a = 0; a < NUMALPHAS; a++) {
							for(int k = 0; k < alignlen; k++) {
								fprintf(stderr, "%5.0lf ", fgmat->mat[k+fpos][a]);
							}
							fprintf(stderr, "\n");
						}
						fprintf(stderr, "\n");
						for(int a = 0; a < NUMALPHAS; a++) {
							for(int k = 0; k < alignlen; k++) {
								fprintf(stderr, "%5.0lf ", bgmat->mat[k+bpos][a]);
							}
							fprintf(stderr, "\n");
						}
					}

					double lpval = computeAlignedSimilarityLPval(distMat, fpos, bpos, alignlen, multiMotifs);
					if(lpval > lpvalMostSimilar) {
						lpvalMostSimilar = lpval;
					}
				}
			}
		}
	}
	return lpvalMostSimilar;

}

static
double computeMotifSimilarityLPval(Profile *fgmat, Profile *bgmat, MultiMotifs *multiMotifs) {
	Profile *fgmatRc = initProfile(fgmat->maxspan, fgmat->maxspan);
	copyProfile(fgmatRc, fgmat);
	revcomplProfile(fgmatRc);

	double lpvalF = computeMotifSimilarityLPval_forward(fgmat, bgmat, multiMotifs);
	double lpvalRc= computeMotifSimilarityLPval_forward(fgmatRc, bgmat, multiMotifs);

	double lpval = NAN;
	if(lpvalF > lpvalRc) { //take max()
		lpval = lpvalF;
	}
	else {
		lpval = lpvalRc;
	}

	nilProfile(fgmatRc);
	return lpval;
}
//-----------------------------------------------------------------------------
// Core functions
//-----------------------------------------------------------------------------

double computeMotifDistance(Profile *countmat1, Profile *countmat2, MultiMotifs *multiMotifs) {
	if(DEBUG1) {
		printCountmat(stderr, countmat1);
		printCountmat(stderr, countmat2);
	}

	//H0: the two motifs are the same
	double lpval1 = computeMotifSimilarityLPval(countmat1, countmat2, multiMotifs);
	//double lpval2 = computeMotifSimilarityLPval(countmat2, countmat1, multiMotifs);
	double lpval2 = lpval1; //temporary hack

	double mean = (lpval1 + lpval2) / 2.0; //geometric mean of pval is the arithmetic mean in log()
	//lower p-value => motifs are different => larger distance
	double dist = (1.0-exp(mean));

	if(DEBUG1) {
		printCountmat(stderr, countmat1);
		printCountmat(stderr, countmat2);
		fprintf(stderr, "lpval1: %.10lg, lpval2: %.10lg\n", lpval1, lpval2);
		fprintf(stderr, "Mean of two lpvals: %.10lg\n", mean);
		fprintf(stderr, "Distance: %.15lf\n", dist);
		fprintf(stderr, "==========================\n");
	}

	return dist;
}
