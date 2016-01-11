#include "wndbgm.h"

TransProb* chooseModelForWndbgm(Wndbgm *wndbgm, int countA, int countC, int countG, int countT) {
	if(wndbgm->maxorder == 0) {
		fprintf(stderr, "Error: chooseModel should not be called when maxorder == 0\n");
		exit(1);
	}

	int total = countA + countC + countG + countT;
	double gcPercent = ((double)countC + countG) / total * 100.0;

	int modelIndex = wndbgm->numModels - 1; //case when it is the last bin (between last location and 100%)
	for(int i = 1; i < wndbgm->numModels; i++) {
		if(gcPercent < wndbgm->gcBinLocations[i]) {
			modelIndex = i - 1;
			break;
		}
	}

	if(DEBUG1) {
		fprintf(stderr, "GC percent: %.2lf, model-index: %d\n", gcPercent, modelIndex);
	}

	return wndbgm->models[modelIndex];
}

void nilWndbgm(Wndbgm *wndbgm) {
	if(wndbgm->numModels > 0) {
		free(wndbgm->gcBinLocations);
		for(int i = 0; i < wndbgm->numModels; i++) {
			free(wndbgm->models[i]);
		}
		free(wndbgm->models);
	}
	free(wndbgm);
}

//--------------------------------------------------------------------------------------------
// Initialization and opening file that contains Markov-chain
//
//--------------------------------------------------------------------------------------------

static
int minInt(int int1, int int2) {
	return (int1 < int2 ? int1 : int2);
}

static
int findMinValInTransCount(TransCount *transcount){
	int minval = INT_MAX;

	for(int a = 0; a < NUMALPHAS && transcount->maxorder >= 0; a++) {
		minval = minInt(minval, transcount->count0[a]);
		for(int b = 0; b < NUMALPHAS && transcount->maxorder >= 1; b++) {
			minval = minInt(minval, transcount->count1[a][b]);
			for(int c = 0; c < NUMALPHAS && transcount->maxorder >= 2; c++) {
				minval = minInt(minval, transcount->count2[a][b][c]);
				for(int d = 0; d < NUMALPHAS && transcount->maxorder >= 3; d++) {
					minval = minInt(minval, transcount->count3[a][b][c][d]);
					for(int e = 0; e < NUMALPHAS && transcount->maxorder >= 4; e++) {
						minval = minInt(minval, transcount->count4[a][b][c][d][e]);
						for(int f = 0; f < NUMALPHAS && transcount->maxorder >= 5; f++) {
							minval = minInt(minval, transcount->count5[a][b][c][d][e][f]);
						}
					}
				}
			}
		}
	}
	return minval;
}


static
FILE *openFile(char *filename) {
	FILE *fptr;
	if((fptr = fopen(filename,"r")) == NULL) {
		fprintf(stderr,"Could not open file \"%s\"\n",filename);
		exit(1);
	}
	return fptr;
}

static 
int readNumOfModels(FILE *fptr) {
	rewind(fptr);
	char c;
	int numModels = 0;
	while ((c = fgetc(fptr)) != EOF) {
		if(c == '>') {
			numModels++;
		}
	}
	return numModels;
}

static
void readBinLocations(FILE *fptr, int *binLocations, int numModels) {
	rewind(fptr);
	char c;
	int loc = -1;
	int index = 0;
	bool readLoc = false;
	while ((c = fgetc(fptr)) != EOF) {
		if(c == '>') {
			readLoc = true;
			loc = 0;
		}
		else if(readLoc) {
			if(c == '\n' || c== '\r') {
				if(loc < 0) {
					fprintf(stderr, "Error: invalid loc found at index %d: %d\n", index, loc);
					exit(1);
				}
				readLoc = false;
				binLocations[index++] = loc;
				loc = -1;
			}
			else if('0' <= c && c <= '9') {
				loc = loc * 10 + (c - '0');
			}
			else {
				fprintf(stderr, "Error: invalid char found at index %d: %c\n", index, c);
				exit(1);				
			}
		}
	}
	if(index != numModels) {
		fprintf(stderr, "Error: invalid number of models found: %d and %d\n", numModels, index);
		exit(1);				
	}

	//assert the bin locations are ascending
	for(int i = 1; i < numModels; i++) {
		if(binLocations[i-1] >= binLocations[i]) {
			fprintf(stderr, "Error: bin locations must be strictly ascending\n");
			exit(1);
		}
	}
}

static
void setTransCount(TransCount *transcount, int nts_read, int nts[MARKOV_ORDER_BOUND+1], int value) {
	if(DEBUG0) {
		assert(nts_read > 0);
	}

	int order = nts_read - 1;
	if(DEBUG0) {
		assert(order <= transcount->maxorder);
	}
	if(order == 0) {
		transcount->count0[nts[0]] = value;
	}
	else if(order == 1) {
		transcount->count1[nts[0]][nts[1]] = value;
	}
	else if(order == 2) {
		transcount->count2[nts[0]][nts[1]][nts[2]] = value;
	}
	else if(order == 3) {
		transcount->count3[nts[0]][nts[1]][nts[2]][nts[3]] = value;
	}
	else if(order == 4) {
		transcount->count4[nts[0]][nts[1]][nts[2]][nts[3]][nts[4]] = value;
	}
	else if(order == 5) {
		transcount->count5[nts[0]][nts[1]][nts[2]][nts[3]][nts[4]][nts[5]] = value;
	}
	else {
		fprintf(stderr, "Error: invalid order found: %d\n", order);
		exit(1);
	}
}


static
void readMarkovModels(FILE *fptr, TransCount **rawModels, int numModels, int maxorder) {
	for(int i = 0; i < numModels; i++) {
		rawModels[i]->maxorder = maxorder;
		resetTransCount(rawModels[i], 0);
	}

	rewind(fptr);
	char c;
	int index = -1;
	bool inLoc = false;
	//bool readModels = false;

	enum ReadModeEnum {NONE, KMER, KMER_VAL};
	enum ReadModeEnum readMode = NONE;

	int nts[MARKOV_ORDER_BOUND + 1];
	int nts_read = 0; // number of nucleotides read

	int value = 0;
	bool isReadingValue = false;

	while ((c = fgetc(fptr)) != EOF) {
		if(c == '>') {
			inLoc = true;
			readMode = NONE;
			index++;
		}
		else if(inLoc) {
			if(c == '\n' || c== '\r') {
				inLoc = false;
				readMode = KMER;
				nts_read = 0;
				value = 0;
			}
		}
		else if(readMode != NONE) {
			if(readMode == KMER) {
				if(c == ' ' || c == '\n' || c == '\r') {
					if(nts_read == 0) { //has not start reading in kmer yet
						continue;
					}
					else { //start reading value (count of the kmer)
						if(nts_read > maxorder + 1) {
							fprintf(stderr, "Error: number of nts read is inconsistent with markov-order: %d\n", nts_read);
							exit(1);
						}
						readMode = KMER_VAL;
						value = 0;
						isReadingValue = false;
					}
				}
				else if(isChar(c)) { //it should be "is IUPAC character"
					int num = charToNum(c);
					if(!isNucleotide(num)) {
						fprintf(stderr, "Error: invalid character found when reading KMER: %c\n", c);
						exit(1);
					}
					nts[nts_read++] = num;
				}
				else {
					fprintf(stderr, "Error: invalid character found when reading KMER: %c\n", c);
					exit(1);
				}
			}
			else { //(readMode == KMER_VAL)
				if(c == ' ' || c == '\n' || c == '\r') {
					if(!isReadingValue) { //has not start reading value
						continue;
					}
					else { //finish reading value
						if(value <= 0) {
							fprintf(stderr, "Error: value must be >= 0: %d\n", value);
							exit(1);
						}
						setTransCount(rawModels[index], nts_read, nts, value);
						readMode = KMER;
						nts_read = 0;
						value = 0;
					}
				}
				else if( '0' <= c && c <= '9') {
					isReadingValue = true;
					value = 10 * value + (c - '0');
				}
				else {
					fprintf(stderr, "Error: invalid character found when reading KMER_VAL: %c\n", c);
				}
			}
		}
	}

	//Possibly no space characters at the end of value
	if(readMode != NONE && nts_read > 0 && value > 0) {
		setTransCount(rawModels[index], nts_read, nts, value);
		if(DEBUG1) {
			fprintf(stderr, "WARNING: MarkovFile does not end with space character.\n");
		}
	}

	//sanity to check all models have non-zero transition count
	for(int i = 0; i < numModels; i++) {
		for(int a = 0; a < NUMALPHAS && maxorder >= 0; a++) {
			assert(rawModels[index]->count0[a] > 0);
			for(int b = 0; b < NUMALPHAS && maxorder >= 1; b++) {
				assert(rawModels[index]->count1[a][b] > 0);
				for(int c = 0; c < NUMALPHAS && maxorder >= 2; c++) {
					assert(rawModels[index]->count2[a][b][c] > 0);
					for(int d = 0; d < NUMALPHAS && maxorder >= 3; d++) {
						assert(rawModels[index]->count3[a][b][c][d] > 0);
						for(int e = 0; e < NUMALPHAS && maxorder >= 4; e++) {
							assert(rawModels[index]->count4[a][b][c][d][e] > 0);
							for(int f = 0; f < NUMALPHAS && maxorder >= 5; f++) {
								assert(rawModels[index]->count5[a][b][c][d][e][f] > 0);
							}
						}
					}
				}
			}
		}
	}
}



Wndbgm* initWndbgm(int markovOrder, char* markovModelsFilename, int bgWndSize) {
	Wndbgm *wndbgm = (Wndbgm*) malloc(sizeof(Wndbgm));
	wndbgm->bgWndSize = bgWndSize;
	wndbgm->maxorder = markovOrder;

	//conditions
	if(markovModelsFilename == NULL && markovOrder != 0) {
		fprintf(stderr, "Error: higher-order wndspecbgm must provide a markovfile.\n");
		exit(1);
	}
	if(markovModelsFilename != NULL && markovOrder == 0) {
		fprintf(stderr, "Error:  0-th order wndspecbgm does not use a markovfile.\n");
		exit(1);
	}

	if(wndbgm->maxorder > 0) {
		FILE *fptr = openFile(markovModelsFilename);

		wndbgm->numModels = readNumOfModels(fptr);

		if(DEBUG1) {
			fprintf(stderr, "Number of models in MarkovFile: %d\n", wndbgm->numModels);
		}

		wndbgm->gcBinLocations = (int*) calloc(wndbgm->numModels, sizeof(int));
		TransCount **rawModels = (TransCount**) malloc(wndbgm->numModels * sizeof(TransCount*));
		for(int i = 0; i < wndbgm->numModels; i++) {
			rawModels[i] = (TransCount*) malloc(sizeof(TransCount));
		}

		readBinLocations(fptr, wndbgm->gcBinLocations, wndbgm->numModels);
		readMarkovModels(fptr, rawModels, wndbgm->numModels, wndbgm->maxorder);


		//accumulate transcount for lower-order based on max-order
		//for(int i = 0; i < wndbgm->numModels; i++) {
		//	marginalizeTransCount(rawModels[i]);
		//}

		wndbgm->models = (TransProb**) calloc(wndbgm->numModels, sizeof(TransProb*));
		for(int i = 0; i < wndbgm->numModels; i++) {
			wndbgm->models[i] = normalizeTransCount(rawModels[i]);
		}

		if(DEBUG1) {
			for(int i = 0; i < wndbgm->numModels; i++) {
				fprintf(stderr, ">%d<\n", wndbgm->gcBinLocations[i]);
				printTransCount(stderr, rawModels[i], true);
				printTransProb(stderr, wndbgm->models[i], true); //modify printTransProb
			}
		}
		//assert TransCount are greater than 0
		for(int i = 0; i < wndbgm->numModels; i++) {
			int minval = findMinValInTransCount(rawModels[i]);
			if(minval <= 0) {
				fprintf(stderr, "Error: invalid count of %d in rawModels for model %d\n", minval, i);
			}
		}

		//deallocate transcount
		for(int i = 0; i < wndbgm->numModels; i++) {
			free(rawModels[i]);
		}
		free(rawModels);
		
		fclose(fptr);
	}
	else {
		wndbgm->numModels = 0;
		wndbgm->gcBinLocations = NULL;
		wndbgm->models = NULL;
	}

	return wndbgm;
}


