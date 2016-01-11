#ifndef _WNDBGM_H
#define _WNDBGM_H

#include "transprob.h"
#include "symbols.h"

typedef struct {
	//models[i] is the Markov-model for GC-content between gcBinLocations[i] and gcBinLocations[i+1].
	//The last model is simply the model between gcBinLocations[n-1] and 100%.
	TransProb **models;
	int *gcBinLocations; 
	int numModels;

	int bgWndSize;
	int maxorder; //max order counted/calculated
} Wndbgm; 

extern TransProb* chooseModelForWndbgm(Wndbgm *wndbgm, int countA, int countC, int countG, int countT);
extern Wndbgm* initWndbgm(int markovOrder, char* markovModelsFilename, int bgWndSize);
extern void nilWndbgm(Wndbgm *wndbgm);

#endif
