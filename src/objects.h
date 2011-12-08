/*
 *  $Id: objects.h,v 1.20 2008/07/06 03:09:13 goswami Exp $
 *  
 *  File:    objects.h
 *  Package: EMC
 *
 *  Created by Gopi Goswami on Tue Apr 12 2006.
 *  Copyright (C) 2006 Gopika R. Goswami
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  For a copy of the GNU General Public License please write to the
 *  Free Software Foundation, Inc.
 *  59 Temple Place, Suite 330.
 *  Boston, MA  02111-1307 USA.
 *
 *  For bugs in the code please contact:
 *  goswami@stat.harvard.edu
 *
 *
 *  SYNOPSIS
 *
 *  This file is the declaration file for objects.c; it also contains
 *  some useful typedefs and enums.
 * 
 */

#ifndef OBJECTS_H
#define OBJECTS_H

#include <R.h>
#include <Rinternals.h>


/*
 * The forward declaration of all the typedefs
 */
typedef struct TimeDetails TimeDetails;
typedef struct SamplerMoveObject SamplerMoveObject;
typedef struct SampleLevelContext SampleLevelContext;
typedef struct ProposalCounter ProposalCounter;
typedef struct ARContext ARContext;
typedef struct ArgsList1 ArgsList1;
typedef struct ArgsList2 ArgsList2;
typedef struct Sampler Sampler;

typedef int    (*SamplerMove) (Sampler *);
typedef int    (*SamplerUtil) (Sampler *, SEXP);
typedef double (*FuncPtr1) (Sampler *, SEXP);
typedef SEXP   (*FuncPtr2) (Sampler *, SEXP);

typedef enum MoveType {
        MH   = 0,
        RC   = 1,
        SC   = 2,
        RE   = 3,
        BCE  = 4,
        BIRE = 5,
        BSE  = 6,
        CE   = 7
} MoveType;

typedef enum SelectionCode {
        RANDOM = 1,
        BEST,
        WORST
} SelectionCode;

#define N_IMPLEMENTED_MOVES 8
#define N_SELECTION 2
#define MAX_SAMPLER_MOVE_NAME 10

struct TimeDetails {
        double usr, sys;
};

struct SamplerMoveObject {
        char name[MAX_SAMPLER_MOVE_NAME];
        SamplerMove func;
};

struct SampleLevelContext {
        double *logWeights, *adjWeights, *partialSum;
        int samp, endLevel;
        double prob, numerator, sum, maxLogWeights;
};

struct ProposalCounter {
    double accepted;
    double proposed;
};

struct ARContext {
        int isSet, minLevels[2], maxLevels[2];
        double min, max;
        char minLevelsLabel[MAX_WORD_LENGTH], maxLevelsLabel[MAX_WORD_LENGTH];
        char minLabel[MAX_WORD_LENGTH], maxLabel[MAX_WORD_LENGTH];
};

struct ArgsList1 {
        int posDraw;
        SEXP argsList;
};

struct ArgsList2 {
        int posTemperature, posBlock, posDraw, posLogDens, posLevel, posIter;
        SEXP argsList;
};

struct Sampler {
        int nIters, thisIter;

        double timeInSecs;
        int nItersActual;

        int verboseLevel, printEstTimeAt, printEstTimeNTimes;
        int printInitialDotsWhen, printDotAt, eachDotWorth, nDotsPerLine;

        int nLevels, thisLevel;
        double *temperLadder, *invTemperLadder;
        double *logDensities, *logDensitiesStore; 
        SampleLevelContext *scratch_SLC;
        ARContext *scratch_ARC;
        
        int sampDim, MHNBlocks, *MHBlockNTimes, thisBlock, thisStep;

        SamplerMoveObject *moveObjs[N_IMPLEMENTED_MOVES];

        double moveProbs[N_IMPLEMENTED_MOVES], cumsumProbs[N_IMPLEMENTED_MOVES];
        Rboolean moveProbsPositive[N_IMPLEMENTED_MOVES];
        int moveNTimes[N_IMPLEMENTED_MOVES], moveNTimesAtIter[N_IMPLEMENTED_MOVES];
        SelectionCode moveSelectionCodes[N_IMPLEMENTED_MOVES][N_SELECTION];
        double moveSelectionTempers[N_IMPLEMENTED_MOVES][N_SELECTION];
        int SCRWMNTimes;
        double SCRWMPropSD;
        
        /*
         * MH: nLevels \times sampDim
         * all other moves: nLevels \times nLevels
         */
        ProposalCounter ***movePropCtrs[N_IMPLEMENTED_MOVES];

        int nLevelsSaveSampFor, *levelsSaveSampFor, nThin, nSave, saveIter;
        Rboolean saveFitness;
        SamplerUtil samplerOneIter, registerThisDraw;
        
        int nProtected;
        
        SEXP logTarDensFunc;
        ArgsList1 *logTarDensArgsList;
        FuncPtr1 logTarDens;

        SEXP oneIterMHFunc;
        ArgsList2 *oneIterMHArgsList;
        FuncPtr2 oneIterMH;
        
        SEXP doCallFuncCall, doCallFuncEnv;
        SEXP procTimeFuncCall, procTimeFuncEnv;
        TimeDetails *timeDetails;
        /*
         * nLevels \times sampDim
         */        
        SEXP SEXPCurrDraws, *SEXPCurrDrawsStore;
        SEXP SEXPPropDraws;
        
        SEXP argDraw, argCurrentDraw, argProposalDraw;
        SEXP argTemperature, argBlock, argLogDens, argLevel, argIter;
        SEXP dotsList;
        
        /*
         * nLevelsSaveSampFor \times sampDim \times nItersActual 
         * NOTE: this might need reallocation
         */
        double ***draws;
};

extern Sampler *
sampler_new (SEXP opts);

extern SEXP
sampler_run (Sampler *ss);

extern int
sampler_print (Sampler *ss);

extern int
sampler_move_RWM (Sampler *ss);

extern int
sampler_move_MH (Sampler *ss);

extern int
sampler_move_RC (Sampler *ss);

extern int
sampler_move_SC (Sampler *ss);

extern int
sampler_move_RE (Sampler *ss);

extern int
sampler_move_BCE (Sampler *ss);

extern int
sampler_move_BIRE (Sampler *ss);

extern int
sampler_move_BSE (Sampler *ss);

extern int
sampler_move_CE (Sampler *ss);

extern int
sampler_init (Sampler *ss);

extern SEXP
sampler_run (Sampler *ss);

extern SEXP
TOEMCMain (SEXP argsList);

#endif /* OBJECTS_H */
