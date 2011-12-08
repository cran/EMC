/*
 *  $Id: objects.c,v 1.28 2008/07/06 03:09:13 goswami Exp $
 *  
 *  File:    objects.c
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
 *  This file contains all the exported and the static functions
 *  related to a sampler object.
 *
 */

#include <assert.h>
#include <R.h>           
#include <Rinternals.h>
#include <Rmath.h>
#include "utils.h"
#include "objects.h"

#ifndef DEBUG_OBJECTS
#define DEBUG_OBJECTS 0
#endif

#if (!DEBUG_OBJECTS)
#define DEBUG(_x) ((void) 0)
#endif

/*
 * TODO: add checks for acceptance prob computations using R_FINITE,
 * now it's only done for MH and RWM
 */

/* The static functions of this file */
static SEXP
getListElement (SEXP list, char *str);

static int
gather_time_details (Sampler *ss);

static ProposalCounter ***
PC_mat_new (int dim1, int dim2);

static ARContext *
ARContext_init (ARContext *arc);

static SelectionCode
selection_codes (const char *codeName);

static int
args_list1_init (ArgsList1 *al);

static int
args_list2_init (ArgsList2 *al);

static double
log_target_density_func_user_Rfunc (Sampler *ss, SEXP draw);

static SEXP
one_iter_MH_user_Rfunc (Sampler *ss, SEXP draw);

static int
copy_draw (Sampler *ss, SEXP dest, SEXP src);

static int
sample_with_details (Sampler *ss);

static int
sample_with_details_excluding_one (Sampler *ss, int sl);

static int
exchange_given_prob (Sampler *ss, int *sl, ProposalCounter *pc, char *move,
                     double prob);

static int
exchange (Sampler *ss, int *sl, ProposalCounter *pc, char *move);

static int
sampler_move_RC_RANDOM_RANDOM (Sampler *ss);

static int
RC_BEST_BEST_prob_and_levels (Sampler *ss, double *prob, int *sl);

static int
RC_BEST_BEST_prob_given_levels (Sampler *ss, double *prob, int *sl,
                                double *logPropDens);

static int
sampler_move_RC_BEST_BEST (Sampler *ss);

static int
sampler_move_SC_given_selected_levels (Sampler *ss, int *sl);

static int
sampler_move_SC_with_two_levels (Sampler *ss);

static int
sampler_move_SC_with_many_levels (Sampler *ss);

static int
register_this_draw_fixed_iter_with_fitness (Sampler *ss, SEXP SEXPDraws);

static int
register_this_draw_fixed_iter (Sampler *ss, SEXP SEXPDraws);

static int
register_this_draw_fixed_time_with_fitness (Sampler *ss, SEXP SEXPDraws);

static int
register_this_draw_fixed_time (Sampler *ss, SEXP SEXPDraws);

static int
sampler_one_iter_MH (Sampler *ss);

static int
sampler_one_iter_with_one_level (Sampler *ss, SEXP notRequired);

static int
sampler_move_n_times_at_iter (Sampler *ss);

static int
sampler_one_iter_with_two_levels (Sampler *ss, SEXP notRequired);

static int
sampler_one_iter_with_many_levels (Sampler *ss, SEXP notRequired);

static int
sampler_make_draws_fixed_iter (Sampler *ss, SEXP SEXPDraws);

static int
sampler_make_draws_fixed_time (Sampler *ss, SEXP SEXPDraws);

static SEXP
make_draws (Sampler *ss);

static SEXP
make_accept_ratios (Sampler *ss);

static SEXP
make_accept_ratios_list_MH (Sampler *ss);

static SEXP
make_accept_ratios_list_move (Sampler *ss, int moveComp);

static SEXP
make_accept_ratios_list_CE (Sampler *ss);

static SEXP
make_accept_ratios_list (Sampler *ss);


/*
 * The following has been taken from the "Manual" section of R
 * website, it gets the list element named str, or returns NULL, if
 * not found. It has been modified a little bit.
 */
static SEXP
getListElement (SEXP list, char *str)
{
        SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
        int ii, found = 0, nn = length(list);
     
        for (ii = 0; ii < nn; ii++) {
                if (strcmp(CHAR(STRING_ELT(names, ii)), str) == 0) {
                        elmt = VECTOR_ELT(list, ii); ++found; break;
                }
        }
        if (found == 0) {
                char *errMsg1, errMsg2[MAX_LINE_LENGTH];

                errMsg1 = (char *) R_alloc(MAX_LINE_LENGTH, sizeof(errMsg1));
                for (ii = 0; ii < nn; ii++) {
                        errMsg1 = strcat(errMsg1, CHAR(STRING_ELT(names, ii)));
                        errMsg1 = strcat(errMsg1,
                                         ((ii == (nn - 1)) ? "" : ", "));
                }
                        
                sprintf(errMsg2,
                        "No element called \"%s\" found in the SEXP list, " \
                        "check your C code; given elements are:\n%s\n",
                        str, errMsg1);
                error(errMsg2);
        }
        return elmt;
}


static int
gather_time_details (Sampler *ss)
{
        SEXP SEXPTmp;
        double *doublesTmp;
        
        PROTECT(SEXPTmp      = eval(ss->procTimeFuncCall, ss->procTimeFuncEnv));
        doublesTmp           = REAL(SEXPTmp);
        ss->timeDetails->usr = doublesTmp[0];
        ss->timeDetails->sys = doublesTmp[1];
        UNPROTECT(1);
        return 0;
}


static ProposalCounter ***
PC_mat_new (int dim1, int dim2)
{
        ProposalCounter ***pppc;
        int ii, jj;
        
        pppc = (ProposalCounter ***) R_alloc(dim1, sizeof(ProposalCounter **));
        for (ii = 0; ii < dim1; ++ii) {
                pppc[ii] = (ProposalCounter **) R_alloc(dim2, sizeof(ProposalCounter *));
                for (jj = 0; jj < dim2; ++jj) {
                        pppc[ii][jj] = (ProposalCounter *) R_alloc(1, sizeof(struct ProposalCounter));
                        (pppc[ii][jj])->accepted = 0.0;
                        (pppc[ii][jj])->proposed = 0.0;
                }
        }
        return pppc;
}


static ARContext *
ARContext_init (ARContext *arc)
{
        arc->isSet = 0; arc->min = 1; arc->max = 0;
        arc->minLevels[0] = -1; arc->minLevels[1] = -1;
        arc->maxLevels[0] = -1; arc->maxLevels[1] = -1;

        return arc;
}


static SelectionCode
selection_codes (const char *codeName)
{
        SelectionCode code;
        
        if (strcmp(codeName, "random") == 0) code = RANDOM;
        else if (strcmp(codeName, "best") == 0) code = BEST;
        else code = WORST;
        return code;
}


static int
args_list1_init (ArgsList1 *al)
{
        int nComps = 0, comp = 0, nProtected = 0, nProtectedCaller = 0;
        SEXP names;

        al->posDraw = nComps++;

        PROTECT(al->argsList = allocVector(VECSXP, nComps)); ++nProtectedCaller;
        PROTECT(names        = allocVector(STRSXP, nComps)); ++nProtected;

        for (comp = 0; comp < nComps; ++comp)
                SET_VECTOR_ELT(al->argsList, comp, R_NilValue);
        comp = 0;
        SET_STRING_ELT(names, comp, mkChar("draw")); ++comp;

        setAttrib(al->argsList, R_NamesSymbol, names);
        UNPROTECT(nProtected);
        return nProtectedCaller;        
}
        

static int
args_list2_init (ArgsList2 *al)
{
        int nComps = 0, comp = 0, nProtected = 0, nProtectedCaller = 0;
        SEXP names;

        al->posTemperature = nComps++;
        al->posBlock       = nComps++;
        al->posDraw        = nComps++;
        al->posLogDens     = nComps++;
        al->posLevel       = nComps++;
        al->posIter        = nComps++;
        
        PROTECT(al->argsList = allocVector(VECSXP, nComps)); ++nProtectedCaller;
        PROTECT(names        = allocVector(STRSXP, nComps)); ++nProtected;

        for (comp = 0; comp < nComps; ++comp)
                SET_VECTOR_ELT(al->argsList, comp, R_NilValue);
        comp = 0;
        SET_STRING_ELT(names, comp, mkChar("temperature")); ++comp;
        SET_STRING_ELT(names, comp, mkChar("block")); ++comp;
        SET_STRING_ELT(names, comp, mkChar("draw")); ++comp;
        SET_STRING_ELT(names, comp, mkChar("logDens")); ++comp;
        SET_STRING_ELT(names, comp, mkChar("level")); ++comp;
        SET_STRING_ELT(names, comp, mkChar("iter")); ++comp;

        setAttrib(al->argsList, R_NamesSymbol, names);
        UNPROTECT(nProtected);
        return nProtectedCaller;        
}
        

Sampler *
sampler_new (SEXP opts)
{
        Sampler *ss;
        SEXP SEXPTmp, SEXPCurr, SEXPProp, SEXPCodes;
        int ii, jj, count, *intsTmp;
        double *doublesTmp, *startingVals;        
        char *name;
        
        ss             = (Sampler *) R_alloc(1, sizeof(struct Sampler));
        ss->nIters     = INTEGER(getListElement(opts, "nIters"))[0];
        ss->timeInSecs = REAL(getListElement(opts, "timeInSecs"))[0];

        ss->verboseLevel         = INTEGER(getListElement(opts, "verboseLevel"))[0];
        ss->printEstTimeAt       = 100;
        ss->printEstTimeNTimes   = 10;
        /* FIXME: The setting for ss->printDotAt, is it all right? */
        ss->printInitialDotsWhen = ss->printEstTimeAt / 10;
        ss->printDotAt           = 0;
        ss->nDotsPerLine         = 20;
        ss->eachDotWorth         = (int) ceil((ss->nIters - ss->printEstTimeAt + 1.0) / \
                                              (ss->printEstTimeNTimes * ss->nDotsPerLine));

        SEXPTmp             = getListElement(opts, "temperLadder");
        ss->nLevels         = length(SEXPTmp);
        doublesTmp          = REAL(SEXPTmp);
        ss->temperLadder    = (double *) R_alloc(ss->nLevels, sizeof(double));
        ss->invTemperLadder = (double *) R_alloc(ss->nLevels, sizeof(double));
        for (ii = 0; ii < ss->nLevels; ++ii) {
                ss->temperLadder[ii]    = doublesTmp[ii];
                ss->invTemperLadder[ii] = 1.0 / doublesTmp[ii];
        }

        ss->logDensities      = (double *) R_alloc(ss->nLevels, sizeof(double));
        ss->logDensitiesStore = (double *) R_alloc(ss->nLevels, sizeof(double));

        ss->scratch_SLC             = (SampleLevelContext *) R_alloc(1, sizeof(struct SampleLevelContext));
        ss->scratch_SLC->logWeights = (double *) R_alloc(ss->nLevels, sizeof(double));
        ss->scratch_SLC->adjWeights = (double *) R_alloc(ss->nLevels, sizeof(double));
        ss->scratch_SLC->partialSum = (double *) R_alloc(ss->nLevels, sizeof(double));

        ss->scratch_ARC = (ARContext *) R_alloc(1, sizeof(struct ARContext));
        
        SEXPTmp      = getListElement(opts, "startingVals");
        startingVals = REAL(SEXPTmp);        
        ss->sampDim  = ncols(SEXPTmp);

        SEXPTmp           = getListElement(opts, "MHBlocks");
        ss->MHNBlocks     = length(SEXPTmp);
        SEXPTmp           = getListElement(opts, "MHBlockNTimes");
        ss->MHBlockNTimes = (int *) R_alloc(ss->MHNBlocks, sizeof(int));
        intsTmp           = INTEGER(SEXPTmp);
        for (ii = 0; ii < ss->MHNBlocks; ++ii)
                ss->MHBlockNTimes[ii] = intsTmp[ii];

        for (ii = 0; ii < N_IMPLEMENTED_MOVES; ++ii)
                ss->moveObjs[ii] = (SamplerMoveObject *) \
                        R_alloc(1, sizeof(struct SamplerMoveObject));
        strcpy(ss->moveObjs[MH]->name,   "MH");
        strcpy(ss->moveObjs[RC]->name,   "RC");
        strcpy(ss->moveObjs[SC]->name,   "SC");
        strcpy(ss->moveObjs[RE]->name,   "RE");
        strcpy(ss->moveObjs[BCE]->name,  "BCE");
        strcpy(ss->moveObjs[BIRE]->name, "BIRE");
        strcpy(ss->moveObjs[BSE]->name,  "BSE");
        strcpy(ss->moveObjs[CE]->name,   "CE");
        
        SEXPTmp = getListElement(opts, "moveProbsList");
        for (ii = MH; ii <= CE; ++ii) {
                name = ss->moveObjs[ii]->name;
                ss->moveProbs[ii] = REAL(getListElement(SEXPTmp, name))[0];
                ss->moveProbsPositive[ii] = FALSE;
                if (ss->moveProbs[ii] > 0.0) {
                        ss->moveProbsPositive[ii] = TRUE;
                }
        }
        ss->cumsumProbs[0] = ss->moveProbs[MH];
        for (ii = RC; ii <= CE; ++ii)
                ss->cumsumProbs[ii] = ss->cumsumProbs[ii - 1] + ss->moveProbs[ii];

        SEXPTmp = getListElement(opts, "moveNTimesList");
        for (ii = MH; ii <= CE; ++ii) {
                name = ss->moveObjs[ii]->name;
                ss->moveNTimes[ii] = INTEGER(getListElement(SEXPTmp, name))[0];
        }

        SEXPTmp = getListElement(opts, "moveSelectionCodesList");
        SEXPCodes = getListElement(SEXPTmp, "RC");
        SEXPTmp = getListElement(opts, "moveSelectionTempersList");
        doublesTmp = REAL(getListElement(SEXPTmp, "RC"));
        for (ii = 0; ii < 2; ++ii) {
                ss->moveSelectionCodes[RC][ii] = \
                        selection_codes(CHAR(STRING_ELT(SEXPCodes, ii)));
                ss->moveSelectionTempers[RC][ii] = doublesTmp[ii];
        }
        
        SEXPTmp = getListElement(opts, "moveSelectionCodesList");
        SEXPCodes = getListElement(SEXPTmp, "SC");
        SEXPTmp = getListElement(opts, "moveSelectionTempersList");
        doublesTmp = REAL(getListElement(SEXPTmp, "SC"));
        ii = 0;        
        ss->moveSelectionCodes[SC][ii] = \
                selection_codes(CHAR(STRING_ELT(SEXPCodes, ii)));
        ss->moveSelectionTempers[SC][ii] = doublesTmp[ii];
        
        SEXPTmp = getListElement(opts, "moveSelectionTempersList");
        doublesTmp = REAL(getListElement(SEXPTmp, "BCE"));
        ii = 0;        
        ss->moveSelectionTempers[BCE][ii] = doublesTmp[ii];
        
        ss->movePropCtrs[MH] = PC_mat_new(ss->nLevels, ss->sampDim);
        for (ii = RC; ii <= CE; ++ii)
                ss->movePropCtrs[ii] = PC_mat_new(ss->nLevels, ss->nLevels);

        ss->SCRWMPropSD = REAL(getListElement(opts, "SCRWMPropSD"))[0];
        ss->SCRWMNTimes = INTEGER(getListElement(opts, "SCRWMNTimes"))[0];
        
        SEXPTmp = getListElement(opts, "levelsSaveSampFor");
        ss->nLevelsSaveSampFor = length(SEXPTmp);
        intsTmp = INTEGER(SEXPTmp);
        ss->levelsSaveSampFor = (int *) R_alloc(ss->nLevelsSaveSampFor, sizeof(int));
        for (ii = 0; ii < ss->nLevelsSaveSampFor; ++ii)
                ss->levelsSaveSampFor[ii] = intsTmp[ii] - 1;
        ss->nThin       = INTEGER(getListElement(opts, "nThin"))[0];
        ss->nSave       = INTEGER(getListElement(opts, "nSave"))[0];
        ss->saveFitness = LOGICAL(getListElement(opts, "saveFitness"))[0];
        
        ss->nProtected = 0;
        /* The user provided functions */
        ss->logTarDensFunc     = getListElement(opts, "logTarDensFunc");
        ss->logTarDensArgsList = (ArgsList1 *) R_alloc(1, sizeof(struct ArgsList1));
        ss->nProtected        += args_list1_init(ss->logTarDensArgsList);

        ss->oneIterMHFunc     = getListElement(opts, "oneIterMHFunc");
        ss->oneIterMHArgsList = (ArgsList2 *) R_alloc(1, sizeof(struct ArgsList2));
        ss->nProtected       += args_list2_init(ss->oneIterMHArgsList);

        PROTECT(ss->argTemperature = allocVector(REALSXP, 1)); ++(ss->nProtected);
        PROTECT(ss->argBlock       = allocVector(INTSXP, 1)); ++(ss->nProtected);
        PROTECT(ss->argLogDens     = allocVector(REALSXP, 1)); ++(ss->nProtected);
        PROTECT(ss->argLevel       = allocVector(INTSXP, 1)); ++(ss->nProtected);
        PROTECT(ss->argIter        = allocVector(INTSXP, 1)); ++(ss->nProtected);
        
        SEXPTmp                    = getListElement(opts, "doCallFunc");
        PROTECT(ss->doCallFuncCall = lang4(SEXPTmp, R_NilValue,
                                           R_NilValue, R_NilValue));
        ++(ss->nProtected);
        ss->doCallFuncEnv          = getListElement(opts, "doCallFuncEnv");
        
        SEXPTmp                      = getListElement(opts, "procTimeFunc");
        PROTECT(ss->procTimeFuncCall = lang1(SEXPTmp));
        ++(ss->nProtected);
        ss->procTimeFuncEnv          = getListElement(opts, "procTimeFuncEnv");

        ss->timeDetails = (TimeDetails *) R_alloc(1, sizeof(struct TimeDetails));

        PROTECT(ss->SEXPCurrDraws = allocVector(VECSXP, ss->nLevels));
        ss->SEXPCurrDrawsStore = (SEXP *) R_alloc(ss->nLevels, sizeof(SEXP));
        PROTECT(ss->SEXPPropDraws = allocVector(VECSXP, ss->nLevels));
        ss->nProtected += 2;
        for (ii = 0; ii < ss->nLevels; ++ii) {
                PROTECT(SEXPTmp = allocVector(REALSXP, ss->sampDim)); 
                SET_VECTOR_ELT(ss->SEXPCurrDraws, ii, SEXPTmp);
                PROTECT(SEXPTmp = allocVector(REALSXP, ss->sampDim)); 
                SET_VECTOR_ELT(ss->SEXPPropDraws, ii, SEXPTmp);
                UNPROTECT(2);
        }

        /* filling the starting values in */
        for (ii = 0; ii < ss->nLevels; ++ii) {
                SEXPCurr = VECTOR_ELT(ss->SEXPCurrDraws, ii);
                SEXPProp = VECTOR_ELT(ss->SEXPPropDraws, ii);
                for (jj = 0; jj < ss->sampDim; ++jj) {
                        count              = jj * ss->nLevels + ii;
                        REAL(SEXPCurr)[jj] = startingVals[count];
                        REAL(SEXPProp)[jj] = 0.0;
                }                
                PHONY(utils_SEXP_darray_print(SEXPCurr, ", "););
        }

        ss->dotsList = getListElement(opts, "dotsList");

        /* fixed iter case */
        if (ss->timeInSecs <= 0) { ss->draws = NULL; return ss; }

        if (ss->timeInSecs <= 0) ss->nItersActual = ss->nSave;
        else                     ss->nItersActual = -1;
        
        /* fixed time case */
        ss->draws = (double ***) R_alloc(ss->nLevelsSaveSampFor, sizeof(double **));
        for (ii = 0; ii < ss->nLevelsSaveSampFor; ++ii) {
                ss->draws[ii] = (double **) R_alloc(ss->sampDim, sizeof(double *));
                for (jj = 0; jj < ss->sampDim; ++jj) 
                        ss->draws[ii][jj] = (double *) R_alloc(ss->nItersActual, sizeof(double));
        }
        return ss;
}


int
sampler_print (Sampler *ss)
{
        PRINT_STUB_INT(ss->nIters);
        PRINT_STUB_INT(ss->nItersActual);
        PRINT_STUB_DOUBLE(ss->timeInSecs);
        PRINT_STUB_INT(ss->verboseLevel);
        PRINT_STUB_INT(ss->printDotAt);
        PRINT_STUB_INT(ss->nDotsPerLine);
        PRINT_STUB_INT(ss->eachDotWorth);
        PRINT_STUB_INT(ss->nLevels);
        Rprintf("The temperature ladder:\n");
        utils_darray_print(ss->temperLadder, ss->nLevels, ", ");
        Rprintf("The log densities:\n");
        utils_darray_print(ss->logDensities, ss->nLevels, ", ");
        PRINT_STUB_INT(ss->sampDim);
        PRINT_STUB_INT(ss->MHNBlocks);
        Rprintf("The MHBlockNTimes:\n");
        utils_iarray_print(ss->MHBlockNTimes, ss->MHNBlocks, ", ");
        PRINT_STUB_INT(ss->nLevelsSaveSampFor);        
        Rprintf("The levels to save samples for:\n");
        utils_iarray_print(ss->levelsSaveSampFor, ss->nLevelsSaveSampFor, ", ");
        PRINT_STUB_INT(ss->saveFitness);


        PRINT_STUB_DOUBLE(ss->moveProbs[MH]);
        PRINT_STUB_DOUBLE(ss->moveProbs[RC]);
        PRINT_STUB_DOUBLE(ss->moveProbs[SC]);
        utils_darray_print(ss->cumsumProbs, N_IMPLEMENTED_MOVES, ", ");
        
        PRINT_STUB_INT(ss->moveNTimes[MH]);
        PRINT_STUB_INT(ss->moveNTimes[RC]);
        PRINT_STUB_INT(ss->moveNTimes[SC]);
        PRINT_STUB_INT(ss->moveNTimes[RE]);
        PRINT_STUB_INT(ss->moveNTimes[BCE]);
        PRINT_STUB_INT(ss->moveNTimes[BIRE]);
        PRINT_STUB_INT(ss->moveNTimes[BSE]);
        PRINT_STUB_INT(ss->moveNTimes[CE]);

        PRINT_STUB_INT(ss->moveSelectionCodes[RC][0]);
        PRINT_STUB_INT(ss->moveSelectionCodes[RC][1]);
        PRINT_STUB_INT(ss->moveSelectionCodes[SC][0]);
        
        PRINT_STUB_DOUBLE(ss->moveSelectionTempers[RC][0]);
        PRINT_STUB_DOUBLE(ss->moveSelectionTempers[RC][1]);
        PRINT_STUB_DOUBLE(ss->moveSelectionTempers[SC][0]);
        PRINT_STUB_DOUBLE(ss->moveSelectionTempers[BCE][0]);
        return 0;
}


static double
log_target_density_func_user_Rfunc (Sampler *ss, SEXP draw)
{
        ArgsList1 *al = ss->logTarDensArgsList;
        SEXP SEXPTmp;
        double res;
        
        SET_VECTOR_ELT(al->argsList, al->posDraw, draw);
        
        SETCADR(ss->doCallFuncCall, ss->logTarDensFunc);
        SETCADDR(ss->doCallFuncCall, al->argsList);
        SETCADDDR(ss->doCallFuncCall, ss->dotsList);
        PROTECT(SEXPTmp = eval(ss->doCallFuncCall, ss->doCallFuncEnv));
        res = REAL(SEXPTmp)[0];
        UNPROTECT(1);
        return res;
}


static SEXP
one_iter_MH_user_Rfunc (Sampler *ss, SEXP draw)
{
        ArgsList2 *al = ss->oneIterMHArgsList;
        int ll = ss->thisLevel;
        
        REAL(ss->argTemperature)[0] = ss->temperLadder[ll];
        INTEGER(ss->argBlock)[0]    = ss->thisBlock + 1;
        REAL(ss->argLogDens)[0]     = ss->logDensities[ll];
        INTEGER(ss->argLevel)[0]    = ll + 1;
        INTEGER(ss->argIter)[0]     = ss->thisIter + 1;
        
        SET_VECTOR_ELT(al->argsList, al->posTemperature, ss->argTemperature);
        SET_VECTOR_ELT(al->argsList, al->posBlock, ss->argBlock);
        SET_VECTOR_ELT(al->argsList, al->posDraw, draw);
        SET_VECTOR_ELT(al->argsList, al->posLogDens, ss->argLogDens);
        SET_VECTOR_ELT(al->argsList, al->posLevel, ss->argLevel);
        SET_VECTOR_ELT(al->argsList, al->posIter, ss->argIter);

        SETCADR(ss->doCallFuncCall, ss->oneIterMHFunc);
        SETCADDR(ss->doCallFuncCall, al->argsList);
        SETCADDDR(ss->doCallFuncCall, ss->dotsList);
        return eval(ss->doCallFuncCall, ss->doCallFuncEnv);
}


static int
copy_draw (Sampler *ss, SEXP dest, SEXP src)
{
        int ii, nn = ss->sampDim;
        double *destTmp = REAL(dest), *srcTmp = REAL(src);

        for (ii = 0; ii < nn; ++ii) 
                destTmp[ii] = srcTmp[ii];
        return 0;
}


int
sampler_move_MH (Sampler *ss)
{
        int ll = ss->thisLevel, bb = ss->thisBlock;
        double logPropDens, alpha;
        ProposalCounter **ppc = ss->movePropCtrs[MH][ll];
        SEXP curr, prop, SEXPTmp;

        PHONY(PRINT_STUB_INT(ll); PRINT_STUB_INT(bb););      
        curr = VECTOR_ELT(ss->SEXPCurrDraws, ll);
        PROTECT(SEXPTmp = ss->oneIterMH(ss, curr));
        alpha = REAL(getListElement(SEXPTmp, "alpha"))[0];
        (ppc[bb])->proposed += 1;
        if (ss->verboseLevel >= 100) 
                Rprintf("MH: level: %d | block: %d | iter: %d | " \
                        "alpha: %5.4g\n", ll, bb, ss->thisIter, alpha);
        
        /* MH acceptance rejection step */
        if (runif(0, 1) <= alpha) {
                if (ss->verboseLevel >= 10) 
                        Rprintf("MH: level: %d | block: %d | iter: %d | " \
                                "alpha: %5.4g [*** accepted]\n",
                                ll, bb, ss->thisIter, alpha);
                prop        = getListElement(SEXPTmp, "prop");
                logPropDens = REAL(getListElement(SEXPTmp, "logPropDens"))[0];
                copy_draw(ss, curr, prop);
                ss->logDensities[ll] = logPropDens;
                (ppc[bb])->accepted += 1;
        }
        UNPROTECT(1);
        return 0;        
}


static int
sample_with_details (Sampler *ss)
{
        SampleLevelContext *slc = ss->scratch_SLC;
        double *lw = slc->logWeights, *aw = slc->adjWeights, *ps = slc->partialSum;
        double mlw = slc->maxLogWeights, uu, sum;
        int ll, nn = slc->endLevel;

        aw[0] = exp(lw[0] - mlw); ps[0] = aw[0];
        for (ll = 1; ll < nn; ++ll) {
                aw[ll] = exp(lw[ll] - mlw); ps[ll] = ps[ll - 1] + aw[ll];
        }
        sum = ps[nn - 1]; uu = runif(0, sum);
        for (ll = 0; ll < nn; ++ll) 
                if (uu <= ps[ll]) {
                        slc->samp = ll; slc->prob = aw[ll] / sum;
                        slc->numerator = aw[ll]; slc->sum = sum;
                        break;
                }
        return 0;
}


static int
sample_with_details_excluding_one (Sampler *ss, int sl)
{
        SampleLevelContext *slc = ss->scratch_SLC;
        double *lw = slc->logWeights, *aw = slc->adjWeights, *ps = slc->partialSum;
        double mlw = slc->maxLogWeights, uu, sum;
        int ll, nn = slc->endLevel;

        if (sl == 0) ps[0] = 0.0;
        else         aw[0] = exp(lw[0] - mlw); ps[0] = aw[0];
        for (ll = 1; ll < nn; ++ll) {
                if (sl == ll) { ps[ll] = ps[ll - 1]; continue; }
                aw[ll] = exp(lw[ll] - mlw); ps[ll] = ps[ll - 1] + aw[ll];
        }
        sum = ps[nn - 1]; uu = runif(0, sum);
        for (ll = 0; ll < nn; ++ll) 
                if ((sl != ll) && (uu <= ps[ll])) {
                        slc->samp = ll; slc->prob = aw[ll] / sum;
                        slc->numerator = aw[ll]; slc->sum = sum;
                        break;
                }
        return 0;
}


static int
exchange_given_prob (Sampler *ss, int *sl, ProposalCounter *pc, char *move,
                     double prob)
{        
        if (runif(0, 1) <= prob) {
                double *logDens = ss->logDensities;
                SEXP SEXPTmp;
                
                pc->accepted += 1; 
                if (ss->verboseLevel >= 10) 
                        Rprintf("%s: levels: %d, %d | iter: %d | " \
                                "alpha: %5.4g [*** accepted]\n",
                                move, sl[0], sl[1], ss->thisIter, prob);
                SEXPTmp = VECTOR_ELT(ss->SEXPCurrDraws, sl[0]);
                SET_VECTOR_ELT(ss->SEXPCurrDraws, sl[0], VECTOR_ELT(ss->SEXPCurrDraws, sl[1]));
                SET_VECTOR_ELT(ss->SEXPCurrDraws, sl[1], SEXPTmp);
                SWAP(double, logDens[sl[0]], logDens[sl[1]]);
        }
        return 0;
}


static int
exchange (Sampler *ss, int *sl, ProposalCounter *pc, char *move)
{
        double *ld = ss->logDensities, *itl = ss->invTemperLadder, alpha;
        
        pc->proposed += 1;
        alpha = exp((ld[sl[0]] - ld[sl[1]]) * (itl[sl[1]] - itl[sl[0]]));
        alpha = MIN(1.0, alpha);
        if (ss->verboseLevel >= 100) 
                Rprintf("%s: levels: %d, %d | iter: %d | " \
                        "alpha: %5.4g\n", move, sl[0], sl[1], ss->thisIter, alpha);
        exchange_given_prob(ss, sl, pc, move, alpha);
        return 0;
}


static int
sampler_move_RC_RANDOM_RANDOM (Sampler *ss)
{
        int ii, jj, sl[2], nLevels = ss->nLevels, pos;
        int sampDim = ss->sampDim, crosspoint;
        double *ld = ss->logDensities, *tl = ss->temperLadder;
        double logPropDens[2], sum, alpha;
        SEXP parents[2], children[2];
        double *parentVals[2], *childVals[2];
        ProposalCounter *pc;
        
        sl[0] = floor(runif(0, nLevels));
        pos = floor(runif(0, nLevels - 1));
        for (ii = 0; ii < nLevels; ++ii) 
                if (ii == pos) { sl[1] = ii + 1; break; }

        for (ii = 0; ii < 2; ++ii) {
                parents[ii] = VECTOR_ELT(ss->SEXPCurrDraws, sl[ii]);
                parentVals[ii] = REAL(parents[ii]);
                children[ii] = VECTOR_ELT(ss->SEXPPropDraws, sl[ii]);
                childVals[ii] = REAL(children[ii]);
        }
        /*
         * Here we choose a crosspoint in \{ 2:d \}, making it a
         * proper crosssover as opposed to exchange, as is the case
         * when coordinate 1 and d are chosen.
         */
        crosspoint = floor(runif(1, sampDim));
        for (jj = 0; jj < crosspoint; ++jj) {
                childVals[0][jj] = parentVals[0][jj];
                childVals[1][jj] = parentVals[1][jj];
        }
        for (jj = crosspoint; jj < sampDim; ++jj) {
                childVals[0][jj] = parentVals[1][jj];
                childVals[1][jj] = parentVals[0][jj];
        }
        for (sum = 0.0, ii = 0; ii < 2; ++ii) {
                logPropDens[ii] = (*(ss->logTarDens))(ss, children[ii]);
                sum += (logPropDens[ii] - ld[sl[ii]]) / tl[sl[ii]];
        }
        alpha = exp(sum); alpha = MIN(1.0, alpha);
        pc = ss->movePropCtrs[RC][sl[0]][sl[1]];
        pc->proposed += 1;
        if (ss->verboseLevel >= 100)
                Rprintf("RC: levels: %d, %d | iter: %d | " \
                        "alpha: %5.4g\n", sl[0], sl[1], ss->thisIter, alpha);
        /* the acceptance rejection step */
        if (runif(0, 1) <= alpha) {
                SEXP SEXPTmp;
                
                pc->accepted += 1;
                if (ss->verboseLevel >= 10)
                        Rprintf("RC: levels: %d, %d | iter: %d | " \
                                "alpha: %5.4g [*** accepted]\n",
                                sl[0], sl[1], ss->thisIter, alpha);
                for (ii = 0; ii < 2; ++ii) {
                        SEXPTmp = VECTOR_ELT(ss->SEXPCurrDraws, sl[ii]);
                        SET_VECTOR_ELT(ss->SEXPCurrDraws, sl[ii], VECTOR_ELT(ss->SEXPPropDraws, sl[ii]));
                        SET_VECTOR_ELT(ss->SEXPPropDraws, sl[ii], SEXPTmp);
                        ld[sl[ii]] = logPropDens[ii];
                }
        }
        return 0;
}


static int
RC_BEST_BEST_prob_and_levels (Sampler *ss, double *prob, int *sl)
{
        int ll, nLevels = ss->nLevels;
        SampleLevelContext *slc = ss->scratch_SLC;
        double *ld = ss->logDensities, selTemper = ss->moveSelectionTempers[RC][0];
        double *lw = slc->logWeights, *aw = slc->adjWeights, *ps = slc->partialSum;
        /* A: All, N: Numerator, S: Sum, B: But, O: One */
        double mlwA, parentNA[2], parentSA, mlwBO[2], parentNBO[2], parentSBO[2], uu;

        /* choose the first parent */
        for (mlwA = R_NegInf, ll = 0; ll < nLevels; ++ll) {
                lw[ll] = ld[ll] / selTemper; mlwA = MAX(mlwA, lw[ll]);
        }
        aw[0] = exp(lw[0] - mlwA); ps[0] = aw[0];
        for (ll = 1; ll < nLevels; ++ll) {
                aw[ll] = exp(lw[ll] - mlwA); ps[ll] = ps[ll - 1] + aw[ll];
        }
        parentSA = ps[nLevels - 1]; uu = runif(0, parentSA);
        for (ll = 0; ll < nLevels; ++ll)
                if (uu <= ps[ll]) { sl[0] = ll; break; }
        parentNA[0] = aw[sl[0]];

        /* choose the second parent */
        for (mlwBO[1] = R_NegInf, ll = 0; ll < nLevels; ++ll)
                if (ll != sl[0]) mlwBO[1] = MAX(mlwBO[1], lw[ll]);
        if (sl[0] == 0) ps[0] = 0.0;
        else            ps[0] = aw[0];
        for (ll = 1; ll < nLevels; ++ll) 
                if (ll == sl[0]) ps[ll] = ps[ll - 1];
                else             ps[ll] = ps[ll - 1] + exp(lw[ll] - mlwBO[1]);
        parentSBO[1] = ps[nLevels - 1]; uu = runif(0, parentSBO[1]);
        for (ll = 0; ll < nLevels; ++ll)
                if (uu <= ps[ll]) { sl[1] = ll; break; }
        parentNA[1] = aw[sl[1]];
        parentNBO[1] = exp(lw[sl[1]] - mlwBO[1]);
        
        /* do the rest of the probability computation */
        for (mlwBO[0] = R_NegInf, ll = 0; ll < nLevels; ++ll)
                if (ll != sl[1]) mlwBO[0] = MAX(mlwBO[0], lw[ll]);
        parentNBO[0] = exp(lw[sl[0]] - mlwBO[0]);
        for (parentSBO[0] = 0.0, ll = 0; ll < nLevels; ++ll) 
                if (ll != sl[1]) parentSBO[0] += exp(lw[ll] - mlwBO[0]);
        *prob = (parentNA[0] / parentSA) * (parentNBO[1] / parentSBO[1]) + \
                (parentNA[1] / parentSA) * (parentNBO[0] / parentSBO[0]);
        return 0;
}


static int
RC_BEST_BEST_prob_given_levels (Sampler *ss, double *prob, int *sl,
                                double *logPropDens)
{
        int ll, nLevels = ss->nLevels;
        SampleLevelContext *slc = ss->scratch_SLC;
        double *ld = ss->logDensities, selTemper = ss->moveSelectionTempers[RC][0];
        double *lw = slc->logWeights;
        /* A: All, N: Numerator, S: Sum, B: But, O: One */
        double mlwA, mlwBO[2] = { 0.0, 0.0 };
        double childNA[2] = { 0.0, 0.0 }, childSA;
        double childNBO[2] = { 0.0, 0.0 }, childSBO[2] = { 0.0, 0.0 };

        mlwA = R_NegInf; mlwBO[0] = R_NegInf; mlwBO[1] = R_NegInf;
        for (ll = 0; ll < nLevels; ++ll) {
                if (ll == sl[0]) {
                        lw[ll] = logPropDens[0] / selTemper;
                        mlwBO[0] = MAX(mlwBO[0], lw[ll]);
                } else if (ll == sl[1]) {
                        lw[ll] = logPropDens[1] / selTemper;
                        mlwBO[1] = MAX(mlwBO[1], lw[ll]);
                } else {
                        lw[ll] = ld[ll] / selTemper;
                        mlwBO[0] = MAX(mlwBO[0], lw[ll]);
                        mlwBO[1] = MAX(mlwBO[1], lw[ll]);
                }
                mlwA = MAX(mlwA, lw[ll]);                
        }
        childSA = 0.0; childSBO[0] = 0.0; childSBO[1] = 0.0;
        for (ll = 0; ll < nLevels; ++ll) {
                if (ll == sl[0]) {
                        childNA[0] = exp(lw[ll] - mlwA); childSA += childNA[0];
                        childNBO[0] = exp(lw[ll] - mlwBO[0]); childSBO[0] += childNBO[0];
                } else if (ll == sl[1]) {
                        childNA[1] = exp(lw[ll] - mlwA); childSA += childNA[1];
                        childNBO[1] = exp(lw[ll] - mlwBO[1]); childSBO[1] += childNBO[1];
                } else {
                        childSA += exp(lw[ll] - mlwA);
                        childNBO[0] += exp(lw[ll] - mlwBO[0]);
                        childNBO[1] += exp(lw[ll] - mlwBO[1]);
                }
        }
        *prob = (childNA[0] / childSA) * (childNBO[1] / childSBO[1]) + \
                (childNA[1] / childSA) * (childNBO[0] / childSBO[0]);
        return 0;
}


static int
sampler_move_RC_BEST_BEST (Sampler *ss)
{
        int ii, jj, sl[2], sampDim = ss->sampDim, crosspoint;
        double *ld = ss->logDensities, *tl = ss->temperLadder;
        SEXP parents[2], children[2];
        double *parentVals[2], *childVals[2];
        ProposalCounter *pc;
        double parentsToChildrenProb, childrenToParentsProb, logPropDens[2];
        double sum, alpha;

        RC_BEST_BEST_prob_and_levels(ss, &parentsToChildrenProb, sl);
        for (ii = 0; ii < 2; ++ii) {
                parents[ii] = VECTOR_ELT(ss->SEXPCurrDraws, sl[ii]);
                parentVals[ii] = REAL(parents[ii]);
                children[ii] = VECTOR_ELT(ss->SEXPPropDraws, sl[ii]);
                childVals[ii] = REAL(children[ii]);
        }
        /*
         * Here we choose a crosspoint in \{ 2:d \}, making it a
         * proper crosssover as opposed to exchange, as is the case
         * when coordinate 1 and d are chosen.
         */
        crosspoint = floor(runif(1, sampDim));
        for (jj = 0; jj < crosspoint; ++jj) {
                childVals[0][jj] = parentVals[0][jj];
                childVals[1][jj] = parentVals[1][jj];
        }
        for (jj = crosspoint; jj < sampDim; ++jj) {
                childVals[0][jj] = parentVals[1][jj];
                childVals[1][jj] = parentVals[0][jj];
        }
        for (sum = 0.0, ii = 0; ii < 2; ++ii) {
                logPropDens[ii] = (*(ss->logTarDens))(ss, children[ii]);
                sum += (logPropDens[ii] - ld[sl[ii]]) / tl[sl[ii]];
        }
        RC_BEST_BEST_prob_given_levels(ss, &childrenToParentsProb, sl, logPropDens);
        alpha = exp(sum) * (childrenToParentsProb / parentsToChildrenProb);
        alpha = MIN(1.0, alpha);
        pc = ss->movePropCtrs[RC][sl[0]][sl[1]];
        pc->proposed += 1;
        if (ss->verboseLevel >= 100)
                Rprintf("RC: levels: %d, %d | iter: %d | " \
                        "alpha: %5.4g\n", sl[0], sl[1], ss->thisIter, alpha);
        /* the acceptance rejection step */
        if (runif(0, 1) <= alpha) {
                SEXP SEXPTmp;
                
                pc->accepted += 1;
                if (ss->verboseLevel >= 10)
                        Rprintf("RC: levels: %d, %d | iter: %d | " \
                                "alpha: %5.4g [*** accepted]\n",
                                sl[0], sl[1], ss->thisIter, alpha);
                for (ii = 0; ii < 2; ++ii) {
                        SEXPTmp = VECTOR_ELT(ss->SEXPCurrDraws, sl[ii]);
                        SET_VECTOR_ELT(ss->SEXPCurrDraws, sl[ii], VECTOR_ELT(ss->SEXPPropDraws, sl[ii]));
                        SET_VECTOR_ELT(ss->SEXPPropDraws, sl[ii], SEXPTmp);
                        ld[sl[ii]] = logPropDens[ii];
                }
        }
        return 0;
}


static int
sampler_move_SC_given_selected_levels (Sampler *ss, int *sl)
{
        int ii, jj;
        int sampDim = ss->sampDim, sampDimM1 = sampDim - 1, accepted;
        double *ld = ss->logDensities, *tl = ss->temperLadder;
        double *currVals, *propVals, *anchorVals, *direcVals;
        double sum, direcNorm, rrCurr, rrProp = 0, rrPropSD = ss->SCRWMPropSD;
        double logCurrDens, logPropDens, rrLogCurrDens, rrLogPropDens, rrAlpha;
        ProposalCounter *pc;
        SEXP prop;

        currVals = REAL(VECTOR_ELT(ss->SEXPCurrDraws, sl[0]));
        prop = VECTOR_ELT(ss->SEXPPropDraws, sl[0]); propVals = REAL(prop);
        anchorVals = REAL(VECTOR_ELT(ss->SEXPCurrDraws, sl[1]));
        direcVals = REAL(VECTOR_ELT(ss->SEXPPropDraws, sl[1]));
        /* form the direction vector */
        for (sum = 0.0, ii = 0; ii < sampDim; ++ii) {
                direcVals[ii] = anchorVals[ii] - currVals[ii];
                sum += SQR(direcVals[ii]);
        }
        for (direcNorm = sqrt(sum), ii = 0; ii < sampDim; ++ii) 
                direcVals[ii] /= direcNorm;
        /* initialize the MH iteration for r */
        logCurrDens = ld[sl[0]];
        rrCurr = -direcNorm; accepted = 0;
        rrLogCurrDens = logCurrDens / tl[sl[0]] + sampDimM1 * log(fabs(rrCurr));
        /* MH iterations for r */
        for (jj = 0; jj < ss->SCRWMNTimes; ++jj) {
                rrProp = rnorm(rrCurr, rrPropSD);
                for (ii = 0; ii < sampDim; ++ii) 
                        propVals[ii] = anchorVals[ii] + rrProp * direcVals[ii];
                logPropDens = (*(ss->logTarDens))(ss, prop);
                rrLogPropDens = logPropDens / tl[sl[0]] + sampDimM1 * log(fabs(rrProp));
                rrAlpha = exp(rrLogPropDens - rrLogCurrDens);
                rrAlpha = MIN(1.0, rrAlpha);

                if (ss->verboseLevel >= 100)                
                        Rprintf("SC: levels: %d, %d | iter: %d | SCRWMiter: %d | " \
                                "rrCurr: %g | rrProp: %g | rrAlpha: %5.4g\n",
                                sl[0], sl[1], ss->thisIter, jj, rrCurr, rrProp, rrAlpha);
                
                if (runif(0, 1) <= rrAlpha) {
                        ++accepted; rrCurr = rrProp;
                        logCurrDens = logPropDens; rrLogCurrDens = rrLogPropDens;
                }
        }
        pc = ss->movePropCtrs[SC][sl[0]][sl[1]];
        pc->proposed += 1;
        if (accepted > 0) {
                SEXP SEXPTmp;

                DEBUG(Rprintf("SC: levels: %d, %d | iter: %d | rrCurr: %g | rrProp: %g" \
                              "[*** accepted]\n",
                              sl[0], sl[1], ss->thisIter, rrCurr, rrProp););
                
                if (ss->verboseLevel >= 10)
                        Rprintf("SC: levels: %d, %d | iter: %d | rrCurr: %g | rrProp: %g" \
                                "[*** accepted]\n",
                                sl[0], sl[1], ss->thisIter, rrCurr, rrProp);
                        
                pc->accepted += 1;
                for (ii = 0; ii < sampDim; ++ii) 
                        propVals[ii] = anchorVals[ii] + rrCurr * direcVals[ii];
                ld[sl[0]] = logCurrDens;
                SEXPTmp = VECTOR_ELT(ss->SEXPCurrDraws, sl[0]);
                SET_VECTOR_ELT(ss->SEXPCurrDraws, sl[0], VECTOR_ELT(ss->SEXPPropDraws, sl[0]));
                SET_VECTOR_ELT(ss->SEXPPropDraws, sl[0], SEXPTmp);
        }
        return 0;
}


static int
sampler_move_SC_with_two_levels (Sampler *ss)
{
        int sl[2];

        if (runif(0, 1) <= 0.5) { sl[0] = 0; sl[1] = 1; }
        else                    { sl[0] = 1; sl[1] = 0; }
        return sampler_move_SC_given_selected_levels(ss, sl);
}


static int
sampler_move_SC_with_many_levels (Sampler *ss)
{
        int ll, sl[2], nLevels = ss->nLevels;
        SampleLevelContext *slc = ss->scratch_SLC;        
        double *ld = ss->logDensities;
        double *lw = slc->logWeights, mlw;        
        double selTemper = ss->moveSelectionTempers[SC][0];
        
        sl[0] = floor(nLevels * runif(0, 1));
        for (mlw = R_NegInf, ll = 0; ll < nLevels; ++ll) 
                if (ll != sl[0]) {
                        lw[ll] = ld[ll] / selTemper; mlw = MAX(mlw, lw[ll]);
                }
        slc->maxLogWeights = mlw; slc->endLevel = nLevels;
        sample_with_details_excluding_one(ss, sl[0]); sl[1] = slc->samp;
        return sampler_move_SC_given_selected_levels(ss, sl);
}


int
sampler_move_RE (Sampler *ss)
{
        int sl[2], nLevels = ss->nLevels;
        ProposalCounter *pc;
        
        /* choose the levels to exchange */
        sl[0] = floor(nLevels * runif(0, 1));
        if (sl[0] == 0)                  sl[1] = sl[0] + 1;
        else if (sl[0] == (nLevels - 1)) sl[1] = sl[0] - 1;
        else {
                if (runif(0, 1) <= 0.5) sl[1] = sl[0] + 1;
                else                    sl[1] = sl[0] - 1;
        }
        pc = ss->movePropCtrs[RE][sl[0]][sl[1]];
        exchange(ss, sl, pc, "RE");        
        return 0;
}


int
sampler_move_BCE (Sampler *ss)
{
        int ll, sl[2];
        ProposalCounter *pc;
        SampleLevelContext *slc = ss->scratch_SLC;        
        double *ld = ss->logDensities;
        double *itl = ss->invTemperLadder, *lw = slc->logWeights;
        double selTemper = ss->moveSelectionTempers[BCE][0];
        double mlw, denom, CC, logR1, logR2, diff, alpha;

        sl[0] = ss->thisStep + 2;
        for (mlw = R_NegInf, ll = 0; ll < sl[0]; ++ll) {
                lw[ll] = ld[ll] / selTemper; mlw = MAX(mlw, lw[ll]);
        }
        slc->maxLogWeights = mlw; slc->endLevel = sl[0];
        sample_with_details(ss); sl[1] = slc->samp;
        denom = exp(ld[sl[0]] / selTemper - mlw);
        CC = slc->sum - slc->numerator;
        logR2 = log(slc->sum / (denom + CC));
        diff = (ld[sl[0]] - ld[sl[1]]);
        logR2 +=  diff / selTemper;
        logR1 = diff * (itl[sl[0]] - itl[sl[1]]);
        alpha = exp(logR1 + logR2); alpha = MIN(1.0, alpha);
        pc = ss->movePropCtrs[BCE][sl[0]][sl[1]];
        pc->proposed += 1;
        exchange_given_prob(ss, sl, pc, "BCE", alpha);
        return 0;
}


int
sampler_move_BIRE (Sampler *ss)
{
        int ll, sl[2];
        ProposalCounter *pc;
        SampleLevelContext *slc = ss->scratch_SLC;        
        double *ld = ss->logDensities, *tl = ss->temperLadder;
        double *lw = slc->logWeights, mlw, denom, CC, alpha;

        sl[0] = ss->thisStep + 2;
        for (mlw = R_NegInf, ll = 0; ll < sl[0]; ++ll) {
                lw[ll] = (ld[ll] / tl[sl[0]]) - (ld[ll] / tl[ll]);
                mlw = MAX(mlw, lw[ll]);
        }
        slc->maxLogWeights = mlw; slc->endLevel = sl[0];
        sample_with_details(ss); sl[1] = slc->samp;
        denom = exp((ld[sl[0]] / tl[sl[0]]) - (ld[sl[0]] / tl[sl[1]]) - mlw);
        CC = slc->sum - slc->numerator;
        alpha = slc->sum / (denom + CC); alpha = MIN(1.0, alpha);
        pc = ss->movePropCtrs[BIRE][sl[0]][sl[1]];
        pc->proposed += 1;
        exchange_given_prob(ss, sl, pc, "BIRE", alpha);
        return 0;
}


int
sampler_move_BSE (Sampler *ss)
{
        int ll, sl[2];
        ProposalCounter *pc;
        SampleLevelContext *slc = ss->scratch_SLC;        
        double *ld = ss->logDensities, *tl = ss->temperLadder;
        double *lw = slc->logWeights, mlw, denom, CC, alpha;

        sl[0] = ss->thisStep + 2;
        for (mlw = R_NegInf, ll = 0; ll < sl[0]; ++ll) {
                lw[ll] = (ld[ll] / tl[sl[0]]) + (ld[sl[0]] / tl[ll]);
                mlw = MAX(mlw, lw[ll]);
        }
        slc->maxLogWeights = mlw; slc->endLevel = sl[0];
        sample_with_details(ss); sl[1] = slc->samp;
        denom = exp((ld[sl[0]] / tl[sl[0]]) + (ld[sl[1]] / tl[sl[1]]) - mlw);
        CC = slc->sum - slc->numerator;
        alpha = slc->sum / (denom + CC); alpha = MIN(1.0, alpha);
        pc = ss->movePropCtrs[BSE][sl[0]][sl[1]];
        pc->proposed += 1;
        exchange_given_prob(ss, sl, pc, "BSE", alpha);
        return 0;
}


static int
cyclic_shift_draws (Sampler *ss, int ladderLength, int selLength)
{
        int jj, shiftBy = selLength + 1, rem;
        SEXP *cds = ss->SEXPCurrDrawsStore, cd = ss->SEXPCurrDraws;
        double *lds = ss->logDensitiesStore, *ld = ss->logDensities;

        for (jj = 0; jj < ladderLength; ++jj) {
                cds[jj] = VECTOR_ELT(cd, jj); lds[jj] = ld[jj];
        }
        for (jj = 0; jj < ladderLength; ++jj) {
                rem = (jj + shiftBy) % ladderLength;
                SET_VECTOR_ELT(cd, jj, cds[rem]); ld[jj] = lds[rem];
        }
        return 0;
}


int
sampler_move_CE (Sampler *ss)
{
        int ll, jj, ladderLength, selLength, shiftBy, rem;
        ProposalCounter *pc;
        SampleLevelContext *slc = ss->scratch_SLC;        
        double *ld = ss->logDensities, *tl = ss->temperLadder;
        double *lw = slc->logWeights, mlw, sum, logDenom, CC, alpha;

        ladderLength = ss->thisStep + 3; 
        for (mlw = R_NegInf, ll = 0; ll < (ladderLength - 1); ++ll) {
                shiftBy = ll + 1; lw[ll] = 0.0;
                for (jj = 0; jj < ladderLength; ++jj) {
                        rem = (jj + shiftBy) % ladderLength;
                        lw[ll] += ld[rem] / tl[jj];
                }
                mlw = MAX(mlw, lw[ll]);
        }
        slc->maxLogWeights = mlw; slc->endLevel = ladderLength - 1;
        sample_with_details(ss); selLength = slc->samp;
        sum = slc->sum; CC = sum - slc->numerator;
        for (logDenom = 0.0, jj = 0; jj < ladderLength; ++jj)
                logDenom += ld[jj] / tl[jj];
        alpha = sum / (exp(logDenom - mlw) + CC); alpha = MIN(alpha, 1.0);
        pc = ss->movePropCtrs[CE][ladderLength - 1][selLength];
        pc->proposed += 1;
        if (runif(0, 1) <= alpha) {
                pc->accepted += 1;
                cyclic_shift_draws(ss, ladderLength, selLength);
        }
        return 0;
}


static int
register_this_draw_fixed_iter_with_fitness (Sampler *ss, SEXP SEXPDraws)
{
        int ii, jj, kk, ll, mm;
        double *dest, *src;
        
        dest = REAL(SEXPDraws);
        for (kk = 0; kk < ss->nLevelsSaveSampFor; ++kk) {
                ll  = ss->levelsSaveSampFor[kk];
                src = REAL(VECTOR_ELT(ss->SEXPCurrDraws, ll));
                mm  = (kk * ss->nSave * (ss->sampDim + 1)); 
                for (jj = 0; jj < ss->sampDim; ++jj) {
                        ii       = mm + ss->saveIter + jj * ss->nSave;
                        dest[ii] = src[jj];
                }
                ii       = mm + ss->saveIter + jj * ss->nSave;
                dest[ii] = -ss->logDensities[ll];
                PHONY(utils_darray_print(src, ss->sampDim, ", "););
        }
        return 0;
}


static int
register_this_draw_fixed_iter (Sampler *ss, SEXP SEXPDraws)
{
        int ii, jj, kk, ll, mm;
        double *dest, *src;
        
        dest = REAL(SEXPDraws);
        for (kk = 0; kk < ss->nLevelsSaveSampFor; ++kk) {
                ll = ss->levelsSaveSampFor[kk];
                src = REAL(VECTOR_ELT(ss->SEXPCurrDraws, ll));
                mm = (kk * ss->nSave * ss->sampDim); 
                for (jj = 0; jj < ss->sampDim; ++jj) {
                        ii = mm + ss->saveIter + jj * ss->nSave;
                        dest[ii] = src[jj];
                }
                PHONY(utils_darray_print(src, ss->sampDim, ", "););
        }
        return 0;
}


static int
register_this_draw_fixed_time_with_fitness (Sampler *ss, SEXP SEXPDraws)
{
        /* WRITME */
        return 0;
}


static int
register_this_draw_fixed_time (Sampler *ss, SEXP SEXPDraws)
{        
        /* WRITME */
        return 0;
}


static int
sampler_one_iter_MH (Sampler *ss)
{
        int bb, tt, nTimes;

        for (bb = 0; bb < ss->MHNBlocks; ++bb) {
                ss->thisBlock = bb;
                nTimes = ss->MHBlockNTimes[bb];
                for (tt = 0; tt < nTimes; ++tt)
                        (*(ss->moveObjs[MH]->func))(ss);
        }
        return 0;
}


static int
sampler_one_iter_with_one_level (Sampler *ss, SEXP notRequired)
{
        int ii;
        
        ss->thisLevel = 0;
        for (ii = 0; ii < ss->moveNTimes[MH]; ++ii)
                sampler_one_iter_MH(ss);
        return 0;
}


static int
sampler_move_n_times_at_iter (Sampler *ss)
{
        int jj;
        double uu = runif(0, 1);
        Rboolean found = FALSE;
        
        for (jj = MH; jj <= CE; ++jj) {
                if (ss->moveProbsPositive[jj] == TRUE) {
                        if ((found == FALSE) && (uu <= ss->cumsumProbs[jj])) {
                                ss->moveNTimesAtIter[jj] = ss->moveNTimes[jj];
                                found                    = TRUE;
                        }
                        else { ss->moveNTimesAtIter[jj] = 0; }
                }
                else { ss->moveNTimesAtIter[jj] = ss->moveNTimes[jj]; }
        }
        return 0;
}


static int
sampler_one_iter_with_two_levels (Sampler *ss, SEXP notRequired)
{
        int ii, jj, sl[2], ll;
        ProposalCounter *pc;
        char *moveName;

        sampler_move_n_times_at_iter(ss);
        for (ii = 0; ii < ss->moveNTimesAtIter[MH]; ++ii)  {
                for (ll = 0; ll < 2; ++ll) {
                        ss->thisLevel = ll;
                        sampler_one_iter_MH(ss);
                }
        }
        for (jj = RC; jj <= SC; ++jj) 
                for (ii = 0; ii < ss->moveNTimesAtIter[jj]; ++ii)  
                        (*(ss->moveObjs[jj]->func))(ss);
        for (jj = RE; jj <= CE; ++jj) {
                for (ii = 0; ii < ss->moveNTimesAtIter[jj]; ++ii) {
                        ss->thisStep = ii; sl[0] = 0; sl[1] = 1;
                        pc = ss->movePropCtrs[jj][sl[0]][sl[1]];
                        moveName = ss->moveObjs[jj]->name;
                        exchange(ss, sl, pc, moveName);
                }
        }
        return 0;
}


static int
sampler_one_iter_with_many_levels (Sampler *ss, SEXP notRequired)
{
        int ii, jj, ll;
        
        sampler_move_n_times_at_iter(ss);
        for (ii = 0; ii < ss->moveNTimesAtIter[MH]; ++ii)  {
                for (ll = 0; ll < ss->nLevels; ++ll) {
                        ss->thisLevel = ll;
                        sampler_one_iter_MH(ss);
                }
        }
        for (jj = RC; jj <= CE; ++jj) {
                for (ii = 0; ii < ss->moveNTimesAtIter[jj]; ++ii) {
                        ss->thisStep = ii; (*(ss->moveObjs[jj]->func))(ss);
                }
        }
        return 0;
}


static int
sampler_make_draws_fixed_iter (Sampler *ss, SEXP SEXPDraws)
{
        int ii, iter = 0, nIters = ss->nIters, nThin = ss->nThin;
        Rboolean initialEstDone = FALSE;
        double timeStartUsr, timeEndUsr, timeToFinish;
        double timeStartSys, timeEndSys;    

        gather_time_details(ss);
        timeStartUsr = ss->timeDetails->usr;
        timeStartSys = ss->timeDetails->sys;
        GetRNGstate( );
        ss->saveIter = 0;
        for (ii = 0; ii < nIters; ++ii) {
                ss->thisIter = ii; (*(ss->samplerOneIter))(ss, R_NilValue);
                if ((ii % nThin) == 0) {
                        (*(ss->registerThisDraw))(ss, SEXPDraws);
                        ++(ss->saveIter);
                }
                if (ss->verboseLevel >= 1) {
                        iter = ii + 1;
                        if (initialEstDone == FALSE) {
                                if ((iter % ss->printInitialDotsWhen) == 0)
                                        Rprintf(".");
                        }
                        else {
                                if (iter == ss->printDotAt) {
                                        Rprintf(".");
                                        ss->printDotAt += ss->eachDotWorth;
                                }
                        }
                        if (iter == ss->printEstTimeAt) {
                                gather_time_details(ss);
                                timeEndUsr = ss->timeDetails->usr;
                                timeToFinish = ((timeEndUsr - timeStartUsr) / iter * \
                                                (nIters - iter));
                                Rprintf("\n[Time to finish (est): %8.0f secs, " \
                                        "this iter: %10d]", ceil(timeToFinish), iter);
                                if (initialEstDone == FALSE) initialEstDone = TRUE;
                                ss->printEstTimeAt += ss->eachDotWorth * ss->nDotsPerLine;
                                ss->printDotAt = iter + 1;
                        }                        
                }
        }
        PutRNGstate( );
        if (ss->verboseLevel >= 1) {
                gather_time_details(ss);
                timeEndUsr = ss->timeDetails->usr;
                timeEndSys = ss->timeDetails->sys;
                Rprintf("\n[Total time: %10.0f secs (usr), %5.0f secs (sys), " \
                        "this iter: %10d]\n", (timeEndUsr - timeStartUsr),
                        (timeEndSys - timeStartSys), iter);
        }
        PHONY(utils_SEXP_darray_print(SEXPDraws, ", "););
        return 0;
}


static int
sampler_make_draws_fixed_time (Sampler *ss, SEXP SEXPDraws)
{
        int ii, iter = 0, nIters = ss->nIters, nThin = ss->nThin;
        Rboolean initialEstDone = FALSE;
        double timeStartUsr, timeEndUsr, timeToFinish;
        double timeStartSys, timeEndSys;    

        error("fixed time case: not yet implemented");
        gather_time_details(ss);
        timeStartUsr = ss->timeDetails->usr;
        timeStartSys = ss->timeDetails->sys;
        GetRNGstate( );
        ss->saveIter = 0;
        /* infinite loop, break out of it if ss->timeInSecs is exceeded  */
        for (ii = 0; ; ++ii) {
                ss->thisIter = ii; (*(ss->samplerOneIter))(ss, R_NilValue);
                if ((ii % nThin) == 0) {
                        (*(ss->registerThisDraw))(ss, SEXPDraws);
                        ++(ss->saveIter);
                }
                if (ss->verboseLevel >= 1) {
                        iter = ii + 1;
                        if (initialEstDone == FALSE) {
                                if ((iter % ss->printInitialDotsWhen) == 0) Rprintf(".");
                        }
                        else {
                                if (iter == ss->printDotAt) {
                                        Rprintf(".");
                                        ss->printDotAt += ss->eachDotWorth;
                                }
                        }
                        if (iter == ss->printEstTimeAt) {
                                gather_time_details(ss);
                                timeEndUsr = ss->timeDetails->usr;
                                timeToFinish = ((timeEndUsr - timeStartUsr) / iter * \
                                                (nIters - iter));
                                Rprintf("\n[Time to finish (est): %8.0f secs, " \
                                        "this iter: %10d]", ceil(timeToFinish), iter);
                                if (initialEstDone == FALSE) initialEstDone = TRUE;
                                ss->printEstTimeAt += ss->eachDotWorth * ss->nDotsPerLine;
                                ss->printDotAt = iter + 1;
                        }                        
                }
                gather_time_details(ss);
                timeEndSys = ss->timeDetails->sys;
                if ((timeEndSys - timeStartSys) > ss->timeInSecs) {
                        ss->nItersActual = ii; break;
                }
        } 
        PutRNGstate( );
        if (ss->verboseLevel >= 1) {
                gather_time_details(ss);
                timeEndUsr = ss->timeDetails->usr;
                timeEndSys = ss->timeDetails->sys;
                Rprintf("\n[Total time: %10.0f secs (usr), %5.0f secs (sys), " \
                        "this iter: %10d]\n", (timeEndUsr - timeStartUsr),
                        (timeEndSys - timeStartSys), iter);
        }
        return 0;
}


int
sampler_init (Sampler *ss)
{
        int ii, nLevels = ss->nLevels;
        SEXP curr;
        
        /* installing the functions */
        if (ss->logTarDensFunc != R_NilValue)
                ss->logTarDens = log_target_density_func_user_Rfunc;
        if (ss->oneIterMHFunc != R_NilValue) {
                ss->oneIterMH          = one_iter_MH_user_Rfunc;
                ss->moveObjs[MH]->func = sampler_move_MH;
        }
        
        if (nLevels == 1) {
                ss->samplerOneIter = sampler_one_iter_with_one_level;
        } else if (nLevels == 2) {
                ss->samplerOneIter = sampler_one_iter_with_two_levels;
                /*
                 * In this case the BEST-BEST case yields to proposal
                 * probability cacellations and hence we just use the
                 * RANDOM-RANDOM case a proxy for the BEST-BEST case.
                 */
                ss->moveObjs[RC]->func = sampler_move_RC_RANDOM_RANDOM;
                ss->moveObjs[SC]->func = sampler_move_SC_with_two_levels;
        } else {
                ss->samplerOneIter = sampler_one_iter_with_many_levels;
                if (ss->moveSelectionCodes[RC][0] == RANDOM)
                        ss->moveObjs[RC]->func = sampler_move_RC_RANDOM_RANDOM;
                else if (ss->moveSelectionCodes[RC][0] == BEST)
                        ss->moveObjs[RC]->func = sampler_move_RC_BEST_BEST;
                ss->moveObjs[SC]->func = sampler_move_SC_with_many_levels;
        }

        ss->moveObjs[RE]->func   = sampler_move_RE;
        ss->moveObjs[BCE]->func  = sampler_move_BCE;        
        ss->moveObjs[BIRE]->func = sampler_move_BIRE;
        ss->moveObjs[BSE]->func  = sampler_move_BSE;
        ss->moveObjs[CE]->func   = sampler_move_CE;                
        
        for (ii = 0; ii < nLevels; ++ii) {
                curr = VECTOR_ELT(ss->SEXPCurrDraws, ii);                
                ss->logDensities[ii] = (*(ss->logTarDens))(ss, curr);
                if (R_FINITE(ss->logDensities[ii]) == FALSE) {
                        char errMsg[MAX_LINE_LENGTH];

                        sprintf(errMsg,
                                "logTarDens evaluation for level [%d] " \
                                "gives [%f]", ii, ss->logDensities[ii]);
                        error(errMsg);
                }
        }
        if (ss->timeInSecs <= 0) {
                if (ss->saveFitness == TRUE)
                        ss->registerThisDraw = register_this_draw_fixed_iter_with_fitness;
                else
                        ss->registerThisDraw = register_this_draw_fixed_iter;
        }
        else {
                if (ss->saveFitness == TRUE)
                        ss->registerThisDraw = register_this_draw_fixed_time_with_fitness;
                else
                        ss->registerThisDraw = register_this_draw_fixed_time;
        }
        return 0;
}


static SEXP
make_draws (Sampler *ss)
{
        /*
         * draws is an array of dimension:
         * c(nSave, sampDim (+ 1), nLevelsSaveSampFor)
         */
        SEXP draws, drawsDim;
        int ii, nn, *intsTmp, dims[3];

        dims[0] = ss->nSave; dims[2] = ss->nLevelsSaveSampFor;
        if (ss->saveFitness == TRUE) dims[1] = ss->sampDim + 1;
        else                         dims[1] = ss->sampDim;
        nn = dims[0] * dims[1] * dims[2];
        PROTECT(draws = allocVector(REALSXP, nn)); ++(ss->nProtected);        
        PROTECT(drawsDim = allocVector(INTSXP, 3)); 
        intsTmp = INTEGER(drawsDim);
        for (ii = 0; ii < 3; ++ii) intsTmp[ii] = dims[ii];
        setAttrib(draws, R_DimSymbol, drawsDim);
        UNPROTECT(1);
        if (ss->timeInSecs <= 0) sampler_make_draws_fixed_iter(ss, draws);
        else                     sampler_make_draws_fixed_time(ss, draws);
        return draws;
}


static SEXP
make_accept_ratios (Sampler *ss)
{
        int ii, bb, ll1, ll2, nLevels = ss->nLevels, nVarious;
        int row, nProtected = 0;
        double accepted, proposed, *doublesTmp;
        ProposalCounter ***pppc;
        SEXP acceptRatios, aRRowNames, aRColNames, aRDimNames;

        for (nVarious = 0, ii = MH; ii <= CE; ++ii) 
                if (ss->moveNTimes[ii] > 0) ++nVarious;        
        PROTECT(acceptRatios = allocMatrix(REALSXP, nVarious, 4));
        ++(ss->nProtected); 
        doublesTmp = REAL(acceptRatios); 
        /* The acceptRatios for MH */
        accepted = 0; proposed = 0; pppc = ss->movePropCtrs[MH];
        for (ll1 = 0; ll1 < nLevels; ++ll1) {
                for (bb = 0; bb < ss->MHNBlocks; ++bb) {
                        accepted += (pppc[ll1][bb])->accepted;
                        proposed += (pppc[ll1][bb])->proposed;
                }
        }
        if (fabs(proposed - 0) <= 0) doublesTmp[MH] = R_NaReal;
        else                         doublesTmp[MH] = accepted / proposed;
        doublesTmp[MH + nVarious] = accepted;
        doublesTmp[MH + 2 * nVarious] = proposed;        
        doublesTmp[MH + 3 * nVarious] = ss->moveNTimes[MH];        
        /* The acceptRatios for RC through CE */
        row = 1;        
        for (ii = RC; ii <= CE; ++ii) {
                /* no row to be filled in if the move wasn't used */
                if (ss->moveNTimes[ii] == 0) continue;                
                accepted = 0; proposed = 0; pppc = ss->movePropCtrs[ii];
                for (ll1 = 0; ll1 < nLevels; ++ll1) {
                        for (ll2 = 0; ll2 < nLevels; ++ll2) {
                                accepted += (pppc[ll1][ll2])->accepted;
                                proposed += (pppc[ll1][ll2])->proposed;
                        }
                }
                if (fabs(proposed - 0) <= 0) doublesTmp[row] = R_NaReal;
                else                         doublesTmp[row] = accepted / proposed;
                doublesTmp[row + nVarious] = accepted;
                doublesTmp[row + 2 * nVarious] = proposed;
                doublesTmp[row + 3 * nVarious] = ss->moveNTimes[ii];        
                ++row;
        }

        PROTECT(aRRowNames = allocVector(STRSXP, nVarious)); ++nProtected; 
        row = 0;
        for (ii = MH; ii <= CE; ++ii) {
                if (ss->moveNTimes[ii] > 0) {                
                        SET_STRING_ELT(aRRowNames, row, mkChar(ss->moveObjs[ii]->name));
                        ++row;
                }
        }
        PROTECT(aRColNames = allocVector(STRSXP, 4)); ++nProtected;  
        SET_STRING_ELT(aRColNames, 0, mkChar("ratio"));
        SET_STRING_ELT(aRColNames, 1, mkChar("accepted"));
        SET_STRING_ELT(aRColNames, 2, mkChar("proposed"));
        SET_STRING_ELT(aRColNames, 3, mkChar("moveNTimes"));
        PROTECT(aRDimNames = allocVector(VECSXP, 2)); ++nProtected; 
        SET_VECTOR_ELT(aRDimNames, 0, aRRowNames);
        SET_VECTOR_ELT(aRDimNames, 1, aRColNames);
        setAttrib(acceptRatios, R_DimNamesSymbol, aRDimNames);
        UNPROTECT(nProtected);
        return acceptRatios;
}


static SEXP
make_accept_ratios_list_MH (Sampler *ss)
{
        int bb, ll, jj, nLevels = ss->nLevels, nBlocks = ss->MHNBlocks;
        int nProtected = 0;
        double *doublesTmp, aR;
        ProposalCounter **ppc;
        SEXP aRMat, aRRowNames, aRColNames, aRDimNames;
        char charsTmp[MAX_LINE_LENGTH];

        PROTECT(aRMat = allocMatrix(REALSXP, nLevels, nBlocks)); ++nProtected;
        doublesTmp = REAL(aRMat); 
        for (ll = 0; ll < nLevels; ++ll) {
                ppc = ss->movePropCtrs[MH][ll];
                for (bb = 0; bb < nBlocks; ++bb) {
                        aR = R_NaReal;
                        if (ppc[bb]->proposed > 0)
                                aR = ppc[bb]->accepted / ppc[bb]->proposed;
                        
                        jj             = (nLevels * bb + ll);
                        doublesTmp[jj] = aR;
                }
        }
        
        PROTECT(aRRowNames = allocVector(STRSXP, nLevels)); ++nProtected;
        for (ll = 0; ll < nLevels; ++ll) {
                sprintf(charsTmp, "level%d", ll + 1);
                SET_STRING_ELT(aRRowNames, ll, mkChar(charsTmp));
        }
        
        PROTECT(aRColNames = allocVector(STRSXP, nBlocks)); ++nProtected;
        for (bb = 0; bb < nBlocks; ++bb) {
                sprintf(charsTmp, "block%d", bb + 1);
                SET_STRING_ELT(aRColNames, bb, mkChar(charsTmp));
        }

        PROTECT(aRDimNames = allocVector(VECSXP, 2)); ++nProtected;
        SET_VECTOR_ELT(aRDimNames, 0, aRRowNames);
        SET_VECTOR_ELT(aRDimNames, 1, aRColNames);
        setAttrib(aRMat, R_DimNamesSymbol, aRDimNames);

        UNPROTECT(nProtected);
        return aRMat;
}


static SEXP
make_accept_ratios_list_move (Sampler *ss, int moveComp)
{
        int ll0, ll1, jj, nLevels = ss->nLevels;
        int nProtected = 0;
        double *doublesTmp, proposed, aR;
        ProposalCounter ***pppc = ss->movePropCtrs[moveComp];
        ProposalCounter *pc01, *pc10;
        SEXP aRMat, aRRowNames, aRColNames, aRDimNames;
        char charsTmp[MAX_LINE_LENGTH];

        PROTECT(aRMat = allocMatrix(REALSXP, nLevels, nLevels)); ++nProtected;
        doublesTmp = REAL(aRMat); 
        for (ll0 = 0; ll0 < nLevels; ++ll0) {
                for (ll1 = 0; ll1 < nLevels; ++ll1) {
                        aR       = R_NaReal;
                        pc01     = pppc[ll0][ll1];
                        pc10     = pppc[ll1][ll0]; 
                        proposed = pc01->proposed + pc10->proposed;
                        if (proposed > 0)
                                aR = (pc01->accepted + pc10->accepted) / proposed;

                        jj             = (nLevels * ll1 + ll0);
                        doublesTmp[jj] = aR;
                }
        }

        PROTECT(aRRowNames = allocVector(STRSXP, nLevels)); ++nProtected;
        PROTECT(aRColNames = allocVector(STRSXP, nLevels)); ++nProtected;
        for (ll0 = 0; ll0 < nLevels; ++ll0) {
                sprintf(charsTmp, "level%d", ll0 + 1);
                SET_STRING_ELT(aRRowNames, ll0, mkChar(charsTmp));
                SET_STRING_ELT(aRColNames, ll0, mkChar(charsTmp));
        }
        
        PROTECT(aRDimNames = allocVector(VECSXP, 2)); ++nProtected;
        SET_VECTOR_ELT(aRDimNames, 0, aRRowNames);
        SET_VECTOR_ELT(aRDimNames, 1, aRColNames);
        setAttrib(aRMat, R_DimNamesSymbol, aRDimNames);

        UNPROTECT(nProtected);
        return aRMat;        
}


static SEXP
make_accept_ratios_list_CE (Sampler *ss)
{
        int ii, ll, sl, jj, nLevels = ss->nLevels, nLadderLengths, nProtected = 0;
        double *doublesTmp, aR;
        ProposalCounter **ppc;
        SEXP aRMat, aRRowNames, aRColNames, aRDimNames;
        char charsTmp[MAX_LINE_LENGTH];

        nLadderLengths = ss->moveNTimes[CE];
        PROTECT(aRMat = allocMatrix(REALSXP, nLadderLengths, nLevels)); ++nProtected;
        doublesTmp = REAL(aRMat); 
        for (ii = 0; ii < nLadderLengths; ++ii) {
                /* initialize the accept ratio matrix with NAs */
                for (ll = 0; ll < nLevels; ++ll) {
                        jj             = (nLadderLengths * ll + ii);
                        doublesTmp[jj] = R_NaReal;
                }
                
                ll = ii + 3; ppc = ss->movePropCtrs[CE][ll - 1];
                for (sl = 0; sl < (ll - 1); ++sl) {
                        aR = R_NaReal;
                        if (ppc[sl]->proposed > 0) 
                                aR = ppc[sl]->accepted / ppc[sl]->proposed;

                        jj             = (nLadderLengths * sl + ii);
                        doublesTmp[jj] = aR;
                }
        }
        
        PROTECT(aRRowNames = allocVector(STRSXP, nLadderLengths)); ++nProtected;
        for (ii = 0; ii < nLadderLengths; ++ii) {
                sprintf(charsTmp, "ladderLength%d", ii + 3);
                SET_STRING_ELT(aRRowNames, ii, mkChar(charsTmp));
        }
        
        PROTECT(aRColNames = allocVector(STRSXP, nLevels)); ++nProtected;
        for (ll = 0; ll < nLevels; ++ll) {
                sprintf(charsTmp, "level%d", ll + 1);
                SET_STRING_ELT(aRColNames, ll, mkChar(charsTmp));
        }

        PROTECT(aRDimNames = allocVector(VECSXP, 2)); ++nProtected;
        SET_VECTOR_ELT(aRDimNames, 0, aRRowNames);
        SET_VECTOR_ELT(aRDimNames, 1, aRColNames);
        setAttrib(aRMat, R_DimNamesSymbol, aRDimNames);

        UNPROTECT(nProtected);
        return aRMat;
}


static SEXP
make_accept_ratios_list (Sampler *ss)
{
        int ii, nComps, comp = 0, nProtected = 0;
        SEXP aRMat, aRList, aRListNames;

        /* MH will always be used */
        nComps = 1;
        for (ii = RC; ii <= CE; ++ii)
                if (ss->moveNTimes[ii] > 0) ++nComps;

        PROTECT(aRList      = allocVector(VECSXP, nComps)); ++(ss->nProtected);
        PROTECT(aRListNames = allocVector(STRSXP, nComps)); ++nProtected;

        /* the MH component: should always be included */
        PROTECT(aRMat = make_accept_ratios_list_MH(ss)); ++nProtected;
        SET_VECTOR_ELT(aRList, comp, aRMat); 
        SET_STRING_ELT(aRListNames, comp, mkChar("MH")); ++comp;

        for (ii = RC; ii <= BSE; ++ii) {
                /* no component to be filled in if the move wasn't used */
                if (ss->moveNTimes[ii] == 0) continue;                
                
                PROTECT(aRMat = make_accept_ratios_list_move(ss, ii));
                ++nProtected;
                SET_VECTOR_ELT(aRList, comp, aRMat); 
                SET_STRING_ELT(aRListNames, comp, mkChar(ss->moveObjs[ii]->name));
                ++comp;
        }

        /* the CE component: if the move was used */
        if (ss->moveNTimes[CE] > 0) {        
                PROTECT(aRMat = make_accept_ratios_list_CE(ss)); ++nProtected;
                SET_VECTOR_ELT(aRList, comp, aRMat); 
                SET_STRING_ELT(aRListNames, comp, mkChar("CE")); ++comp;
        }
        
        setAttrib(aRList, R_NamesSymbol, aRListNames);
        UNPROTECT(nProtected);
        return aRList;
}


SEXP
sampler_run (Sampler *ss)
{
        int nComps = 0, comp = 0;
        SEXP draws, acceptRatios, acceptRatiosList;
        SEXP samplerObj, names;

        DEBUG(sampler_print(ss););
        draws                = make_draws(ss); ++nComps;
        acceptRatios         = make_accept_ratios(ss); ++nComps;
        acceptRatiosList     = make_accept_ratios_list(ss); ++nComps;
        PROTECT(samplerObj   = allocVector(VECSXP, nComps)); ++(ss->nProtected);
        PROTECT(names        = allocVector(STRSXP, nComps)); 

        /* fill in the samplerObj */
        SET_VECTOR_ELT(samplerObj, comp, draws);
        SET_STRING_ELT(names, comp, mkChar("draws")); ++comp;
        SET_VECTOR_ELT(samplerObj, comp, acceptRatios);
        SET_STRING_ELT(names, comp, mkChar("acceptRatios")); ++comp;
        SET_VECTOR_ELT(samplerObj, comp, acceptRatiosList);
        SET_STRING_ELT(names, comp, mkChar("acceptRatiosList")); ++comp;
        setAttrib(samplerObj, R_NamesSymbol, names);
        UNPROTECT(1 + ss->nProtected);
        return samplerObj;
}


SEXP
TOEMCMainC (SEXP argsList)
{
        Sampler *ss;
        SEXP samplerObj;
        
        ss = sampler_new(argsList);
        sampler_init(ss);
        PROTECT(samplerObj = sampler_run(ss));        
        UNPROTECT(1);
        return samplerObj;
}
