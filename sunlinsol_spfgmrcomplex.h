/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 * Based on code sundials_spfgmr.h by: Daniel R. Reynolds and
 *    Hilari C. Tiedeman @ SMU
 * Edited by Sylvia Amihere @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for the SPFGMR implementation of the
 * SUNLINSOL module, SUNLINSOL_SPFGMR.  The SPFGMR algorithm is based
 * on the Scaled Preconditioned FGMRES (Flexible Generalized Minimal
 * Residual) method [Y. Saad, SIAM J. Sci. Comput., 1993].
 *
 * Note:
 *   - The definition of the generic SUNLinearSolver structure can
 *     be found in the header file sundials_linearsolver.h.
 * -----------------------------------------------------------------
 */

#ifndef _SUNLINSOL_SPFGMRComplex_H
#define _SUNLINSOL_SPFGMRComplex_H

#include <stdio.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
#include "nvector_serialcomplex.h"
#include <sundials/sundials_math.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Default SPFGMR solver parameters */
#define SUNSPFGMRComplex_MAXL_DEFAULT   5
#define SUNSPFGMRComplex_MAXRS_DEFAULT  0
#define SUNSPFGMRComplex_GSTYPE_DEFAULT SUN_MODIFIED_GS

/* -----------------------------------------
 * SPFGMR Implementation of SUNLinearSolver
 * ----------------------------------------- */

struct _SUNLinearSolverContent_SPFGMRComplex
{
  int maxl;
  int pretype;
  int gstype;
  int max_restarts;
  sunbooleantype zeroguess;
  int numiters;
  sunrealtype resnorm;
  int last_flag;

  SUNATimesFn ATimes;
  void* ATData;
  SUNPSetupFn Psetup;
  SUNPSolveFn Psolve;
  void* PData;

  N_Vector s1;
  N_Vector s2;
  N_Vector* V;
  N_Vector* Z;
  //sunrealtype** Hes;
  //sunrealtype* givens;
  suncomplextype** Hes;
  suncomplextype* givens;
  N_Vector xcor;
  //sunrealtype* yg;
  suncomplextype* yg;
  N_Vector vtemp;

  //sunrealtype* cv;
  suncomplextype* cv;
  N_Vector* Xv;
};

typedef struct _SUNLinearSolverContent_SPFGMRComplex* SUNLinearSolverContent_SPFGMRComplex;

/* ----------------------------------------
 * Exported Functions for SUNLINSOL_SPFGMR
 * ---------------------------------------- */

SUNDIALS_EXPORT SUNLinearSolver SUNLinSol_SPFGMRComplex(N_Vector y, int pretype,
                                                 int maxl, SUNContext sunctx);
SUNDIALS_EXPORT SUNErrCode SUNLinSol_SPFGMRComplex_SetPrecType(SUNLinearSolver S,
                                                       int pretype);
SUNDIALS_EXPORT SUNErrCode SUNLinSol_SPFGMRComplex_SetGSType(SUNLinearSolver S,
                                                     int gstype);
SUNDIALS_EXPORT SUNErrCode SUNLinSol_SPFGMRComplex_SetMaxRestarts(SUNLinearSolver S,
                                                          int maxrs);
SUNDIALS_EXPORT SUNLinearSolver_Type SUNLinSolGetType_SPFGMRComplex(SUNLinearSolver S);
SUNDIALS_EXPORT SUNLinearSolver_ID SUNLinSolGetID_SPFGMRComplex(SUNLinearSolver S);
SUNDIALS_EXPORT SUNErrCode SUNLinSolInitialize_SPFGMRComplex(SUNLinearSolver S);
SUNDIALS_EXPORT SUNErrCode SUNLinSolSetATimes_SPFGMRComplex(SUNLinearSolver S,
                                                     void* A_data,
                                                     SUNATimesFn ATimes);
SUNDIALS_EXPORT SUNErrCode SUNLinSolSetPreconditioner_SPFGMRComplex(SUNLinearSolver S,
                                                             void* P_data,
                                                             SUNPSetupFn Pset,
                                                             SUNPSolveFn Psol);
SUNDIALS_EXPORT SUNErrCode SUNLinSolSetScalingVectors_SPFGMRComplex(SUNLinearSolver S,
                                                             N_Vector s1,
                                                             N_Vector s2);
SUNDIALS_EXPORT SUNErrCode SUNLinSolSetZeroGuess_SPFGMRComplex(SUNLinearSolver S,
                                                        sunbooleantype onoff);
SUNDIALS_EXPORT int SUNLinSolSetup_SPFGMRComplex(SUNLinearSolver S, SUNMatrix A);
SUNDIALS_EXPORT int SUNLinSolSolve_SPFGMRComplex(SUNLinearSolver S, SUNMatrix A,
                                          N_Vector x, N_Vector b,
                                          sunrealtype tol);
SUNDIALS_EXPORT int SUNLinSolNumIters_SPFGMRComplex(SUNLinearSolver S);
SUNDIALS_EXPORT sunrealtype SUNLinSolResNorm_SPFGMRComplex(SUNLinearSolver S);
SUNDIALS_EXPORT N_Vector SUNLinSolResid_SPFGMRComplex(SUNLinearSolver S);
SUNDIALS_EXPORT sunindextype SUNLinSolLastFlag_SPFGMRComplex(SUNLinearSolver S);
SUNDIALS_EXPORT SUNErrCode SUNLinSolSpace_SPFGMRComplex(SUNLinearSolver S,
                                                 long int* lenrwLS,
                                                 long int* leniwLS);
SUNDIALS_EXPORT SUNErrCode SUNLinSolFree_SPFGMRComplex(SUNLinearSolver S);

#ifdef __cplusplus
}
#endif

#endif
