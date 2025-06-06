/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 * Based on code sundials_spgmr.h by: Scott D. Cohen,
 *      Alan C. Hindmarsh and Radu Serban @ LLNL
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
 * This is the header file for the SPGMR implementation of the
 * SUNLINSOL module, SUNLINSOL_SPGMR.  The SPGMR algorithm is based
 * on the Scaled Preconditioned GMRES (Generalized Minimal Residual)
 * method.
 *
 * Note:
 *   - The definition of the generic SUNLinearSolver structure can
 *     be found in the header file sundials_linearsolver.h.
 * -----------------------------------------------------------------
 */

#ifndef _SUNLINSOL_SPGMRComplex_H
#define _SUNLINSOL_SPGMRComplex_H

#include <stdio.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
#include "nvector_serialcomplex.h"
#include <sundials/sundials_math.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Default SPGMR solver parameters */
#define SUNSPGMRComplex_MAXL_DEFAULT   5
#define SUNSPGMRComplex_MAXRS_DEFAULT  0
#define SUNSPGMRComplex_GSTYPE_DEFAULT SUN_MODIFIED_GS

/* ----------------------------------------
 * SPGMR Implementation of SUNLinearSolver
 * ---------------------------------------- */

struct _SUNLinearSolverContent_SPGMRComplex
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
  // sunrealtype** Hes;
  // sunrealtype* givens;
  suncomplextype** Hes;
  suncomplextype* givens;
  N_Vector xcor;
  // sunrealtype* yg;
  suncomplextype* yg;
  N_Vector vtemp;

  // sunrealtype* cv;
  suncomplextype* cv;
  N_Vector* Xv;
};

typedef struct _SUNLinearSolverContent_SPGMRComplex* SUNLinearSolverContent_SPGMRComplex;

/* ---------------------------------------
 * Exported Functions for SUNLINSOL_SPGMR
 * --------------------------------------- */

SUNDIALS_EXPORT SUNLinearSolver SUNLinSol_SPGMRComplex(N_Vector y, int pretype,
                                                int maxl, SUNContext sunctx);
SUNDIALS_EXPORT SUNErrCode SUNLinSol_SPGMRComplex_SetPrecType(SUNLinearSolver S,
                                                      int pretype);
SUNDIALS_EXPORT SUNErrCode SUNLinSol_SPGMRComplex_SetGSType(SUNLinearSolver S,
                                                    int gstype);
SUNDIALS_EXPORT SUNErrCode SUNLinSol_SPGMRComplex_SetMaxRestarts(SUNLinearSolver S,
                                                         int maxrs);
SUNDIALS_EXPORT SUNLinearSolver_Type SUNLinSolGetType_SPGMRComplex(SUNLinearSolver S);
SUNDIALS_EXPORT SUNLinearSolver_ID SUNLinSolGetID_SPGMRComplex(SUNLinearSolver S);
SUNDIALS_EXPORT SUNErrCode SUNLinSolInitialize_SPGMRComplex(SUNLinearSolver S);
SUNDIALS_EXPORT SUNErrCode SUNLinSolSetATimes_SPGMRComplex(SUNLinearSolver S,
                                                    void* A_data,
                                                    SUNATimesFn ATimes);
SUNDIALS_EXPORT SUNErrCode SUNLinSolSetPreconditioner_SPGMRComplex(SUNLinearSolver S,
                                                            void* P_data,
                                                            SUNPSetupFn Pset,
                                                            SUNPSolveFn Psol);
SUNDIALS_EXPORT SUNErrCode SUNLinSolSetScalingVectors_SPGMRComplex(SUNLinearSolver S,
                                                            N_Vector s1,
                                                            N_Vector s2);
SUNDIALS_EXPORT SUNErrCode SUNLinSolSetZeroGuess_SPGMRComplex(SUNLinearSolver S,
                                                       sunbooleantype onff);
SUNDIALS_EXPORT int SUNLinSolSetup_SPGMRComplex(SUNLinearSolver S, SUNMatrix A);
SUNDIALS_EXPORT int SUNLinSolSolve_SPGMRComplex(SUNLinearSolver S, SUNMatrix A,
                                         N_Vector x, N_Vector b, sunrealtype tol);
SUNDIALS_EXPORT int SUNLinSolNumIters_SPGMRComplex(SUNLinearSolver S);
SUNDIALS_EXPORT sunrealtype SUNLinSolResNorm_SPGMRComplex(SUNLinearSolver S);
SUNDIALS_EXPORT N_Vector SUNLinSolResid_SPGMRComplex(SUNLinearSolver S);
SUNDIALS_EXPORT sunindextype SUNLinSolLastFlag_SPGMRComplex(SUNLinearSolver S);
SUNDIALS_EXPORT SUNErrCode SUNLinSolSpace_SPGMRComplex(SUNLinearSolver S,
                                                long int* lenrwLS,
                                                long int* leniwLS);
SUNDIALS_EXPORT SUNErrCode SUNLinSolFree_SPGMRComplex(SUNLinearSolver S);

#ifdef __cplusplus
}
#endif

#endif
