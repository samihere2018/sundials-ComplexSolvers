/* -----------------------------------------------------------------
 * Programmer(s): Sylvia Amihere @ SMU
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
 * ------------------------------------------------------------------------
 * These test functions check some components of a complex-valued
 * SUNLINEARSOLVER module implementation (for more thorough tests,
 * see the main SUNDIALS repository, inside examples/sunlinsol/).

 * The solvers tested are the Classical Gram-Schmidt, 
 * Modified Gram-Schmidt and QR factorization (using Givens Rotation).
 * ------------------------------------------------------------------------
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nvector_serialcomplex.h"
#include "sundials_iterativecomplex.h"
// #include "sundials_iterativecomplex.h" //Amihere
#include "sundials_iterativecomplex_impl.h" //Amihere


#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

SUNErrCode SUNClassicalGSComplex(N_Vector*, suncomplextype**, int k, int p,
                          suncomplextype* new_vk_norm, suncomplextype* stemp,
                          N_Vector* vtemp);

SUNErrCode SUNModifiedGSComplex(N_Vector* v, suncomplextype** h, int k, int p,
                         suncomplextype* new_vk_norm);

int SUNQRfactComplex(int n, suncomplextype** h, suncomplextype* q, int job);

int main(int argc, char* argv[])
{
  N_Vector* V;
  N_Vector x;
  N_Vector* vtemp;
  suncomplextype* stemp;
  suncomplextype* vdata;
  suncomplextype* givens;
  SUNContext sunctx;
  suncomplextype** H;
  int k, l, job, n, krydim;
  int comp1, comp2;
  suncomplextype vnorm, ssnorm;

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx))
  {
    printf("ERROR: SUNContext_Create failed\n");
    return (-1);
  }

  /* Create vectors */
  x = N_VNew_SComplex(5, sunctx);
  V = N_VCloneVectorArray(5, x);
  N_Vector CheckV =  N_VClone(x);

  H = (suncomplextype**)malloc((5+1) * sizeof(suncomplextype*)); //Amihere
  for (k = 0; k <= 5; k++)
  {
    H[k] = NULL;
    H[k] = (suncomplextype*)malloc((5) * sizeof(suncomplextype)); //Amihere
  }

  vtemp = (N_Vector*)malloc((5) * sizeof(N_Vector)); //Amihere
  stemp = (suncomplextype*)malloc((5) * sizeof(suncomplextype)); //Amihere

  /* set up matrix */
  vdata = N_VGetArrayPointer_SComplex(V[0]);
  vdata[0] = SUN_RCONST(1.0) + 1.0*I; 
  vdata[1] = 2.0 - 2.0*I;
  vdata[2] = 3.0 + 3.0*I;
  vdata[3] = 2.0 - 3.0*I;
  vdata[4] = 6.0 + 1.0*I;
  vdata = N_VGetArrayPointer_SComplex(V[1]);
  vdata[0] = 2.0 - 1.0*I;
  vdata[1] = 1.0 + 1.0*I;
  vdata[2] = 4.0 - 3.0*I;
  vdata[3] = 1.0 - 2.0*I;
  vdata[4] = 3.0 + 5.0*I; 
  vdata = N_VGetArrayPointer_SComplex(V[2]);
  vdata[0] = 3.0 + 2.0*I;
  vdata[1] = 1.0 - 1.0*I;
  vdata[2] = 5.0 + 4.0*I;
  vdata[3] = 4.0 + 3.0*I;
  vdata[4] = 5.0 - 1.0*I;
  vdata = N_VGetArrayPointer_SComplex(V[3]);
  vdata[0] = 1.0 - 4.0*I;
  vdata[1] = 5.0 + 2.0*I;
  vdata[2] = 3.0 - 4.0*I;
  vdata[3] = 1.0 - 3.0*I;
  vdata[4] = 1.0 + 1.0*I;
  vdata = N_VGetArrayPointer_SComplex(V[4]);
  vdata[0] = 4.0 + 5.0*I;
  vdata[1] = 6.0 + 3.0*I;
  vdata[2] = 2.0 - 1.0*I;
  vdata[3] = 2.0 + 5.0*I;
  vdata[4] = 3.0 - 5.0*I;

  // /*Amihere: for debugging*/
  // suncomplextype* dfdata = N_VGetArrayPointer_SComplex(CheckV); 
  // dfdata[0] = 0.0;
  // dfdata[1] = 0.0;
  // dfdata[2] = 0.0;
  // dfdata[3] = 0.0;
  // dfdata[4] = 0.0;

  /* perform Gram-Schmidt process for all vectors in V */
  char functionName[100] = "ClassicalGS"; // default function name

  if (argc > 1) {
    strcpy(functionName, argv[1]); //if a function name (2nd argument) is provided after executable name
  }

  if (strcmp(functionName, "ClassicalGS") == 0) {
    printf("Using Classical Gram Schmidt \n");
    for (k=0; k<5; k++){
      SUNClassicalGSComplex(V, H, k, 5, &vnorm, stemp, vtemp);
      N_VScale_SComplex(1.0/vnorm, V[k], V[k]);
    }
  }
  else if (strcmp(functionName, "ModifiedGS") == 0) {
    printf("Using Modified Gram Schmidt \n");
    for (k=0; k<5; k++){
      SUNModifiedGSComplex(V, H, k, 5, &vnorm);
      N_VScale_SComplex(1.0/vnorm, V[k], V[k]);
    }
  }
  else {
    printf("Incorrect function name, use: ClassicalGS or ModifiedGS. \nUsing default: ClassicalGS\n");
    for (k=0; k<5; k++){
      SUNClassicalGSComplex(V, H, k, 5, &vnorm, stemp, vtemp);// Default function
      N_VScale_SComplex(1.0/vnorm, V[k], V[k]);
    }
  }

   /* Threshold for orthogonal vectors in matrix Q and imaginary component for the norm of a column vector in Q */
   sunrealtype tolerance = 1e-14;

   /* check dot product results */
   int unit_vectorsReal = 0;
   int unit_vectorsImag = 0;
   int orthogonalReal = 0;
   int orthogonalImag = 0;

  /* check dot product results */
  for (k=0; k<5; k++)
  {
    for (l=0; l<5; l++)
    {
        float vnorm = N_VDotProd_SComplex(V[k],V[l]);
        // printf("<V[%i],V[%i]> = %e + %e I\n", l, k, creal(vnorm), cimag(vnorm));
        if ((k==l) && (fabs(fabs(creal(vnorm))-SUN_RCONST(1.0)))>tolerance){unit_vectorsReal = 1;} //unit vectors
        if ((k==l) && (fabs(cimag(vnorm))>tolerance)){unit_vectorsImag = 1;}
        if ((k!=l) && (fabs(creal(vnorm))>tolerance)) {orthogonalReal = 1;}//orthogonal vectors
        if ((k!=l) && (fabs(cimag(vnorm))>tolerance)) {orthogonalImag = 1;}
    }
  }

  /* Check if the columns of Q are unit vectors. */
  if ((unit_vectorsReal==0) && (unit_vectorsImag==0)) {
    printf("The columns of Q are unit vectors!\n");
  } 
  else if ((unit_vectorsReal==1) && (unit_vectorsImag==0)){
    printf("The columns of Q are not unit vectors!\n");
    // return 1;
  } 
  else if ((unit_vectorsReal==0) && (unit_vectorsImag==1)){
    printf("The columns of Q are not unit vectors!\n");
    // return 1;
  } 
  else if ((unit_vectorsReal==1) && (unit_vectorsImag==1)){
    printf("The columns of Q are not unit vectors!\n");
    // return 1;
  } 

  /* Check if the columns of Q are orthogonal. */
  if ((orthogonalReal==0) && (orthogonalImag==0)) {
    printf("The columns of Q are orthogonal!\n");
  } 
  else if ((orthogonalReal==1) && (orthogonalImag==0)){
    printf("The columns of Q are not orthogonal!\n");
    // return 1;
  } 
  else if ((orthogonalReal==0) && (orthogonalImag==1)){
    printf("The columns of Q are not orthogonal!\n");
    // return 1;
  } 
  else if ((orthogonalReal==1) && (orthogonalImag==1)){
    printf("The columns of Q are not orthogonal!\n");
    // return 1;
  } 

  /* Check if the columns of Q are orthonormal. */
  if ((orthogonalReal==0) && (orthogonalImag==0) && (unit_vectorsReal==0) && (unit_vectorsImag==0)){
    printf("The columns of Q are orthonormal!\nTest Passed!\n");
  } 
  else {
    printf("The columns of Q are not orthonormal!\nTest failed!\n");
    // return 1;
  } 

  /* return with success */
  return 0;
}
