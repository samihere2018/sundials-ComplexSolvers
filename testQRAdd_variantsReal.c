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

 * The solvers tested are the different QR decomposition variants
   used in Anderson Acceleration.
 * ------------------------------------------------------------------------
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nvector_serialcomplex.h"
#include "sundials_iterativecomplex.h"
#include "sundials_iterativecomplex_impl.h" 


#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

SUNErrCode SUNQRAdd_MGSComplex(N_Vector* Q, suncomplextype* R, N_Vector df, int m,
                               int mMax, void* QRdata);

SUNErrCode SUNQRAdd_ICWYComplex(N_Vector* Q, suncomplextype* R, N_Vector df, int m,
                                int mMax, void* QRdata);

SUNErrCode SUNQRAdd_CGS2Complex(N_Vector* Q, suncomplextype* R, N_Vector df, int m,
                                int mMax, void* QRdata);

SUNErrCode SUNQRAdd_DCGS2Complex(N_Vector* Q, suncomplextype* R, N_Vector df, int m,
                                int mMax, void* QRdata);


int main(int argc, char* argv[])
{
  N_Vector* Q;
  N_Vector x;
  SUNContext sunctx;
  int k, l, m, mMax;
  suncomplextype vnorm;
  suncomplextype* temp_array;

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx))
  {
    printf("ERROR: SUNContext_Create failed\n");
    return (-1);
  }
  SUNQRData qrdata;
  qrdata = (SUNQRData)malloc(sizeof(*qrdata));

  /* Create vectors */
  x = N_VNew_SComplex(3, sunctx);
  Q = N_VCloneVectorArray(3, x);
  N_Vector vtemp = N_VClone(x);
  N_Vector vtemp2 = N_VClone(x);
  N_Vector df_True = N_VClone(x);
  N_Vector R_Approx =  N_VClone(x);
  N_Vector df_Approx =  N_VClone(x);

  qrdata->vtemp = vtemp;
  qrdata->vtemp2 = vtemp2;//together with temp_array, they are used in all the other QRAdd variants except QRAdd_MGS

  /* this stores the elements of the correction matrix (square matrix) as a one column vector by stacking the columns together starting with the first column */
  qrdata->temp_array = (suncomplextype*)malloc((9) * sizeof(suncomplextype));

  m = 2;  //number of vectors already orthogonalised (and are othornormal)
  mMax = 3; //number of rows = number of columns of the matrix Q


  /* the vector to orthogonalise */
  suncomplextype* dfdata = N_VGetArrayPointer_SComplex(df_True);
  dfdata[0] = 4.0   + 0.0*I;
  dfdata[1] = -68.0 + 0.0*I;
  dfdata[2] = -41.0 + 0.0*I;

  /* set up matrix, last column is obtained from any of the QRAdd functions */
  suncomplextype* vdata = N_VGetArrayPointer_SComplex(Q[0]);
  vdata[0] = SUN_RCONST(6.0/7.0) + 0.0*I;
  vdata[1] = (3.0/7.0)  + 0.0*I;
  vdata[2] = -(2.0/7.0) + 0.0*I;

  vdata = N_VGetArrayPointer_SComplex(Q[1]);
  vdata[0] = -(69.0/175.0) + 0.0*I;
  vdata[1] = (158.0/175.0) + 0.0*I;
  vdata[2] = (6.0/35.0)  + 0.0*I;

  vdata = N_VGetArrayPointer_SComplex(Q[2]);
  vdata[0] = 0.0 + 0.0*I;
  vdata[1] = 0.0 + 0.0*I;
  vdata[2] = 0.0 + 0.0*I;

  /* upper trinagular matrix R, the last column is obtained from any of the QRAdd functions*/
  suncomplextype R[9] = {14.0 + 0.0*I, 0.0, 0.0,
                         21.0 + 0.0*I, 175.0 + 0.0*I, 0.0,
                         0.0, 0.0, 0.0 };//R matrix stored as a column vector by stacking columns together starting with the first column


  /* perform QR decomposition using Gram-Schmidt process for df, enter any QRAdd function to use */                       
  // char functionName[100] ;  
  // printf("Enter the function name to call: \n'sunqradd_mgscomplex' or 'sunqradd_icwycomplex' or 'sunqradd_cgs2complex' or 'sunqradd_dcgs2complex': \n ");
  // scanf("%s", functionName);

  // if(strcmp(functionName,"sunqradd_mgscomplex")==0){
  //   SUNQRAdd_MGSComplex(Q, R, df_True, m, mMax, qrdata);
  //   printf("Using the function: SUNQRAdd_MGSComplex \n");
  // } 
  // else if(strcmp(functionName,"sunqradd_icwycomplex")==0){
  //   SUNQRAdd_ICWYComplex(Q, R, df_True, m, mMax, qrdata);
  //   printf("Using the function: SUNQRAdd_ICWYComplex \n");
  // } 
  // else if (strcmp(functionName,"sunqradd_cgs2complex")==0){
  //   SUNQRAdd_CGS2Complex(Q, R, df_True, m, mMax, qrdata);
  //   printf("Using the function: SUNQRAdd_CGS2Complex \n");
  // }
  // else if (strcmp(functionName,"sunqradd_dcgs2complex")==0){
  //   SUNQRAdd_DCGS2Complex(Q, R, df_True, m, mMax, qrdata);
  //   printf("Using the function: SUNQRAdd_DCGS2Complex \n");
  // }
  // else{
  //   printf("Incorrect function name! \n");
  //   return 1;
  // }

  char functionName[100] = "sunqradd_mgscomplex"; // default function name

  if (argc > 1) {
    strcpy(functionName, argv[1]); //if a function name (2nd argument) is provided after executable name
  }

  if (strcmp(functionName, "sunqradd_mgscomplex") == 0) {
    SUNQRAdd_MGSComplex(Q, R, df_True, m, mMax, qrdata);
    printf("Using SUNQRAdd_MGSComplex \n");
  }
  else if (strcmp(functionName, "sunqradd_icwycomplex") == 0) {
    SUNQRAdd_ICWYComplex(Q, R, df_True, m, mMax, qrdata);
    printf("Using SUNQRAdd_ICWYComplex \n");
  }
  else if (strcmp(functionName, "sunqradd_cgs2complex") == 0) {
    SUNQRAdd_CGS2Complex(Q, R, df_True, m, mMax, qrdata);
    printf("Using SUNQRAdd_CGS2Complex \n");
  }
  else if (strcmp(functionName, "sunqradd_dcgs2complex") == 0) {
    SUNQRAdd_DCGS2Complex(Q, R, df_True, m, mMax, qrdata);
    printf("Using SUNQRAdd_DCGS2Complex \n");
  }
  else {
    printf("Incorrect function name, use: sunqradd_mgscomplex or sunqradd_icwycomplex or sunqradd_cgs2complex or sunqradd_dcgs2complex. \nUsing default: sunqradd_mgscomplex\n");
    SUNQRAdd_MGSComplex(Q, R, df_True, m, mMax, qrdata); // Default function
  }


  /* Threshold for orthogonal vectors in matrix Q and imaginary component for the norm of a column vector in Q */
  sunrealtype tolerance = 1e-14;

  /* check dot product results */
  int unit_vectorsReal = 0;
  int unit_vectorsImag = 0;
  int orthogonalReal = 0;
  int orthogonalImag = 0;
  int solnCheckReal = 0;
  int solnCheckImag = 0;


  for (k=0; k<3; k++)
  {
    for (l=0; l<3; l++)
    {
        float vnorm = N_VDotProd_SComplex(Q[k],Q[l]);
        // printf("<Q[%i],Q[%i]> = %e + %e I\n", l, k, creal(vnorm), cimag(vnorm)); //orthonormal vectors
        if ((k==l) && (fabs(fabs(creal(vnorm))-SUN_RCONST(1.0)))>tolerance){unit_vectorsReal = 1;} //unit vectors
        if ((k==l) && (fabs(cimag(vnorm))>tolerance)){unit_vectorsImag = 1;}
        if ((k!=l) && (fabs(creal(vnorm))>tolerance)) {orthogonalReal = 1;}//orthogonal vectors
        if ((k!=l) && (fabs(cimag(vnorm))>tolerance)) {orthogonalImag = 1;}
    }
  }

  // for (k=0; k<3; k++)
  // {
  //   printf("Q[%i] = \n",k);
  //   N_VPrint(Q[k]);
  // }

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
    printf("The columns of Q are orthonormal!\n");
  } 
  else {
    printf("The columns of Q are not orthonormal!\n");
    // return 1;
  } 

  /* the last column in R */
  suncomplextype* Rdata = N_VGetArrayPointer_SComplex(R_Approx);
  Rdata[0] = R[6];
  Rdata[1] = R[7];
  Rdata[2] = R[8];

  /* use the last column in R to check if the product of the last column of Q and R gives df_True */
  suncomplextype* finalR = N_VGetArrayPointer_SComplex(df_Approx);
  finalR[0] = 0.0;
  finalR[1] = 0.0;
  finalR[2] = 0.0;

  /* multiply Q by the last column of R (the result) and the final answer should be df */
  N_VLinearCombination_SComplex(3, Rdata, Q, df_Approx);
  for (l=0;l<3;l++){
    // printf("df_Approx[%i] = %e + %e I\n", k, creal(finalR[l]), cimag(finalR[l]));
    // printf("df_True[%i] = %e + %e I\n", k, creal(dfdata[l]), cimag(dfdata[l]));
    if (fabs(creal(dfdata[l]) - creal(finalR[l]))>tolerance ){solnCheckReal = 1;}
    if (fabs(cimag(dfdata[l]) - cimag(finalR[l]))>tolerance ){solnCheckImag = 1;}
  }

  // for (l=0;l<3;l++){
  //   printf("df_Approx[%i] = %e + %e I\n", k, creal(Rdata[l]), cimag(Rdata[l]));
  // }

  /* Check if the computed last columns of Q and R are correct. */
  if ((solnCheckReal==0) && (solnCheckImag==0) && (orthogonalReal==0) && (orthogonalImag==0) && (unit_vectorsReal==0) && (unit_vectorsImag==0)) {
    printf("Test passed!\n");
  } 
  else {
    printf("Test failed!\n");
    return 1;
  }
  // else if ((solnCheckReal==1) && (solnCheckImag==0)){
  //   printf("Test failed!\n");
  // return 1;
  // } 
  // else if ((solnCheckReal==0) && (solnCheckImag==1)){
  //   printf("Test failed!\n");
  // return 1;
  // } 
  // else if ((solnCheckReal==1) && (solnCheckImag==1)){
  //   printf("Test failed!\n");
  // return 1;
  // } 


  /* return with success */
  return 0;
}
