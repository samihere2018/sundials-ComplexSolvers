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
  x = N_VNew_SComplex(5, sunctx);
  Q = N_VCloneVectorArray(5, x);
  N_Vector vtemp = N_VClone(x);
  N_Vector vtemp2 = N_VClone(x);
  N_Vector df_True = N_VClone(x);
  N_Vector R_Approx =  N_VClone(x);
  N_Vector df_Approx =  N_VClone(x);

  qrdata->vtemp = vtemp;
  qrdata->vtemp2 = vtemp2;//together with temp_array, they are used in all the other QRAdd variants except QRAdd_MGS

  /* this stores the elements of the correction matrix (square matrix) as a one column vector by stacking the columns together starting with the first column */
  qrdata->temp_array = (suncomplextype*)malloc((16) * sizeof(suncomplextype));

  m = 4;  //number of vectors already orthogonalised (and are othornormal)
  mMax = 5; //number of rows = number of columns of the matrix Q


  /* the vector to orthogonalise */
  suncomplextype* dfdata = N_VGetArrayPointer_SComplex(df_True);
  dfdata[0] = 4.0+5.0*I;
  dfdata[1] = 6.0+3.0*I;
  dfdata[2] = 2.0-1.0*I;
  dfdata[3] = 2.0+5.0*I;
  dfdata[4] = 3.0-5.0*I;

  /* set up matrix, last column is obtained from any of the QRAdd functions */
  suncomplextype* vdata = N_VGetArrayPointer_SComplex(Q[0]);
  vdata[0] = SUN_RCONST(0.113227703414460) + 0.226455406828919*I;
  vdata[1] = 0.226455406828919 + 0.339683110243379*I;
  vdata[2] = 0.339683110243379 + 0.113227703414460*I;
  vdata[3] = 0.226455406828919 - 0.339683110243379*I;
  vdata[4] = 0.679366220486758 + 0.113227703414460*I;

  vdata = N_VGetArrayPointer_SComplex(Q[1]);
  vdata[0] = 0.358047898868247 - 0.209079065032553*I;
  vdata[1] = 0.449519989819989 + 0.164649763713135*I;
  vdata[2] = 0.352820922242433 + 0.044429301319417*I;
  vdata[3] = -0.334526504052085 - 0.007840464938720*I;
  vdata[4] = -0.376342317058596 + 0.467814408010337*I;

  vdata = N_VGetArrayPointer_SComplex(Q[2]);
  vdata[0] = 0.368417696619559 + 0.108720463349240*I;
  vdata[1] = 0.382885110056019 - 0.076920802132466*I;
  vdata[2] = -0.108648842490643 + 0.349438169091529*I;
  vdata[3] = 0.326877598633683 + 0.632698664840043*I;
  vdata[4] = 0.056007511422335 - 0.236062349933527*I;

  vdata = N_VGetArrayPointer_SComplex(Q[3]);
  vdata[0] = -0.173120531596438 - 0.317326783017719*I;
  vdata[1] = 0.305340355271806 + 0.559154947299423*I;
  vdata[2] = -0.270428892207880 - 0.452932177178935*I;
  vdata[3] = 0.395829946721830 + 0.018686126433033*I;
  vdata[4] = -0.144400733195965 - 0.085349715976197*I;

  vdata = N_VGetArrayPointer_SComplex(Q[4]);
  vdata[0] = 0.0 + 0.0*I;
  vdata[1] = 0.0 + 0.0*I;
  vdata[2] = 0.0 + 0.0*I;
  vdata[3] = 0.0 + 0.0*I;
  vdata[4] = 0.0 + 0.0*I;

  /* upper trinagular matrix R, the last column is obtained from any of the QRAdd functions*/
  suncomplextype R[25] = {8.831760866327848+0.0*I, 0.0, 0.0, 0.0, 0.0,
                         7.586256128768794+2.717464881947030*I, 4.905517563326268-0.000000000000002*I, 0.0, 0.0, 0.0,
                         5.887840577551898+0.113227703414459*I, -0.606329288594409-0.786659982184980*I, 7.438685615532769+0.000000000000000*I, 0.0, 0.0,
                         3.623286509262707 - 3.396831102433787*I, 4.432476178690121 - 2.524629710268074*I, -1.780852648997921 - 2.366997755750347*I,  4.851661421774982 + 0.000000000000001*I, 0.0,
                         0.0, 0.0, 0.0, 0.0, 0.0 };//R matrix stored as a column vector by stacking columns together starting with the first column


  /* perform QR decomposition using Gram-Schmidt process for df, enter any QRAdd function to use */                       
  char functionName[100] = "sunqradd_mgscomplex"; // default function name

  if (argc > 1) {
    strcpy(functionName, argv[1]); //if a function name (2nd argument) is provided after executable name
  }

  if (strcmp(functionName, "sunqradd_mgscomplex") == 0) {
    printf("Using SUNQRAdd_MGSComplex \n");
    SUNQRAdd_MGSComplex(Q, R, df_True, m, mMax, qrdata);
  }
  else if (strcmp(functionName, "sunqradd_icwycomplex") == 0) {
    printf("Using SUNQRAdd_ICWYComplex \n");
    SUNQRAdd_ICWYComplex(Q, R, df_True, m, mMax, qrdata);
  }
  else if (strcmp(functionName, "sunqradd_cgs2complex") == 0) {
    printf("Using SUNQRAdd_CGS2Complex \n");
    SUNQRAdd_CGS2Complex(Q, R, df_True, m, mMax, qrdata);
  }
  else if (strcmp(functionName, "sunqradd_dcgs2complex") == 0) {
    printf("Using SUNQRAdd_DCGS2Complex \n");
    SUNQRAdd_DCGS2Complex(Q, R, df_True, m, mMax, qrdata);
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


  for (k=0; k<5; k++)
  {
    for (l=0; l<5; l++)
    {
        float vnorm = N_VDotProd_SComplex(Q[k],Q[l]);
        // printf("<Q[%i],Q[%i]> = %e + %e I\n", l, k, creal(vnorm), cimag(vnorm)); //orthonormal vectors
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
    printf("The columns of Q are orthonormal!\n");
  } 
  else {
    printf("The columns of Q are not orthonormal!\n");
    // return 1;
  } 

  /* the last column in R */
  suncomplextype* Rdata = N_VGetArrayPointer_SComplex(R_Approx);
  Rdata[0] = R[20];
  Rdata[1] = R[21];
  Rdata[2] = R[22];
  Rdata[3] = R[23];
  Rdata[4] = R[24];

  /* use the last column in R to check if the product of the last column of Q and R gives df_True */
  suncomplextype* finalR = N_VGetArrayPointer_SComplex(df_Approx);
  finalR[0] = 0.0;
  finalR[1] = 0.0;
  finalR[2] = 0.0;
  finalR[3] = 0.0;
  finalR[4] = 0.0;

  /* multiply Q by the last column of R (the result) and the final answer should be df */
  N_VLinearCombination_SComplex(5, Rdata, Q, df_Approx);
  for (l=0;l<5;l++){
    // printf("df_Approx[%i] = %e + %e I\n", k, creal(finalR[l]), cimag(finalR[l]));
    // printf("df_True[%i] = %e + %e I\n", k, creal(dfdata[l]), cimag(dfdata[l]));
    if (fabs(creal(dfdata[l]) - creal(finalR[l]))>tolerance ){solnCheckReal = 1;}
    if (fabs(cimag(dfdata[l]) - cimag(finalR[l]))>tolerance ){solnCheckImag = 1;}
  }

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
