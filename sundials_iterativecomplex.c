/* -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 *                Shelby Lockhart @ LLNL
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
 * This is the implementation file for the iterative.h header
 * file. It contains the implementation of functions that may be
 * useful for many different iterative solvers of A x = b.
 * -----------------------------------------------------------------*/

#include <stdio.h>
// #include <sundials/priv/sundials_errors_impl.h>
// #include <sundials/sundials_core.h>

// #include "sundials/sundials_errors.h"
#include "sundials_iterativecomplex_impl.h"

#define FACTOR SUN_RCONST(1000.0)
#define ZERO   SUN_RCONST(0.0)
#define ONE    SUN_RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * Function : SUNModifiedGS
 * -----------------------------------------------------------------
 * This implementation of SUNModifiedGS is a slight modification of
 * a previous modified Gram-Schmidt routine (called mgs) written by
 * Milo Dorr.
 * -----------------------------------------------------------------
 */

// SUNErrCode SUNModifiedGS(N_Vector* v, sunrealtype** h, int k, int p,
//                          sunrealtype* new_vk_norm)
SUNErrCode SUNModifiedGSComplex(N_Vector* v, suncomplextype** h, int k, int p,
                         suncomplextype* new_vk_norm)
{
  // SUNFunctionBegin(v[0]->sunctx);

  int i, k_minus_1, i0;
  // sunrealtype new_norm_2, new_product, vk_norm, temp;
  suncomplextype new_product, temp;
  sunrealtype new_norm_2, vk_norm;

  vk_norm = N_VDotProd_SComplex(v[k], v[k]);
  // SUNCheckLastErr();
  vk_norm   = SUNRsqrt(vk_norm);
  k_minus_1 = k - 1;
  i0        = SUNMAX(k - p, 0);

  /* Perform modified Gram-Schmidt */

  for (i = i0; i < k; i++)
  {
    // h[i][k_minus_1] = N_VDotProd(v[i], v[k]);
    // h[i][k_minus_1] = N_VDotProd_SComplex(v[k], v[i]); //Amihere: original
    h[i][k_minus_1] = N_VDotProd_SComplex(v[i], v[k]); //Amihere: we use this instead of original because we need to take the complex conjuagte of the bases we are projecting onto
    // printf("h[%i][%i]=<v[%i],v[%i]> = %e + %e I\n", i, k_minus_1, k, i, creal(h[i][k_minus_1]), cimag(h[i][k_minus_1]));//Amihere:debugging
    // SUNCheckLastErr();
    N_VLinearSum_SComplex(ONE, v[k], -h[i][k_minus_1], v[i], v[k]);
    // printf("V[%i] = \n",k); //Amihere: debugging
    // N_VPrint(v[k]); //Amihere: debugging
    // SUNCheckLastErr();
  }

  /* Compute the norm of the new vector at v[k] */

  *new_vk_norm = N_VDotProd_SComplex(v[k], v[k]);
  // SUNCheckLastErr();
  *new_vk_norm = SUNRsqrt((sunrealtype)*new_vk_norm);

  /* If the norm of the new vector at v[k] is less than
     FACTOR (== 1000) times unit roundoff times the norm of the
     input vector v[k], then the vector will be reorthogonalized
     in order to ensure that nonorthogonality is not being masked
     by a very small vector length. */

  temp = FACTOR * vk_norm;
  if ((temp + (*new_vk_norm)) != temp) { return SUN_SUCCESS; }

  new_norm_2 = ZERO;

  for (i = i0; i < k; i++)
  {
    // new_product = N_VDotProd(v[i], v[k]);
    // new_product = N_VDotProd_SComplex(v[k], v[i]); //Amihere: original
    new_product = N_VDotProd_SComplex(v[i], v[k]); //Amihere: we use this instead of original because we need to take the complex conjuagte of the bases we are projecting onto
    // printf("newProd = <v[%i],v[%i]> = %e + %e I\n", k, i, creal(new_product), cimag(new_product));//Amihere:debugging
    // SUNCheckLastErr();
    temp = FACTOR * h[i][k_minus_1];
    // if ((temp + new_product) == temp) { continue; }
    if ((creal(temp) + creal(new_product) == creal(temp)) && (cimag(temp) + cimag(new_product) == cimag(temp))) { continue; }
    h[i][k_minus_1] += new_product;
    N_VLinearSum_SComplex(ONE, v[k], -new_product, v[i], v[k]);
    // SUNCheckLastErr();
    // new_norm_2 += SUNSQR(new_product);
    new_norm_2 += new_product * conj(new_product);
  }

  if (new_norm_2 != ZERO)
  {
    new_product  = SUNSQR(*new_vk_norm) - new_norm_2;
    *new_vk_norm = ((sunrealtype) new_product > ZERO) ? SUNRsqrt((sunrealtype) new_product) : ZERO;
  }

  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * Function : SUNClassicalGS
 * -----------------------------------------------------------------
 * This implementation of SUNClassicalGS was contributed by Homer
 * Walker and Peter Brown.
 * -----------------------------------------------------------------
 */

// SUNErrCode SUNClassicalGS(N_Vector* v, sunrealtype** h, int k, int p,
//                           sunrealtype* new_vk_norm, sunrealtype* stemp,
//                           N_Vector* vtemp)
SUNErrCode SUNClassicalGSComplex(N_Vector* v, suncomplextype** h, int k, int p,
                          suncomplextype* new_vk_norm, suncomplextype* stemp,
                          N_Vector* vtemp)
{
  // SUNFunctionBegin(v[0]->sunctx);
  int i, i0, k_minus_1;
  sunrealtype vk_norm;

  k_minus_1 = k - 1;
  i0        = SUNMAX(k - p, 0);

  /* Perform Classical Gram-Schmidt */

  N_VDotProdMulti_SComplex(k - i0 + 1, v[k], v + i0, stemp);
  // for (int ij=0; ij<k - i0 + 1; ij++) //Amihere:debugging
  // {
  //   printf("stemp[%i] =  %e + %e I \n",ij, creal(stemp[ij]), cimag(stemp[ij]));
  // }
  vk_norm = SUNRsqrt((sunrealtype)(stemp[k - i0] * conj(stemp[k - i0])));
  for (i = k - i0 - 1; i >= 0; i--)
  {
    h[i][k_minus_1] = stemp[i];
    stemp[i + 1]    = -stemp[i];
    vtemp[i + 1]    = v[i];
  }
  stemp[0] = ONE;
  vtemp[0] = v[k];

  // SUNCheckCall(N_VLinearCombination_SComplex(k - i0 + 1, stemp, vtemp, v[k]));
  N_VLinearCombination_SComplex(k - i0 + 1, stemp, vtemp, v[k]);
  

  /* Compute the norm of the new vector at v[k] */

  *new_vk_norm = SUNRsqrt((sunrealtype) N_VDotProd_SComplex(v[k], v[k]));
  // N_VScale_SComplex(1.0/(*new_vk_norm), v[k], v[k]);//Amihere: use when running unit test for Gram Schmidt
  // SUNCheckLastErr();

  /* Reorthogonalize if necessary */

  if ((FACTOR * ((sunrealtype)*new_vk_norm)) < vk_norm)
  {
    // SUNCheckCall(N_VDotProdMulti_SComplex(k - i0, v[k], v + i0, stemp + 1));
    //Sylvia: this line is supposed to be "N_VDotProdMulti_SComplex(k - i0, v + i0, v[k], stemp + 1);" because v+io is Av and v[k] is v.
    //Sylvia: the algorithm is dotprod of Av and v not dotprod of v and Av.
    //Sylvia: because of the structure of dotprodmulti, we have to do dotprod of v and Av. Since we are dealing with complex values, we just have to take the complex conjugate of the second entry Av and multiply by the second entry.
    N_VDotProdMulti_SComplex(k - i0, v[k], v + i0, stemp + 1);

    stemp[0] = ONE;
    vtemp[0] = v[k];
    for (i = i0; i < k; i++)
    {
      h[i][k_minus_1] += stemp[i - i0 + 1];
      stemp[i - i0 + 1] = -stemp[i - i0 + 1];
      vtemp[i - i0 + 1] = v[i - i0];
    }

    // SUNCheckCall(N_VLinearCombination_SComplex(k + 1, stemp, vtemp, v[k]));
    N_VLinearCombination_SComplex(k + 1, stemp, vtemp, v[k]);

    *new_vk_norm = SUNRsqrt((sunrealtype)N_VDotProd_SComplex(v[k], v[k]));
    // SUNCheckLastErr();
  }

  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * Function : SUNQRfact
 * -----------------------------------------------------------------
 * This implementation of SUNQRfact is a slight modification of a
 * previous routine (called qrfact) written by Milo Dorr.
 * -----------------------------------------------------------------
 */

// int SUNQRfact(int n, sunrealtype** h, sunrealtype* q, int job)
int SUNQRfactComplex(int n, suncomplextype** h, suncomplextype* q, int job)
{
  // sunrealtype c, s, temp1, temp2, temp3;
  suncomplextype c, s, temp1, temp2, temp3;
  suncomplextype RotationColumn1, RotationColumn2, RotationColumn, signf, denomCS, sumCS;
  suncomplextype absTemp1, absTemp2;
  int i, j, k, q_ptr, n_minus_1, code = 0;

  switch (job)
  {
  case 0:

    /* Compute a new factorization of H */

    code = 0;
    for (k = 0; k < n; k++)
    {
      /* Multiply column k by the previous k-1 Givens rotations */

      for (j = 0; j < k - 1; j++)
      {
        i           = 2 * j;
        temp1       = h[j][k];
        temp2       = h[j + 1][k];
        c           = q[i];
        s           = q[i + 1];
        h[j][k]     = c * temp1 - s * temp2;
        // h[j + 1][k] = s * temp1 + c * temp2;
        h[j + 1][k] = conj(s) * temp1 + conj(c) * temp2;
      }
      // // Amihere: debugging
      // printf ("rotation c: %f + i %f\n", creal(c), cimag(c));
      // printf ("rotation s: %f + i %f\n", creal(s), cimag(s));
      // printf ("rotation q[i]: %f + i %f\n", creal(q[1]), cimag(q[1]));
      // printf ("rotation s: %f + i %f\n", creal(s), cimag(s));
      // RotationColumn1 = c * -s;
      // RotationColumn2 = conj(s) * conj(c); //((c * -s) + (conj(s) * c));
      // RotationColumn = RotationColumn1 + RotationColumn2;
      // printf("Rotation value: %f + i %f\n", creal(RotationColumn), cimag(RotationColumn));

      /* Compute the Givens rotation components c and s */

      q_ptr = 2 * k;
      temp1 = h[k][k];
      temp2 = h[k + 1][k];
      if (temp2 == ZERO)
      {
        c = ONE;
        s = ZERO;
      }
      // else if(temp1 == ZERO)
      // {
      //   c = ZERO;
      //   s = ONE;
      // }
      // else if (SUNRabs(temp2) >= SUNRabs(temp1))
      else if (SUNRsqrt((sunrealtype)(conj(temp2) * temp2)) >= SUNRsqrt((sunrealtype)(conj(temp1) * temp1)))
      {
        // temp3 = temp1 / temp2;
        // temp3 = (conj(temp1)*temp1) / (conj(temp2)*temp2); //Amihere (use this and next four lines)
        // absTemp2 = SUNRsqrt((sunrealtype)(conj(temp2)*temp2));
        // absTemp1 = SUNRsqrt((sunrealtype)(conj(temp1)*temp1));
        // c = ONE / ((absTemp2/absTemp1) * SUNRsqrt(ONE + (sunrealtype)temp3)) ;
        // s = -c * (temp1 * conj(temp2)) / ((sunrealtype)(conj(temp1)*temp1));

        //New - 30th March 2025
        temp3 = (conj(temp1)*temp1) / (conj(temp2)*temp2); //Amihere (use this and next three lines)
        absTemp2 = SUNRsqrt((sunrealtype)(conj(temp2)*temp2));
        s =  -ONE / ((temp2/absTemp2) * SUNRsqrt(ONE + (sunrealtype)temp3)) ;
        c = -s * (conj(temp1)/conj(temp2));


        // s     = -ONE / SUNRsqrt(ONE + SUNSQR((sunrealtype)temp3));
        // // s     = -ONE / SUNRsqrt(ONE + ((sunrealtype)(conj(temp3) * temp3)));
        // c     = -s * temp3; 

        // denomCS = SUNRsqrt((sunrealtype)((conj(temp1)*temp1)+(conj(temp2)*temp2))); //Amihere (use this and next two lines)
        // c = conj(temp1)/denomCS;
        // s = -conj(temp2)/denomCS;

        // c  = SUNRsqrt((sunrealtype)(conj(temp1)*temp1))/denomCS; 
        // signf = temp1/SUNRsqrt((sunrealtype)(conj(temp1)*temp1));
        // s  = -signf * conj(temp2)/denomCS; 
        // s = c * temp3;
      }
      else 
      {
        // temp3 = temp2 / temp1;
        // temp3 = (conj(temp2)*temp2) / (conj(temp1)*temp1); //Amihere (use this and next two lines)
        // c  = ONE / SUNRsqrt(ONE + (sunrealtype)temp3);
        // s = -c * (temp1 * conj(temp2)) / ((sunrealtype)(conj(temp1)*temp1));

        //New - 30th March 2025
        temp3 = (conj(temp2)*temp2) / (conj(temp1)*temp1); //Amihere (use this and next three lines)
        absTemp1 = SUNRsqrt((sunrealtype)(conj(temp1)*temp1));
        c  = ONE / ((temp1/absTemp1)* SUNRsqrt(ONE + (sunrealtype)temp3));
        s = -c * (conj(temp2)/conj(temp1));

        // c     = ONE / SUNRsqrt(ONE + SUNSQR((sunrealtype)temp3));
        // // c     = ONE / SUNRsqrt(ONE + ((sunrealtype)(conj(temp3) * temp3)));
        // s     = -c * temp3;

        // denomCS = SUNRsqrt((sunrealtype)((conj(temp1)*temp1)+(conj(temp2)*temp2))); //Amihere (use this and next two lines)
        // c = conj(temp1)/denomCS;
        // s = -conj(temp2)/denomCS;

        // signf = temp1/SUNRsqrt((sunrealtype)(conj(temp1)*temp1));
        // s  = -signf * conj(temp2)/denomCS; 
        // s = c * temp3;
        
      }
      q[q_ptr]     = c;
      q[q_ptr + 1] = s;
      if ((h[k][k] = c * temp1 - s * temp2) == ZERO) { code = k + 1; }
    }

    // // Amihere: debugging
    //   sumCS = (c * c) + (s * s);
    //   printf ("sum of c and s: %f + i %f\n", creal(sumCS), cimag(sumCS));
    //   printf ("rotation c: %f + i %f\n", creal(c), cimag(c));
    //   printf ("rotation s: %f + i %f\n", creal(s), cimag(s));
    //   printf ("rotation q[i]: %f + i %f\n", creal(q[1]), cimag(q[1]));
    //   printf ("rotation s: %f + i %f\n", creal(s), cimag(s));
    //   RotationColumn1 = c * -s;
    //   RotationColumn2 = conj(s) * c; //((c * -s) + (conj(s) * c));
    //   RotationColumn = RotationColumn1 + RotationColumn2;
    //   printf("Rotation value: %f + i %f\n", creal(RotationColumn), cimag(RotationColumn));
    break;

  default:

    /* Update the factored H to which a new column has been added */

    n_minus_1 = n - 1;
    code      = 0;

    /* Multiply the new column by the previous n-1 Givens rotations */

    for (k = 0; k < n_minus_1; k++)
    {
      i                   = 2 * k;
      temp1               = h[k][n_minus_1];
      temp2               = h[k + 1][n_minus_1];
      c                   = q[i];
      s                   = q[i + 1];
      h[k][n_minus_1]     = c * temp1 - s * temp2;
      h[k + 1][n_minus_1] = conj(s) * temp1 + conj(c) * temp2;
    }

    // // Amihere: debugging
    //   // sumCS = (c * c) + (s * s);
    //   // printf ("sum of c and s: %f + i %f\n", creal(sumCS), cimag(sumCS));
    //   printf ("rotation c: %f + i %f\n", creal(c), cimag(c));
    //   printf ("rotation s: %f + i %f\n", creal(s), cimag(s));
    //   printf ("rotation q[i]: %f + i %f\n", creal(q[1]), cimag(q[1]));
    //   printf ("rotation s: %f + i %f\n", creal(s), cimag(s));
    //   RotationColumn1 = c * -s;
    //   RotationColumn2 = conj(s) * conj(c); //((c * -s) + (conj(s) * c));
    //   RotationColumn = RotationColumn1 + RotationColumn2;
    //   printf("Rotation value: %f + i %f\n", creal(RotationColumn), cimag(RotationColumn));

    /* Compute new Givens rotation and multiply it times the last two
       entries in the new column of H.  Note that the second entry of
       this product will be 0, so it is not necessary to compute it. */

    temp1 = h[n_minus_1][n_minus_1];
    temp2 = h[n][n_minus_1];
    if (temp2 == ZERO)
    {
      c = ONE;
      s = ZERO;
    }
    // else if (temp1 == ZERO)
    // {
    //   c = ZERO;
    //   s = ONE;
    // }
    // else if (SUNRabs(temp2) >= SUNRabs(temp1))
    else if (SUNRsqrt((sunrealtype)(conj(temp2) * temp2)) >= SUNRsqrt((sunrealtype)(conj(temp1) * temp1)))
    {
      // // temp3 = temp1 / temp2;
      // temp3 = (conj(temp1)*temp1) / (conj(temp2)*temp2); //Amihere (use this and next four lines)
      // absTemp2 = SUNRsqrt((sunrealtype)(conj(temp2)*temp2));
      // absTemp1 = SUNRsqrt((sunrealtype)(conj(temp1)*temp1));
      // c = ONE / ((absTemp2/absTemp1) * SUNRsqrt(ONE + (sunrealtype)temp3)) ;
      // s = -c * (temp1 * conj(temp2)) / ((sunrealtype)(conj(temp1)*temp1));

       //New - 30th March 2025
       temp3 = (conj(temp1)*temp1) / (conj(temp2)*temp2); //Amihere (use this and next three lines)
       absTemp2 = SUNRsqrt((sunrealtype)(conj(temp2)*temp2));
       s =  -ONE / ((temp2/absTemp2) * SUNRsqrt(ONE + (sunrealtype)temp3)) ;
       c = -s * (conj(temp1)/conj(temp2));

      // s     = -ONE / SUNRsqrt(ONE + SUNSQR((sunrealtype)temp3));
      // // s     = -ONE / SUNRsqrt(ONE + ((sunrealtype)(conj(temp3) * temp3)));
      // c     = -s * temp3; 

      // denomCS = SUNRsqrt((sunrealtype)((conj(temp1)*temp1)+(conj(temp2)*temp2))); //Amihere (use this and next two lines)
      // c = conj(temp1)/denomCS;
      // s = -conj(temp2)/denomCS;
     
      // c  = SUNRsqrt((sunrealtype)(conj(temp1)*temp1))/denomCS; 
      // signf = temp1/SUNRsqrt((sunrealtype)(conj(temp1)*temp1));
      // s  = -signf * conj(temp2)/denomCS; 
      // s = c * temp3;
    }
    else 
    {
      // // temp3 = temp2 / temp1;
      // temp3 = (conj(temp2)*temp2) / (conj(temp1)*temp1); //Amihere (use this and next two lines)
      // c  = ONE / SUNRsqrt(ONE + (sunrealtype)temp3);
      // s = -c * (temp1 * conj(temp2)) / ((sunrealtype)(conj(temp1)*temp1));

      //New - 30th March 2025
      temp3 = (conj(temp2)*temp2) / (conj(temp1)*temp1); //Amihere (use this and next three lines)
      absTemp1 = SUNRsqrt((sunrealtype)(conj(temp1)*temp1));
      c  = ONE / ((temp1/absTemp1)* SUNRsqrt(ONE + (sunrealtype)temp3));
      s = -c * (conj(temp2)/conj(temp1));

      // c     = ONE / SUNRsqrt(ONE + SUNSQR((sunrealtype)temp3));
      // // c     = ONE / SUNRsqrt(ONE + ((sunrealtype)(conj(temp3) * temp3)));
      // s     = -c * temp3;

      // denomCS = SUNRsqrt((sunrealtype)((conj(temp1)*temp1)+(conj(temp2)*temp2))); //Amihere (use this and next two lines)
      // c = conj(temp1)/denomCS;
      // s = -conj(temp2)/denomCS;

      // signf = temp1/SUNRsqrt((sunrealtype)(conj(temp1)*temp1));
      // s  = -signf * conj(temp2)/denomCS; 
      // s = c * temp3;
    }
    // // Amihere: debugging
    //   sumCS = (c * c) + (s * s);
    //   printf ("sum of c and s: %f + i %f\n", creal(sumCS), cimag(sumCS));
    //   printf ("rotation c: %f + i %f\n", creal(c), cimag(c));
    //   printf ("rotation s: %f + i %f\n", creal(s), cimag(s));
    //   printf ("rotation q[i]: %f + i %f\n", creal(q[1]), cimag(q[1]));
    //   printf ("rotation s: %f + i %f\n", creal(s), cimag(s));
    //   RotationColumn1 = c * -s;
    //   RotationColumn2 = conj(s) * conj(c); //((c * -s) + (conj(s) * c));
    //   RotationColumn = RotationColumn1 + RotationColumn2;
    //   printf("Rotation value: %f + i %f\n", creal(RotationColumn), cimag(RotationColumn));

    q_ptr        = 2 * n_minus_1;
    q[q_ptr]     = c;
    q[q_ptr + 1] = s;
    if ((h[n_minus_1][n_minus_1] = c * temp1 - s * temp2) == ZERO) { code = n; }
  }
  // // Amihere: debugging
  //     sumCS = (c * c) + (s * s);
  //     printf ("sum of c and s: %f + i %f\n", creal(sumCS), cimag(sumCS));
  //     printf ("rotation c: %f + i %f\n", creal(c), cimag(c));
  //     printf ("rotation s: %f + i %f\n", creal(s), cimag(s));
  //     printf ("rotation q[i]: %f + i %f\n", creal(q[1]), cimag(q[1]));
  //     printf ("rotation s: %f + i %f\n", creal(s), cimag(s));
  //     RotationColumn1 = c * -s;
  //     RotationColumn2 = conj(s) * conj(c); //((c * -s) + (conj(s) * c));
  //     RotationColumn = RotationColumn1 + RotationColumn2;
  //     printf("Rotation value: %f + i %f\n", creal(RotationColumn), cimag(RotationColumn));

  return (code);
}

/*
 * -----------------------------------------------------------------
 * Function : SUNQRsol
 * -----------------------------------------------------------------
 * This implementation of SUNQRsol is a slight modification of a
 * previous routine (called qrsol) written by Milo Dorr.
 * -----------------------------------------------------------------
 */

// int SUNQRsolComplex(int n, sunrealtype** h, sunrealtype* q, sunrealtype* b)
int SUNQRsolComplex(int n, suncomplextype** h, suncomplextype* q, suncomplextype* b)
{
  // sunrealtype c, s, temp1, temp2;
  suncomplextype c, s, temp1, temp2;
  int i, k, q_ptr, code = 0;

  /* Compute Q*b */

  for (k = 0; k < n; k++)
  {
    q_ptr    = 2 * k;
    c        = q[q_ptr];
    s        = q[q_ptr + 1];
    temp1    = b[k];
    temp2    = b[k + 1];
    b[k]     = c * temp1 - s * temp2;
    b[k + 1] = conj(s) * temp1 + conj(c) * temp2;
  }

  /* Solve  R*x = Q*b */

  for (k = n - 1; k >= 0; k--)
  {
    if (h[k][k] == ZERO)
    {
      code = k + 1;
      break;
    }
    b[k] /= h[k][k];
    for (i = 0; i < k; i++) { b[i] -= b[k] * h[i][k]; }
  }

  return (code);
}

/*
 * -----------------------------------------------------------------
 * Function : SUNQRAdd_MGS
 * -----------------------------------------------------------------
 * Implementation of QRAdd to be called in Anderson Acceleration
 * -----------------------------------------------------------------
 */

SUNErrCode SUNQRAdd_MGSComplex(N_Vector* Q, suncomplextype* R, N_Vector df, int m,
                        int mMax, void* QRdata)
{
  // SUNFunctionBegin(Q[0]->sunctx);

  sunindextype j;
  SUNQRData qrdata = (SUNQRData)QRdata;

  N_VScale_SComplex(ONE, df, qrdata->vtemp);
  // SUNCheckLastErr();
  for (j = 0; j < m; j++)
  {
    R[m * mMax + j] = N_VDotProd_SComplex(Q[j], qrdata->vtemp);
    // SUNCheckLastErr();
    N_VLinearSum_SComplex(ONE, qrdata->vtemp, -R[m * mMax + j], Q[j], qrdata->vtemp);
    // SUNCheckLastErr();
  }
  R[m * mMax + m] = (sunrealtype)N_VDotProd_SComplex(qrdata->vtemp, qrdata->vtemp);
  // SUNCheckLastErr();
  R[m * mMax + m] = SUNRsqrt((sunrealtype)R[m * mMax + m]);
  N_VScale_SComplex((1 / R[m * mMax + m]), qrdata->vtemp, Q[m]);
  // SUNCheckLastErr();

  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * Function : SUNQRAdd_ICWY
 * -----------------------------------------------------------------
 * Low synchronous implementation of QRAdd to be called in
 * Anderson Acceleration.
 * -----------------------------------------------------------------
 */

SUNErrCode SUNQRAdd_ICWYComplex(N_Vector* Q, suncomplextype* R, N_Vector df, int m,
                         int mMax, void* QRdata)
{
  // SUNFunctionBegin(Q[0]->sunctx);
  sunindextype j, k;
  SUNQRData qrdata = (SUNQRData)QRdata;

  N_VScale_SComplex(ONE, df, qrdata->vtemp);
  // SUNCheckLastErr(); /* stores d_fi in temp */

  if (m > 0)
  {
    /* T(1:k-1,k-1)^T = Q(:,1:k-1)^T * Q(:,k-1) */
    // SUNCheckCall(
    //   N_VDotProdMulti(m, Q[m - 1], Q, qrdata->temp_array + (m - 1) * mMax));
    N_VDotProdMulti_SComplex(m, Q[m - 1], Q, qrdata->temp_array + (m - 1) * mMax);

    /* T(k-1,k-1) = 1.0 */
    qrdata->temp_array[(m - 1) * mMax + (m - 1)] = ONE;

    /* R(1:k-1,k) = Q_k-1^T * df */
    // SUNCheckCall(N_VDotProdMulti(m, qrdata->vtemp, Q, R + m * mMax));
    N_VDotProdMulti_SComplex(m, qrdata->vtemp, Q, R + m * mMax);

    /* Solve T^T * R(1:k-1,k) = R(1:k-1,k) */
    for (k = 0; k < m; k++)
    {
      /* Skip setting the diagonal element because it doesn't change */
      for (j = k + 1; j < m; j++)
      {
        R[m * mMax + j] -= R[m * mMax + k] * qrdata->temp_array[j * mMax + k];
      }
    }
    /* end */

    /* Q(:,k-1) = df - Q_k-1 R(1:k-1,k) */
    // SUNCheckCall(N_VLinearCombination(m, R + m * mMax, Q, qrdata->vtemp2));
    N_VLinearCombination_SComplex(m, R + m * mMax, Q, qrdata->vtemp2);
    N_VLinearSum_SComplex(ONE, qrdata->vtemp, -ONE, qrdata->vtemp2, qrdata->vtemp);
    // SUNCheckLastErr();
  }

  /* R(k,k) = \| df \| */
  R[m * mMax + m] = (sunrealtype)N_VDotProd_SComplex(qrdata->vtemp, qrdata->vtemp);
  // SUNCheckLastErr();
  R[m * mMax + m] = SUNRsqrt((sunrealtype)R[m * mMax + m]);
  /* Q(:,k) = df / \| df \| */
  N_VScale_SComplex((1 / R[m * mMax + m]), qrdata->vtemp, Q[m]);
  // SUNCheckLastErr();

  return SUN_SUCCESS;
}

// /*
//  * -----------------------------------------------------------------
//  * Function : SUNQRAdd_ICWY_SB
//  * -----------------------------------------------------------------
//  * Low synchronous implementation of QRAdd to be called in
//  * Anderson Acceleration which utilizes a single buffer reduction.
//  * -----------------------------------------------------------------
//  */

// SUNErrCode SUNQRAdd_ICWY_SB(N_Vector* Q, sunrealtype* R, N_Vector df, int m,
//                             int mMax, void* QRdata)
// {
//   SUNFunctionBegin(Q[0]->sunctx);
//   sunindextype j, k;
//   SUNQRData qrdata = (SUNQRData)QRdata;

//   N_VScale(ONE, df, qrdata->vtemp);
//   SUNCheckLastErr(); /* stores d_fi in temp */

//   if (m > 0)
//   {
//     /* T(1:k-1,k-1)^T = Q(:,1:k-1)^T * Q(:,k-1) */
//     SUNCheckCall(N_VDotProdMultiLocal(m, Q[m - 1], Q,
//                                       qrdata->temp_array + (m - 1) * mMax));

//     /* R(1:k-1,k) = Q_k-1^T * df */
//     /* Put R values at end of temp_array */
//     SUNCheckCall(N_VDotProdMultiLocal(m, qrdata->vtemp, Q,
//                                       qrdata->temp_array + (m - 1) * mMax + m));

//     SUNCheckCall(N_VDotProdMultiAllReduce(m + m, qrdata->vtemp,
//                                           qrdata->temp_array + (m - 1) * mMax));

//     /* Move the last values from temp array into R */
//     for (k = 0; k < m; k++)
//     {
//       R[m * mMax + k] = qrdata->temp_array[(m - 1) * mMax + m + k];
//     }

//     /* T(k-1,k-1) = 1.0 */
//     qrdata->temp_array[(m - 1) * mMax + (m - 1)] = ONE;

//     /* Solve T^T * R(1:k-1,k) = R(1:k-1,k) */
//     for (k = 0; k < m; k++)
//     {
//       /* Skip setting the diagonal element because it doesn't change */
//       for (j = k + 1; j < m; j++)
//       {
//         R[m * mMax + j] -= R[m * mMax + k] * qrdata->temp_array[j * mMax + k];
//       }
//     }
//     /* end */

//     /* Q(:,k-1) = df - Q_k-1 R(1:k-1,k) */
//     SUNCheckCall(N_VLinearCombination(m, R + m * mMax, Q, qrdata->vtemp2));
//     N_VLinearSum(ONE, qrdata->vtemp, -ONE, qrdata->vtemp2, qrdata->vtemp);
//     SUNCheckLastErr();
//   }

//   /* R(k,k) = \| df \| */
//   R[m * mMax + m] = N_VDotProd(qrdata->vtemp, qrdata->vtemp);
//   SUNCheckLastErr();
//   R[m * mMax + m] = SUNRsqrt(R[m * mMax + m]);
//   /* Q(:,k) = df / \| df \| */
//   N_VScale((1 / R[m * mMax + m]), qrdata->vtemp, Q[m]);
//   SUNCheckLastErr();

//   return SUN_SUCCESS;
// }

/*
 * -----------------------------------------------------------------
 * Function : SUNQRAdd_CGS2
 * -----------------------------------------------------------------
 * Low synchronous Implementation of QRAdd to be called in
 * Anderson Acceleration.
 * -----------------------------------------------------------------
 */

SUNErrCode SUNQRAdd_CGS2Complex(N_Vector* Q, suncomplextype* R, N_Vector df, int m,
                         int mMax, void* QRdata)
{
  // SUNFunctionBegin(Q[0]->sunctx);
  sunindextype j;
  SUNQRData qrdata = (SUNQRData)QRdata;

  N_VScale_SComplex(ONE, df, qrdata->vtemp);
  // SUNCheckLastErr(); /* temp = df */

  if (m > 0)
  {
    /* s_k = Q_k-1^T df_aa -- update with sdata as a sunrealtype* array */
    // SUNCheckCall(N_VDotProdMulti(m, qrdata->vtemp, Q, R + m * mMax));
    N_VDotProdMulti_SComplex(m, qrdata->vtemp, Q, R + m * mMax);

    /* y = df - Q_k-1 s_k */
    // SUNCheckCall(N_VLinearCombination(m, R + m * mMax, Q, qrdata->vtemp2));
    N_VLinearCombination_SComplex(m, R + m * mMax, Q, qrdata->vtemp2);
    N_VLinearSum_SComplex(ONE, qrdata->vtemp, -ONE, qrdata->vtemp2, qrdata->vtemp2);
    // SUNCheckLastErr();

    /* z_k = Q_k-1^T y */
    // SUNCheckCall(N_VDotProdMulti(m, qrdata->vtemp2, Q, qrdata->temp_array));
    N_VDotProdMulti_SComplex(m, qrdata->vtemp2, Q, qrdata->temp_array);

    /* df = y - Q_k-1 z_k */
    // SUNCheckCall(N_VLinearCombination(m, qrdata->temp_array, Q, Q[m]));
    N_VLinearCombination_SComplex(m, qrdata->temp_array, Q, Q[m]);
    N_VLinearSum_SComplex(ONE, qrdata->vtemp2, -ONE, Q[m], qrdata->vtemp);
    // SUNCheckLastErr();

    /* R(1:k-1,k) = s_k + z_k */
    for (j = 0; j < m; j++)
    {
      R[m * mMax + j] = R[m * mMax + j] + qrdata->temp_array[j];
    }
  }

  /* R(k,k) = \| df \| */
  R[m * mMax + m] = (sunrealtype)N_VDotProd_SComplex(qrdata->vtemp, qrdata->vtemp);
  // SUNCheckLastErr();
  R[m * mMax + m] = SUNRsqrt((sunrealtype)R[m * mMax + m]);
  /* Q(:,k) = df / R(k,k) */
  N_VScale_SComplex((1 / R[m * mMax + m]), qrdata->vtemp, Q[m]);
  // SUNCheckLastErr();

  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * Function : SUNQRAdd_DCGS2
 * -----------------------------------------------------------------
 * Low synchronous Implementation of QRAdd to be called in
 * Anderson Acceleration.
 * -----------------------------------------------------------------
 */

SUNErrCode SUNQRAdd_DCGS2Complex(N_Vector* Q, suncomplextype* R, N_Vector df, int m,
                          int mMax, void* QRdata)
{
  // SUNFunctionBegin(Q[0]->sunctx);
  sunindextype j;
  SUNQRData qrdata = (SUNQRData)QRdata;

  N_VScale_SComplex(ONE, df, qrdata->vtemp);
  // SUNCheckLastErr(); /* temp = df */

  if (m > 0)
  {
    /* R(1:k-1,k) = Q_k-1^T df_aa */
    // SUNCheckCall(N_VDotProdMulti(m, qrdata->vtemp, Q, R + m * mMax));
    N_VDotProdMulti_SComplex(m, qrdata->vtemp, Q, R + m * mMax);
    /* Delayed reorthogonalization */
    if (m > 1)
    {
      /* s = Q_k-2^T Q(:,k-1) */
      // SUNCheckCall(N_VDotProdMulti(m - 1, Q[m - 1], Q, qrdata->temp_array));
      N_VDotProdMulti_SComplex(m - 1, Q[m - 1], Q, qrdata->temp_array);

      /* Q(:,k-1) = Q(:,k-1) - Q_k-2 s */
      // SUNCheckCall(
      //   N_VLinearCombination(m - 1, qrdata->temp_array, Q, qrdata->vtemp2));
      N_VLinearCombination_SComplex(m - 1, qrdata->temp_array, Q, qrdata->vtemp2);
      N_VLinearSum_SComplex(ONE, Q[m - 1], -ONE, qrdata->vtemp2, Q[m - 1]);
      // SUNCheckLastErr();

      /* R(1:k-2,k-1) = R(1:k-2,k-1) + s */
      for (j = 0; j < m - 1; j++)
      {
        R[(m - 1) * mMax + j] = R[(m - 1) * mMax + j] + qrdata->temp_array[j];
      }
    }

    /* df = df - Q(:,k-1) R(1:k-1,k) */
    // SUNCheckCall(N_VLinearCombination(m, R + m * mMax, Q, qrdata->vtemp2));
    N_VLinearCombination_SComplex(m, R + m * mMax, Q, qrdata->vtemp2);
    N_VLinearSum_SComplex(ONE, qrdata->vtemp, -ONE, qrdata->vtemp2, qrdata->vtemp);
    // SUNCheckLastErr();
  }

  /* R(k,k) = \| df \| */
  R[m * mMax + m] = (sunrealtype)N_VDotProd_SComplex(qrdata->vtemp, qrdata->vtemp);
  // SUNCheckLastErr();
  R[m * mMax + m] = SUNRsqrt((sunrealtype)R[m * mMax + m]);
  /* Q(:,k) = df / R(k,k) */
  N_VScale_SComplex((1 / R[m * mMax + m]), qrdata->vtemp, Q[m]);
  // SUNCheckLastErr();

  return SUN_SUCCESS;
}

// /*
//  * -----------------------------------------------------------------
//  * Function : SUNQRAdd_DCGS2_SB
//  * -----------------------------------------------------------------
//  * Low synchronous Implementation of QRAdd to be called in
//  * Anderson Acceleration which utilizes a single buffer reduction.
//  * -----------------------------------------------------------------
//  */

// SUNErrCode SUNQRAdd_DCGS2_SB(N_Vector* Q, sunrealtype* R, N_Vector df, int m,
//                              int mMax, void* QRdata)
// {
//   SUNFunctionBegin(Q[0]->sunctx);
//   sunindextype j;
//   SUNQRData qrdata = (SUNQRData)QRdata;

//   N_VScale(ONE, df, qrdata->vtemp);
//   SUNCheckLastErr(); /* temp = df */

//   if (m > 0)
//   {
//     if (m == 1)
//     {
//       /* R(1:k-1,k) = Q_k-1^T df_aa */
//       SUNCheckCall(N_VDotProdMulti(m, qrdata->vtemp, Q, R + m * mMax));
//     }
//     /* Delayed reorthogonalization */
//     else if (m > 1)
//     {
//       /* R(1:k-1,k) = Q_k-1^T df_aa */
//       /* Put R values at beginning of temp array */
//       SUNCheckCall(N_VDotProdMultiLocal(m, qrdata->vtemp, Q, qrdata->temp_array));

//       /* s = Q_k-2^T Q(:,k-1) */
//       SUNCheckCall(
//         N_VDotProdMultiLocal(m - 1, Q[m - 1], Q, qrdata->temp_array + m));
//       SUNCheckCall(
//         N_VDotProdMultiAllReduce(m + m - 1, qrdata->vtemp, qrdata->temp_array));

//       /* Move R values to R */
//       for (j = 0; j < m; j++) { R[m * mMax + j] = qrdata->temp_array[j]; }

//       /* Q(:,k-1) = Q(:,k-1) - Q_k-2 s */
//       SUNCheckCall(
//         N_VLinearCombination(m - 1, qrdata->temp_array + m, Q, qrdata->vtemp2));
//       N_VLinearSum(ONE, Q[m - 1], -ONE, qrdata->vtemp2, Q[m - 1]);
//       SUNCheckLastErr();

//       /* R(1:k-2,k-1) = R(1:k-2,k-1) + s */
//       for (j = 0; j < m - 1; j++)
//       {
//         R[(m - 1) * mMax + j] = R[(m - 1) * mMax + j] + qrdata->temp_array[m + j];
//       }
//     }

//     /* df = df - Q(:,k-1) R(1:k-1,k) */
//     SUNCheckCall(N_VLinearCombination(m, R + m * mMax, Q, qrdata->vtemp2));
//     N_VLinearSum(ONE, qrdata->vtemp, -ONE, qrdata->vtemp2, qrdata->vtemp);
//     SUNCheckLastErr();
//   }

//   /* R(k,k) = \| df \| */
//   R[m * mMax + m] = N_VDotProd(qrdata->vtemp, qrdata->vtemp);
//   SUNCheckLastErr();
//   R[m * mMax + m] = SUNRsqrt(R[m * mMax + m]);
//   /* Q(:,k) = df / R(k,k) */
//   N_VScale((1 / R[m * mMax + m]), qrdata->vtemp, Q[m]);
//   SUNCheckLastErr();

//   return SUN_SUCCESS;
// }
