#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "common.h"
#include "macdrp.h"
#include "params.h"

/*
 **************************************************************************
 *  Free surface condtion for Velocity                                    *
 *  Transform parallel derivatives to normal, i.e. Vx2Vz                  *
 **************************************************************************
 */
int coef_surface(float *M, float *matVx2Vz, float *matVy2Vz, float *pml){
  int i, j, k, pos_m;
  int ij;
  k = nk2-1;
  float e11, e12, e13, e21, e22, e23, e31, e32, e33;
  float lam2mu, lam, mu;
  float A[3][3], B[3][3], C[3][3];
  float AB[3][3], AC[3][3];
  float rb1, rb2, rb3;
  float *p0, *p1, *p2;

  for (i = ni1; i < ni2; i++)
    for (j = nj1; j < nj2; j++){

#ifdef usePML
  p0 = pml;
  p1 = pml + ni * 3;
  p2 = pml + ni * 3 + nj * 3;
  rb1 = 1.0f / p0[ni + (i-ni1)];
  rb2 = 1.0f / p1[nj + (j-nj1)];
  rb3 = 1.0f / p2[nk + (k-nk1)];
#else
  rb1 = 1.0f;
  rb2 = 1.0f;
  rb3 = 1.0f;
#endif
 //     pos   = (i * ny * lnz + j * lnz + k) * WSIZE;
      pos_m = (i * ny * lnz + j * lnz + k) * MSIZE;
      e11 = M[pos_m + 0];
      e12 = M[pos_m + 1];
      e13 = M[pos_m + 2];
      e21 = M[pos_m + 3];
      e22 = M[pos_m + 4];
      e23 = M[pos_m + 5];
      e31 = M[pos_m + 6];
      e32 = M[pos_m + 7];
      e33 = M[pos_m + 8];

      lam = M[pos_m + 10];
      mu  = M[pos_m + 11];
      lam2mu = lam + 2.0f * mu;

      A[0][0] = rb3*(lam2mu*e31*e31 + mu*(e32*e32+e33*e33));
      A[0][1] = rb3*(lam*e31*e32 + mu*e32*e31);
      A[0][2] = rb3*(lam*e31*e33 + mu*e33*e31);
      A[1][0] = rb3*(lam*e32*e31 + mu*e31*e32);
      A[1][1] = rb3*(lam2mu*e32*e32 + mu*(e31*e31+e33*e33));
      A[1][2] = rb3*(lam*e32*e33 + mu*e33*e32);
      A[2][0] = rb3*(lam*e33*e31 + mu*e31*e33);
      A[2][1] = rb3*(lam*e33*e32 + mu*e32*e33);
      A[2][2] = rb3*(lam2mu*e33*e33 + mu*(e31*e31+e32*e32));
      invert3x3(A);

      B[0][0] = rb1*(-lam2mu*e31*e11 - mu*(e32*e12+e33*e13));
      B[0][1] = rb1*(-lam*e31*e12 - mu*e32*e11);
      B[0][2] = rb1*(-lam*e31*e13 - mu*e33*e11);
      B[1][0] = rb1*(-lam*e32*e11 - mu*e31*e12);
      B[1][1] = rb1*(-lam2mu*e32*e12 - mu*(e31*e11+e33*e13));
      B[1][2] = rb1*(-lam*e32*e13 - mu*e33*e12);
      B[2][0] = rb1*(-lam*e33*e11 - mu*e31*e13);
      B[2][1] = rb1*(-lam*e33*e12 - mu*e32*e13);
      B[2][2] = rb1*(-lam2mu*e33*e13 - mu*(e31*e11+e32*e12));

      C[0][0] = rb2*(-lam2mu*e31*e21 - mu*(e32*e22+e33*e23));
      C[0][1] = rb2*(-lam*e31*e22 - mu*e32*e21);
      C[0][2] = rb2*(-lam*e31*e23 - mu*e33*e21);
      C[1][0] = rb2*(-lam*e32*e21 - mu*e31*e22);
      C[1][1] = rb2*(-lam2mu*e32*e22 - mu*(e31*e21+e33*e23));
      C[1][2] = rb2*(-lam*e32*e23 - mu*e33*e22);
      C[2][0] = rb2*(-lam*e33*e21 - mu*e31*e23);
      C[2][1] = rb2*(-lam*e33*e22 - mu*e32*e23);
      C[2][2] = rb2*(-lam2mu*e33*e23 - mu*(e31*e21+e32*e22));

      // B[0][0] = -rb1*(-lam2mu*e31*e11 - mu*(e32*e12+e33*e13));
      // B[0][1] = -rb1*(-lam*e31*e12 - mu*e32*e11);
      // B[0][2] = -rb1*(-lam*e31*e13 - mu*e33*e11);
      // B[1][0] = -rb1*(-lam*e32*e11 - mu*e31*e12);
      // B[1][1] = -rb1*(-lam2mu*e32*e12 - mu*(e31*e11+e33*e13));
      // B[1][2] = -rb1*(-lam*e32*e13 - mu*e33*e12);
      // B[2][0] = -rb1*(-lam*e33*e11 - mu*e31*e13);
      // B[2][1] = -rb1*(-lam*e33*e12 - mu*e32*e13);
      // B[2][2] = -rb1*(-lam2mu*e33*e13 - mu*(e31*e11+e32*e12));

      // C[0][0] = -rb2*(-lam2mu*e31*e21 - mu*(e32*e22+e33*e23));
      // C[0][1] = -rb2*(-lam*e31*e22 - mu*e32*e21);
      // C[0][2] = -rb2*(-lam*e31*e23 - mu*e33*e21);
      // C[1][0] = -rb2*(-lam*e32*e21 - mu*e31*e22);
      // C[1][1] = -rb2*(-lam2mu*e32*e22 - mu*(e31*e21+e33*e23));
      // C[1][2] = -rb2*(-lam*e32*e23 - mu*e33*e22);
      // C[2][0] = -rb2*(-lam*e33*e21 - mu*e31*e23);
      // C[2][1] = -rb2*(-lam*e33*e22 - mu*e32*e23);
      // C[2][2] = -rb2*(-lam2mu*e33*e23 - mu*(e31*e21+e32*e22));

      matmul3x3(A, B, AB);
      matmul3x3(A, C, AC);

      ij = (i*ny+j)*9;

      int m, n;
      for(m = 0; m < 3; m++)
        for(n = 0; n < 3; n++){
          matVx2Vz[ij + m*3 + n] = AB[m][n];
          matVy2Vz[ij + m*3 + n] = AC[m][n];
          //matF2Vz [ij + m*3 + n] =  A[m][n];
        }
  }

  return 0 ;
}
