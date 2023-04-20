/*
********************************************************************************
* Curve grid metric calculation using MacCormack scheme                        *
********************************************************************************
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "common.h"
#include "params.h"
#include "macdrp.h"

//void cal_metric(struct Coord *C, float *M) {
void cal_metric(float *C, float *M) {

  int i, j, k, pos_c, pos_m, slice, segment;
  float vp, vs, rho, lam, mu;
//  int n;

  float x_xi, x_et, x_zt;
  float y_xi, y_et, y_zt;
  float z_xi, z_et, z_zt;
  float jac;
  float vec1[3], vec2[3], vec3[3], vecg[3];


  for (i = ni1; i < ni2; i++)
    for (j = nj1; j < nj2; j++)
      for (k = nk1; k < nk2; k++){

        pos_c = (i * ny * lnz + j * lnz + k) * CSIZE;
        pos_m = (i * ny * lnz + j * lnz + k) * MSIZE;
        slice = ny * lnz * CSIZE; segment = lnz * CSIZE;

        x_xi = 0.0; x_et = 0.0; x_zt = 0.0;
        y_xi = 0.0; y_et = 0.0; y_zt = 0.0;
        z_xi = 0.0; z_et = 0.0; z_zt = 0.0;

        // MacCormack scheme 
//        for (n = 0; n < 5; n++){
//            x_xi += ((C->x[i+n-1][j][k] * coeF[n] + C->x[i+n-3][j][k] * coeB[n]) / (2.0*DH));
//            y_xi += ((C->y[i+n-1][j][k] * coeF[n] + C->y[i+n-3][j][k] * coeB[n]) / (2.0*DH));
//            z_xi += ((C->z[i+n-1][j][k] * coeF[n] + C->z[i+n-3][j][k] * coeB[n]) / (2.0*DH));
//
//            x_et += ((C->x[i][j+n-1][k] * coeF[n] + C->x[i][j+n-3][k] * coeB[n]) / (2.0*DH));
//            y_et += ((C->y[i][j+n-1][k] * coeF[n] + C->y[i][j+n-3][k] * coeB[n]) / (2.0*DH));
//            z_et += ((C->z[i][j+n-1][k] * coeF[n] + C->z[i][j+n-3][k] * coeB[n]) / (2.0*DH));
//
//            x_zt += ((C->x[i][j][k+n-1] * coeF[n] + C->x[i][j][k+n-3] * coeB[n]) / (2.0*DH));
//            y_zt += ((C->y[i][j][k+n-1] * coeF[n] + C->y[i][j][k+n-3] * coeB[n]) / (2.0*DH));
//            z_zt += ((C->z[i][j][k+n-1] * coeF[n] + C->z[i][j][k+n-3] * coeB[n]) / (2.0*DH));
//        }

        x_xi = 0.5f * ( L(C, (pos_c + 0), slice, FWD) + L(C, (pos_c + 0), slice, BWD) );
        y_xi = 0.5f * ( L(C, (pos_c + 1), slice, FWD) + L(C, (pos_c + 1), slice, BWD) );
        z_xi = 0.5f * ( L(C, (pos_c + 2), slice, FWD) + L(C, (pos_c + 2), slice, BWD) );

        x_et = 0.5f * ( L(C, (pos_c + 0), segment, FWD) + L(C, (pos_c + 0), segment, BWD) );
        y_et = 0.5f * ( L(C, (pos_c + 1), segment, FWD) + L(C, (pos_c + 1), segment, BWD) );
        z_et = 0.5f * ( L(C, (pos_c + 2), segment, FWD) + L(C, (pos_c + 2), segment, BWD) );

        x_zt = 0.5f * ( L(C, (pos_c + 0), CSIZE, FWD) + L(C, (pos_c + 0), CSIZE, BWD) );
        y_zt = 0.5f * ( L(C, (pos_c + 1), CSIZE, FWD) + L(C, (pos_c + 1), CSIZE, BWD) );
        z_zt = 0.5f * ( L(C, (pos_c + 2), CSIZE, FWD) + L(C, (pos_c + 2), CSIZE, BWD) );

        vec1[0] = x_xi; vec1[1] = y_xi; vec1[2] = z_xi;
        vec2[0] = x_et; vec2[1] = y_et; vec2[2] = z_et;
        vec3[0] = x_zt; vec3[1] = y_zt; vec3[2] = z_zt;

        cross_product(vec1, vec2, vecg);
        jac = dot_product(vecg, vec3);
        M[pos_m + 9]  = jac;

        cross_product(vec2, vec3, vecg);
        M[pos_m + 0] = vecg[0] / jac;
        M[pos_m + 1] = vecg[1] / jac;
        M[pos_m + 2] = vecg[2] / jac;

        cross_product(vec3, vec1, vecg);
        M[pos_m + 3] = vecg[0] / jac;
        M[pos_m + 4] = vecg[1] / jac;
        M[pos_m + 5] = vecg[2] / jac;

        cross_product(vec1, vec2, vecg);
        M[pos_m + 6] = vecg[0] / jac;
        M[pos_m + 7] = vecg[1] / jac;
        M[pos_m + 8] = vecg[2] / jac;
        
//        if( tmpz > -3e3){
//          M[pos_m + 10] = M[pos_m + 10] / 1.5f;
//          M[pos_m + 11] = M[pos_m + 11] / 1.5f ;
//          M[pos_m + 12] = M[pos_m + 12] / 1.5f;
//        }
        //if (masternode) {
        	//printf("--->%s,i:%d,j:%d,k:%d,rho:%f\n",__FUNCTION__,i,j,k,M[pos_m + 12]);
        //}

      }

  return;
}
