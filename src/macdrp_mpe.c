#ifndef ARCH_SW
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "common.h"
#include "macdrp.h"
#include "params.h"

int RK_init(float *W, float *mW){
  int blocksize = nx * ny * lnz * WSIZE * sizeof(float);
  int pos = 0; //int len = (nx * ny * nz) * WSIZE;
  memcpy(&mW[pos + 0], &W[pos + 0], blocksize);

  return 0 ;
}

int RK_init2(float *W, float *mW, int nnx, int nny, int nnz){
  int i = 0, j = 0, k = 0, iw = 0;
  for(i = 0 ; i < nnx ; i ++) {
    for(j = 0 ; j < nny ; j ++) {
      int pos = (i * nny * nnz + j * nnz) * WSIZE;
      memcpy(mW + pos, W + pos, sizeof(float) * WSIZE * nnz);
    }
  }

  return 0 ;
}


int Calculate(float *W, float *mW, float *hW, float *tW, float *M, float *matVx2Vz, float *matVy2Vz, int Flag[3], float rka, float rkb,int flag, int idx0, float *pml, struct aux *Aux, int isx1, int isx2, int isy1, int isy2, int isz1, int isz2) {
  static int cnt = -1;
  cnt ++;
  //fprintf(stderr, "cnt = %d\n", cnt);
  float rrho;
  //rrho = 1.0/rho;
  float lam, mu, lam2mu;
  float d1;//, d2, d3;
  float ad1;//, ad2, ad3;

  float b1;//, b2, b3;

  int FlagX = Flag[0];
  int FlagY = Flag[1];
  int FlagZ = Flag[2];

  int n, l, i, ii, j, jj, k, kk, iii, iiii;
  int pos, pos_m, pos_d, pos_a, slice, segment;
  float rrhojac ;
  //float vecTx[7], vecTy[7], vecTz[7];
  float vecT11[7],vecT12[7],vecT13[7];
  float vecT21[7],vecT22[7],vecT23[7];
  float vecT31[7],vecT32[7],vecT33[7];
  //float DxTx, DyTy, DzTz;
  float DxT11, DxT12, DxT13;
  float DyT21, DyT22, DyT23;
  float DzT31, DzT32, DzT33;

  float DxTxx, DxTyy, DxTzz, DxTxy, DxTxz, DxTyz;
  float DyTxx, DyTyy, DyTzz, DyTxy, DyTxz, DyTyz;
  float DzTxx, DzTyy, DzTzz, DzTxy, DzTxz, DzTyz;
  float DxVx, DxVy, DxVz;
  float DyVx, DyVy, DyVz;
  float DzVx, DzVy, DzVz;
  float DzVx1, DzVy1, DzVz1;

  float xix, xiy, xiz;
  float etx, ety, etz;
  float ztx, zty, ztz;

  //Grid3D matVx2Vz, matVy2Vz, matF2Vz;

  //matVx2Vz = Alloc3D(nx*ny, 3, 3);
  //matVy2Vz = Alloc3D(nx*ny, 3, 3);
  //matF2Vz  = Alloc3D(nx*ny, 3, 3);
  //coef_fdxy2fdz(M,matVx2Vz,matVy2Vz,matF2Vz);

  float *Aux_Wx, *Aux_hWx, *Aux_Wy, *Aux_hWy, *Aux_Wz, *Aux_hWz;
  float *oW = W, *omW = mW, *ohW = hW, *otW = tW;
  rka = rka * DT; rkb = rkb * DT;

  for (i = ni1; i < ni2; i++) {
    for (j = nj1; j < nj2; j++) {

      int ij = (i*ny+j)*9;

      for (k = nk1; k < nk2; k++) {

        //====================== calcu_deriv ================
        pos   = (i * ny * lnz + j * lnz + k) * WSIZE;
        pos_m = (i * ny * lnz + j * lnz + k) * MSIZE;
        slice = ny * lnz * WSIZE; segment = lnz * WSIZE;
        if(isnan(W[pos + 0])) {
          fprintf(stderr, "0.1 pos: %d %d %d\n", i, j , k);
          exit(1);
        }

        DxVx  = L(W,(pos + 0), slice, FlagX); DyVx  = L(W,(pos + 0), segment, FlagY); DzVx  = L(W, (pos + 0), WSIZE, FlagZ);
        DxVy  = L(W,(pos + 1), slice, FlagX); DyVy  = L(W,(pos + 1), segment, FlagY); DzVy  = L(W, (pos + 1), WSIZE, FlagZ);
        DxVz  = L(W,(pos + 2), slice, FlagX); DyVz  = L(W,(pos + 2), segment, FlagY); DzVz  = L(W, (pos + 2), WSIZE, FlagZ);
        DxTxx = L(W,(pos + 3), slice, FlagX); DyTxx = L(W,(pos + 3), segment, FlagY); DzTxx = L(W, (pos + 3), WSIZE, FlagZ);
        DxTyy = L(W,(pos + 4), slice, FlagX); DyTyy = L(W,(pos + 4), segment, FlagY); DzTyy = L(W, (pos + 4), WSIZE, FlagZ);
        DxTzz = L(W,(pos + 5), slice, FlagX); DyTzz = L(W,(pos + 5), segment, FlagY); DzTzz = L(W, (pos + 5), WSIZE, FlagZ);
        DxTxy = L(W,(pos + 6), slice, FlagX); DyTxy = L(W,(pos + 6), segment, FlagY); DzTxy = L(W, (pos + 6), WSIZE, FlagZ);
        DxTxz = L(W,(pos + 7), slice, FlagX); DyTxz = L(W,(pos + 7), segment, FlagY); DzTxz = L(W, (pos + 7), WSIZE, FlagZ);
        DxTyz = L(W,(pos + 8), slice, FlagX); DyTyz = L(W,(pos + 8), segment, FlagY); DzTyz = L(W, (pos + 8), WSIZE, FlagZ);

        xix = M[pos_m + 0]; xiy = M[pos_m + 1]; xiz = M[pos_m + 2];
        etx = M[pos_m + 3]; ety = M[pos_m + 4]; etz = M[pos_m + 5];
        ztx = M[pos_m + 6]; zty = M[pos_m + 7]; ztz = M[pos_m + 8];

        lam = M[pos_m + 10];
        mu  = M[pos_m + 11];
        lam2mu = lam + 2.0f*mu;
        rrho = 1.0f/M[pos_m + 12];

        float *p0, *p1, *p2;
        float rb1, rb2, rb3;
#ifdef usePML
        p0 = pml;
        p1 = pml + ni * 3;
        p2 = pml + ni * 3 + nj * 3;
        rb1 = 1.0f/p0[ni + (i - ni1)];
        rb2 = 1.0f/p1[nj + (j - nj1)];
        rb3 = 1.0f/p2[nk + (k - nk1)];
#else
        rb1 = 1.0f;
        rb2 = 1.0f;
        rb3 = 1.0f;
#endif

        // Moment equation
        hW[pos + 0] = ( (DxTxx*xix + DxTxy*xiy + DxTxz*xiz) * rb1
                      + (DyTxx*etx + DyTxy*ety + DyTxz*etz) * rb2
                      + (DzTxx*ztx + DzTxy*zty + DzTxz*ztz) * rb3 ) * rrho;
        if(isnan(DxTxx) || isnan(DxTxy) || isnan(DxTxz) || isnan(DyTyy) || isnan(DyTzz) || isnan(DzTyz)) {
          fprintf(stderr, "%f %f %f %f %f %f\n", DxTxx, DxTyy, DxTzz, DxTxy, DxTxz, DxTyz);
          fprintf(stderr, "0.2 pos: %d %d %d\n", i, j , k);
          exit(1);
        }
        if(isnan(hW[pos + 0])) {
          fprintf(stderr, "0.3 pos: %d %d %d\n", i, j , k);
          exit(1);
        }

        hW[pos + 1] = ( (DxTxy*xix + DxTyy*xiy + DxTyz*xiz) * rb1
                      + (DyTxy*etx + DyTyy*ety + DyTyz*etz) * rb2
                      + (DzTxy*ztx + DzTyy*zty + DzTyz*ztz) * rb3 ) * rrho;

        hW[pos + 2] = ( (DxTxz*xix + DxTyz*xiy + DxTzz*xiz) * rb1
                      + (DyTxz*etx + DyTyz*ety + DyTzz*etz) * rb2
                      + (DzTxz*ztx + DzTyz*zty + DzTzz*ztz) * rb3 ) * rrho;

        // Hooke's law
        hW[pos + 3] = (lam2mu*DxVx*xix + lam*DxVy*xiy + lam*DxVz*xiz) * rb1
                    + (lam2mu*DyVx*etx + lam*DyVy*ety + lam*DyVz*etz) * rb2
                    + (lam2mu*DzVx*ztx + lam*DzVy*zty + lam*DzVz*ztz) * rb3;

        hW[pos + 4] = (lam*DxVx*xix + lam2mu*DxVy*xiy + lam*DxVz*xiz) * rb1
                    + (lam*DyVx*etx + lam2mu*DyVy*ety + lam*DyVz*etz) * rb2
                    + (lam*DzVx*ztx + lam2mu*DzVy*zty + lam*DzVz*ztz) * rb3;

        hW[pos + 5] = (lam*DxVx*xix + lam*DxVy*xiy + lam2mu*DxVz*xiz) * rb1
                    + (lam*DyVx*etx + lam*DyVy*ety + lam2mu*DyVz*etz) * rb2
                    + (lam*DzVx*ztx + lam*DzVy*zty + lam2mu*DzVz*ztz) * rb3;

        hW[pos + 6] = mu * ( (DxVx*xiy + DxVy*xix) * rb1
                           + (DyVx*ety + DyVy*etx) * rb2
                           + (DzVx*zty + DzVy*ztx) * rb3 );

        hW[pos + 7] = mu * ( (DxVx*xiz + DxVz*xix) * rb1
                           + (DyVx*etz + DyVz*etx) * rb2
                           + (DzVx*ztz + DzVz*ztx) * rb3 );

        hW[pos + 8] = mu * ( (DxVy*xiz + DxVz*xiy) * rb1
                           + (DyVy*etz + DyVz*ety) * rb2
                           + (DzVy*ztz + DzVz*zty) * rb3 );

        // ================end of Calcu_deriv==============
        // =================free surface====================
#ifdef FreeSurface
        if(freenode && (k >= nk2 -3) && (k <= nk2 - 1)) {          
          //  if not free surface then do nothing!
          //
          //  use Traction Image method to calculate fd of stress compoents,
          //  then assemble the right hand side to update velocities
          //
          n = nk2 - 1 - k;
          rrhojac = rrho/M[pos_m + 9];

          for (l = 0; l < 7; l++) {
            ii = i + l - 3; 
            pos   = (ii * ny * lnz + j * lnz + k) * WSIZE;
            pos_m = (ii * ny * lnz + j * lnz + k) * MSIZE;
            // Txx 3 Tyy 4 Tzz 5 Txy 6 Txz 7 Tyz 8
            vecT11[l] = M[pos_m + 9] * ( M[pos_m + 0] * W[pos + 3]
                                       + M[pos_m + 1] * W[pos + 6]
                                       + M[pos_m + 2] * W[pos + 7] );
            vecT12[l] = M[pos_m + 9] * ( M[pos_m + 0] * W[pos + 6]
                                       + M[pos_m + 1] * W[pos + 4]
                                       + M[pos_m + 2] * W[pos + 8] );
            vecT13[l] = M[pos_m + 9] * ( M[pos_m + 0] * W[pos + 7]
                                       + M[pos_m + 1] * W[pos + 8]
                                       + M[pos_m + 2] * W[pos + 5] );
            jj = j + l - 3; 
            pos   = (i * ny * lnz + jj * lnz + k) * WSIZE;
            pos_m = (i * ny * lnz + jj * lnz + k) * MSIZE;
            // Txx 3 Tyy 4 Tzz 5 Txy 6 Txz 7 Tyz 8
            vecT21[l] = M[pos_m + 9] * ( M[pos_m + 3] * W[pos + 3]
                                       + M[pos_m + 4] * W[pos + 6]
                                       + M[pos_m + 5] * W[pos + 7] );
            vecT22[l] = M[pos_m + 9] * ( M[pos_m + 3] * W[pos + 6]
                                       + M[pos_m + 4] * W[pos + 4]
                                       + M[pos_m + 5] * W[pos + 8] );
            vecT23[l] = M[pos_m + 9] * ( M[pos_m + 3] * W[pos + 7]
                                       + M[pos_m + 4] * W[pos + 8]
                                       + M[pos_m + 5] * W[pos + 5] );
            kk = k + l - 3; 
            pos   = (i * ny * lnz + j * lnz + kk) * WSIZE;
            pos_m = (i * ny * lnz + j * lnz + kk) * MSIZE;
            // Txx 3 Tyy 4 Tzz 5 Txy 6 Txz 7 Tyz 8
            vecT31[l] = M[pos_m + 9] * ( M[pos_m + 6] * W[pos + 3]
                                       + M[pos_m + 7] * W[pos + 6]
                                       + M[pos_m + 8] * W[pos + 7] );
            vecT32[l] = M[pos_m + 9] * ( M[pos_m + 6] * W[pos + 6]
                                       + M[pos_m + 7] * W[pos + 4]
                                       + M[pos_m + 8] * W[pos + 8] );
            vecT33[l] = M[pos_m + 9] * ( M[pos_m + 6] * W[pos + 7]
                                       + M[pos_m + 7] * W[pos + 8]
                                       + M[pos_m + 8] * W[pos + 5] );
          }

          vecT31[n+3] = 0.0f;               // FreeSurface set to zero
          vecT32[n+3] = 0.0f;               // FreeSurface set to zero
          vecT33[n+3] = 0.0f;               // FreeSurface set to zero
          for (l = 0; l < 3 - n; l++) {
            vecT31[6-l] = -vecT31[l+2*n];   // Traction Image
            vecT32[6-l] = -vecT32[l+2*n];   // Traction Image
            vecT33[6-l] = -vecT33[l+2*n];   // Traction Image
          }

          DxT11 = vec_L(vecT11,3,FlagX);
          DxT12 = vec_L(vecT12,3,FlagX);
          DxT13 = vec_L(vecT13,3,FlagX);
          DyT21 = vec_L(vecT21,3,FlagY);
          DyT22 = vec_L(vecT22,3,FlagY);
          DyT23 = vec_L(vecT23,3,FlagY);
          DzT31 = vec_L(vecT31,3,FlagZ);
          DzT32 = vec_L(vecT32,3,FlagZ);
          DzT33 = vec_L(vecT33,3,FlagZ);

          pos = (i * ny * lnz + j * lnz + k)  *  WSIZE; // reset the position of index !!!!!
          hW[pos + 0] = ( DxT11*rb1 + DyT21*rb2 + DzT31*rb3 ) * rrhojac;
        if(isnan(hW[pos + 0])) {
          fprintf(stderr, "0.4 pos: %d %d %d\n", i, j , k);
          exit(1);
        }

          hW[pos + 1] = ( DxT12*rb1 + DyT22*rb2 + DzT32*rb3 ) * rrhojac;
          hW[pos + 2] = ( DxT13*rb1 + DyT23*rb2 + DzT33*rb3 ) * rrhojac;

          // ============================================================================
          //  use Velocity Low scheme to calculate velocities fd with
          //  respect to eta and assemble the right hand side to update stresses
          // ============================================================================

          pos   = (i * ny * lnz + j * lnz + k) * WSIZE;
          pos_m = (i * ny * lnz + j * lnz + k) * MSIZE;

          xix = M[pos_m + 0];
          xiy = M[pos_m + 1];
          xiz = M[pos_m + 2];
          etx = M[pos_m + 3];
          ety = M[pos_m + 4];
          etz = M[pos_m + 5];
          ztx = M[pos_m + 6];
          zty = M[pos_m + 7];
          ztz = M[pos_m + 8];

          // Lx
          DxVx  = L(W,(pos + 0), slice, FlagX);
          DxVy  = L(W,(pos + 1), slice, FlagX);
          DxVz  = L(W,(pos + 2), slice, FlagX);
          // Ly
          DyVx  = L(W,(pos + 0), segment, FlagY);
          DyVy  = L(W,(pos + 1), segment, FlagY);
          DyVz  = L(W,(pos + 2), segment, FlagY);

          // Get Lz by Lx and Ly
          if(k == nk2-1) {
            DzVx = matVx2Vz[ij+3*0+0] * DxVx
                 + matVx2Vz[ij+3*0+1] * DxVy
                 + matVx2Vz[ij+3*0+2] * DxVz
                 + matVy2Vz[ij+3*0+0] * DyVx
                 + matVy2Vz[ij+3*0+1] * DyVy
                 + matVy2Vz[ij+3*0+2] * DyVz;

            DzVy = matVx2Vz[ij+3*1+0] * DxVx
                 + matVx2Vz[ij+3*1+1] * DxVy
                 + matVx2Vz[ij+3*1+2] * DxVz
                 + matVy2Vz[ij+3*1+0] * DyVx
                 + matVy2Vz[ij+3*1+1] * DyVy
                 + matVy2Vz[ij+3*1+2] * DyVz;

            DzVz = matVx2Vz[ij+3*2+0] * DxVx
                 + matVx2Vz[ij+3*2+1] * DxVy
                 + matVx2Vz[ij+3*2+2] * DxVz
                 + matVy2Vz[ij+3*2+0] * DyVx
                 + matVy2Vz[ij+3*2+1] * DyVy
                 + matVy2Vz[ij+3*2+2] * DyVz;
          } else if(k == nk2-2) {
            // Get Lz directly
            DzVx = L22(W,(pos + 0), WSIZE, FlagZ);
            DzVy = L22(W,(pos + 1), WSIZE, FlagZ);
            DzVz = L22(W,(pos + 2), WSIZE, FlagZ);
          }  else {
            // Get Lz directly
            DzVx = L24(W,(pos + 0), WSIZE, FlagZ);
            DzVy = L24(W,(pos + 1), WSIZE, FlagZ);
            DzVz = L24(W,(pos + 2), WSIZE, FlagZ);
          }

          hW[pos + 3] = (lam2mu * xix * DxVx + lam * xiy * DxVy + lam * xiz * DxVz) * rb1
                      + (lam2mu * etx * DyVx + lam * ety * DyVy + lam * etz * DyVz) * rb2
                      + (lam2mu * ztx * DzVx + lam * zty * DzVy + lam * ztz * DzVz) * rb3;

          hW[pos + 4] = (lam * xix * DxVx + lam2mu * xiy * DxVy + lam * xiz * DxVz) * rb1
                      + (lam * etx * DyVx + lam2mu * ety * DyVy + lam * etz * DyVz) * rb2
                      + (lam * ztx * DzVx + lam2mu * zty * DzVy + lam * ztz * DzVz) * rb3;

          hW[pos + 5] = (lam * xix * DxVx + lam * xiy * DxVy + lam2mu * xiz * DxVz) * rb1
                      + (lam * etx * DyVx + lam * ety * DyVy + lam2mu * etz * DyVz) * rb2
                      + (lam * ztx * DzVx + lam * zty * DzVy + lam2mu * ztz * DzVz) * rb3;

          hW[pos + 6] = mu * ( (xiy * DxVx + xix * DxVy) * rb1
                             + (ety * DyVx + etx * DyVy) * rb2
                             + (zty * DzVx + ztx * DzVy) * rb3 );

          hW[pos + 7] = mu * ( (xiz * DxVx + xix * DxVz) * rb1
                             + (etz * DyVx + etx * DyVz) * rb2
                             + (ztz * DzVx + ztx * DzVz) * rb3 );

          hW[pos + 8] = mu * ( (xiz * DxVy + xiy * DxVz) * rb1
                             + (etz * DyVy + ety * DyVz) * rb2
                             + (ztz * DzVy + zty * DzVz) * rb3 );


        } //end of freenode
#endif

#ifdef usePML
        float hW_add[WSIZE];
        memset(hW_add, 0, sizeof(float) * WSIZE);
        for(iii = 0 ; iii < 6 ; iii ++) {
          //x direction
          if((iii == 0 && isx1 && i < PML_ND + ni1 ) || (iii == 1 && isx2 && i >= ni2 - PML_ND)) {
            int i1 = i;
            int j1 = j;
            int k1 = k;

            int pni = PML_ND;
            int pnj = ny;
            int pnk = lnz;

            if(iii == 0) {
              idx0 = ni1;
              Aux_Wx = Aux->Wx1;
              Aux_hWx = Aux->hWx1;
            }
            else if(iii == 1) {
              idx0 = ni2 - PML_ND;
              Aux_Wx = Aux->Wx2;
              Aux_hWx = Aux->hWx2;
            }

            float Deriv[WSIZE];
            if ((i < PML_ND + ni1 || i >= ni2 - PML_ND) && j < pnj && k < pnk) {
              // idx0 = 0 or ni-ND or nj-ND or nk-ND
              pos = (k1 + j1 * lnz + i1 * ny * lnz) * WSIZE;

              DxVx  = L(W, (pos + 0), slice, FlagX);
              DxVy  = L(W, (pos + 1), slice, FlagX);
              DxVz  = L(W, (pos + 2), slice, FlagX);
              DxTxx = L(W, (pos + 3), slice, FlagX);
              DxTyy = L(W, (pos + 4), slice, FlagX);
              DxTzz = L(W, (pos + 5), slice, FlagX);
              DxTxy = L(W, (pos + 6), slice, FlagX);
              DxTxz = L(W, (pos + 7), slice, FlagX);
              DxTyz = L(W, (pos + 8), slice, FlagX);
              int pos_m = (i * ny * lnz + j * lnz + k) * MSIZE;

              float *p0 = pml;
              b1 = p0[ni + (i - ni1)];
              d1 = p0[ni * 2 + (i - ni1)];
              ad1 = -(p0[(i - ni1)] + d1);

              lam = M[pos_m + 10]; mu = M[pos_m + 11];
              rrho = 1.0f/M[pos_m + 12];
              lam2mu = lam + 2.0f*mu;
              rb1 = 1.0f/b1;

              Deriv[0] = (DxTxx*xix + DxTxy*xiy + DxTxz*xiz)*rrho;
              Deriv[1] = (DxTxy*xix + DxTyy*xiy + DxTyz*xiz)*rrho;
              Deriv[2] = (DxTxz*xix + DxTyz*xiy + DxTzz*xiz)*rrho;
              Deriv[3] = (DxVx*xix*lam2mu + DxVy*xiy*lam    + DxVz*xiz*lam   );
              Deriv[4] = (DxVx*xix*lam    + DxVy*xiy*lam2mu + DxVz*xiz*lam   );
              Deriv[5] = (DxVx*xix*lam    + DxVy*xiy*lam    + DxVz*xiz*lam2mu);
              Deriv[6] = (DxVx*xiy + DxVy*xix)*mu;
              Deriv[7] = (DxVx*xiz + DxVz*xix)*mu;
              Deriv[8] = (DxVy*xiz + DxVz*xiy)*mu;

              // Auxiliary Equations
              pos_a = (k + j * pnk + (i - idx0) * pnj * pnk) * WSIZE;

              hW_add[0] += -Aux_Wx[pos_a + 0] * rb1;
              hW_add[1] += -Aux_Wx[pos_a + 1] * rb1;
              hW_add[2] += -Aux_Wx[pos_a + 2] * rb1;
              hW_add[3] += -Aux_Wx[pos_a + 3] * rb1;
              hW_add[4] += -Aux_Wx[pos_a + 4] * rb1;
              hW_add[5] += -Aux_Wx[pos_a + 5] * rb1;
              hW_add[6] += -Aux_Wx[pos_a + 6] * rb1;
              hW_add[7] += -Aux_Wx[pos_a + 7] * rb1;
              hW_add[8] += -Aux_Wx[pos_a + 8] * rb1;

              Aux_hWx[pos_a + 0] = ad1 * Aux_Wx[pos_a + 0] + d1 * Deriv[0];
              Aux_hWx[pos_a + 1] = ad1 * Aux_Wx[pos_a + 1] + d1 * Deriv[1];
              Aux_hWx[pos_a + 2] = ad1 * Aux_Wx[pos_a + 2] + d1 * Deriv[2];
              Aux_hWx[pos_a + 3] = ad1 * Aux_Wx[pos_a + 3] + d1 * Deriv[3];
              Aux_hWx[pos_a + 4] = ad1 * Aux_Wx[pos_a + 4] + d1 * Deriv[4];
              Aux_hWx[pos_a + 5] = ad1 * Aux_Wx[pos_a + 5] + d1 * Deriv[5];
              Aux_hWx[pos_a + 6] = ad1 * Aux_Wx[pos_a + 6] + d1 * Deriv[6];
              Aux_hWx[pos_a + 7] = ad1 * Aux_Wx[pos_a + 7] + d1 * Deriv[7];
              Aux_hWx[pos_a + 8] = ad1 * Aux_Wx[pos_a + 8] + d1 * Deriv[8];

#ifdef FreeSurface
              if(freenode && k == nk2-1) {
                // Get Lz by Lx and Ly
                int ij = (j1 + i1 * ny) * 9;
                DzVx1 = matVx2Vz[ij+3*0+0] * DxVx
                  + matVx2Vz[ij+3*0+1] * DxVy
                  + matVx2Vz[ij+3*0+2] * DxVz;

                DzVy1 = matVx2Vz[ij+3*1+0] * DxVx
                  + matVx2Vz[ij+3*1+1] * DxVy
                  + matVx2Vz[ij+3*1+2] * DxVz;

                DzVz1 = matVx2Vz[ij+3*2+0] * DxVx
                  + matVx2Vz[ij+3*2+1] * DxVy
                  + matVx2Vz[ij+3*2+2] * DxVz;

                ztx = M[pos_m + 6]; zty = M[pos_m + 7]; ztz = M[pos_m + 8];

                Aux_hWx[pos_a + 3] += d1*b1*(lam2mu*ztx*DzVx1 + lam   *zty*DzVy1 + lam   *ztz*DzVz1);
                Aux_hWx[pos_a + 4] += d1*b1*(lam   *ztx*DzVx1 + lam2mu*zty*DzVy1 + lam   *ztz*DzVz1);
                Aux_hWx[pos_a + 5] += d1*b1*(lam   *ztx*DzVx1 + lam   *zty*DzVy1 + lam2mu*ztz*DzVz1);
                Aux_hWx[pos_a + 6] += d1*b1*mu*(zty*DzVx1 + ztx*DzVy1);
                Aux_hWx[pos_a + 7] += d1*b1*mu*(ztz*DzVx1 + ztx*DzVz1);
                Aux_hWx[pos_a + 8] += d1*b1*mu*(ztz*DzVy1 + zty*DzVz1);
              }
#endif
            }
          }
          //y direction
          if((iii == 2 && isy1 && j < PML_ND + nj1) || (iii == 3 && isy2 && j >= nj2 - PML_ND)) {
            int i1 = i;
            int j1 = j;
            int k1 = k;

            int pni = nx;
            int pnj = PML_ND;
            int pnk = lnz;

            if(iii == 2) {
              idx0 = nj1;
              Aux_Wy = Aux->Wy1;
              Aux_hWy = Aux->hWy1;
            }
            else if(iii == 3) {
              idx0 = nj2 - PML_ND;
              Aux_Wy = Aux->Wy2;
              Aux_hWy = Aux->hWy2;
            }

            float Deriv[WSIZE];
            if ((j < PML_ND + nj1 || j >= nj2 - PML_ND) && i < pni && k < pnk) {
              pos = (k1 + j1 * lnz + i1 * ny * lnz) * WSIZE;
              DxVx  = L(W, (pos + 0), segment, FlagY);
              DxVy  = L(W, (pos + 1), segment, FlagY);
              DxVz  = L(W, (pos + 2), segment, FlagY);
              DxTxx = L(W, (pos + 3), segment, FlagY);
              DxTyy = L(W, (pos + 4), segment, FlagY);
              DxTzz = L(W, (pos + 5), segment, FlagY);
              DxTxy = L(W, (pos + 6), segment, FlagY);
              DxTxz = L(W, (pos + 7), segment, FlagY);
              DxTyz = L(W, (pos + 8), segment, FlagY);
              int pos_m = (i * ny * lnz + j * lnz + k) * MSIZE;

              float *p1 = pml + ni * 3;
              b1 = p1[nj + (j - nj1)];
              d1 = p1[nj * 2 + (j - nj1)];
              ad1 = -(p1[(j - nj1)] + d1);

              lam = M[pos_m + 10]; mu = M[pos_m + 11];
              rrho = 1.0f/M[pos_m + 12];
              lam2mu = lam + 2.0f*mu;
              rb1 = 1.0f/b1;

              //add
              xix = M[pos_m + 3]; 
              xiy = M[pos_m + 4]; 
              xiz = M[pos_m + 5];

              Deriv[0] = (DxTxx*xix + DxTxy*xiy + DxTxz*xiz)*rrho;
              Deriv[1] = (DxTxy*xix + DxTyy*xiy + DxTyz*xiz)*rrho;
              Deriv[2] = (DxTxz*xix + DxTyz*xiy + DxTzz*xiz)*rrho;
              Deriv[3] = (DxVx*xix*lam2mu + DxVy*xiy*lam    + DxVz*xiz*lam   );
              Deriv[4] = (DxVx*xix*lam    + DxVy*xiy*lam2mu + DxVz*xiz*lam   );
              Deriv[5] = (DxVx*xix*lam    + DxVy*xiy*lam    + DxVz*xiz*lam2mu);
              Deriv[6] = (DxVx*xiy + DxVy*xix)*mu;
              Deriv[7] = (DxVx*xiz + DxVz*xix)*mu;
              Deriv[8] = (DxVy*xiz + DxVz*xiy)*mu;

              pos_a = (k + (j - idx0) * pnk + i * pnj * pnk) * WSIZE;
              hW_add[0] += -Aux_Wy[pos_a + 0] * rb1;
              hW_add[1] += -Aux_Wy[pos_a + 1] * rb1;
              hW_add[2] += -Aux_Wy[pos_a + 2] * rb1;
              hW_add[3] += -Aux_Wy[pos_a + 3] * rb1;
              hW_add[4] += -Aux_Wy[pos_a + 4] * rb1;
              hW_add[5] += -Aux_Wy[pos_a + 5] * rb1;
              hW_add[6] += -Aux_Wy[pos_a + 6] * rb1;
              hW_add[7] += -Aux_Wy[pos_a + 7] * rb1;
              hW_add[8] += -Aux_Wy[pos_a + 8] * rb1;

              Aux_hWy[pos_a + 0] = ad1 * Aux_Wy[pos_a + 0] + d1 * Deriv[0];
              Aux_hWy[pos_a + 1] = ad1 * Aux_Wy[pos_a + 1] + d1 * Deriv[1];
              Aux_hWy[pos_a + 2] = ad1 * Aux_Wy[pos_a + 2] + d1 * Deriv[2];
              Aux_hWy[pos_a + 3] = ad1 * Aux_Wy[pos_a + 3] + d1 * Deriv[3];
              Aux_hWy[pos_a + 4] = ad1 * Aux_Wy[pos_a + 4] + d1 * Deriv[4];
              Aux_hWy[pos_a + 5] = ad1 * Aux_Wy[pos_a + 5] + d1 * Deriv[5];
              Aux_hWy[pos_a + 6] = ad1 * Aux_Wy[pos_a + 6] + d1 * Deriv[6];
              Aux_hWy[pos_a + 7] = ad1 * Aux_Wy[pos_a + 7] + d1 * Deriv[7];
              Aux_hWy[pos_a + 8] = ad1 * Aux_Wy[pos_a + 8] + d1 * Deriv[8];

#ifdef FreeSurface
              if(freenode && k == nk2-1) {
                // Get Lz by Lx and Ly
                int ij = (j1 + i1 * ny) * 9;
                DzVx1 = matVy2Vz[ij+3*0+0] * DxVx
                  + matVy2Vz[ij+3*0+1] * DxVy
                  + matVy2Vz[ij+3*0+2] * DxVz;

                DzVy1 = matVy2Vz[ij+3*1+0] * DxVx
                  + matVy2Vz[ij+3*1+1] * DxVy
                  + matVy2Vz[ij+3*1+2] * DxVz;

                DzVz1 = matVy2Vz[ij+3*2+0] * DxVx
                  + matVy2Vz[ij+3*2+1] * DxVy
                  + matVy2Vz[ij+3*2+2] * DxVz;

                ztx = M[pos_m + 6]; zty = M[pos_m + 7]; ztz = M[pos_m + 8];

                Aux_hWy[pos_a + 3] += d1*b1*(lam2mu*ztx*DzVx1 + lam   *zty*DzVy1 + lam   *ztz*DzVz1);
                Aux_hWy[pos_a + 4] += d1*b1*(lam   *ztx*DzVx1 + lam2mu*zty*DzVy1 + lam   *ztz*DzVz1);
                Aux_hWy[pos_a + 5] += d1*b1*(lam   *ztx*DzVx1 + lam   *zty*DzVy1 + lam2mu*ztz*DzVz1);
                Aux_hWy[pos_a + 6] += d1*b1*mu*(zty*DzVx1 + ztx*DzVy1);
                Aux_hWy[pos_a + 7] += d1*b1*mu*(ztz*DzVx1 + ztx*DzVz1);
                Aux_hWy[pos_a + 8] += d1*b1*mu*(ztz*DzVy1 + zty*DzVz1);
              }
#endif
            }
          }
          //z direction
          if((iii == 4 && isz1 && k < PML_ND + nk1 ) || (iii == 5 && isz2 && k >= nk2 - PML_ND)) {
            int i1 = i;
            int j1 = j;
            int k1 = k;

            int pni = nx;
            int pnj = ny;
            int pnk = PML_ND;

            if(iii == 4) {
              idx0 = nk1;
              Aux_Wz = Aux->Wz1;
              Aux_hWz = Aux->hWz1;
            }
            else if(iii == 5) {
              idx0 = nk2 - PML_ND;
              Aux_Wz = Aux->Wz2;
              Aux_hWz = Aux->hWz2;
            }

            float Deriv[WSIZE];
            if ((k < PML_ND + nk1 || k >= nk2 - PML_ND) && i < pni && j < pnj) {
              // idx0 = 0 or ni-ND or nj-ND or nk-ND
              if (k == 19 || k == 25)
                printf("debug PML z k is not normal! i:%d, j:%d, k:%d, isz1:%d, isz2:%d\n",i, j, k, isz1, isz2);
              // printf("debug PML z i:%d, j:%d, k:%d, isz1:%d, isz2:%d\n",i, j, k, isz1, isz2);
              pos = (k1 + j1 * lnz + i1 * ny * lnz) * WSIZE;

              DxVx  = L(W, (pos + 0), WSIZE, FlagZ);
              DxVy  = L(W, (pos + 1), WSIZE, FlagZ);
              DxVz  = L(W, (pos + 2), WSIZE, FlagZ);
              DxTxx = L(W, (pos + 3), WSIZE, FlagZ);
              DxTyy = L(W, (pos + 4), WSIZE, FlagZ);
              DxTzz = L(W, (pos + 5), WSIZE, FlagZ);
              DxTxy = L(W, (pos + 6), WSIZE, FlagZ);
              DxTxz = L(W, (pos + 7), WSIZE, FlagZ);
              DxTyz = L(W, (pos + 8), WSIZE, FlagZ);
              int pos_m = (i * ny * lnz + j * lnz + k) * MSIZE;

              float *p2 = pml + ni * 3 + nj * 3;
              b1 = p2[nk + (k - nk1)];
              d1 = p2[nk * 2 + (k - nk1)];
              ad1 = -(p2[(k - nk1)] + d1);

              lam = M[pos_m + 10]; mu = M[pos_m + 11];
              rrho = 1.0f/M[pos_m + 12];
              lam2mu = lam + 2.0f*mu;
              rb1 = 1.0f/b1;

              //add
              xix = M[pos_m + 6]; 
              xiy = M[pos_m + 7]; 
              xiz = M[pos_m + 8];

              Deriv[0] = (DxTxx*xix + DxTxy*xiy + DxTxz*xiz)*rrho;
              Deriv[1] = (DxTxy*xix + DxTyy*xiy + DxTyz*xiz)*rrho;
              Deriv[2] = (DxTxz*xix + DxTyz*xiy + DxTzz*xiz)*rrho;
              Deriv[3] = (DxVx*xix*lam2mu + DxVy*xiy*lam    + DxVz*xiz*lam   );
              Deriv[4] = (DxVx*xix*lam    + DxVy*xiy*lam2mu + DxVz*xiz*lam   );
              Deriv[5] = (DxVx*xix*lam    + DxVy*xiy*lam    + DxVz*xiz*lam2mu);
              Deriv[6] = (DxVx*xiy + DxVy*xix)*mu;
              Deriv[7] = (DxVx*xiz + DxVz*xix)*mu;
              Deriv[8] = (DxVy*xiz + DxVz*xiy)*mu;

              // Auxiliary Equations
              pos_a = ((k - idx0) + j * pnk + i * pnj * pnk) * WSIZE;

              hW_add[0] += -Aux_Wz[pos_a + 0] * rb1;
              hW_add[1] += -Aux_Wz[pos_a + 1] * rb1;
              hW_add[2] += -Aux_Wz[pos_a + 2] * rb1;
              hW_add[3] += -Aux_Wz[pos_a + 3] * rb1;
              hW_add[4] += -Aux_Wz[pos_a + 4] * rb1;
              hW_add[5] += -Aux_Wz[pos_a + 5] * rb1;
              hW_add[6] += -Aux_Wz[pos_a + 6] * rb1;
              hW_add[7] += -Aux_Wz[pos_a + 7] * rb1;
              hW_add[8] += -Aux_Wz[pos_a + 8] * rb1;

              Aux_hWz[pos_a + 0] = ad1 * Aux_Wz[pos_a + 0] + d1 * Deriv[0];
              Aux_hWz[pos_a + 1] = ad1 * Aux_Wz[pos_a + 1] + d1 * Deriv[1];
              Aux_hWz[pos_a + 2] = ad1 * Aux_Wz[pos_a + 2] + d1 * Deriv[2];
              Aux_hWz[pos_a + 3] = ad1 * Aux_Wz[pos_a + 3] + d1 * Deriv[3];
              Aux_hWz[pos_a + 4] = ad1 * Aux_Wz[pos_a + 4] + d1 * Deriv[4];
              Aux_hWz[pos_a + 5] = ad1 * Aux_Wz[pos_a + 5] + d1 * Deriv[5];
              Aux_hWz[pos_a + 6] = ad1 * Aux_Wz[pos_a + 6] + d1 * Deriv[6];
              Aux_hWz[pos_a + 7] = ad1 * Aux_Wz[pos_a + 7] + d1 * Deriv[7];
              Aux_hWz[pos_a + 8] = ad1 * Aux_Wz[pos_a + 8] + d1 * Deriv[8];
            }
          }
        }
        for(iiii = 0 ; iiii < WSIZE ; iiii ++) {
          hW[pos + iiii] += hW_add[iiii];
        }
        fflush(stdout);
#endif

        pos = (i * ny * lnz + j * lnz + k)  *  WSIZE;

          if(flag == 1){ //RK_begin
            tW[pos + 0] = mW[pos + 0] + rkb * hW[pos + 0];
            tW[pos + 1] = mW[pos + 1] + rkb * hW[pos + 1];
            tW[pos + 2] = mW[pos + 2] + rkb * hW[pos + 2];
            tW[pos + 3] = mW[pos + 3] + rkb * hW[pos + 3];
            tW[pos + 4] = mW[pos + 4] + rkb * hW[pos + 4];
            tW[pos + 5] = mW[pos + 5] + rkb * hW[pos + 5];
            tW[pos + 6] = mW[pos + 6] + rkb * hW[pos + 6];
            tW[pos + 7] = mW[pos + 7] + rkb * hW[pos + 7];
            tW[pos + 8] = mW[pos + 8] + rkb * hW[pos + 8];
          }

          if(flag == 2) { // RK_inner
            tW[pos + 0] += rkb * hW[pos + 0];
            tW[pos + 1] += rkb * hW[pos + 1];
            tW[pos + 2] += rkb * hW[pos + 2];
            tW[pos + 3] += rkb * hW[pos + 3];
            tW[pos + 4] += rkb * hW[pos + 4];
            tW[pos + 5] += rkb * hW[pos + 5];
            tW[pos + 6] += rkb * hW[pos + 6];
            tW[pos + 7] += rkb * hW[pos + 7];
            tW[pos + 8] += rkb * hW[pos + 8];
          }

          hW[pos + 0] = mW[pos + 0] + rka * hW[pos + 0];
          hW[pos + 1] = mW[pos + 1] + rka * hW[pos + 1];
          hW[pos + 2] = mW[pos + 2] + rka * hW[pos + 2];
          hW[pos + 3] = mW[pos + 3] + rka * hW[pos + 3];
          hW[pos + 4] = mW[pos + 4] + rka * hW[pos + 4];
          hW[pos + 5] = mW[pos + 5] + rka * hW[pos + 5];
          hW[pos + 6] = mW[pos + 6] + rka * hW[pos + 6];
          hW[pos + 7] = mW[pos + 7] + rka * hW[pos + 7];
          hW[pos + 8] = mW[pos + 8] + rka * hW[pos + 8];

          //deriv_abs rk
#ifdef usePML
          for(iii = 0 ; iii < 6 ; iii ++) {
            if((iii == 0 && isx1 && i < PML_ND + ni1) || (iii == 1 && isx2 && i >= ni2 - PML_ND) || (iii == 2 && isy1 && j < PML_ND + nj1) || (iii == 3 && isy2 && j >= nj2 - PML_ND) || (iii == 4 && isz1 && k < PML_ND + nk1) || (iii == 5 && isz2 && k >= nk2 - PML_ND)) {
              if(iii == 0) {
                idx0 = ni1;
                W = Aux->Wx1;
                hW = Aux->hWx1;
                if(flag != 3) {
                  tW = Aux->tWx1;
                  mW = Aux->mWx1;
                }
                else {
                  tW = Aux->mWx1;
                  mW = Aux->tWx1;
                }
              }
              else if(iii == 1) {
                idx0 = ni2 - PML_ND;
                W = Aux->Wx2;
                hW = Aux->hWx2;
                if(flag != 3) {
                  tW = Aux->tWx2;
                  mW = Aux->mWx2;
                }
                else {
                  tW = Aux->mWx2;
                  mW = Aux->tWx2;
                }
              }
              else if(iii == 2) {
                idx0 = nj1;
                W = Aux->Wy1;
                hW = Aux->hWy1;
                if(flag != 3) {
                  tW = Aux->tWy1;
                  mW = Aux->mWy1;
                }
                else {
                  tW = Aux->mWy1;
                  mW = Aux->tWy1;
                }
              }
              else if(iii == 3) {
                idx0 = nj2 - PML_ND;
                W = Aux->Wy2;
                hW = Aux->hWy2;
                if(flag != 3) {
                  tW = Aux->tWy2;
                  mW = Aux->mWy2;
                }
                else {
                  tW = Aux->mWy2;
                  mW = Aux->tWy2;
                }
              }
              else if(iii == 4) {
                idx0 = nk1;
                W = Aux->Wz1;
                hW = Aux->hWz1;
                if(flag != 3) {
                  tW = Aux->tWz1;
                  mW = Aux->mWz1;
                }
                else {
                  tW = Aux->mWz1;
                  mW = Aux->tWz1;
                }
              }
              else if(iii == 5) {
                idx0 = nk2 - PML_ND;
                W = Aux->Wz2;
                hW = Aux->hWz2;
                if(flag != 3) {
                  tW = Aux->tWz2;
                  mW = Aux->mWz2;
                }
                else {
                  tW = Aux->mWz2;
                  mW = Aux->tWz2;
                }
              }
              if(iii == 0 || iii == 1)
                pos_a = (k + j * lnz + (i - idx0) * ny * lnz) * WSIZE;
              if(iii == 2 || iii == 3)
                pos_a = (k + (j - idx0) * lnz + i * PML_ND * lnz) * WSIZE;
              if(iii == 4 || iii == 5)
                pos_a = ((k - idx0) + j * PML_ND + i * ny * PML_ND) * WSIZE;

              if(flag == 1){ //RK_begin
                tW[pos_a + 0] = mW[pos_a + 0] + rkb * hW[pos_a + 0];
                tW[pos_a + 1] = mW[pos_a + 1] + rkb * hW[pos_a + 1];
                tW[pos_a + 2] = mW[pos_a + 2] + rkb * hW[pos_a + 2];
                tW[pos_a + 3] = mW[pos_a + 3] + rkb * hW[pos_a + 3];
                tW[pos_a + 4] = mW[pos_a + 4] + rkb * hW[pos_a + 4];
                tW[pos_a + 5] = mW[pos_a + 5] + rkb * hW[pos_a + 5];
                tW[pos_a + 6] = mW[pos_a + 6] + rkb * hW[pos_a + 6];
                tW[pos_a + 7] = mW[pos_a + 7] + rkb * hW[pos_a + 7];
                tW[pos_a + 8] = mW[pos_a + 8] + rkb * hW[pos_a + 8];
              }

              if(flag == 2) { // RK_inner
                tW[pos_a + 0] += rkb * hW[pos_a + 0];
                tW[pos_a + 1] += rkb * hW[pos_a + 1];
                tW[pos_a + 2] += rkb * hW[pos_a + 2];
                tW[pos_a + 3] += rkb * hW[pos_a + 3];
                tW[pos_a + 4] += rkb * hW[pos_a + 4];
                tW[pos_a + 5] += rkb * hW[pos_a + 5];
                tW[pos_a + 6] += rkb * hW[pos_a + 6];
                tW[pos_a + 7] += rkb * hW[pos_a + 7];
                tW[pos_a + 8] += rkb * hW[pos_a + 8];
              }

              hW[pos_a + 0] = mW[pos_a + 0] + rka * hW[pos_a + 0];
              hW[pos_a + 1] = mW[pos_a + 1] + rka * hW[pos_a + 1];
              hW[pos_a + 2] = mW[pos_a + 2] + rka * hW[pos_a + 2];
              hW[pos_a + 3] = mW[pos_a + 3] + rka * hW[pos_a + 3];
              hW[pos_a + 4] = mW[pos_a + 4] + rka * hW[pos_a + 4];
              hW[pos_a + 5] = mW[pos_a + 5] + rka * hW[pos_a + 5];
              hW[pos_a + 6] = mW[pos_a + 6] + rka * hW[pos_a + 6];
              hW[pos_a + 7] = mW[pos_a + 7] + rka * hW[pos_a + 7];
              hW[pos_a + 8] = mW[pos_a + 8] + rka * hW[pos_a + 8];
              W = oW;
              hW = ohW;
              tW = otW;
              mW = omW;
            }
          }
#endif

        // ================end of RK=================

      } // end of for i j k
    }  // end of for i j k
  }   // end of for i j k

//  Delloc3D(matVx2Vz);
//  Delloc3D(matVy2Vz);
//  Delloc3D(matF2Vz);

return 0;
}


int RK_Syn(float *W, float *M, float *matVx2Vz, float *matVy2Vz, int istep, float *pml, int isx1, int isx2, int isy1, int isy2, int isz1, int isz2, struct aux *Aux, float DH, int freenode){
  // step 1: BBB FFF BBB FFF
  // step 2: FFB BBF FFB BBF
  // step 3: FFF BBB FFF BBB
  // step 4: BBF FFB BBF FFB
  // step 5: BFB FBF BFB FBF
  // step 6: FBB BFF FBB BFF
  // step 7: FBF BFB FBF BFB
  // step 8: BFF FBB BFF FBB
  int Flag[3]; 
  int FlagMinus[3];

  switch (istep){
    default:
      break;
    case 1:
      Flag[0] = BWD; Flag[1] = BWD; Flag[2] = BWD;
      break;
    case 2:
      Flag[0] = FWD; Flag[1] = FWD; Flag[2] = BWD;
      break;
    case 3:
      Flag[0] = FWD; Flag[1] = FWD; Flag[2] = FWD;
      break;
    case 4:
      Flag[0] = BWD; Flag[1] = BWD; Flag[2] = FWD;
      break;
    case 5:
      Flag[0] = BWD; Flag[1] = FWD; Flag[2] = BWD;
      break;
    case 6:
      Flag[0] = FWD; Flag[1] = BWD; Flag[2] = BWD;
      break;
    case 7:
      Flag[0] = FWD; Flag[1] = BWD; Flag[2] = FWD;
      break;
    case 8:
      Flag[0] = BWD; Flag[1] = FWD; Flag[2] = FWD;
      break;
  }

  FlagMinus[0] = -Flag[0];
  FlagMinus[1] = -Flag[1];
  FlagMinus[2] = -Flag[2];

  float *hW = Alloc_wave();
  float *mW = Alloc_wave();
  float *tW = Alloc_wave();
  float *T;
  RK_init(W,hW);
  exchange_Wave(W);
  RK_init(W,mW);
#ifdef usePML
  if(isx1) {
    RK_init2(Aux->Wx1,Aux->hWx1, PML_ND, ny, lnz);
    RK_init2(Aux->Wx1,Aux->mWx1, PML_ND, ny, lnz);
  }
  if(isx2) {
    RK_init2(Aux->Wx2,Aux->hWx2, PML_ND, ny, lnz);
    RK_init2(Aux->Wx2,Aux->mWx2, PML_ND, ny, lnz);
  }
  if(isy1) {
    RK_init2(Aux->Wy1,Aux->hWy1, nx, PML_ND, lnz);
    RK_init2(Aux->Wy1,Aux->mWy1, nx, PML_ND, lnz);
  }
  if(isy2) {
    RK_init2(Aux->Wy2,Aux->hWy2, nx, PML_ND, lnz);
    RK_init2(Aux->Wy2,Aux->mWy2, nx, PML_ND, lnz);
  }
  if(isz1) {
    RK_init2(Aux->Wz1,Aux->hWz1, nx, ny, PML_ND);
    RK_init2(Aux->Wz1,Aux->mWz1, nx, ny, PML_ND);
  }
  if(isz2) {
    RK_init2(Aux->Wz2,Aux->hWz2, nx, ny, PML_ND);
    RK_init2(Aux->Wz2,Aux->mWz2, nx, ny, PML_ND);
  }
#endif
  Calculate(W,mW,hW,tW,M,matVx2Vz,matVy2Vz,Flag,RK4a[0],RK4b[0],1, 0, pml, Aux, isx1, isx2, isy1, isy2, isz1, isz2);
#ifdef usePML
  /*
  if(isx1) 
    abs_deriv_x(W, hW, Aux, M, pml, istep, 0, Flag[0], 0, DH, freenode, matVx2Vz);
  if(isx2)
    abs_deriv_x(W, hW, Aux, M, pml, istep, 0, Flag[0], ni - PML_ND, DH, freenode, matVx2Vz);
  if(isy1)
    abs_deriv_y(W, hW, Aux, M, pml, istep, 0, Flag[1], 0, DH, freenode, matVy2Vz);
  if(isy2)
    abs_deriv_y(W, hW, Aux, M, pml, istep, 0, Flag[1], nj - PML_ND, DH, freenode, matVy2Vz);
  if(isz1)
    abs_deriv_z(W, hW, Aux, M, pml, istep, 0, Flag[2], 0, DH, freenode);
  if(isz2)
    abs_deriv_z(W, hW, Aux, M, pml, istep, 0, Flag[2], nk - PML_ND, DH, freenode);
    */
#endif
  T = hW; hW = W; W = T;
#ifdef usePML
  if(isx1) {
    T = Aux->hWx1; Aux->hWx1 = Aux->Wx1; Aux->Wx1 = T;
  }
  if(isx2) {
    T = Aux->hWx2; Aux->hWx2 = Aux->Wx2; Aux->Wx2 = T;
  }
  if(isy1) {
    T = Aux->hWy1; Aux->hWy1 = Aux->Wy1; Aux->Wy1 = T;
  }
  if(isy2) {
    T = Aux->hWy2; Aux->hWy2 = Aux->Wy2; Aux->Wy2 = T;
  }
  if(isz1) {
    T = Aux->hWz1; Aux->hWz1 = Aux->Wz1; Aux->Wz1 = T;
  }
  if(isz2) {
    T = Aux->hWz2; Aux->hWz2 = Aux->Wz2; Aux->Wz2 = T;
  }
#endif

  exchange_Wave(W);
  Calculate(W,mW,hW,tW,M,matVx2Vz,matVy2Vz,FlagMinus,RK4a[1],RK4b[1],2, 0, pml, Aux, isx1, isx2, isy1, isy2, isz1, isz2);
#ifdef usePML
  /*
  if(isx1) 
    abs_deriv_x(W, hW, Aux, M, pml, istep, 1, Flag[0], 0, DH, freenode, matVx2Vz);
  if(isx2)
    abs_deriv_x(W, hW, Aux, M, pml, istep, 1, Flag[0], ni - PML_ND, DH, freenode, matVx2Vz);
  if(isy1)
    abs_deriv_y(W, hW, Aux, M, pml, istep, 1, Flag[1], 0, DH, freenode, matVy2Vz);
  if(isy2)
    abs_deriv_y(W, hW, Aux, M, pml, istep, 1, Flag[1], nj - PML_ND, DH, freenode, matVy2Vz);
  if(isz1)
    abs_deriv_z(W, hW, Aux, M, pml, istep, 1, Flag[2], 0, DH, freenode);
  if(isz2)
    abs_deriv_z(W, hW, Aux, M, pml, istep, 1, Flag[2], nk - PML_ND, DH, freenode);
    */
#endif
  T = hW; hW = W; W = T;
#ifdef usePML
  if(isx1) {
    T = Aux->hWx1; Aux->hWx1 = Aux->Wx1; Aux->Wx1 = T;
  }
  if(isx2) {
    T = Aux->hWx2; Aux->hWx2 = Aux->Wx2; Aux->Wx2 = T;
  }
  if(isy1) {
    T = Aux->hWy1; Aux->hWy1 = Aux->Wy1; Aux->Wy1 = T;
  }
  if(isy2) {
    T = Aux->hWy2; Aux->hWy2 = Aux->Wy2; Aux->Wy2 = T;
  }
  if(isz1) {
    T = Aux->hWz1; Aux->hWz1 = Aux->Wz1; Aux->Wz1 = T;
  }
  if(isz2) {
    T = Aux->hWz2; Aux->hWz2 = Aux->Wz2; Aux->Wz2 = T;
  }
#endif

  exchange_Wave(W);
  Calculate(W,mW,hW,tW,M,matVx2Vz,matVy2Vz,Flag,RK4a[2],RK4b[2],2, 0, pml, Aux, isx1, isx2, isy1, isy2, isz1, isz2);
#ifdef usePML
  /*
  if(isx1) 
    abs_deriv_x(W, hW, Aux, M, pml, istep, 2, Flag[0], 0, DH, freenode, matVx2Vz);
  if(isx2)
    abs_deriv_x(W, hW, Aux, M, pml, istep, 2, Flag[0], ni - PML_ND, DH, freenode, matVx2Vz);
  if(isy1)
    abs_deriv_y(W, hW, Aux, M, pml, istep, 2, Flag[1], 0, DH, freenode, matVy2Vz);
  if(isy2)
    abs_deriv_y(W, hW, Aux, M, pml, istep, 2, Flag[1], nj - PML_ND, DH, freenode, matVy2Vz);
  if(isz1)
    abs_deriv_z(W, hW, Aux, M, pml, istep, 2, Flag[2], 0, DH, freenode);
  if(isz2)
    abs_deriv_z(W, hW, Aux, M, pml, istep, 2, Flag[2], nk - PML_ND, DH, freenode);
    */
#endif
  T = hW; hW = W; W = T;
#ifdef usePML
  if(isx1) {
    T = Aux->hWx1; Aux->hWx1 = Aux->Wx1; Aux->Wx1 = T;
  }
  if(isx2) {
    T = Aux->hWx2; Aux->hWx2 = Aux->Wx2; Aux->Wx2 = T;
  }
  if(isy1) {
    T = Aux->hWy1; Aux->hWy1 = Aux->Wy1; Aux->Wy1 = T;
  }
  if(isy2) {
    T = Aux->hWy2; Aux->hWy2 = Aux->Wy2; Aux->Wy2 = T;
  }
  if(isz1) {
    T = Aux->hWz1; Aux->hWz1 = Aux->Wz1; Aux->Wz1 = T;
  }
  if(isz2) {
    T = Aux->hWz2; Aux->hWz2 = Aux->Wz2; Aux->Wz2 = T;
  }
#endif

  exchange_Wave(W);
  Calculate(W,tW,hW,mW,M,matVx2Vz,matVy2Vz,FlagMinus,RK4b[3],RK4b[3],3, 0, pml, Aux, isx1, isx2, isy1, isy2, isz1, isz2);
#ifdef usePML
  /*
  if(isx1) 
    abs_deriv_x(W, hW, Aux, M, pml, istep, 2, Flag[0], 0, DH, freenode, matVx2Vz);
  if(isx2)
    abs_deriv_x(W, hW, Aux, M, pml, istep, 2, Flag[0], ni - PML_ND, DH, freenode, matVx2Vz);
  if(isy1)
    abs_deriv_y(W, hW, Aux, M, pml, istep, 2, Flag[1], 0, DH, freenode, matVy2Vz);
  if(isy2)
    abs_deriv_y(W, hW, Aux, M, pml, istep, 2, Flag[1], nj - PML_ND, DH, freenode, matVy2Vz);
  if(isz1)
    abs_deriv_z(W, hW, Aux, M, pml, istep, 2, Flag[2], 0, DH, freenode);
  if(isz2)
    abs_deriv_z(W, hW, Aux, M, pml, istep, 2, Flag[2], nk - PML_ND, DH, freenode);
    */
#endif
  T = hW; hW = W; W = T;
#ifdef usePML
  if(isx1) {
    T = Aux->hWx1; Aux->hWx1 = Aux->Wx1; Aux->Wx1 = T;
  }
  if(isx2) {
    T = Aux->hWx2; Aux->hWx2 = Aux->Wx2; Aux->Wx2 = T;
  }
  if(isy1) {
    T = Aux->hWy1; Aux->hWy1 = Aux->Wy1; Aux->Wy1 = T;
  }
  if(isy2) {
    T = Aux->hWy2; Aux->hWy2 = Aux->Wy2; Aux->Wy2 = T;
  }
  if(isz1) {
    T = Aux->hWz1; Aux->hWz1 = Aux->Wz1; Aux->Wz1 = T;
  }
  if(isz2) {
    T = Aux->hWz2; Aux->hWz2 = Aux->Wz2; Aux->Wz2 = T;
  }
#endif

  free(hW);
  free(mW);
  free(tW);

  return 0 ;
}
#endif
