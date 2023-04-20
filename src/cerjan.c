#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include "common.h"
//#include "macdrp.h"
#include "params.h"

int abs_exp(Grid3D damp) {

  //int ND = 16 ;
  int ND = 10 ;
  int i, j, k;

  Grid3D tmp;
  tmp = Alloc3D(nx, ny, lnz);

  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      for (k = 0; k < nz; k++){
        damp[i][j][k] = 1.0;
      }

  if(!absnode) return 0;

  // x1
  if(neigxid[0] == MPI_PROC_NULL){
    for (i = 0; i < ND; i++)
      for (j = nj1; j < nj2; j++)
        for (k = nk1; k < nk2; k++){
          damp[ni1+i][j][k] = exp( -pow( 0.015 * (ND - i), 2) );
        }
  }
  // x2
  if(neigxid[1] == MPI_PROC_NULL){
    for (i = 0; i < ND; i++)
      for (j = nj1; j < nj2; j++)
        for (k = nk1; k < nk2; k++){
          damp[ni2-i-1][j][k] = exp( -pow( 0.015 * (ND - i), 2) );
        }
  }

  memcpy(&tmp[0][0][0], &damp[0][0][0], nx*ny*lnz*sizeof(float));

  // y1
  if(neigyid[0] == MPI_PROC_NULL){
    for (j = 0; j < ND; j++)
      for (i = ni1; i < ni2; i++)
        for (k = nk1; k < nk2; k++){
          damp[i][nj1+j][k] = exp( -pow( 0.015 * (ND - j), 2) );
          damp[i][nj1+j][k] = _min( damp[i][nj1+j][k],
                                    tmp[i][nj1+j][k] );
        }
  }

  // y2
  if(neigyid[1] == MPI_PROC_NULL){
    for (j = 0; j < ND; j++)
      for (i = ni1; i < ni2; i++)
        for (k = nk1; k < nk2; k++){
          damp[i][nj2-j-1][k] = exp( -pow( 0.015 * (ND - j), 2) );
          damp[i][nj2-j-1][k] = _min( damp[i][nj2-j-1][k],
                                      tmp[i][nj2-j-1][k] );
        }
  }

  memcpy(&tmp[0][0][0], &damp[0][0][0], nx*ny*lnz*sizeof(float));

  // z1
  if(neigzid[0] == MPI_PROC_NULL){
    for (k = 0; k < ND; k++)
      for (i = ni1; i < ni2; i++)
        for (j = nj1; j < nj2; j++){
          damp[i][j][nk1+k] = exp( -pow( 0.015 * (ND - k), 2) );
          damp[i][j][nk1+k] = _min( damp[i][j][nk1+k],
                                    tmp[i][j][nk1+k] );
        }
  }
  // z2
  //if(neigzid[1] == MPI_PROC_NULL){
  //  for (k = 0; k < ND; k++)
  //    for (i = ni1; i < ni2; i++)
  //      for (j = nj1; j < nj2; j++){
  //        damp[i][j][nk2-k-1] = exp( -pow( 0.015 * (ND - k), 2) );
  //        damp[i][j][nk2-k-1] = min( damp[i][j][nk2-k-1],
  //                                    tmp[i][j][nk2-k-1] );
  //      }
  //}

  Delloc3D(tmp);
  return 0;
}

void Apply_abs(float *W, Grid3D damp){

  int i, j, k, pos;
  for (i = ni1; i < ni2; i++)
    for (j = nj1; j < nj2; j++)
      for (k = nk1; k < nk2; k++){
    	pos = (i * ny * lnz + j * lnz + k) * WSIZE;
        W[pos + 0] *= damp[i][j][k];
        W[pos + 1] *= damp[i][j][k];
        W[pos + 2] *= damp[i][j][k];
        W[pos + 3] *= damp[i][j][k];
        W[pos + 4] *= damp[i][j][k];
        W[pos + 5] *= damp[i][j][k];
        W[pos + 6] *= damp[i][j][k];
        W[pos + 7] *= damp[i][j][k];
        W[pos + 8] *= damp[i][j][k];
      }

  return ;
}

void Apply_attenuation(float *W, float *Qs){

  int i, j, k, q, pos;
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      for (k = 0; k < lnz; k++){
    	  pos = i * ny * lnz + j * lnz + k;
        for (q = 0; q < WSIZE; q++)
        W[ pos*WSIZE + q] *= exp( (- PI * Qsf0 * DT) / Qs[pos]); 
      }
  return ;
}

void cal_attenu(float *M, float *Qs){
  
  int i, j, k, pos;
  for(i = 0; i < nx; i++)
    for(j = 0; j < ny; j++)
      for(k = 0; k < lnz; k++){
    	  pos = i * ny * lnz + j * lnz + k;
        Qs[pos] = 1.4 * M[pos*MSIZE+11] / M[pos*MSIZE+12];// vs
        //Qs[pos] = 1.4*( M[pos*MSIZE+10] + 2*M[pos*MSIZE+11]) / M[pos*MSIZE+12];// vp
      }

  return ;  
}


