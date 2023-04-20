//float *Ax, *Bx, *Dx;
//float *Ay, *By, *Dy;
//float *Az, *Bz, *Dz;
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "common.h"
#include "macdrp.h"
#include "params.h"
#include "pml.h"

//#define DEBUG
#define CONSPD 2.0f // power for d
#define CONSPB 2.0f // power for beta
#define CONSPA 1.0f // power for alpha
#define AbsVzero

inline float cal_pml_R(int N){
  return (float) (pow(10, -( (log10((double)N)-1.0)/log10(2.0) + 3.0)));
}

inline float cal_pml_dmax(float L, float Vp, float Rpp){
  return (float) (-Vp / (2.0 * L) * log(Rpp) * (CONSPD + 1.0));
}

inline float cal_pml_amax(float fc){return PI*fc;}

//inline float cal_pml_bmax(){return 3.0;}

inline float cal_pml_d(float x, float L, float dmax){
  return (x<0) ? 0.0f : (float) (dmax * pow(x/L, CONSPD));
}

inline float cal_pml_a(float x, float L, float amax){
  return (x<0) ? 0.0f : (float) (amax * (1.0 - pow(x/L, CONSPA)));
}

inline float cal_pml_b(float x, float L, float bmax){
  return (x<0) ? 1.0f : (float) (1.0 + (bmax-1.0) * pow(x/L, CONSPB));
}

void abs_init(float *pml){
  int i = 0;
  float pml_fc = 1.6;
  float Vs = 5716.0;
  float Rpp = cal_pml_R(PML_ND);

  float amax = cal_pml_amax(pml_fc);
  float bmax = 3.0;
  float L0 = DH*(PML_ND-1);
  float dmax = cal_pml_dmax(L0, Vs, Rpp);

  float *p0 = pml;
  for (i = 0; i < ni; i++){
    p0[ni * 0 + i] = 0.0f;
    p0[ni * 1 + i] = 1.0f;
    p0[ni * 2 + i] = 0.0f;
  }

  float *p1 = pml + ni * 3;
  for (i = 0; i < nj; i++){
    p1[nj * 0 + i] = 0.0f;
    p1[nj * 1 + i] = 1.0f;
    p1[nj * 2 + i] = 0.0f;
  }

  float *p2 = pml + ni * 3 + nj * 3;
  for (i = 0; i < nk; i++){
    p2[nk * 0 + i] = 0.0f;
    p2[nk * 1 + i] = 1.0f;
    p2[nk * 2 + i] = 0.0f;
  }

  float *Lx = (float *) malloc(sizeof(float)*ni);
  float *Ly = (float *) malloc(sizeof(float)*nj);
  float *Lz = (float *) malloc(sizeof(float)*nk);

  for (i = 0; i < ni; i++) Lx[i] = -1.0f;
  for (i = 0; i < nj; i++) Ly[i] = -1.0f;
  for (i = 0; i < nk; i++) Lz[i] = -1.0f;
#ifdef usePML
  if(isx1){
    for (i = 0; i < PML_ND; i++)
    Lx[i+0] = (PML_ND-1-i) * DH;
  }
  if(isx2){
    for (i = 0; i < PML_ND; i++)
    Lx[ni-PML_ND+i-0] = i * DH;
  }
  if(isy1){
    for (i = 0; i < PML_ND; i++)
    Ly[i+0] = (PML_ND-1-i) * DH;
  }
  // BUG
  if(isy2){
    for (i = 0; i < PML_ND; i++)
    Ly[nj-PML_ND+i-0] = i * DH;
  }
  if(isz1){
    for (i = 0; i < PML_ND; i++)
    Lz[i+0] = (PML_ND-1-i) * DH;
  }
  if(isz2){
    for (i = 0; i < PML_ND; i++)
    Lz[nk-PML_ND+i-0] = i * DH;
  }
#endif

  for (i = 0; i < ni; i++){
    p0[ni * 0 + i] = cal_pml_a(Lx[i], L0, amax);
    p0[ni * 1 + i] = cal_pml_b(Lx[i], L0, bmax);
    p0[ni * 2 + i] = cal_pml_d(Lx[i], L0, dmax);
  }
  for (i = 0; i < nj; i++){
    p1[nj * 0 + i] = cal_pml_a(Ly[i], L0, amax);
    p1[nj * 1 + i] = cal_pml_b(Ly[i], L0, bmax);
    p1[nj * 2 + i] = cal_pml_d(Ly[i], L0, dmax);
  }
  for (i = 0; i < nk; i++){
    p2[nk * 0 + i] = cal_pml_a(Lz[i], L0, amax);
    p2[nk * 1 + i] = cal_pml_b(Lz[i], L0, bmax);
    p2[nk * 2 + i] = cal_pml_d(Lz[i], L0, dmax);
  }

  // convert d_x to d_x/beta_x since only d_x/beta_x needed
  for (i = 0; i < ni; i++) p0[ni * 2 + i] /= p0[ni + i];
  for (i = 0; i < nj; i++) p1[nj * 2 + i] /= p1[nj + i];
  for (i = 0; i < nk; i++) p2[nk * 2 + i] /= p2[nk + i];
  //// reset to normal: A = 0, B = 1, D = 0
  //for (i = 0; i < ni; i++) {
  //  p0[ni * 0 + i] = 0;
  //  p0[ni * 1 + i] = 1;
  //  p0[ni * 2 + i] = 0;
  //}
  //for (i = 0; i < nj; i++) {
  //  p1[nj * 0 + i] = 0;
  //  p1[nj * 1 + i] = 1;
  //  p1[nj * 2 + i] = 0;
  //}
  //for (i = 0; i < nk; i++) {
  //  p2[nk * 0 + i] = 0;
  //  p2[nk * 1 + i] = 1;
  //  p2[nk * 2 + i] = 0;
  //}

#ifdef DEBUG_XX
  FILE *fp;
  char fnm[500];
  sprintf(fnm, "output/ABD%03d%03d%03d.txt", thisid[0], thisid[1], thisid[2]);
  //sprintf(fnm, "output/ABD.txt");
  fp = fopen(fnm, "w");
  for (i = 0; i < ni; i++)
  fprintf(fp, "ABDx[%05d] = %f %f %f\n", i, p0[i], p0[ni + i], p0[ni * 2 + i]);
  for (i = 0; i < nj; i++)
  fprintf(fp, "ABDy[%05d] = %f %f %f\n", i, p1[i], p1[nj + i], p1[nj * 2 + i]);
  for (i = 0; i < nk; i++)
  fprintf(fp, "ABDz[%05d] = %f %f %f\n", i, p2[i], p2[nk + i], p2[nk * 2 + i]);
  fclose(fp);
#endif


  free(Lx);
  free(Ly);
  free(Lz);
  return;
}

