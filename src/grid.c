/*
********************************************************************************
* Grid3D.c                                                                     *
* programming in C language                                                    *
* 3D data structure                                                            *
********************************************************************************
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "common.h"
#include "params.h"
#include "macdrp.h"
#include "pml.h"

// float * Alloc_wave() {
// 	int i;
// 	float *W = (float *) malloc(sizeof(float) * nx * ny * lnz * WSIZE);
// 	for (i = 0; i < nx * ny * lnz * WSIZE; i ++)
// 		W[i] = 0.0f;

// 	return W;
// }

float * Alloc_wave_8() {
	int i;
	float *W_8 = (float *) malloc(sizeof(float) * nx * ny * lnz * WSIZE_V);
	for (i = 0; i < nx * ny * lnz * WSIZE_V; i ++)
		W_8[i] = 0.0f;

	return W_8;
}

float * Alloc_wave_1() {
	int i;
	float *W_1 = (float *) malloc(sizeof(float) * nx * ny * lnz);
	for (i = 0; i < nx * ny * lnz; i ++)
		W_1[i] = 0.0f;

	return W_1;
}


float * Alloc_metric() {
	int i;
	float *M = (float *) malloc(sizeof(float) * nx * ny * lnz * MSIZE);
	for (i = 0; i < nx * ny * lnz * MSIZE; i ++)
		M[i] = 0.0f;

	return M;
}

float *Alloc_media() {
  int i;
  float *D = (float *) malloc(sizeof(float) * nx * ny * lnz * DSIZE);
  for (i = 0; i < nx * ny * lnz * DSIZE; i ++) D[i] = 0.0f;
  return D;
}

float *Alloc_coord() {
  int i;
  float *C = (float *) malloc(sizeof(float) * nx * ny * lnz * CSIZE);
  for (i = 0; i < nx * ny * lnz * CSIZE; i ++) C[i] = 0.0f;
  return C;
}

float *Alloc_attenu() {
  int i;
  float *Qs = (float *) malloc(sizeof(float) * nx * ny * lnz);
  for (i = 0; i < nx * ny * lnz; i ++) Qs[i] = 0.0f;
  return Qs;
}

float *Alloc_pml() {
  float *p = NULL; 
  p = (float*)malloc(sizeof(float)*(ni + nj + nk) * 3);
  return p;
}

void Alloc_aux(struct aux *Aux, int isx1, int isx2, int isy1, int isy2, int isz1, int isz2, int this_rank) {
  int i = 0;
  if(isx1) {
    size_t ibytes = PML_ND*ny*lnz*WSIZE*sizeof(float);
    Aux-> Wx1 = (float *) malloc(ibytes);
    Aux->hWx1 = (float *) malloc(ibytes);
    Aux->mWx1 = (float *) malloc(ibytes);
    Aux->tWx1 = (float *) malloc(ibytes);
    memset(Aux->Wx1, 0, ibytes);
  }
  if(isx2){
    size_t ibytes = PML_ND*ny*lnz*WSIZE*sizeof(float);
    Aux-> Wx2 = (float *) malloc(ibytes);
    Aux->hWx2 = (float *) malloc(ibytes);
    Aux->mWx2 = (float *) malloc(ibytes);
    Aux->tWx2 = (float *) malloc(ibytes);
    memset(Aux->Wx2, 0, ibytes);
  }
  if(isy1) {
    size_t ibytes = PML_ND*nx*lnz*WSIZE*sizeof(float);
    Aux-> Wy1 = (float *) malloc(ibytes);
    Aux->hWy1 = (float *) malloc(ibytes);
    Aux->mWy1 = (float *) malloc(ibytes);
    Aux->tWy1 = (float *) malloc(ibytes);
    memset(Aux->Wy1, 0, ibytes);
  }
  if(isy2){
    size_t ibytes = PML_ND*nx*lnz*WSIZE*sizeof(float);
    Aux-> Wy2 = (float *) malloc(ibytes);
    Aux->hWy2 = (float *) malloc(ibytes);
    Aux->mWy2 = (float *) malloc(ibytes);
    Aux->tWy2 = (float *) malloc(ibytes);
    memset(Aux->Wy2, 0, ibytes);
  }
  if(isz1) {
    size_t ibytes = PML_ND*nx*ny*WSIZE*sizeof(float);
    Aux-> Wz1 = (float *) malloc(ibytes);
    Aux->hWz1 = (float *) malloc(ibytes);
    Aux->mWz1 = (float *) malloc(ibytes);
    Aux->tWz1 = (float *) malloc(ibytes);
    memset(Aux->Wz1, 0, ibytes);
  }
  if(isz2){
    size_t ibytes = PML_ND*nx*ny*WSIZE*sizeof(float);
    Aux-> Wz2 = (float *) malloc(ibytes);
    Aux->hWz2 = (float *) malloc(ibytes);
    Aux->mWz2 = (float *) malloc(ibytes);
    Aux->tWz2 = (float *) malloc(ibytes);
    memset(Aux->Wz2, 0, ibytes);
  }
}

Grid3D Alloc3D(int nx, int ny, int nz) {
  int i, j, k;
  Grid3D U = (Grid3D)malloc(sizeof(float **) * nx + sizeof(float *) * nx * ny + sizeof(float) * nx * ny * nz);

  if (!U) {
    printf("Cannot allocate 3D float array\n");
    exit(-1);
  }
  for (i = 0; i < nx; i++) {
    U[i] = ((float **)U) + nx + i * ny;
  }

  float *Ustart = (float *)(U[nx - 1] + ny);
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++) U[i][j] = Ustart + i * ny * nz + j * nz;

  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      for (k = 0; k < nz; k++) U[i][j][k] = 0.0f;

  return U;
}

Grid2D Alloc2D(int nx, int ny){
  int i, j;
  Grid2D U = malloc(sizeof(float *) * nx);
  U[0] = malloc(sizeof(float) * nx * ny);

  if(!U) {
    printf("Cannot allocate 2D float array\n");
    exit(-1);
  }

  for (i = 1; i < nx; i++) 
    U[i] = U[i-1] + ny ;

  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      U[i][j] = 0.0f;

  return U;
}

Grid1D Alloc1D(int nx) {
  int i;
  Grid1D U = (Grid1D)malloc(sizeof(float) * nx);

  if (!U) {
    printf("Cannot allocate 1D float array\n");
    exit(-1);
  }

  for (i = 0; i < nx; i++) U[i] = 0.0f;

  return U;
}

int* Alloc1P(int nx) {
  int i;
  int* U = (int *)malloc(sizeof(int) * nx);

  if (!U) {
    printf("Cannot allocate 2D integer array\n");
    exit(-1);
  }

  for (i = 0; i < nx; i++) U[i] = 0;

  return U;
}

unsigned int *Alloc1PU(unsigned int nx) {
  int i;
  unsigned int *U = (unsigned int *)malloc(sizeof(unsigned int) * nx);

  if (!U) {
    printf("Cannot allocate 2D integer array\n");
    exit(-1);
  }

  for (i = 0; i < nx; i++) U[i] = 0;

  return U;
}

void Delloc3D(Grid3D U) {
  if (U) {
    free((void*)U);
    U = NULL;
  }

  return;
}

void Delloc2D(Grid2D U) {
  if (U) {
    free((void*)U[0]);
    free((void*)U);
  }

  return;
}

void Delloc1D(Grid1D U) {
  if (U) {
    free(U);
    U = NULL;
  }

  return;
}
/*
void Alloc_Coord(struct Coord *C) {
  C->x = Alloc3D( nx, ny, nz);
  C->y = Alloc3D( nx, ny, nz);
  C->z = Alloc3D( nx, ny, nz);
  return ;
}

void coord_extend(Grid3D U) {
  int i, j, k, n;
  int i1 = ni1; int i2 = ni2 - 1;
  int j1 = nj1; int j2 = nj2 - 1;
  int k1 = nk1; int k2 = nk2 - 1;
  // extend x direction
  for (j = ny1; j < ny2; j++)
    for (k = nz1; k < nz2; k++) 
      for (n = 1; n <= 3; n++) {
        U[i1-n][j][k] = 2.0 * U[i1][j][k] - U[i1+n][j][k];
        U[i2+n][j][k] = 2.0 * U[i2][j][k] - U[i2-n][j][k];
      }
  
  // extend y direction
  for (i = nx1; i < nx2; i++)
    for (k = nz1; k < nz2; k++) 
      for (n = 1; n <= 3; n++) {
        U[i][j1-n][k] = 2.0 * U[i][j1][k] - U[i][j1+n][k];
        U[i][j2+n][k] = 2.0 * U[i][j2][k] - U[i][j2-n][k];
      }
  
  // extend z direction
  for (i = nx1; i < nx2; i++) 
    for (j = ny1; j < ny2; j++) 
      for (n = 1; n <= 3; n++) {
        U[i][j][k1-n] = 2.0 * U[i][j][k1] - U[i][j][k1+n];
        U[i][j][k2+n] = 2.0 * U[i][j][k2] - U[i][j][k2-n];
      }

  return;
}


void extend_Symm(Grid3D U) {
  int i, j, k, n;
  int i1 = ni1; int i2 = ni2 - 1;
  int j1 = ni1; int j2 = nj2 - 1;
  int k1 = ni1; int k2 = nk2 - 1;
  // extend x direction
  for (j = 0; j < ny; j++) 
    for (k = 0; k < nz; k++) 
      for (n = 1; n <= 3; n++) {
        U[i1-n][j][k] = U[i1+n][j][k];
        U[i2+n][j][k] = U[i2-n][j][k];
      }
  
  // extend y direction
  for (i = 0; i < nx; i++) 
    for (k = 0; k < nz; k++) 
      for (n = 1; n <= 3; n++) {
        U[i][j1-n][k] = U[i][j1+n][k];
        U[i][j2+n][k] = U[i][j2-n][k];
      }

  // extend z direction
  for (i = 0; i < nx; i++) 
    for (j = 0; j < ny; j++) 
      for (n = 1; n <= 3; n++) {
        U[i][j][k1-n] = U[i][j][k1+n];
        U[i][j][k2+n] = U[i][j][k2-n];
      }

  return;
}

*/

void extend_coord(float *C, int SIZE){
  int i, j, k, n, l;
  int i1 = ni1; int i2 = ni2 - 1;
  int j1 = nj1; int j2 = nj2 - 1;
  int k1 = nk1; int k2 = nk2 - 1;
  int pos_s_1, pos_s_2;
  int pos_0_1, pos_0_2;
  int pos_d_1, pos_d_2;
  // extend x direction
  for (j = 0; j < ny; j++)
    for (k = 0; k < nz; k++)
      for (n = 1; n <= 3; n++) {
        pos_s_1 = ((i1-n) * ny * lnz + j * lnz + k) * SIZE;
        pos_0_1 = ((i1  ) * ny * lnz + j * lnz + k) * SIZE;
        pos_d_1 = ((i1+n) * ny * lnz + j * lnz + k) * SIZE;
        pos_s_2 = ((i2+n) * ny * lnz + j * lnz + k) * SIZE;
        pos_0_2 = ((i2  ) * ny * lnz + j * lnz + k) * SIZE;
        pos_d_2 = ((i2-n) * ny * lnz + j * lnz + k) * SIZE;
        for (l = 0; l < SIZE; l++) {
          C[pos_s_1 + l] = 2.0f * C[pos_0_1 + l] - C[pos_d_1 + l];
          C[pos_s_2 + l] = 2.0f * C[pos_0_2 + l] - C[pos_d_2 + l];
        }
      }
  // extend y direction
  for (i = 0; i < nx; i++)
    for (k = 0; k < nz; k++)
      for (n = 1; n <= 3; n++) {
        pos_s_1 = (i * ny * lnz + (j1-n) * lnz + k) * SIZE;
        pos_0_1 = (i * ny * lnz + (j1  ) * lnz + k) * SIZE;
        pos_d_1 = (i * ny * lnz + (j1+n) * lnz + k) * SIZE;
        pos_s_2 = (i * ny * lnz + (j2+n) * lnz + k) * SIZE;
        pos_0_2 = (i * ny * lnz + (j2  ) * lnz + k) * SIZE;
        pos_d_2 = (i * ny * lnz + (j2-n) * lnz + k) * SIZE;
        for (l = 0; l < SIZE; l++) {
          C[pos_s_1 + l] = 2.0f * C[pos_0_1 + l] - C[pos_d_1 + l];
          C[pos_s_2 + l] = 2.0f * C[pos_0_2 + l] - C[pos_d_2 + l];
        }
      }
  // extend z direction
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      for (n = 1; n <= 3; n++) {
        pos_s_1 = (i * ny * lnz + j * lnz + (k1-n)) * SIZE;
        pos_0_1 = (i * ny * lnz + j * lnz + (k1  )) * SIZE;
        pos_d_1 = (i * ny * lnz + j * lnz + (k1+n)) * SIZE;
        pos_s_2 = (i * ny * lnz + j * lnz + (k2+n)) * SIZE;
        pos_0_2 = (i * ny * lnz + j * lnz + (k2  )) * SIZE;
        pos_d_2 = (i * ny * lnz + j * lnz + (k2-n)) * SIZE;
        for (l = 0; l < SIZE; l++) {
          C[pos_s_1 + l] = 2.0f * C[pos_0_1 + l] - C[pos_d_1 + l];
          C[pos_s_2 + l] = 2.0f * C[pos_0_2 + l] - C[pos_d_2 + l];
        }
      }
  return;
}


void extend_Symm_array(float *W, int SIZE) {
  int i, j, k, n, l, pos_s_1, pos_d_1, pos_s_2, pos_d_2;
  int i1 = ni1; int i2 = ni2 - 1;
  int j1 = nj1; int j2 = nj2 - 1;
  int k1 = nk1; int k2 = nk2 - 1;
  // extend x direction
  for (j = 0; j < ny; j++)
    for (k = 0; k < nz; k++)
      for (n = 1; n <= 3; n++) {
    	int pos_s_1 = ((i1-n) * ny * lnz + j * lnz + k) * SIZE;
    	int pos_d_1 = ((i1+n) * ny * lnz + j * lnz + k) * SIZE;
    	int pos_s_2 = ((i2+n) * ny * lnz + j * lnz + k) * SIZE;
    	int pos_d_2 = ((i2-n) * ny * lnz + j * lnz + k) * SIZE;
    	for (l = 0; l < SIZE; l++) {
    		W[pos_s_1 + l] = W[pos_d_1 + l];
			W[pos_s_2 + l] = W[pos_d_2 + l];
    	}
      }

  // extend y direction
  for (i = 0; i < nx; i++)
    for (k = 0; k < nz; k++)
      for (n = 1; n <= 3; n++) {
      	int pos_s_1 = (i * ny * lnz + (j1-n) * lnz + k) * SIZE;
      	int pos_d_1 = (i * ny * lnz + (j1+n) * lnz + k) * SIZE;
      	int pos_s_2 = (i * ny * lnz + (j2+n) * lnz + k) * SIZE;
      	int pos_d_2 = (i * ny * lnz + (j2-n) * lnz + k) * SIZE;
    	for (l = 0; l < SIZE; l++) {
    		W[pos_s_1 + l] = W[pos_d_1 + l];
			W[pos_s_2 + l] = W[pos_d_2 + l];
    	}
      }

  // extend z direction
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      for (n = 1; n <= 3; n++) {
        int pos_s_1 = (i * ny * lnz + j * lnz + (k1-n)) * SIZE;
        int pos_d_1 = (i * ny * lnz + j * lnz + (k1+n)) * SIZE;
        int pos_s_2 = (i * ny * lnz + j * lnz + (k2+n)) * SIZE;
        int pos_d_2 = (i * ny * lnz + j * lnz + (k2-n)) * SIZE;
    	for (l = 0; l < SIZE; l++) {
    		W[pos_s_1 + l] = W[pos_d_1 + l];
			W[pos_s_2 + l] = W[pos_d_2 + l];
    	}
      }

  return;
}


int Malloc( void ** mem, long long  size  )
{
	*mem = malloc( size );
	if ( *mem == NULL )
	{
		printf( "can not malloc, Error: %s:%d\n", __FILE__, __LINE__ );
	}
	return 0;

}

/*
int Gauss3D(Grid3D U, struct Coord *C) {
    int i, j, k;
    float a = 4 * DH; 
    a = a*a;
    float r;
    float xsrc = 0.0;
    float ysrc = 0.0;
    float zsrc = -4.0e3;
    
    for (i = 0; i < nx; i++)
      for (j= 0; j < ny; j++)
        for (k= 0; k < nz; k++){
          r = pow(C->x[i][j][k] - xsrc, 2) 
            + pow(C->y[i][j][k] - ysrc, 2) 
            + pow(C->z[i][j][k] - zsrc, 2) ;

          U[i][j][k] = exp( -r/a );
        }

   return 0;
}

*/
