/*
********************************************************************************
* modules used by seismic wave mpi                                             *
* programming in C language                                                    *
********************************************************************************
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "common.h"
#include "params.h"

int swmpi_datatype(void) {
  // LenFDS = 1; LenFDL = 3;
  // count, blocklen, stride
  // x plane
  MPI_Type_vector(LenFDS, ny*lnz, ny*lnz, MPI_FLOAT, &DTypeXS);
  MPI_Type_vector(LenFDL, ny*lnz, ny*lnz, MPI_FLOAT, &DTypeXL);
  // y plane
  MPI_Type_vector(nx, lnz*LenFDS, ny*lnz, MPI_FLOAT, &DTypeYS);
  MPI_Type_vector(nx, lnz*LenFDL, ny*lnz, MPI_FLOAT, &DTypeYL);
  // z plane
  MPI_Type_vector(nx*ny, LenFDS, lnz, MPI_FLOAT, &DTypeZS);
  MPI_Type_vector(nx*ny, LenFDL, lnz, MPI_FLOAT, &DTypeZL);

  MPI_Type_commit(&DTypeXS);
  MPI_Type_commit(&DTypeXL);
  MPI_Type_commit(&DTypeYS);
  MPI_Type_commit(&DTypeYL);
  MPI_Type_commit(&DTypeZS);
  MPI_Type_commit(&DTypeZL);

  // x plane
  MPI_Type_vector(LenFDS, ny*lnz*WSIZE, ny*lnz*WSIZE, MPI_FLOAT, &WDTypeXS);
  MPI_Type_vector(LenFDL, ny*lnz*WSIZE, ny*lnz*WSIZE, MPI_FLOAT, &WDTypeXL);
  // y plane
  MPI_Type_vector(nx, lnz*LenFDS*WSIZE, ny*lnz*WSIZE, MPI_FLOAT, &WDTypeYS);
  MPI_Type_vector(nx, lnz*LenFDL*WSIZE, ny*lnz*WSIZE, MPI_FLOAT, &WDTypeYL);
  // z plane
  MPI_Type_vector(nx*ny, LenFDS*WSIZE, lnz*WSIZE, MPI_FLOAT, &WDTypeZS);
  MPI_Type_vector(nx*ny, LenFDL*WSIZE, lnz*WSIZE, MPI_FLOAT, &WDTypeZL);

  MPI_Type_commit(&WDTypeXS);
  MPI_Type_commit(&WDTypeXL);
  MPI_Type_commit(&WDTypeYS);
  MPI_Type_commit(&WDTypeYL);
  MPI_Type_commit(&WDTypeZS);
  MPI_Type_commit(&WDTypeZL);

  // x plane
  MPI_Type_vector(LenFDS, ny*lnz*MSIZE, ny*lnz*MSIZE, MPI_FLOAT, &MDTypeXS);
  MPI_Type_vector(LenFDL, ny*lnz*MSIZE, ny*lnz*MSIZE, MPI_FLOAT, &MDTypeXL);
  // y plane
  MPI_Type_vector(nx, lnz*LenFDS*MSIZE, ny*lnz*MSIZE, MPI_FLOAT, &MDTypeYS);
  MPI_Type_vector(nx, lnz*LenFDL*MSIZE, ny*lnz*MSIZE, MPI_FLOAT, &MDTypeYL);
  // z plane
  MPI_Type_vector(nx*ny, LenFDS*MSIZE, lnz*MSIZE, MPI_FLOAT, &MDTypeZS);
  MPI_Type_vector(nx*ny, LenFDL*MSIZE, lnz*MSIZE, MPI_FLOAT, &MDTypeZL);

  MPI_Type_commit(&MDTypeXS);
  MPI_Type_commit(&MDTypeXL);
  MPI_Type_commit(&MDTypeYS);
  MPI_Type_commit(&MDTypeYL);
  MPI_Type_commit(&MDTypeZS);
  MPI_Type_commit(&MDTypeZL);

  // x plane
  MPI_Type_vector(LenFDS, ny*lnz*CSIZE, ny*lnz*CSIZE, MPI_FLOAT, &CDTypeXS);
  MPI_Type_vector(LenFDL, ny*lnz*CSIZE, ny*lnz*CSIZE, MPI_FLOAT, &CDTypeXL);
  // y plane
  MPI_Type_vector(nx, lnz*LenFDS*CSIZE, ny*lnz*CSIZE, MPI_FLOAT, &CDTypeYS);
  MPI_Type_vector(nx, lnz*LenFDL*CSIZE, ny*lnz*CSIZE, MPI_FLOAT, &CDTypeYL);
  // z plane
  MPI_Type_vector(nx*ny, LenFDS*CSIZE, lnz*CSIZE, MPI_FLOAT, &CDTypeZS);
  MPI_Type_vector(nx*ny, LenFDL*CSIZE, lnz*CSIZE, MPI_FLOAT, &CDTypeZL);

  MPI_Type_commit(&CDTypeXS);
  MPI_Type_commit(&CDTypeXL);
  MPI_Type_commit(&CDTypeYS);
  MPI_Type_commit(&CDTypeYL);
  MPI_Type_commit(&CDTypeZS);
  MPI_Type_commit(&CDTypeZL);

  return 0;

}
  
int swmpi_datatype_free(void) {
  MPI_Type_free(&DTypeXS);
  MPI_Type_free(&DTypeXL);
  MPI_Type_free(&DTypeYS);
  MPI_Type_free(&DTypeYL);
  MPI_Type_free(&DTypeZS);
  MPI_Type_free(&DTypeZL);

  MPI_Type_free(&WDTypeXS);
  MPI_Type_free(&WDTypeXL);
  MPI_Type_free(&WDTypeYS);
  MPI_Type_free(&WDTypeYL);
  MPI_Type_free(&WDTypeZS);
  MPI_Type_free(&WDTypeZL);

  MPI_Type_free(&MDTypeXS);
  MPI_Type_free(&MDTypeXL);
  MPI_Type_free(&MDTypeYS);
  MPI_Type_free(&MDTypeYL);
  MPI_Type_free(&MDTypeZS);
  MPI_Type_free(&MDTypeZL);

  MPI_Type_free(&CDTypeXS);
  MPI_Type_free(&CDTypeXL);
  MPI_Type_free(&CDTypeYS);
  MPI_Type_free(&CDTypeYL);
  MPI_Type_free(&CDTypeZS);
  MPI_Type_free(&CDTypeZL);

  return 0;
}

//int coord_exchange(struct Coord *C) {
int coord_exchange(float *C) {

  MPI_Status stat;
  int pos_s, pos_d;
  pos_s = (ni1 * ny * lnz) * CSIZE; pos_d = (ni2 * ny * lnz) * CSIZE;
  MPI_Sendrecv( &C[pos_s], 1, CDTypeXL, neigxid[0], 311,
                &C[pos_d], 1, CDTypeXL, neigxid[1], 311,
                SWMPI_COMM, &stat );

  pos_s = ((ni2-LenFDL) * ny * lnz) * CSIZE; pos_d = (0 * ny * lnz) * CSIZE;
  MPI_Sendrecv( &C[pos_s], 1, CDTypeXL, neigxid[1], 312,
                &C[pos_d], 1, CDTypeXL, neigxid[0], 312,
                SWMPI_COMM, &stat );

  //exchange y plane data
  pos_s = (nj1 * lnz) * CSIZE; pos_d = (nj2 * lnz) * CSIZE;
  MPI_Sendrecv( &C[pos_s], 1, CDTypeYL, neigyid[0], 321,
                &C[pos_d], 1, CDTypeYL, neigyid[1], 321,
                SWMPI_COMM, &stat );
  pos_s = ((nj2-LenFDL) * lnz) * CSIZE; pos_d = (0 * lnz) * CSIZE;
  MPI_Sendrecv( &C[pos_s], 1, CDTypeYL, neigyid[1], 322,
                &C[pos_d], 1, CDTypeYL, neigyid[0], 322,
                SWMPI_COMM, &stat );

  // exchange z plane data
  pos_s = (nk1) * CSIZE; pos_d = (nk2) * CSIZE;
  MPI_Sendrecv( &C[pos_s], 1, CDTypeZL, neigzid[0], 331,
                &C[pos_d], 1, CDTypeZL, neigzid[1], 331,
                SWMPI_COMM, &stat);
  pos_s = (nk2-LenFDL) * CSIZE; pos_d = (0) * CSIZE;
  MPI_Sendrecv( &C[pos_s], 1, CDTypeZL, neigzid[1], 332,
                &C[pos_d], 1, CDTypeZL, neigzid[0], 332,
                SWMPI_COMM, &stat);

  return 0;
}

int metric_exchange(float *M) {
  MPI_Status stat;

  //exchange x plane data
  int pos_m_s, pos_m_d;
  pos_m_s = (ni1 * ny * lnz) * MSIZE; pos_m_d = (ni2 * ny * lnz) * MSIZE;
  MPI_Sendrecv( &M[pos_m_s], 1, MDTypeXL, neigxid[0], 411,
                &M[pos_m_d], 1, MDTypeXL, neigxid[1], 411,
                SWMPI_COMM, &stat );

  pos_m_s = ((ni2-LenFDL) * ny * lnz) * MSIZE; pos_m_d = (0 * ny * lnz) * MSIZE;
  MPI_Sendrecv( &M[pos_m_s], 1, MDTypeXL, neigxid[1], 412,
                &M[pos_m_d], 1, MDTypeXL, neigxid[0], 412,
                SWMPI_COMM, &stat );

  //exchange y plane data
  pos_m_s = (nj1 * lnz) * MSIZE; pos_m_d = (nj2 * lnz) * MSIZE;
  MPI_Sendrecv( &M[pos_m_s], 1, MDTypeYL, neigyid[0], 421,
                &M[pos_m_d], 1, MDTypeYL, neigyid[1], 421,
                SWMPI_COMM, &stat );
  pos_m_s = ((nj2-LenFDL) * lnz) * MSIZE; pos_m_d = (0 * lnz) * MSIZE;
  MPI_Sendrecv( &M[pos_m_s], 1, MDTypeYL, neigyid[1], 422,
                &M[pos_m_d], 1, MDTypeYL, neigyid[0], 422,
                SWMPI_COMM, &stat );

  // exchange z plane data
  pos_m_s = (nk1) * MSIZE; pos_m_d = (nk2) * MSIZE;
  MPI_Sendrecv( &M[pos_m_s], 1, MDTypeZL, neigzid[0], 431,
                &M[pos_m_d], 1, MDTypeZL, neigzid[1], 431,
                SWMPI_COMM, &stat);
  pos_m_s = (nk2-LenFDL) * MSIZE; pos_m_d = (0) * MSIZE;
  MPI_Sendrecv( &M[pos_m_s], 1, MDTypeZL, neigzid[1], 432,
                &M[pos_m_d], 1, MDTypeZL, neigzid[0], 432,
                SWMPI_COMM, &stat);
  return 0;
}

void pack(float *w, float *wp, int count, int blocklen, int stride) {
  int i, j;
  // for(i = 0 ; i < count ; i ++) { 
  //   for(j = 0 ; j < blocklen ; j ++) {
  //     wp[i * blocklen + j] = w[i * stride + j];
  //   }
  // }

  for(i = 0 ; i < count ; i ++) { 
    memcpy(&wp[i * blocklen], &w[i * stride], sizeof(float)*blocklen);
  }
}

void unpack(float *wp, float *w, int count, int blocklen, int stride) {
  int i, j;
  // for(i = 0 ; i < count ; i ++) { 
  //   for(j = 0 ; j < blocklen ; j ++) {
  //     w[i * stride + j] = wp[i * blocklen + j];
  //   }
  // }
  for(i = 0 ; i < count ; i ++) { 
    memcpy(&w[i * stride], &wp[i * blocklen], sizeof(float)*blocklen);
  }
}

// int exchange_Wave2(float *W) {
//   MPI_Status stat;
//   //exchange x plane data
//   int pos_w_s, pos_w_d;
//   pos_w_s = (ni1 * ny * lnz) * WSIZE; pos_w_d = (ni2 * ny * lnz) * WSIZE;
//   //MPI_Sendrecv( &W[pos_w_s], 1, WDTypeXL, neigxid[0], 5111,
//   //              &W[pos_w_d], 1, WDTypeXL, neigxid[1], 5111,
//   //              SWMPI_COMM, &stat );
//   pack(&W[pos_w_s], wp_yzs0, LenFDL, ny * lnz * WSIZE, ny * lnz * WSIZE); 
//   MPI_Sendrecv(wp_yzs0, ny * lnz * WSIZE * LenFDL, MPI_FLOAT, neigxid[0], 5111,
//                 wp_yzr0, ny * lnz * WSIZE * LenFDL, MPI_FLOAT, neigxid[1], 5111,
//                 SWMPI_COMM, &stat );
//   unpack(wp_yzr0, &W[pos_w_d], LenFDL, ny * lnz * WSIZE, ny * lnz * WSIZE); 




//   pos_w_s = ((ni2-LenFDL) * ny * lnz) * WSIZE; pos_w_d = (0 * ny * lnz) * WSIZE;
//   //MPI_Sendrecv( &W[pos_w_s], 1, WDTypeXL, neigxid[1], 5112,
//   //              &W[pos_w_d], 1, WDTypeXL, neigxid[0], 5112,
//   //              SWMPI_COMM, &stat );
//   pack(&W[pos_w_s], wp_yzs1, LenFDL, ny * lnz * WSIZE, ny * lnz * WSIZE); 
//   MPI_Sendrecv(wp_yzs1, ny * lnz * WSIZE * LenFDL, MPI_FLOAT, neigxid[1], 5112,
//                 wp_yzr1, ny * lnz * WSIZE * LenFDL, MPI_FLOAT, neigxid[0], 5112,
//                 SWMPI_COMM, &stat );
//   unpack(wp_yzr1, &W[pos_w_d], LenFDL, ny * lnz * WSIZE, ny * lnz * WSIZE); 

//   // exchange y plane data
//   pos_w_s = (nj1 * lnz) * WSIZE; pos_w_d = (nj2 * lnz) * WSIZE;
//   //MPI_Sendrecv( &W[pos_w_s], 1, WDTypeYL, neigyid[0], 5121,
//   //              &W[pos_w_d], 1, WDTypeYL, neigyid[1], 5121,
//   //              SWMPI_COMM, &stat );
//   pack(&W[pos_w_s], wp_xzs0, nx, lnz * WSIZE * LenFDL, ny * lnz * WSIZE); 
//   MPI_Sendrecv(wp_xzs0, nx * lnz * WSIZE * LenFDL, MPI_FLOAT, neigyid[0], 5121,
//                 wp_xzr0, nx * lnz * WSIZE * LenFDL, MPI_FLOAT, neigyid[1], 5121,
//                 SWMPI_COMM, &stat );
//   unpack(wp_xzr0, &W[pos_w_d], nx, lnz * WSIZE * LenFDL, ny * lnz * WSIZE); 

//   pos_w_s = ((nj2-LenFDL) * lnz) * WSIZE; pos_w_d = (0 * lnz) * WSIZE;
//   //MPI_Sendrecv( &W[pos_w_s], 1, WDTypeYL, neigyid[1], 5122,
//   //              &W[pos_w_d], 1, WDTypeYL, neigyid[0], 5122,
//   //              SWMPI_COMM, &stat );
//   pack(&W[pos_w_s], wp_xzs1, nx, lnz * WSIZE * LenFDL, ny * lnz * WSIZE); 
//   MPI_Sendrecv(wp_xzs1, nx * lnz * WSIZE * LenFDL, MPI_FLOAT, neigyid[1], 5122,
//                 wp_xzr1, nx * lnz * WSIZE * LenFDL, MPI_FLOAT, neigyid[0], 5122,
//                 SWMPI_COMM, &stat );
//   unpack(wp_xzr1, &W[pos_w_d], nx, lnz * WSIZE * LenFDL, ny * lnz * WSIZE); 

//   // exchange z plane data
//   pos_w_s = (nk1) * WSIZE; pos_w_d = (nk2) * WSIZE;
//   //MPI_Sendrecv( &W[pos_w_s], 1, WDTypeZL, neigzid[0], 5131,
//   //              &W[pos_w_d], 1, WDTypeZL, neigzid[1], 5131,
//   //              SWMPI_COMM, &stat);
//   pack(&W[pos_w_s], wp_xys0, nx * ny, WSIZE * LenFDL, lnz * WSIZE); 
//   MPI_Sendrecv(wp_xys0, nx * ny * WSIZE * LenFDL, MPI_FLOAT, neigzid[0], 5131,
//                 wp_xyr0, nx * ny * WSIZE * LenFDL, MPI_FLOAT, neigzid[1], 5131,
//                 SWMPI_COMM, &stat );
//   unpack(wp_xyr0, &W[pos_w_d], nx * ny, WSIZE * LenFDL, lnz * WSIZE); 
  
//   pos_w_s = (nk2-LenFDL) * WSIZE; pos_w_d = (0) * WSIZE;
//   //MPI_Sendrecv( &W[pos_w_s], 1, WDTypeZL, neigzid[1], 5132,
//   //              &W[pos_w_d], 1, WDTypeZL, neigzid[0], 5132,
//   //              SWMPI_COMM, &stat);
//   pack(&W[pos_w_s], wp_xys1, nx * ny, WSIZE * LenFDL, lnz * WSIZE); 
//   MPI_Sendrecv(wp_xys1, nx * ny * WSIZE * LenFDL, MPI_FLOAT, neigzid[1], 5132,
//                 wp_xyr1, nx * ny * WSIZE * LenFDL, MPI_FLOAT, neigzid[0], 5132,
//                 SWMPI_COMM, &stat );
//   unpack(wp_xyr1, &W[pos_w_d], nx * ny, WSIZE * LenFDL, lnz * WSIZE); 
//   return 0;
// }

// int exchange_Wave2(float *my_W, float *my_W_Tyz) {
//   MPI_Status stat;
  
//   // splitW(W, my_W, my_W_Tyz);
//   int WSIZE1 = WSIZE - 1;
//   int seg0 = ny * lnz * (WSIZE-1) * LenFDL;
//   int seg1 = nx * lnz * (WSIZE-1) * LenFDL;
//   int seg2 = nx * ny * (WSIZE-1) * LenFDL;
//   //exchange x plane data
//   int pos_w_s, pos_w_d;
//   pos_w_s = (ni1 * ny * lnz) * WSIZE1; pos_w_d = (ni2 * ny * lnz) * WSIZE1;
//   //MPI_Sendrecv( &W[pos_w_s], 1, WDTypeXL, neigxid[0], 5111,
//   //              &W[pos_w_d], 1, WDTypeXL, neigxid[1], 5111,
//   //              SWMPI_COMM, &stat );
//   pack(&my_W[pos_w_s], wp_yzs0, LenFDL, ny * lnz * WSIZE1, ny * lnz * WSIZE1); 
//   MPI_Sendrecv(wp_yzs0, ny * lnz * WSIZE1 * LenFDL, MPI_FLOAT, neigxid[0], 5111,
//                 wp_yzr0, ny * lnz * WSIZE1 * LenFDL, MPI_FLOAT, neigxid[1], 5111,
//                 SWMPI_COMM, &stat );
//   unpack(wp_yzr0, &my_W[pos_w_d], LenFDL, ny * lnz * WSIZE1, ny * lnz * WSIZE1); 

//   //hd
//   pos_w_s = ni1 * ny * lnz; pos_w_d = ni2 * ny * lnz;
//   pack(&my_W_Tyz[pos_w_s], &wp_yzs0[seg0], LenFDL, ny * lnz, ny * lnz); 
//   MPI_Sendrecv(&wp_yzs0[seg0], ny * lnz * LenFDL, MPI_FLOAT, neigxid[0], 5113,
//                 &wp_yzr0[seg0], ny * lnz * LenFDL, MPI_FLOAT, neigxid[1], 5113,
//                 SWMPI_COMM, &stat );
//   unpack(&wp_yzr0[seg0], &my_W_Tyz[pos_w_d], LenFDL, ny * lnz, ny * lnz);  




//   pos_w_s = ((ni2-LenFDL) * ny * lnz) * WSIZE1; pos_w_d = (0 * ny * lnz) * WSIZE1;
//   //MPI_Sendrecv( &W[pos_w_s], 1, WDTypeXL, neigxid[1], 5112,
//   //              &W[pos_w_d], 1, WDTypeXL, neigxid[0], 5112,
//   //              SWMPI_COMM, &stat );
//   pack(&my_W[pos_w_s], wp_yzs1, LenFDL, ny * lnz * WSIZE1, ny * lnz * WSIZE1); 
//   MPI_Sendrecv(wp_yzs1, ny * lnz * WSIZE1 * LenFDL, MPI_FLOAT, neigxid[1], 5112,
//                 wp_yzr1, ny * lnz * WSIZE1 * LenFDL, MPI_FLOAT, neigxid[0], 5112,
//                 SWMPI_COMM, &stat );
//   unpack(wp_yzr1, &my_W[pos_w_d], LenFDL, ny * lnz * WSIZE1, ny * lnz * WSIZE1); 

//   //hd
//   pos_w_s = (ni2-LenFDL) * ny * lnz; pos_w_d = 0 * ny * lnz;
//   pack(&my_W_Tyz[pos_w_s], &wp_yzs1[seg0], LenFDL, ny * lnz, ny * lnz); 
//   MPI_Sendrecv(&wp_yzs1[seg0], ny * lnz * LenFDL, MPI_FLOAT, neigxid[1], 5114,
//                 &wp_yzr1[seg0], ny * lnz * LenFDL, MPI_FLOAT, neigxid[0], 5114,
//                 SWMPI_COMM, &stat );
//   unpack(&wp_yzr1[seg0], &my_W_Tyz[pos_w_d], LenFDL, ny * lnz, ny * lnz); 

//   // mergeW(W, my_W, my_W_Tyz);

//   // exchange y plane data
//   // splitW(W, my_W, my_W_Tyz);
  
//   pos_w_s = (nj1 * lnz) * WSIZE1; pos_w_d = (nj2 * lnz) * WSIZE1;
//   //MPI_Sendrecv( &W[pos_w_s], 1, WDTypeYL, neigyid[0], 5121,
//   //              &W[pos_w_d], 1, WDTypeYL, neigyid[1], 5121,
//   //              SWMPI_COMM, &stat );
//   pack(&my_W[pos_w_s], wp_xzs0, nx, lnz * WSIZE1 * LenFDL, ny * lnz * WSIZE1); 
//   MPI_Sendrecv(wp_xzs0, nx * lnz * WSIZE1 * LenFDL, MPI_FLOAT, neigyid[0], 5121,
//                 wp_xzr0, nx * lnz * WSIZE1 * LenFDL, MPI_FLOAT, neigyid[1], 5121,
//                 SWMPI_COMM, &stat );
//   unpack(wp_xzr0, &my_W[pos_w_d], nx, lnz * WSIZE1 * LenFDL, ny * lnz * WSIZE1); 

//   //hd
//   pos_w_s = nj1 * lnz; pos_w_d = nj2 * lnz;
//   pack(&my_W_Tyz[pos_w_s], &wp_xzs0[seg1], nx, lnz * LenFDL, ny * lnz); 
//   MPI_Sendrecv(&wp_xzs0[seg1], nx * lnz * LenFDL, MPI_FLOAT, neigyid[0], 5123,
//                 &wp_xzr0[seg1], nx * lnz * LenFDL, MPI_FLOAT, neigyid[1], 5123,
//                 SWMPI_COMM, &stat );
//   unpack(&wp_xzr0[seg1], &my_W_Tyz[pos_w_d], nx, lnz * LenFDL, ny * lnz); 

//   // mergeW(W, my_W, my_W_Tyz);

//   pos_w_s = ((nj2-LenFDL) * lnz) * WSIZE1; pos_w_d = (0 * lnz) * WSIZE1;
//   //MPI_Sendrecv( &W[pos_w_s], 1, WDTypeYL, neigyid[1], 5122,
//   //              &W[pos_w_d], 1, WDTypeYL, neigyid[0], 5122,
//   //              SWMPI_COMM, &stat );
//   pack(&my_W[pos_w_s], wp_xzs1, nx, lnz * WSIZE1 * LenFDL, ny * lnz * WSIZE1); 
//   MPI_Sendrecv(wp_xzs1, nx * lnz * WSIZE1 * LenFDL, MPI_FLOAT, neigyid[1], 5122,
//                 wp_xzr1, nx * lnz * WSIZE1 * LenFDL, MPI_FLOAT, neigyid[0], 5122,
//                 SWMPI_COMM, &stat );
//   unpack(wp_xzr1, &my_W[pos_w_d], nx, lnz * WSIZE1 * LenFDL, ny * lnz * WSIZE1); 

//   //hd
//   pos_w_s = (nj2-LenFDL) * lnz; pos_w_d = 0 * lnz;
//   pack(&my_W_Tyz[pos_w_s], &wp_xzs1[seg1], nx, lnz * LenFDL, ny * lnz); 
//   MPI_Sendrecv(&wp_xzs1[seg1], nx * lnz * LenFDL, MPI_FLOAT, neigyid[1], 5124,
//                 &wp_xzr1[seg1], nx * lnz * LenFDL, MPI_FLOAT, neigyid[0], 5124,
//                 SWMPI_COMM, &stat );
//   unpack(&wp_xzr1[seg1], &my_W_Tyz[pos_w_d], nx, lnz * LenFDL, ny * lnz);
 
//   // mergeW(W, my_W, my_W_Tyz);

//   // exchange z plane data

//   // splitW(W, my_W, my_W_Tyz);
//   pos_w_s = (nk1) * WSIZE1; pos_w_d = (nk2) * WSIZE1;
//   //MPI_Sendrecv( &W[pos_w_s], 1, WDTypeZL, neigzid[0], 5131,
//   //              &W[pos_w_d], 1, WDTypeZL, neigzid[1], 5131,
//   //              SWMPI_COMM, &stat);
//   pack(&my_W[pos_w_s], wp_xys0, nx * ny, WSIZE1 * LenFDL, lnz * WSIZE1); 
//   MPI_Sendrecv(wp_xys0, nx * ny * WSIZE1 * LenFDL, MPI_FLOAT, neigzid[0], 5131,
//                 wp_xyr0, nx * ny * WSIZE1 * LenFDL, MPI_FLOAT, neigzid[1], 5131,
//                 SWMPI_COMM, &stat );
//   unpack(wp_xyr0, &my_W[pos_w_d], nx * ny, WSIZE1 * LenFDL, lnz * WSIZE1); 

//   //hd
//   pos_w_s = nk1; pos_w_d = nk2;
//   pack(&my_W_Tyz[pos_w_s], &wp_xys0[seg2], nx * ny, LenFDL, lnz); 
//   MPI_Sendrecv(&wp_xys0[seg2], nx * ny * LenFDL, MPI_FLOAT, neigzid[0], 5133,
//                 &wp_xyr0[seg2], nx * ny * LenFDL, MPI_FLOAT, neigzid[1], 5133,
//                 SWMPI_COMM, &stat );
//   unpack(&wp_xyr0[seg2], &my_W_Tyz[pos_w_d], nx * ny, LenFDL, lnz); 

//   // mergeW(W, my_W, my_W_Tyz);

  
//   pos_w_s = (nk2-LenFDL) * WSIZE1; pos_w_d = (0) * WSIZE1;
//   //MPI_Sendrecv( &W[pos_w_s], 1, WDTypeZL, neigzid[1], 5132,
//   //              &W[pos_w_d], 1, WDTypeZL, neigzid[0], 5132,
//   //              SWMPI_COMM, &stat);
//   pack(&my_W[pos_w_s], wp_xys1, nx * ny, WSIZE1 * LenFDL, lnz * WSIZE1); 
//   MPI_Sendrecv(wp_xys1, nx * ny * WSIZE1 * LenFDL, MPI_FLOAT, neigzid[1], 5132,
//                 wp_xyr1, nx * ny * WSIZE1 * LenFDL, MPI_FLOAT, neigzid[0], 5132,
//                 SWMPI_COMM, &stat );
//   unpack(wp_xyr1, &my_W[pos_w_d], nx * ny, WSIZE1 * LenFDL, lnz * WSIZE1); 
  
//   // printf("rank: %d x: %d %d y: %d %d z: %d %d \n", this_rank, neigxid[0], neigxid[1], neigyid[0], neigyid[1], neigzid[0], neigzid[1]);

//   //hd
//   pos_w_s = nk2-LenFDL; pos_w_d = 0;
//   pack(&my_W_Tyz[pos_w_s], &wp_xys1[seg2], nx * ny, LenFDL, lnz); 
//   MPI_Sendrecv(&wp_xys1[seg2], nx * ny * LenFDL, MPI_FLOAT, neigzid[1], 5134,
//                 &wp_xyr1[seg2], nx * ny * LenFDL, MPI_FLOAT, neigzid[0], 5134,
//                 SWMPI_COMM, &stat );
//   unpack(&wp_xyr1[seg2], &my_W_Tyz[pos_w_d], nx * ny, LenFDL, lnz);
//   // mergeW(W, my_W, my_W_Tyz);
//   return 0;
// }

int exchange_Wave2(float *my_W, float *my_W_Tyz) {
  MPI_Status stat;
  int i;
  int WSIZE1 = WSIZE - 1;
  int seg0 = ny * lnz * (WSIZE-1) * LenFDL;
  int seg1 = nx * lnz * (WSIZE-1) * LenFDL;
  int seg2 = nx * ny * (WSIZE-1) * LenFDL;
  
  //exchange x plane data
  int pos_w_s, pos_w_d;
  pos_w_s = (ni1 * ny * lnz) * WSIZE1; pos_w_d = (ni2 * ny * lnz) * WSIZE1;
  
  // pack(&my_W[pos_w_s], wp_yzs0, LenFDL, ny * lnz * WSIZE1, ny * lnz * WSIZE1); 
  // MPI_Sendrecv(wp_yzs0, ny * lnz * WSIZE1 * LenFDL, MPI_FLOAT, neigxid[0], 5111,
  //               wp_yzr0, ny * lnz * WSIZE1 * LenFDL, MPI_FLOAT, neigxid[1], 5111,
  //               SWMPI_COMM, &stat );
  // unpack(wp_yzr0, &my_W[pos_w_d], LenFDL, ny * lnz * WSIZE1, ny * lnz * WSIZE1); 

  int stride_8 = ny * lnz * WSIZE1;
  int blocklen_8 = ny * lnz * WSIZE1;
  for(i = 0; i < LenFDL; i++){
      MPI_Sendrecv(&my_W[pos_w_s + i*stride_8], blocklen_8, MPI_FLOAT, neigxid[0], 5111,
                &my_W[pos_w_d + i*stride_8], blocklen_8, MPI_FLOAT, neigxid[1], 5111,
                SWMPI_COMM, &stat );
  }

  //hd
  pos_w_s = ni1 * ny * lnz; pos_w_d = ni2 * ny * lnz;
  // pack(&my_W_Tyz[pos_w_s], &wp_yzs0[seg0], LenFDL, ny * lnz, ny * lnz); 
  // MPI_Sendrecv(&wp_yzs0[seg0], ny * lnz * LenFDL, MPI_FLOAT, neigxid[0], 5113,
  //               &wp_yzr0[seg0], ny * lnz * LenFDL, MPI_FLOAT, neigxid[1], 5113,
  //               SWMPI_COMM, &stat );
  // unpack(&wp_yzr0[seg0], &my_W_Tyz[pos_w_d], LenFDL, ny * lnz, ny * lnz);  

  int stride_1 = ny * lnz;
  int blocklen_1 = ny * lnz;
  for(i = 0; i < LenFDL; i++){
      MPI_Sendrecv(&my_W_Tyz[pos_w_s + i*stride_1], blocklen_1, MPI_FLOAT, neigxid[0], 5113,
                &my_W_Tyz[pos_w_d + i*stride_1], blocklen_1, MPI_FLOAT, neigxid[1], 5113,
                SWMPI_COMM, &stat );
  }




  pos_w_s = ((ni2-LenFDL) * ny * lnz) * WSIZE1; pos_w_d = (0 * ny * lnz) * WSIZE1;
  
  // pack(&my_W[pos_w_s], wp_yzs1, LenFDL, ny * lnz * WSIZE1, ny * lnz * WSIZE1); 
  // MPI_Sendrecv(wp_yzs1, ny * lnz * WSIZE1 * LenFDL, MPI_FLOAT, neigxid[1], 5112,
  //               wp_yzr1, ny * lnz * WSIZE1 * LenFDL, MPI_FLOAT, neigxid[0], 5112,
  //               SWMPI_COMM, &stat );
  // unpack(wp_yzr1, &my_W[pos_w_d], LenFDL, ny * lnz * WSIZE1, ny * lnz * WSIZE1); 

  for(i = 0; i < LenFDL; i++){
      MPI_Sendrecv(&my_W[pos_w_s + i*stride_8], blocklen_8, MPI_FLOAT, neigxid[1], 5112,
                &my_W[pos_w_d + i*stride_8], blocklen_8, MPI_FLOAT, neigxid[0], 5112,
                SWMPI_COMM, &stat );
  }

  //hd
  pos_w_s = (ni2-LenFDL) * ny * lnz; pos_w_d = 0 * ny * lnz;
  // pack(&my_W_Tyz[pos_w_s], &wp_yzs1[seg0], LenFDL, ny * lnz, ny * lnz); 
  // MPI_Sendrecv(&wp_yzs1[seg0], ny * lnz * LenFDL, MPI_FLOAT, neigxid[1], 5114,
  //               &wp_yzr1[seg0], ny * lnz * LenFDL, MPI_FLOAT, neigxid[0], 5114,
  //               SWMPI_COMM, &stat );
  // unpack(&wp_yzr1[seg0], &my_W_Tyz[pos_w_d], LenFDL, ny * lnz, ny * lnz); 

  for(i = 0; i < LenFDL; i++){
      MPI_Sendrecv(&my_W_Tyz[pos_w_s + i*stride_1], blocklen_1, MPI_FLOAT, neigxid[1], 5114,
                &my_W_Tyz[pos_w_d + i*stride_1], blocklen_1, MPI_FLOAT, neigxid[0], 5114,
                SWMPI_COMM, &stat );
  }

  
  pos_w_s = (nj1 * lnz) * WSIZE1; pos_w_d = (nj2 * lnz) * WSIZE1;
  
  pack(&my_W[pos_w_s], wp_xzs0, nx, lnz * WSIZE1 * LenFDL, ny * lnz * WSIZE1); 
  MPI_Sendrecv(wp_xzs0, nx * lnz * WSIZE1 * LenFDL, MPI_FLOAT, neigyid[0], 5121,
                wp_xzr0, nx * lnz * WSIZE1 * LenFDL, MPI_FLOAT, neigyid[1], 5121,
                SWMPI_COMM, &stat );
  unpack(wp_xzr0, &my_W[pos_w_d], nx, lnz * WSIZE1 * LenFDL, ny * lnz * WSIZE1); 

  // stride_8 = ny * lnz * WSIZE1;
  // blocklen_8 = lnz * WSIZE1 * LenFDL;
  // for(i = 0; i < nx; i++){
  //     MPI_Sendrecv(&my_W[pos_w_s + i*stride_8], blocklen_8, MPI_FLOAT, neigxid[0], 5121,
  //               &my_W[pos_w_d + i*stride_8], blocklen_8, MPI_FLOAT, neigxid[1], 5121,
  //               SWMPI_COMM, &stat );
  // }

  //hd
  pos_w_s = nj1 * lnz; pos_w_d = nj2 * lnz;
  pack(&my_W_Tyz[pos_w_s], &wp_xzs0[seg1], nx, lnz * LenFDL, ny * lnz); 
  MPI_Sendrecv(&wp_xzs0[seg1], nx * lnz * LenFDL, MPI_FLOAT, neigyid[0], 5123,
                &wp_xzr0[seg1], nx * lnz * LenFDL, MPI_FLOAT, neigyid[1], 5123,
                SWMPI_COMM, &stat );
  unpack(&wp_xzr0[seg1], &my_W_Tyz[pos_w_d], nx, lnz * LenFDL, ny * lnz); 
  // stride_1 = ny * lnz;
  // blocklen_1 = lnz * LenFDL;
  // for(i = 0; i < nx; i++){
  //     MPI_Sendrecv(&my_W_Tyz[pos_w_s + i*stride_1], blocklen_1, MPI_FLOAT, neigxid[0], 5113,
  //               &my_W_Tyz[pos_w_d + i*stride_1], blocklen_1, MPI_FLOAT, neigxid[1], 5113,
  //               SWMPI_COMM, &stat );
  // }

  pos_w_s = ((nj2-LenFDL) * lnz) * WSIZE1; pos_w_d = (0 * lnz) * WSIZE1;
  
  pack(&my_W[pos_w_s], wp_xzs1, nx, lnz * WSIZE1 * LenFDL, ny * lnz * WSIZE1); 
  MPI_Sendrecv(wp_xzs1, nx * lnz * WSIZE1 * LenFDL, MPI_FLOAT, neigyid[1], 5122,
                wp_xzr1, nx * lnz * WSIZE1 * LenFDL, MPI_FLOAT, neigyid[0], 5122,
                SWMPI_COMM, &stat );
  unpack(wp_xzr1, &my_W[pos_w_d], nx, lnz * WSIZE1 * LenFDL, ny * lnz * WSIZE1); 

  // for(i = 0; i < nx; i++){
  //     MPI_Sendrecv(&my_W[pos_w_s + i*stride_8], blocklen_8, MPI_FLOAT, neigxid[1], 5122,
  //               &my_W[pos_w_d + i*stride_8], blocklen_8, MPI_FLOAT, neigxid[0], 5122,
  //               SWMPI_COMM, &stat );
  // }

  //hd
  pos_w_s = (nj2-LenFDL) * lnz; pos_w_d = 0 * lnz;
  pack(&my_W_Tyz[pos_w_s], &wp_xzs1[seg1], nx, lnz * LenFDL, ny * lnz); 
  MPI_Sendrecv(&wp_xzs1[seg1], nx * lnz * LenFDL, MPI_FLOAT, neigyid[1], 5124,
                &wp_xzr1[seg1], nx * lnz * LenFDL, MPI_FLOAT, neigyid[0], 5124,
                SWMPI_COMM, &stat );
  unpack(&wp_xzr1[seg1], &my_W_Tyz[pos_w_d], nx, lnz * LenFDL, ny * lnz);

  // for(i = 0; i < nx; i++){
  //     MPI_Sendrecv(&my_W_Tyz[pos_w_s + i*stride_1], blocklen_1, MPI_FLOAT, neigxid[0], 5113,
  //               &my_W_Tyz[pos_w_d + i*stride_1], blocklen_1, MPI_FLOAT, neigxid[1], 5113,
  //               SWMPI_COMM, &stat );
  // }

  // splitW(W, my_W, my_W_Tyz);
  pos_w_s = (nk1) * WSIZE1; pos_w_d = (nk2) * WSIZE1;
  
  pack(&my_W[pos_w_s], wp_xys0, nx * ny, WSIZE1 * LenFDL, lnz * WSIZE1); 
  MPI_Sendrecv(wp_xys0, nx * ny * WSIZE1 * LenFDL, MPI_FLOAT, neigzid[0], 5131,
                wp_xyr0, nx * ny * WSIZE1 * LenFDL, MPI_FLOAT, neigzid[1], 5131,
                SWMPI_COMM, &stat );
  unpack(wp_xyr0, &my_W[pos_w_d], nx * ny, WSIZE1 * LenFDL, lnz * WSIZE1); 

  // stride_8 = lnz * WSIZE1;
  // blocklen_8 = WSIZE1 * LenFDL;
  // for(i = 0; i < nx * ny; i++){
  //     MPI_Sendrecv(&my_W[pos_w_s + i*stride_8], blocklen_8, MPI_FLOAT, neigxid[0], 5111,
  //               &my_W[pos_w_d + i*stride_8], blocklen_8, MPI_FLOAT, neigxid[1], 5111,
  //               SWMPI_COMM, &stat );
  // }

  //hd
  pos_w_s = nk1; pos_w_d = nk2;
  pack(&my_W_Tyz[pos_w_s], &wp_xys0[seg2], nx * ny, LenFDL, lnz); 
  MPI_Sendrecv(&wp_xys0[seg2], nx * ny * LenFDL, MPI_FLOAT, neigzid[0], 5133,
                &wp_xyr0[seg2], nx * ny * LenFDL, MPI_FLOAT, neigzid[1], 5133,
                SWMPI_COMM, &stat );
  unpack(&wp_xyr0[seg2], &my_W_Tyz[pos_w_d], nx * ny, LenFDL, lnz); 

  // stride_1 = lnz;
  // blocklen_1 = LenFDL;
  // for(i = 0; i < nx * ny; i++){
  //     MPI_Sendrecv(&my_W_Tyz[pos_w_s + i*stride_1], blocklen_1, MPI_FLOAT, neigxid[0], 5113,
  //               &my_W_Tyz[pos_w_d + i*stride_1], blocklen_1, MPI_FLOAT, neigxid[1], 5113,
  //               SWMPI_COMM, &stat );
  // }

  // mergeW(W, my_W, my_W_Tyz);

  
  pos_w_s = (nk2-LenFDL) * WSIZE1; pos_w_d = (0) * WSIZE1;

  pack(&my_W[pos_w_s], wp_xys1, nx * ny, WSIZE1 * LenFDL, lnz * WSIZE1); 
  MPI_Sendrecv(wp_xys1, nx * ny * WSIZE1 * LenFDL, MPI_FLOAT, neigzid[1], 5132,
                wp_xyr1, nx * ny * WSIZE1 * LenFDL, MPI_FLOAT, neigzid[0], 5132,
                SWMPI_COMM, &stat );
  unpack(wp_xyr1, &my_W[pos_w_d], nx * ny, WSIZE1 * LenFDL, lnz * WSIZE1); 
  // for(i = 0; i < nx * ny; i++){
  //     MPI_Sendrecv(&my_W[pos_w_s + i*stride_8], blocklen_8, MPI_FLOAT, neigxid[0], 5111,
  //               &my_W[pos_w_d + i*stride_8], blocklen_8, MPI_FLOAT, neigxid[1], 5111,
  //               SWMPI_COMM, &stat );
  // }

  //hd
  pos_w_s = nk2-LenFDL; pos_w_d = 0;
  pack(&my_W_Tyz[pos_w_s], &wp_xys1[seg2], nx * ny, LenFDL, lnz); 
  MPI_Sendrecv(&wp_xys1[seg2], nx * ny * LenFDL, MPI_FLOAT, neigzid[1], 5134,
                &wp_xyr1[seg2], nx * ny * LenFDL, MPI_FLOAT, neigzid[0], 5134,
                SWMPI_COMM, &stat );
  unpack(&wp_xyr1[seg2], &my_W_Tyz[pos_w_d], nx * ny, LenFDL, lnz);
  // for(i = 0; i < nx * ny; i++){
  //     MPI_Sendrecv(&my_W_Tyz[pos_w_s + i*stride_1], blocklen_1, MPI_FLOAT, neigxid[0], 5113,
  //               &my_W_Tyz[pos_w_d + i*stride_1], blocklen_1, MPI_FLOAT, neigxid[1], 5113,
  //               SWMPI_COMM, &stat );
  // }

  return 0;
}

int exchange_Wave(float *W) {
  MPI_Status stat;
  //exchange x plane data
  int pos_w_s, pos_w_d;
  pos_w_s = (ni1 * ny * lnz) * WSIZE; pos_w_d = (ni2 * ny * lnz) * WSIZE;
  MPI_Sendrecv( &W[pos_w_s], 1, WDTypeXL, neigxid[0], 5111,
                &W[pos_w_d], 1, WDTypeXL, neigxid[1], 5111,
                SWMPI_COMM, &stat );
  pos_w_s = ((ni2-LenFDL) * ny * lnz) * WSIZE; pos_w_d = (0 * ny * lnz) * WSIZE;
  MPI_Sendrecv( &W[pos_w_s], 1, WDTypeXL, neigxid[1], 5112,
                &W[pos_w_d], 1, WDTypeXL, neigxid[0], 5112,
                SWMPI_COMM, &stat );

  // exchange y plane data
  pos_w_s = (nj1 * lnz) * WSIZE; pos_w_d = (nj2 * lnz) * WSIZE;
  MPI_Sendrecv( &W[pos_w_s], 1, WDTypeYL, neigyid[0], 5121,
                &W[pos_w_d], 1, WDTypeYL, neigyid[1], 5121,
                SWMPI_COMM, &stat );
  pos_w_s = ((nj2-LenFDL) * lnz) * WSIZE; pos_w_d = (0 * lnz) * WSIZE;
  MPI_Sendrecv( &W[pos_w_s], 1, WDTypeYL, neigyid[1], 5122,
                &W[pos_w_d], 1, WDTypeYL, neigyid[0], 5122,
                SWMPI_COMM, &stat );

  // exchange z plane data
  pos_w_s = (nk1) * WSIZE; pos_w_d = (nk2) * WSIZE;
  MPI_Sendrecv( &W[pos_w_s], 1, WDTypeZL, neigzid[0], 5131,
                &W[pos_w_d], 1, WDTypeZL, neigzid[1], 5131,
                SWMPI_COMM, &stat);
  pos_w_s = (nk2-LenFDL) * WSIZE; pos_w_d = (0) * WSIZE;
  MPI_Sendrecv( &W[pos_w_s], 1, WDTypeZL, neigzid[1], 5132,
                &W[pos_w_d], 1, WDTypeZL, neigzid[0], 5132,
                SWMPI_COMM, &stat);
  return 0;
}
