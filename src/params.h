#include <mpi.h>

#ifndef PARAMS_H
#define PARAMS_H

#define NUM_RANK_PER_DIR (5000)
#define NUM_ROUND_OUTPUT (30) /// 16 is suggested by Wang, Xiyang

extern int this_rank;
extern int sub_rank_x;
extern int sub_rank_y;
extern int NX, NY, NZ, PX, PY, PZ; 
extern float DH, DT, TMAX;
extern int NT, RNT;
extern int TSKP0;
extern int TSKP;
extern float rDH;

extern int icoord;
extern int isource;
extern int nrec;//, nrec_local;
//extern int nsrclocal;
extern int SAVE_CKPT;
extern int SAVE_CKPT_START_STEP;
extern int SAVE_CKPT_PER_STEP;
extern int CUR_CKPT_IDX;
extern int NUM_CKPT_DIR;
extern int LOAD_CKPT;

extern int SAMP;
extern int SAMPLE_SIZE;

extern int Gather;
extern int Gather_SIZE_X;
extern int Gather_SIZE_Y;
//snapshot settings
extern int NBGX;
extern int NBGY;
extern int NBGZ;

extern int NEDX;
extern int NEDY;
extern int NEDZ;

extern int NSKPX;
extern int NSKPY;
extern int NSKPZ;

extern int WRITE_STEP;

extern int SAVE_FULL_IMG;
extern int SAVE_FULL_IMG_START_STEP;
extern int SAVE_FULL_IMG_PER_STEP;

extern char OUT[100];

extern MPI_Comm SWMPI_COMM;

extern int thisid[3];
extern int masternode, freenode, absnode;
extern int neigxid[2], neigyid[2], neigzid[2];
extern int WSIZE, MSIZE, DSIZE, CSIZE;
extern int WSIZE_V;

extern int ni1, ni2, nj1, nj2, nk1, nk2;
extern int nx1, nx2, ny1, ny2, nz1, nz2;
extern int ni, nj, nk, nx, ny, nz;
extern int isx1, isx2, isy1, isy2, isz1, isz2;
//for aligned padding
extern int lnz;
extern int ngi1, ngi2, ngj1, ngj2, ngk1, ngk2;
extern int ngx1, ngx2, ngy1, ngy2, ngz1, ngz2;
// i, j : exclude ghost points
// x, z : include ghost points
extern MPI_Datatype DTypeXS, DTypeYS, DTypeZS; // for shorter part of difference
extern MPI_Datatype DTypeXL, DTypeYL, DTypeZL; // for larger part of difference

extern MPI_Datatype WDTypeXS, WDTypeYS, WDTypeZS; // for shorter part of difference
extern MPI_Datatype WDTypeXL, WDTypeYL, WDTypeZL; // for larger part of difference

extern MPI_Datatype MDTypeXS, MDTypeYS, MDTypeZS; // for shorter part of difference
extern MPI_Datatype MDTypeXL, MDTypeYL, MDTypeZL; // for larger part of difference

extern MPI_Datatype CDTypeXS, CDTypeYS, CDTypeZS; // for shorter part of difference
extern MPI_Datatype CDTypeXL, CDTypeYL, CDTypeZL; // for larger part of difference
//extern MPI_Request reqXF[36], reqXB[36];
//extern MPI_Request reqYF[36], reqYB[36];
//extern MPI_Request reqZF[36], reqZB[36];

//extern float vp, vs, rho, lam, lam2mu, mu;
// source
extern int nsrc, src_nt;
extern float src_dt;
extern int *srci, *srcj, *srck;
extern float *area, *strike, *dip, *rake, *rate;
extern float rickerfc, Qsf0;
extern float gauss_height, gauss_width;

extern char INGRD[128];
extern char INVEL[128];
extern char INSRC[128];
extern char INREC[128];

extern int NUM_INTER;
extern char IN_TOPO[512];
extern char IN_VELO[1024];
extern char IN_INTER[1024];
extern char IN_SOURCE[512];
extern char IN_REC[1024];
extern char OUT_DIR[512];

extern float *wp_yzs0, *wp_yzr0;
extern float *wp_yzs1, *wp_yzr1;
extern float *wp_xzs0, *wp_xzr0;
extern float *wp_xzs1, *wp_xzr1;
extern float *wp_xys0, *wp_xyr0;
extern float *wp_xys1, *wp_xyr1;

extern char* buff_addr;
extern float *momentRate;
#endif
