#include "params.h"
#include "macdrp.h"

/*
 *********************************************
 *    Global variables                       *
 *    ONLY ALLOW DEFINE ONCE HERE !          *
 *********************************************
 */
int this_rank;
int sub_rank_x;
int sub_rank_y;
int NX, NY, NZ, PX, PY, PZ;
float DH, DT, TMAX;
int NT, RNT;
int TSKP;
int TSKP0;

int SAVE_CKPT;
int SAVE_CKPT_START_STEP;
int SAVE_CKPT_PER_STEP;
int CUR_CKPT_IDX;
int NUM_CKPT_DIR;
int LOAD_CKPT;

int SAMP;
int SAMPLE_SIZE;

int Gather;
int Gather_SIZE_X;
int Gather_SIZE_Y;
//snapshot settings
int NBGX;
int NBGY;
int NBGZ;

int NEDX;
int NEDY;
int NEDZ;

int NSKPX;
int NSKPY;
int NSKPZ;

int WRITE_STEP;

int SAVE_FULL_IMG;
int SAVE_FULL_IMG_START_STEP;
int SAVE_FULL_IMG_PER_STEP;

char OUT[100];
char OUT_DIR[512];
//int nsrclocal;
int icoord = 0;
int isource = 0;

char INGRD[128], INVEL[128], INSRC[128], INREC[128];

char IN_TOPO[512], IN_SOURCE[512], IN_REC[1024];
char IN_INTER[1024];
char IN_VELO[1024];
int NUM_INTER;

MPI_Comm SWMPI_COMM;

int thisid[3];
int masternode = 0; 
int absnode = 0; 
int freenode =0;
int neigxid[2], neigyid[2], neigzid[2];

int ni1, ni2, nj1, nj2, nk1, nk2;
int nx1, nx2, ny1, ny2, nz1, nz2;
//for aligned padding
int lnz;
int ni, nj, nk, nx, ny, nz;
int isx1, isx2, isy1, isy2, isz1, isz2;
int ngi1, ngi2, ngj1, ngj2, ngk1, ngk2;
int ngx1, ngx2, ngy1, ngy2, ngz1, ngz2;
// i, j : exclude ghost points
// x, z : include ghost points
MPI_Datatype DTypeXS, DTypeYS, DTypeZS; // shorter part FD
MPI_Datatype DTypeXL, DTypeYL, DTypeZL; // longer part FD
MPI_Datatype WDTypeXS, WDTypeYS, WDTypeZS; // shorter part FD
MPI_Datatype WDTypeXL, WDTypeYL, WDTypeZL; // longer part FD
MPI_Datatype MDTypeXS, MDTypeYS, MDTypeZS; // shorter part FD
MPI_Datatype MDTypeXL, MDTypeYL, MDTypeZL; // longer part FD
MPI_Datatype CDTypeXS, CDTypeYS, CDTypeZS; // shorter part FD
MPI_Datatype CDTypeXL, CDTypeYL, CDTypeZL; // longer part FD
//MPI_Request reqXF[36], reqXB[36];
//MPI_Request reqYF[36], reqYB[36];
//MPI_Request reqZF[36], reqZB[36];

int WSIZE = 9; int MSIZE = 14;//13; // stride size times  when transform struct
int WSIZE_V = 8;
int DSIZE = 3; // 4 if include attenuation
int CSIZE = 3; // (x, y, z)
int HALO = 3;
//float rho = 1.2e3;
//float vp = 3.0e3;
//float vs = 1.5e3;
//float lam, lam2mu, mu ;
// source
int nsrc, src_nt;
int nrec;//, nrec_local;
float src_dt;
int *srci, *srcj, *srck;
float *area, *strike, *dip, *rake, *rate;
float rickerfc, Qsf0;
float gauss_height, gauss_width;
//int   vecF[5] = {-1, 0, 1, 2, 3};
//int   vecB[5] = {-3, -2, -1, 0, 1};

float coeF[5] = {-0.30874, -0.6326, 1.2330, -0.3334, 0.04168};
float coeB[5] = {-0.04168, 0.3334, -1.2330, 0.6326, 0.30874};

float coeF_DH[5];
float coeB_DH[5];

float coe24F[3] = {-7.0/6.0, 8.0/6.0, -1.0/6.0};
float coe24B[3] = {1.0/6.0, -8.0/6.0, 7.0/6.0};

float coe24F_DH[3];
float coe24B_DH[3];

float coe22F[2] = {-1.0, 1.0};
float coe22B[2] = {-1.0, 1.0};

float rDH;

float RK4a[3] = {0.5, 0.5, 1.0};
float RK4b[4] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};

float RK6a[5] = {0.353323, 0.999597, 0.152188, 0.534216, 0.603907};
float RK6b[6] = {0.0467621, 0.137286, 0.170975, 0.197572, 0.282263, 0.165142};

float *wp_yzs0, *wp_yzr0;
float *wp_yzs1, *wp_yzr1;
float *wp_xzs0, *wp_xzr0;
float *wp_xzs1, *wp_xzr1;
float *wp_xys0, *wp_xyr0;
float *wp_xys1, *wp_xyr1;
float *wp_xzs0_test;

