/*
********************************************************************************
* common.h                                                                     *
* programming in C language                                                    *
* all data types and functions are defined here                                *
********************************************************************************
*/
#include <mpi.h>
// #include "macdrp.h"
#include "pml_macro.h"
#ifdef __RESTRICT
#define RESTRICT restrict
#else
#define RESTRICT 
#endif

#ifndef _COMMON_H
#define _COMMON_H

void mergeW(float* W, float* my_W, float* my_W_Tyz);
void splitW(float* W, float* my_W, float* my_W_Tyz);

#define _max(a,b)(a>b?a:b)
#define _min(a,b)(a<b?a:b)

#define LenFDL 3
#define LenFDS 1

#define PI 3.141592653589793238463

typedef float *** Grid3D;
typedef float ** Grid2D;
typedef float * Grid1D;

int Malloc( void ** mem, long long  size);

Grid3D Alloc3D(int nx, int ny, int nz);
Grid2D Alloc2D(int nx, int nz);
Grid1D Alloc1D(int nx);

void Delloc3D(Grid3D U);
void Delloc2D(Grid2D U);
void Delloc1D(Grid1D U);

//void Alloc_Coord(struct Coord *);
//void Alloc_Metric(struct Metric *);
//void Alloc_WaveField(struct Wave *);
// float *Alloc_wave();
float *Alloc_wave_8();
float *Alloc_wave_1();
float *Alloc_metric();
float *Alloc_media();
float *Alloc_coord();
float *Alloc_attenu();

int swmpi_datatype(void);
int swmpi_datatype_free(void);
//void coord_extend(Grid3D U);
//int  coord_exchange(struct Coord*);
int coord_exchange(float *c);
int metric_exchange(float *m);
void extend_coord(float *C, int SIZE);
void extend_Symm_array(float *W, int SIZE);
//void cal_metric(struct Coord *C, float *M);
//void extend_Symm(Grid3D U);

/* coord */
int construct_coord(float *C);
int read_coord(float *C);
int read_coord_part(float *C);
int gauss_coord(float *C, float gauss_height, float gauss_width);
int cal_range_steph(float *C, float *range);
/* metric */
void cal_metric(float *C, float *M);
/* media */
int construct_media(float *C, float *M);
int cal_range_media(float *M, float *range);
/* math */
//void invert3x3(float A[][3]);
//void matmul3x3(float A[][3], float B[][3], float C[][3]);
void cross_product(float *A, float *B, float *C);
float dot_product(float *A, float *B);
float dist_point2plane(float x0[3], float x1[3], float x2[3], float x3[3]);

/* source */
//int Init_Source(struct Source *S);
//int Read_SourceHeader(char * filename, struct Source *S);
//int Locate_Source(struct Source*, struct Coord*);
//int Alloc_Source(struct Source *S);
//int Load_Source(char * filename, struct Source *S);
//int Add_Source(float *W, float *M, struct Source *S, int it);
//int Check_Source(struct Source*);
int cal_source(int it, float *S);
int read_source(char *filename, char *);
// int AddSourceRicker(float *W, float *M, int it);
int AddSourceRicker(float *W_8, float *W_1, float *M, int it);
int add_source(float *W, float *M, float *S);
//int add_source(float *W, float *M, float *D, float *S);
//int Add_Source(float *W, float *M, int it);
float ricker_deriv(float rickerfc, float t, float tdelay);
/* macdrp */
//int coef_fdxy2fdz(struct Metric M, Grid3D matVx2Vz, Grid3D matVy2Vz, Grid3D matF2Vz);
//int LxF_LyF_LzF_TIMG(struct Wave W, struct Wave hW, struct Metric M);
//int LxB_LyF_LzB_TIMG(struct Wave W, struct Wave hW, struct Metric M);
//int LxF_LzF_VHOC(struct Wave W, struct Metric M);
//int LxB_LzB_VHOC(struct Wave W, struct Metric M);
//int LxF_LzF_VLOW(struct Wave W, struct Wave hW, struct Metric M);
//int LxB_LzB_VLOW(struct Wave W, struct Wave hW, struct Metric M);
//int macdrp_mesg_init(struct Wave W, 
//    MPI_Request *reqXF, MPI_Request *reqXB,
//    MPI_Request *reqYF, MPI_Request *reqYB,
//    MPI_Request *reqZF, MPI_Request *reqZB );
/*
int exchange_FFF(struct Wave*);
int exchange_BBB(struct Wave*);
int exchange_BBF(struct Wave*);
int exchange_FFB(struct Wave*);
int exchange_BFF(struct Wave*);
int exchange_FBB(struct Wave*);
int exchange_FBF(struct Wave*);
int exchange_BFB(struct Wave*);
int exchange_Wave(float *);

int RK_step1(struct Wave*, struct Metric*); 
int RK_step2(struct Wave*, struct Metric*); 
int RK_step3(struct Wave*, struct Metric*); 
int RK_step4(struct Wave*, struct Metric*); 
int RK_step5(struct Wave*, struct Metric*); 
int RK_step6(struct Wave*, struct Metric*); 
int RK_step7(struct Wave*, struct Metric*); 
int RK_step8(struct Wave*, struct Metric*); 
*/
//int Calculate(float *W, float *mW, float *hW, float *tW, float *M, int Flag[3], float rka, float rkb,int flag);
//int Calculate(float *W, float *mW, float *hW, float *tW, float *M, int Flag[3], float rka, float rkb,int flag, int nii1, int nii2, int njj1, int njj2, int nkk1, int nkk2);


/* Absorbing Boundary Condition */
int abs_exp(Grid3D damp);
void Apply_abs(float *, Grid3D);
void Apply_attenuation(float*, float*);
void cal_attenu(float*, float*);

/* IO */
void paras_msg(void);
//int read_Coord(struct Coord *C);
int write_xplane(char *, Grid3D U, int gi);
int write_yplane(char *, Grid3D U, int gj);
int write_zplane(char *, Grid3D U, int gk);
//int write_coord_xplane(struct Coord*, int );
//int write_coord_yplane(struct Coord*, int );
//int write_coord_zplane(struct Coord*, int );
/*
int write_metric_xplane(struct Metric*, int );
int write_metric_yplane(struct Metric*, int );
int write_metric_zplane(struct Metric*, int );
*/
// IO new
int write_xplane_tmp(char *prefix, float *w, int shift, int Unit, int gi);
int write_yplane_tmp(char *prefix, float *w, int shift, int Unit, int gj);
int write_zplane_tmp(char *prefix, float *w, int shift, int Unit, int gk);

/* receiver */
//int Locate_Receiver(char * filename, struct Receiver *R, struct Coord *C, int rank);
//int keep_seismo(struct Receiver*, float*, int it);
//int export_seismo(struct Receiver*);
void geo2cart(float* , float*,  float, float, float, float, float);
void seek_index(int*, float*, float*, float, float);
void get_nrec();
int Load_Receiver(int*, int);
int keep_seismo(float*, int*, float*, int);
int export_seismo(float*, int*, int);
/*
float * transform_wave(struct Wave * W);
void inverse_wave(float *w, struct Wave *W);
float * transform_metric(struct Metric *M);
void inverse_metric(float *m, struct Metric *M);
*/
//int Gauss3D(Grid3D U, struct Coord C);
void athread_get_local(float *source, float *local, int blocksize);

void athread_put_local(float *local, float *source, int blocksize);





int signal_bind(void);

void calcRecordingPoints(int *rec_nbgx, int *rec_nedx, int *rec_nbgy, int *rec_nedy, 
                         int *rec_nbgz, int *rec_nedz, int *rec_nxt, int *rec_nyt,
                         int *rec_nzt, MPI_Offset *displacement,
                         long int nxt, long int nyt, long int nzt, 
                         int rec_NX, int rec_NY, int rec_NZ,int NBGX, int NEDX, int NSKPX,
                         int NBGY, int NEDY, int NSKPY, int NBGZ, int NEDZ,int NSKPZ,int *coord);
unsigned int *Alloc1PU(unsigned int nx);
int* Alloc1P(int nx);
void write_metadata(const char *OUT, const char *fn, unsigned int *disp_xyzw, 
                    int max_send_size,int RNX, int RNY, int RNZ, int WRITE_STEP, int size);
void write_model(int rank, long int cur_step, char *filenamebase, float *buf, int max_send_size,     float *snapshots, int snapshots_size, MPI_Comm MCW);
//void load_ckpt(int *global_step, int PX, int PY, int PZ, int rank, int n1, int n2, int n3,
//               int WSIZE, int MSIZE, int CUR_CKPT_IDX, int NUM_CKPT_DIR,
//               const float *W, const float *hW, const float *mW,
//               const float *tW, const float *M, MPI_Comm MCW);
//void save_ckpt(int global_step, int PX, int PY, int PZ, int rank, int n1, int n2, int n3,
//               int WSIZE, int MSIZE, int CUR_CKPT_IDX, int NUM_CKPT_DIR, const float *W,
//               const float *hW, const float *mW, const float *tW, const float *M, MPI_Comm MCW);
// void load_ckpt(int *global_step, int PX, int PY, int PZ, int rank, int n1, int n2, int n3,
//                int WSIZE, int MSIZE, int CUR_CKPT_IDX, int NUM_CKPT_DIR,
//                const float *W, MPI_Comm MCW);
// void save_ckpt(int global_step, int PX, int PY, int PZ, int rank, int n1, int n2, int n3,
//                int WSIZE, int MSIZE, int CUR_CKPT_IDX, int NUM_CKPT_DIR, const float *W,
//                MPI_Comm MCW);
void load_ckpt(int *global_step, int PX, int PY, int PZ, int rank, int n1, int n2, int n3,
               int WSIZE, int MSIZE, int CUR_CKPT_IDX, int NUM_CKPT_DIR,
               const float *W_8, const float *W_1, MPI_Comm MCW);
void save_ckpt(int global_step, int PX, int PY, int PZ, int rank, int n1, int n2, int n3,
               int WSIZE, int MSIZE, int CUR_CKPT_IDX, int NUM_CKPT_DIR, const float *W_8, const float *W_1,
               MPI_Comm MCW);


#endif
