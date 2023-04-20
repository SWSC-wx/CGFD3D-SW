#ifndef MACDRP_H
#define MACDRP_H

#define PI 3.141592653589793238463

typedef float *** Grid3D;
typedef float ** Grid2D;
typedef float * Grid1D;

struct aux {
    // 1~ND
    float *Wx1, *hWx1, *mWx1, *tWx1;
    float *Wy1, *hWy1, *mWy1, *tWy1;
    float *Wz1, *hWz1, *mWz1, *tWz1;
    // nx-ND:nx
    float *Wx2, *hWx2, *mWx2, *tWx2;
    float *Wy2, *hWy2, *mWy2, *tWy2;
    float *Wz2, *hWz2, *mWz2, *tWz2;
};

typedef struct {
	// float *wave;
	// float *mwave;
	// float *hwave;
	// float *twave;
	float *wave_8;
	float *mwave_8;
	float *hwave_8;
	float *twave_8;
	float *wave_1;
	float *mwave_1;
	float *hwave_1;
	float *twave_1;
	float *metric;
	float *wp_tmp;
    float *pml;
	float DH;
	float DT;
	float rka;
	float rkb;
	int istep;
	int nk1;
	int dim_x;
	int dim_y;
	int dim_z;
	int ny;
	int nz;
	int lnz;
	int freenode;
	int flag;
	int FlagX;
	int FlagY;
	int FlagZ;
    float Qsa;
    int isx1;
    int isx2;
    int isy1;
    int isy2;
    int isz1;
    int isz2;
    struct aux *Aux;
    float *matVx2Vz;
    float *matVy2Vz;
    int rk_step;
} Param_cal;

typedef struct {
	float *source;
	float *des;
	int len_s_x;
	int len_s_y;
	int len_s_z;
	int len_d_x;
	int len_d_y;
	int len_d_z;
	int wx;
	int wy;
	int wz;
	int SIZE;
} Param_swap;
//void command( int argc, char **argv ); 
void get_params( int argc, char **argv ); 



//extern int   vecF[5], vecB[5];

extern float coeF[5], coeB[5];
extern float coeF_DH[5], coeB_DH[5];
extern float coe24F[3], coe24B[3];
extern float coe24F_DH[3], coe24B_DH[3];
extern float coe22F[2], coe22B[2];

extern float RK4a[3];
extern float RK4b[4];

extern float RK6a[5];
extern float RK6b[6];

int exchange_Wave(float *);
//int RK_Syn(float *, float *, int istep);
//int RK_Syn(float *, float *, float *, int istep);
#ifdef ARCH_SW
// int RK_Syn(float *W, float *M, float *hW, float *mW, float *tW, float *matVx2Vz, float *matVy2Vz, int istep, float *pml, int isx1, int isx2, int isy1, int isy2, int isz1, int isz2, struct aux* Aux);
// int RK_Syn(float *W, float *W_8, float *W_1, float *M, float *hW, float *hW_8, float *hW_1, float *mW, float *mW_8, float *mW_1, float *tW, float *tW_8, float *tW_1, float *matVx2Vz, float *matVy2Vz, int istep, float *pml, int isx1, int isx2, int isy1, int isy2, int isz1, int isz2, struct aux* Aux);
int RK_Syn(float *W_8, float *W_1, float *M, float *hW_8, float *hW_1, float *mW_8, float *mW_1, float *tW_8, float *tW_1, float *matVx2Vz, float *matVy2Vz, int istep, float *pml, int isx1, int isx2, int isy1, int isy2, int isz1, int isz2, struct aux* Aux);
#elif ARCH_x86
int RK_Syn(float *W, float *M, float *matVx2Vz, float *matVy2Vz, int istep, float *pml, int isx1, int isx2, int isy1, int isy2, int isz1, int isz2, struct aux *Aux, float DH, int freenode);
int coef_fdxy2fdz(float *M, Grid3D matVx2Vz, Grid3D matVy2Vz, Grid3D matF2Vz, int nii1, int nii2, int njj1, int njj2);
int coef_surface(float *M, float *matVx2Vz, float *matVy2Vz, float *pml);
void invert3x3(float m[][3]);
void matmul3x3(float A[][3], float B[][3], float C[][3]);
#elif ARCH_MPE
int RK_Syn(float *W, float *M, float *matVx2Vz, float *matVy2Vz, int istep, float *pml, int isx1, int isx2, int isy1, int isy2, int isz1, int isz2, struct aux *Aux, float DH, int freenode);
int coef_fdxy2fdz(float *M, Grid3D matVx2Vz, Grid3D matVy2Vz, Grid3D matF2Vz, int nii1, int nii2, int njj1, int njj2);
void invert3x3(float m[][3]);
void matmul3x3(float A[][3], float B[][3], float C[][3]);
#endif

void calculate_kernel_slave(Param_cal *param);
// void calculate_kernel_slave_xz(Param_cal *param);
// void calculate_kernel_slave_l1(Param_cal *param);
// void calculate_kernel_slave_l2(Param_cal *param);
void swap_kernel_slave(Param_swap *param);

// add FLAG to switch BWD and FWD
// FWD: Forward BWD: Backward
#define BWD -1
#define FWD 1

#define vec_L(var,i,FLAG) ((FLAG==1)?vec_LF(var,i):vec_LB(var,i))
#define vec_LF(var,i) (coeF_DH[0]*var[i-1]+coeF_DH[1]*var[i]+coeF_DH[2]*var[i+1]+coeF_DH[3]*var[i+2]+coeF_DH[4]*var[i+3])
#define vec_LB(var,i) (coeB_DH[0]*var[i-3]+coeB_DH[1]*var[i-2]+coeB_DH[2]*var[i-1]+coeB_DH[3]*var[i]+coeB_DH[4]*var[i+1])


#define L(var,idx,stride,FLAG) ((FLAG==1)?LF(var,idx,stride):LB(var,idx, stride))
#define LF(var,idx,stride) (coeF_DH[0]*var[idx - stride]+coeF_DH[1]*var[idx]+coeF_DH[2]*var[idx + stride]+coeF_DH[3]*var[idx + 2 * stride]+coeF_DH[4]*var[idx + 3 * stride])
#define LB(var,idx,stride) (coeB_DH[0]*var[idx - 3 * stride]+coeB_DH[1]*var[idx - 2 * stride]+coeB_DH[2]*var[idx - stride]+coeB_DH[3]*var[idx]+coeB_DH[4]*var[idx + stride])

#define L_SIMD(var,idx,stride,FLAG) {               \
    if(FLAG == 1){                                  \
        simd_load(va, &var[idx - stride]);          \
        va = va * (-0.30874f);                      \
        simd_load(vb, &var[idx]);                   \
        vb = vb * (-0.6326f);                       \
        simd_load(vc, &var[idx + stride]);          \
        vc = vc * (1.2330f);                        \
        simd_load(vd, &var[idx + 2*stride]);        \
        vd = vd * (-0.3334f);                       \
        simd_load(ve, &var[idx + 3*stride]);        \
        ve = ve * (0.04168f);                       \
    }                                               \
    else{                                           \
        simd_load(va, &var[idx - 3*stride]);        \
        va = va * (-0.04168f);                      \
        simd_load(vb, &var[idx - 2*stride]);        \
        vb = vb * (0.3334f);                        \
        simd_load(vc, &var[idx - stride]);          \
        vc = vc * (-1.2330f);                       \
        simd_load(vd, &var[idx]);                   \
        vd = vd * (0.6326f);                        \
        simd_load(ve, &var[idx + stride]);          \
        ve = ve * (0.30874f);                       \
    }                                               \
    va = va + vb;                                   \
    vc = vc + vd;                                   \
    va = va + vc + ve;                              \
    va = va * rDH;                                  \
    simd_store(va, temp);                           \
}


#define L22(var,idx,stride,FLAG) ((FLAG==1)?L22F(var,idx,stride):L22B(var,idx,stride))
#define L22F(var,idx,stride) (var[idx + stride]-var[idx])*rDH
#define L22B(var,idx,stride) (var[idx]-var[idx - stride])*rDH

#define L24(var,idx,stride,FLAG) ((FLAG==1)?L24F(var,idx,stride):L24B(var,idx,stride))
#define L24F(var,idx,stride) (coe24F_DH[0]*var[idx]+coe24F_DH[1]*var[idx + stride]+coe24F_DH[2]*var[idx + 2 * stride])
#define L24B(var,idx,stride) (coe24B_DH[0]*var[idx - 2 * stride]+coe24B_DH[1]*var[idx - stride]+coe24B_DH[2]*var[idx])


#define LxF(var,i,j,k) (coeF_DH[0]*var[i-1][j][k]+coeF_DH[1]*var[i][j][k]+coeF_DH[2]*var[i+1][j][k]+coeF_DH[3]*var[i+2][j][k]+coeF_DH[4]*var[i+3][j][k])
#define LxB(var,i,j,k) (coeB_DH[0]*var[i-3][j][k]+coeB_DH[1]*var[i-2][j][k]+coeB_DH[2]*var[i-1][j][k]+coeB_DH[3]*var[i][j][k]+coeB_DH[4]*var[i+1][j][k])

#define LyF(var,i,j,k) (coeF_DH[0]*var[i][j-1][k]+coeF_DH[1]*var[i][j][k]+coeF_DH[2]*var[i][j+1][k]+coeF_DH[3]*var[i][j+2][k]+coeF_DH[4]*var[i][j+3][k])
#define LyB(var,i,j,k) (coeB_DH[0]*var[i][j-3][k]+coeB_DH[1]*var[i][j-2][k]+coeB_DH[2]*var[i][j-1][k]+coeB_DH[3]*var[i][j][k]+coeB_DH[4]*var[i][j+1][k])

#define LzF(var,i,j,k) (coeF_DH[0]*var[i][j][k-1]+coeF_DH[1]*var[i][j][k]+coeF_DH[2]*var[i][j][k+1]+coeF_DH[3]*var[i][j][k+2]+coeF_DH[4]*var[i][j][k+3])
#define LzB(var,i,j,k) (coeB_DH[0]*var[i][j][k-3]+coeB_DH[1]*var[i][j][k-2]+coeB_DH[2]*var[i][j][k-1]+coeB_DH[3]*var[i][j][k]+coeB_DH[4]*var[i][j][k+1])

#define vec_LF(var,i) (coeF_DH[0]*var[i-1]+coeF_DH[1]*var[i]+coeF_DH[2]*var[i+1]+coeF_DH[3]*var[i+2]+coeF_DH[4]*var[i+3])
#define vec_LB(var,i) (coeB_DH[0]*var[i-3]+coeB_DH[1]*var[i-2]+coeB_DH[2]*var[i-1]+coeB_DH[3]*var[i]+coeB_DH[4]*var[i+1])

#define L24xF(var,i,j,k) (coe24F_DH[0]*var[i][j][k]+coe24F_DH[1]*var[i+1][j][k]+coe24F_DH[2]*var[i+2][j][k])
#define L24xB(var,i,j,k) (coe24B_DH[0]*var[i-2][j][k]+coe24B_DH[1]*var[i-1][j][k]+coe24B_DH[2]*var[i][j][k])

#define L24yF(var,i,j,k) (coe24F_DH[0]*var[i][j][k]+coe24F_DH[1]*var[i][j+1][k]+coe24F_DH[2]*var[i][j+2][k])
#define L24yB(var,i,j,k) (coe24B_DH[0]*var[i][j-2][k]+coe24B_DH[1]*var[i][j-1][k]+coe24B_DH[2]*var[i][j][k])

#define L24zF(var,i,j,k) (coe24F_DH[0]*var[i][j][k]+coe24F_DH[1]*var[i][j][k+1]+coe24F_DH[2]*var[i][j][k+2])
#define L24zB(var,i,j,k) (coe24B_DH[0]*var[i][j][k-2]+coe24B_DH[1]*var[i][j][k-1]+coe24B_DH[2]*var[i][j][k])

#define L22xF(var,i,j,k) (var[i+1][j][k]-var[i][j][k])*rDH
#define L22xB(var,i,j,k) (var[i][j][k]-var[i-1][j][k])*rDH

#define L22yF(var,i,j,k) (var[i][j+1][k]-var[i][j][k])*rDH
#define L22yB(var,i,j,k) (var[i][j][k]-var[i][j-1][k])*rDH

#define L22zF(var,i,j,k) (var[i][j][k+1]-var[i][j][k])*rDH
#define L22zB(var,i,j,k) (var[i][j][k]-var[i][j][k-1])*rDH

#define Lx(var,i,j,k,FLAG) ((FLAG==1)?LxF(var,i,j,k):LxB(var,i,j,k))
#define Ly(var,i,j,k,FLAG) ((FLAG==1)?LyF(var,i,j,k):LyB(var,i,j,k))
#define Lz(var,i,j,k,FLAG) ((FLAG==1)?LzF(var,i,j,k):LzB(var,i,j,k))

#define vec_L(var,i,FLAG) ((FLAG==1)?vec_LF(var,i):vec_LB(var,i))

#define L24x(var,i,j,k,FLAG) ((FLAG==1)?L24xF(var,i,j,k):L24xB(var,i,j,k))
#define L24y(var,i,j,k,FLAG) ((FLAG==1)?L24yF(var,i,j,k):L24yB(var,i,j,k))
#define L24z(var,i,j,k,FLAG) ((FLAG==1)?L24zF(var,i,j,k):L24zB(var,i,j,k))

#define L22x(var,i,j,k,FLAG) ((FLAG==1)?L22xF(var,i,j,k):L22xB(var,i,j,k))
#define L22y(var,i,j,k,FLAG) ((FLAG==1)?L22yF(var,i,j,k):L22yB(var,i,j,k))
#define L22z(var,i,j,k,FLAG) ((FLAG==1)?L22zF(var,i,j,k):L22zB(var,i,j,k))

#define vec_LL(var,i,FLAG) ((FLAG==1)?vec_LLF(var,i):vec_LLB(var,i))
#define vec_LLF(var,i) ((-0.30874*rDH)*var[i]+(-0.6326*rDH)*var[i+1]+(1.2330*rDH)*var[i+2]+(-0.3334*rDH)*var[i+3]+(0.04168*rDH)*var[i+4])
#define vec_LLB(var,i) ((-0.04168 * rDH)*var[i]+(0.3334*rDH)*var[i+1]+(-1.2330*rDH)*var[i+2]+(0.6326*rDH)*var[i+3]+(0.30874*rDH)*var[i+4])

#endif
