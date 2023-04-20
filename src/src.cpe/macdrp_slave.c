/*
 * macdrp_slave_dma.c
 *
 *  Created on: Feb 5, 2018
 *      Author: bingo
 */

#ifdef ARCH_SW

#include <slave.h>
#include <simd.h>
#include "../macdrp.h"
#include "../pml_macro.h"
#define WSIZE 9
#define WSIZE_V 8
#define MSIZE 14
#define PSIZE 6

#define wx 1
#define wx_c 3     // Here, the wx_c is the number of plane advanced in a slave core
#define wy 5
#define DX 10

#define pe_get(mem, ldm, size, reply) athread_get(PE_MODE, mem, ldm, size, (void*)(reply), 0, 0, 0)                                                       
#define pe_put(mem, ldm, size, reply) athread_put(PE_MODE, ldm, mem, size, (void*)(reply), 0, 0)

static const double ln2 = 0.6931471805599453;
static inline double my_exp_ln2(const double x) {
  double a1 = -0.9999999995;
  double a2 = 0.4999999206;
  double a3 = -0.1666653019;
  double a4 = 0.0416573475;
  double a5 = -0.0083013598;
  double a6 = 0.0013298820;
  double a7 = -0.0001413161;
  double exp_x = 1.0;
  exp_x = a7 * x + a6;
  exp_x = exp_x * x + a5;
  exp_x = exp_x * x + a4;
  exp_x = exp_x * x + a3;
  exp_x = exp_x * x + a2;
  exp_x = exp_x * x + a1;
  exp_x = exp_x * x + 1.0;
  return 1.0 / exp_x;
}

static inline double my_exp_pos(const double x) {
  long long k, twok;
  double x_;
  k = x / ln2;
  twok = 1 << k;
  x_ = x - (double) k * ln2;
  return (double)twok * my_exp_ln2(x_);
}

static inline double myexp(const double x) {
  if(x >= 0.0) {
    return my_exp_pos(x);
  } else {
    return 1.0 / my_exp_pos(-x);
  }
}

inline  addn(int x, int n, int x1) {
	return x + n < x1 ? x + n : x + n - x1;
}


inline void invert3x3(float m[][3]){
  float inv[3][3];
  float det;
  int i, j;

  inv[0][0] = m[1][1]*m[2][2] - m[2][1]*m[1][2];
  inv[0][1] = m[2][1]*m[0][2] - m[0][1]*m[2][2];
  inv[0][2] = m[0][1]*m[1][2] - m[0][2]*m[1][1];
  inv[1][0] = m[1][2]*m[2][0] - m[1][0]*m[2][2];
  inv[1][1] = m[0][0]*m[2][2] - m[2][0]*m[0][2];
  inv[1][2] = m[1][0]*m[0][2] - m[0][0]*m[1][2];
  inv[2][0] = m[1][0]*m[2][1] - m[1][1]*m[2][0];
  inv[2][1] = m[2][0]*m[0][1] - m[0][0]*m[2][1];
  inv[2][2] = m[0][0]*m[1][1] - m[0][1]*m[1][0];

  det = inv[0][0] * m[0][0]
      + inv[0][1] * m[1][0]
      + inv[0][2] * m[2][0];

  det = 1.0f / det;

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      m[i][j] = inv[i][j] * det;
}


inline void matmul3x3(float A[][3], float B[][3], float C[][3]){

  int i, j, k;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++){
      C[i][j] = 0.0;
      for (k = 0; k < 3; k++)
        C[i][j] += A[i][k] * B[k][j];
    }
}
extern int this_rank;
void calculate_kernel_slave(Param_cal *param) {

	int id = athread_get_id(-1);
	volatile unsigned long get_reply, put_reply;
	int plane_comp = 0;

    float *wave_8 = param->wave_8;
	float *mwave_8 = param->mwave_8;
	float *hwave_8 = param->hwave_8;
	float *twave_8 = param->twave_8;
    float *wave_1 = param->wave_1;
	float *mwave_1 = param->mwave_1;
	float *hwave_1 = param->hwave_1;
	float *twave_1 = param->twave_1;
	float *metric = param->metric;
	float *wp_xzs0 = param->wp_tmp;
	float DH = param->DH;
	float DT = param->DT;
	float rka = param->rka;
	float rkb = param->rkb;
	int istep = param->istep;
	int flag = param->flag;
	int nk1 = param->nk1;
	int dim_x = param->dim_x;
	int dim_y = param->dim_y;
	int dim_z = param->dim_z;
	int ny = param->ny;
	int nz = param->nz;
	int lnz = param->lnz;
	int freenode = param->freenode;
	int FlagX = param->FlagX;
	int FlagY = param->FlagY;
	int FlagZ = param->FlagZ;
    float Qsa = param->Qsa;
    float *matVx2Vz_m = param->matVx2Vz;
    float *matVy2Vz_m = param->matVy2Vz;
    int rk_step =  param->rk_step;
#ifdef usePML
    float* pml = param->pml;
    int isx1 = param->isx1;
    int isx2 = param->isx2;
    int isy1 = param->isy1;
    int isy2 = param->isy2;
    int isz1 = param->isz1;
    int isz2 = param->isz2;
    struct aux* Aux = param->Aux;
#endif


	int wz = 25;
	wz = (dim_z <= wz) ? dim_z : wz;                    // consider the deserved data block smaller than obtained block
	const int MX = 1;
	const int MY = 8;  // block of row cores 
	const int MZ = 8;
	int NX = (dim_x + wx_c * MX - 1) / (wx_c * MX);    // Here, the wx_c is the number of plane advanced in a slave core
    int NY = (dim_y + wy * MY - 1) / (wy * MY);
    int	NZ = (dim_z + wz * MZ - 1) / (wz * MZ);

	int iix, iiy, iiz, ix, iy, iz, izbeg, izend, izend_out, iybeg, iyend, ixbeg, ixend, ixn, iyn, izn, izn_final, n, l, i, ii, j, jj, k, kk, n_val;
	int pos, pos_m1, pos_m2, pos_m3, pos_p1, pos_p2, pos_p3, pos_hw, pos_m, pos_mm;
    int tpos8, tpos1, tpos_hw8, tpos_hw1; //hd

	float mu, lam, lam2mu, rrho, rrhojac;
	float vecTx[5], vecTy[5], vecTz_tmp[7];
	int tmp_size = (dim_y%wy == 0) ? 1 : dim_y%wy;
	float vecTz_tmp1[tmp_size], vecTz_tmp2[tmp_size], vecTz_tmp3[tmp_size], vecTz_tmp4[tmp_size], vecTz_tmp5[tmp_size], vecTz_tmp6[tmp_size];
	float *vecTz = vecTz_tmp + 2;
	float DxTx, DyTy, DzTz;

	float DxTxx, DxTyy, DxTzz, DxTxy, DxTxz, DxTyz;
	float DyTxx, DyTyy, DyTzz, DyTxy, DyTxz, DyTyz;
	float DzTxx, DzTyy, DzTzz, DzTxy, DzTxz, DzTyz;
	float DxVx, DxVy, DxVz;
	float DyVx, DyVy, DyVz;
	float DzVx, DzVy, DzVz;
    // float DzVx, DzVy, DzVz;
    float DzVx1, DzVy1, DzVz1;
    float rb1, rb2, rb3;
    float temp[WSIZE_V] __attribute__((aligned(32))); //hd
    // memset(temp, 0, sizeof(float) * WSIZE_V);
    floatv8 va, vb, vc, vd, ve; //hd
    volatile int reply_w_8 = 0, reply_w_1 = 0, reply_m = 0, reply_mw_8 = 0, reply_mw_1 = 0, reply_tw_8 = 0, reply_tw_1 = 0;
    volatile int get_aux_mw = 0, put_aux_mw = 0, get_aux_tw = 0, put_aux_tw = 0, put_aux_hw = 0, get_pml_x = 0, get_pml_y = 0, get_pml_z = 0, get_mat = 0;
    volatile int get_hw = 0, put_hw = 0;
    volatile int iget_reply = 0, iput_reply = 0;
    volatile int iget_reply_pml = 0, iput_reply_pml = 0;
    // volatile int get_hw = 0, get_metric_l5 = 0, get_metric = 0, get_wave = 0, put_hw = 0;

    int iii, jjj, pos_a;
    int idx0;
    float *Aux_W, *Aux_hW, *Aux_tW, *Aux_mW, *Aux_Wx, *Aux_hWx, *Aux_Wy, *Aux_hWy, *Aux_Wz, *Aux_hWz;
    float *Aux_W2, *Aux_hW2, *Aux_tW2, *Aux_mW2;
    float *Aux_W3, *Aux_hW3, *Aux_tW3, *Aux_mW3;
    float b1, d1, ad1;

	rka = rka * DT; rkb = rkb * DT;


	int x0 = (FlagX == 1) ? 1 : 3;
	int x1 = (FlagX == 1) ? 3 : 1;
	int y0 = (FlagY == 1) ? 1 : 3;
	int y1 = (FlagY == 1) ? 3 : 1;
	int z0 = (FlagZ == 1) ? 1 : 3;
	int z1 = (FlagZ == 1) ? 3 : 1;

    int plane_last = x0 + wx + x1; //hd

	float A[3][3];
	float rDH = 1.0 / DH;


    float wave_8_l[(x0 + wx + x1 + 1) * (y0 + wy + y1) * (z0 + wz + z1) * WSIZE_V] __attribute__((aligned(32)));  //hd
    float wave_1_l[(x0 + wx + x1 + 1) * (y0 + wy + y1) * (z0 + wz + z1)];  //hd
    float metric_l5[(x0 + wx + x1 + 1) * (y0 + wy + y1) * (z0 + wz + z1) * MSIZE];  //hd
    float hwave_8_l[wz * WSIZE_V] __attribute__((aligned(32)));
    float hwave_1_l[wz];

    int xstep_l_w_8 = (y0 + wy + y1) * (z0 + wz + z1) * WSIZE_V;
	int ystep_l_w_8 = (z0 + wz + z1) * WSIZE_V;
    int xstep_l_hw_8 = wy * wz * WSIZE_V;
    int ystep_l_hw_8 = wz * WSIZE_V;
    int xstep_l_w_1 = (y0 + wy + y1) * (z0 + wz + z1);
	int ystep_l_w_1 = (z0 + wz + z1);
    int xstep_l_hw_1 = wy * wz;
    int ystep_l_hw_1 = wz;
    int xstep_l_m5 = (y0 + wy + y1) * (z0 + wz + z1) * MSIZE;
    int ystep_l_m5 = (z0 + wz + z1) * MSIZE;
    int xstep_l_m = wy * wz * MSIZE;
    int ystep_l_m = wz * MSIZE;
    int xstep_w_8 = ny * lnz * WSIZE_V;
    int ystep_w_8 = lnz * WSIZE_V;
    int xstep_w_1 = ny * lnz;
    int ystep_w_1 = lnz;
    int xstep_m = ny * lnz * MSIZE;
    int ystep_m = lnz * MSIZE;
    int segment_8 = ystep_l_w_8;
    int segment_1 = ystep_l_w_1;

    float *metric_s, * metric_s_tmp;
    float *metric_l5_tmp, *metric_l_tmp;
    float *wave_8_s, *hwave_8_s, *mwave_8_s, *twave_8_s, *wave_8_s_tmp, *hwave_8_s_tmp, *mwave_8_s_tmp, *twave_8_s_tmp, *wave_1_s, *hwave_1_s, *mwave_1_s, *twave_1_s, *wave_1_s_tmp, *hwave_1_s_tmp, *mwave_1_s_tmp, *twave_1_s_tmp;
    float *wave_8_l_tmp, *wave_1_l_tmp;

    float *wave_8_c = wave_8_l + y0 * ystep_l_w_8 + z0 * WSIZE_V;
    float *wave_1_c = wave_1_l + y0 * ystep_l_w_1 + z0;
    float *metric_c5 = metric_l5 + y0 * ystep_l_m5 + z0 * MSIZE;

    int ni1 = 3;
	int nj1 = 3;
    
	int xid = id / (MY * MZ);
    int yid = (id - xid * (MY * MZ)) / MZ;
    int zid = (id - xid * (MY * MZ)) % MZ;

    float ax[DX], bx[DX], dx[DX];
    float ay[wy], by[wy], dy[wy];
    float az[wz], bz[wz], dz[wz];
    int ni = dim_x;
    int nj = dim_y;
    int nk = dim_z;
    // volatile int ax_reply = 0;

    for (iiz = 0; iiz < NZ; iiz++) {
        izbeg = dim_z - iiz*MZ*wz - wz*(zid + 1);
        izend = dim_z - iiz*MZ*wz - wz*zid;
        izbeg = izbeg < 0 ? 0 : izbeg;
        izn_final = izend - izbeg; // contain the 3 freesurface layer
        izend_out = izend;
        izend = (freenode && (izend >= (nz - nk1 - 6))) ? (dim_z - 3) : izend;
        izn = izend - izbeg; // if in freesurace layer izn is izn_final - 3
#ifdef usePML
        //get pml z 
        iget_reply_pml = 0;

        int *p2 = pml + ni * 3 + nj * 3;
        athread_dma_iget(az, p2 + izbeg, sizeof(float)*wz, &iget_reply_pml);
        athread_dma_iget(bz, p2 + nk + izbeg, sizeof(float)*wz, &iget_reply_pml);
        athread_dma_iget(dz, p2 + nk*2 + izbeg, sizeof(float)*wz, &iget_reply_pml);
        athread_dma_wait_value(&iget_reply_pml, 3);

#endif
		for (iiy = 0; iiy < NY; iiy++) {
			iybeg = wy * MY * iiy + wy * yid;
			iyend = wy * MY * iiy + wy * (yid + 1);
			iyend = iyend < dim_y ? iyend : dim_y;
			iyn = iyend - iybeg;
#ifdef usePML
            //get pml y
            iget_reply_pml = 0;

            int *p1 = pml + ni * 3;
            athread_dma_iget(ay, p1 + iybeg, sizeof(float)*wy, &iget_reply_pml);
            athread_dma_iget(by, p1 + nj + iybeg, sizeof(float)*wy, &iget_reply_pml);
            athread_dma_iget(dy, p1 + nj*2 + iybeg, sizeof(float)*wy, &iget_reply_pml);
            athread_dma_wait_value(&iget_reply_pml, 3);
#endif
			// set computing block position data source main memory
            wave_8_s = wave_8 + iybeg*ystep_w_8 + izbeg*WSIZE_V;
			hwave_8_s = hwave_8 + iybeg*ystep_w_8 + izbeg*WSIZE_V;
            wave_1_s = wave_1 + iybeg*ystep_w_1 + izbeg;
			hwave_1_s = hwave_1 + iybeg*ystep_w_1 + izbeg;
			metric_s = metric + iybeg*ystep_m + izbeg*MSIZE;
            twave_8_s = twave_8 + iybeg*ystep_w_8 + izbeg*WSIZE_V;
			mwave_8_s = mwave_8 + iybeg*ystep_w_8 + izbeg*WSIZE_V;
            twave_1_s = twave_1 + iybeg*ystep_w_1 + izbeg;
			mwave_1_s = mwave_1 + iybeg*ystep_w_1 + izbeg;

	        // set local computing block memory position
            wave_8_l_tmp = wave_8_l;
            wave_1_l_tmp = wave_1_l;

 			for (ix = -x0; ix < wx + x1; ix++) {
                wave_8_s_tmp = wave_8_s + ix * xstep_w_8 - y0 * ystep_w_8 - z0 * WSIZE_V;
                wave_1_s_tmp = wave_1_s + ix * xstep_w_1 - y0 * ystep_w_1 - z0;
				for (iy = -y0; iy < wy + y1; iy++) {
                    athread_dma_get(wave_8_l_tmp, wave_8_s_tmp, sizeof(float)*(wz + z0 + z1)*WSIZE_V);
					wave_8_s_tmp += ystep_w_8;
					wave_8_l_tmp += ystep_l_w_8;

                    athread_dma_get(wave_1_l_tmp, wave_1_s_tmp, sizeof(float)*(wz + z0 + z1));
					wave_1_s_tmp += ystep_w_1;
					wave_1_l_tmp += ystep_l_w_1;
				}
			}

			if (freenode && (izend >= (nz - nk1 - 6))) {
				metric_l5_tmp = metric_l5;
				for (ix = -x0; ix < wx + x1; ix++) {
					metric_s_tmp = metric_s + ix * xstep_m - y0 * ystep_m - z0 * MSIZE;
					for (iy = -y0; iy < wy + y1; iy++) {
                        athread_dma_get(metric_l5_tmp, metric_s_tmp, sizeof(float)*(z0 + wz + z1)*MSIZE);
						metric_s_tmp += ystep_m;
						metric_l5_tmp += ystep_l_m5;
					}
				}
			}
            
			/************************ compute hwave from ix plane to ix++ plane *******************************/
			plane_comp = 0;
            plane_last = x0 + wx + x1; //hd
			for (ix = 0; ix < dim_x; ix++) {
#ifdef usePML
                if(ix % DX == 0) {
                    //get pml x
                    iget_reply_pml = 0;

                    int *p0 = pml;
                    athread_dma_iget(ax, p0 + ix, sizeof(float)*DX, &iget_reply_pml);
                    athread_dma_iget(bx, p0 + ni + ix, sizeof(float)*DX, &iget_reply_pml);
                    athread_dma_iget(dx, p0 + ni*2 + ix, sizeof(float)*DX, &iget_reply_pml);
                    athread_dma_wait_value(&iget_reply_pml, 3);
                }
                rb1 = 1.0f/bx[ix % DX];
#else
                rb1 = 1.0f;
#endif

               /* get next plane wave date */
                wave_8_s_tmp = wave_8_s + (ix + wx + x1) * xstep_w_8 - y0 * ystep_w_8 - z0 * WSIZE_V;
                wave_8_l_tmp = wave_8_l + plane_last * xstep_l_w_8;  //hd
                wave_1_s_tmp = wave_1_s + (ix + wx + x1) * xstep_w_1 - y0 * ystep_w_1 - z0;
                wave_1_l_tmp = wave_1_l + plane_last * xstep_l_w_1;  //hd
                reply_w_8 = 0;
                reply_w_1 = 0;
				for (iy = -y0; iy < wy + y1; iy++) {
                    athread_dma_iget(wave_8_l_tmp, wave_8_s_tmp, sizeof(float)*(z0 + wz + z1)*WSIZE_V, &reply_w_8);
					wave_8_s_tmp += ystep_w_8;
					wave_8_l_tmp += ystep_l_w_8;

                    athread_dma_iget(wave_1_l_tmp, wave_1_s_tmp, sizeof(float)*(z0 + wz + z1), &reply_w_1);
					wave_1_s_tmp += ystep_w_1;
					wave_1_l_tmp += ystep_l_w_1;
				}
                
				/* get next plane metric data */
				if (freenode && (izend >= (nz - nk1 - 6))) {
					metric_s_tmp = metric_s + (ix + wx + x1) * xstep_m - y0 * ystep_m - z0 * MSIZE;
                    metric_l5_tmp = metric_l5 + plane_last * xstep_l_m5;  //hd
                    reply_m = 0;
					for (iy = -y0; iy < wy + y1; iy++) {
                        athread_dma_iget(metric_l5_tmp, metric_s_tmp, sizeof(float)*(z0 + wz + z1)*MSIZE, &reply_m);
						metric_s_tmp += ystep_m;
						metric_l5_tmp += ystep_l_m5;
					}
				}

				// To obtain metric data
				metric_s_tmp = metric_s + ix * xstep_m;

                hwave_8_s_tmp = hwave_8_s + ix*xstep_w_8;
                hwave_1_s_tmp = hwave_1_s + ix*xstep_w_1;

				// To obtain twave data
                twave_8_s_tmp = twave_8_s + ix*xstep_w_8;
                twave_1_s_tmp = twave_1_s + ix*xstep_w_1;

				// To obtain mwave data
                mwave_8_s_tmp = mwave_8_s + ix * xstep_w_8;
                mwave_1_s_tmp = mwave_1_s + ix * xstep_w_1;

               //hd
                int plane1 = plane_comp;
				int plane2 = addn(plane_comp, 1, wx + x0 + x1 + 1);
				int plane3 = addn(plane_comp, 2, wx + x0 + x1 + 1);
				int plane4 = addn(plane_comp, 3, wx + x0 + x1 + 1);
				int plane5 = addn(plane_comp, 4, wx + x0 + x1 + 1);
				int plane = (FlagX == 1) ? plane2 : plane4;

				/********************************************* Calculate step 1  *******************************************/
				for (iy = 0; iy < iyn; iy++) {
#ifdef usePML
                    rb2 = 1.0f/by[iy];
                    int iput_count = 0;
                    iput_reply_pml = 0;
#else
                    rb2 = 1.0f;
#endif

                    float mwave_8_l[wz*WSIZE_V] __attribute__((aligned(32)));
                    float twave_8_l[wz*WSIZE_V] __attribute__((aligned(32)));
                    float mwave_1_l[wz];
				    float twave_1_l[wz];
				    float metric_l[wz*MSIZE];

                    reply_tw_8 = 0; 
                    athread_dma_iget(twave_8_l, twave_8_s_tmp, sizeof(float)*wz*WSIZE_V, &reply_tw_8);
                    reply_tw_1 = 0;
                    athread_dma_iget(twave_1_l, twave_1_s_tmp, sizeof(float)*wz, &reply_tw_1);
                    
                    reply_mw_8 = 0;
                    athread_dma_iget(mwave_8_l, mwave_8_s_tmp, sizeof(float)*wz*WSIZE_V, &reply_mw_8);
                    reply_mw_1 = 0;
                    athread_dma_iget(mwave_1_l, mwave_1_s_tmp, sizeof(float)*wz, &reply_mw_1);

                    athread_dma_get(metric_l, metric_s_tmp, sizeof(float)*wz*MSIZE);

					// get matVx2Vz matVy2Vz
                    float matVx2Vz[WSIZE];
                    float matVy2Vz[WSIZE];
                    int ij = ((ni1+ix)*ny+(nj1+iybeg+iy))*WSIZE;
                    athread_dma_get(matVx2Vz, matVx2Vz_m + ij, sizeof(float)*WSIZE);
                    athread_dma_get(matVy2Vz, matVy2Vz_m + ij, sizeof(float)*WSIZE);

					for (iz = 0; iz < izn_final; iz++) {
#ifdef usePML
                        rb1 = 1.0f/bx[ix % DX];
                        rb3 = 1.0f/bz[iz];
#else
                        rb1 = 1.0f;
                        rb3 = 1.0f;
#endif
                        int pos = plane * xstep_l_w_8 + iy * ystep_l_w_8 + iz * WSIZE_V;
						int pos_hw = iz * WSIZE_V;
						int pos_m = iz * MSIZE;
						int pos_p1 = plane1 * xstep_l_w_8 + iy * ystep_l_w_8 + iz * WSIZE_V;
						int pos_p2 = plane2 * xstep_l_w_8 + iy * ystep_l_w_8 + iz * WSIZE_V;
						int pos_p3 = plane3 * xstep_l_w_8 + iy * ystep_l_w_8 + iz * WSIZE_V;
						int pos_p4 = plane4 * xstep_l_w_8 + iy * ystep_l_w_8 + iz * WSIZE_V;
						int pos_p5 = plane5 * xstep_l_w_8 + iy * ystep_l_w_8 + iz * WSIZE_V;

                        int tpos = plane * xstep_l_w_1 + iy * ystep_l_w_1 + iz;
						int tpos_p1 = plane1 * xstep_l_w_1 + iy * ystep_l_w_1 + iz;
						int tpos_p2 = plane2 * xstep_l_w_1 + iy * ystep_l_w_1 + iz;
						int tpos_p3 = plane3 * xstep_l_w_1 + iy * ystep_l_w_1 + iz;
						int tpos_p4 = plane4 * xstep_l_w_1 + iy * ystep_l_w_1 + iz;
						int tpos_p5 = plane5 * xstep_l_w_1 + iy * ystep_l_w_1 + iz;
                        int tpos_hw = iz;

						lam = metric_l[pos_m + 10];
						mu  = metric_l[pos_m + 11];
						lam2mu = lam + 2.0f*mu;
						rrho = 1.0f/metric_l[pos_m + 12];

                        float old = DxVx;
                        DxTyz = (FlagX == 1) ? ((-0.30874*rDH)*wave_1_c[tpos_p1] + (-0.6326*rDH)*wave_1_c[tpos_p2] + (1.2330*rDH)*wave_1_c[tpos_p3] + (-0.3334*rDH)*wave_1_c[tpos_p4] + (0.04168*rDH)*wave_1_c[tpos_p5])
				        		             : ((-0.04168 * rDH)*wave_1_c[tpos_p1] + (0.3334*rDH)*wave_1_c[tpos_p2] + (-1.2330*rDH)*wave_1_c[tpos_p3] + (0.6326*rDH)*wave_1_c[tpos_p4] + (0.30874*rDH)*wave_1_c[tpos_p5]);
				        DyTyz = L(wave_1_c,tpos, segment_1, FlagY); DzTyz = L(wave_1_c, tpos, 1, FlagZ);
                    
                        if(FlagX == 1){
                            simd_load(va, &wave_8_c[pos_p1 + 0]);
                            va = va * (-0.30874f);
                            simd_load(vb, &wave_8_c[pos_p2 + 0]);
                            vb = vb * (-0.6326f);
                            simd_load(vc, &wave_8_c[pos_p3 + 0]);
                            vc = vc * (1.2330f);
                            simd_load(vd, &wave_8_c[pos_p4 + 0]);
                            vd = vd * (-0.3334f);
                            simd_load(ve, &wave_8_c[pos_p5 + 0]);
                            ve = ve * (0.04168f);
                        }
                        else{
                            simd_load(va, &wave_8_c[pos_p1 + 0]);
                            va = va * (-0.04168f);
                            simd_load(vb, &wave_8_c[pos_p2 + 0]);
                            vb = vb * (0.3334f);
                            simd_load(vc, &wave_8_c[pos_p3 + 0]);
                            vc = vc * (-1.2330f);
                            simd_load(vd, &wave_8_c[pos_p4 + 0]);
                            vd = vd * (0.6326f);
                            simd_load(ve, &wave_8_c[pos_p5 + 0]);
                            ve = ve * (0.30874f);
                        }

                        va = va + vb;
                        vc = vc + vd;
                        va = va + vc + ve;
                        va = va * rDH;
                        simd_store(va, temp);
                        DxVx = temp[0];
                        DxVy = temp[1];  DxVz = temp[2]; DxTxx = temp[3];  DxTyy = temp[4]; DxTzz = temp[5]; DxTxy = temp[6]; DxTxz = temp[7];

                        L_SIMD(wave_8_c,pos,segment_8,FlagY);
                        DyVx = temp[0];  DyVy = temp[1];  DyVz = temp[2]; DyTxx = temp[3];  DyTyy = temp[4]; DyTzz = temp[5]; DyTxy = temp[6]; DyTxz = temp[7];
    
                        L_SIMD(wave_8_c,pos,8,FlagZ);
                        DzVx = temp[0];  DzVy = temp[1];  DzVz = temp[2]; DzTxx = temp[3];  DzTyy = temp[4]; DzTzz = temp[5]; DzTxy = temp[6]; DzTxz = temp[7];

				        // Moment equation
				        hwave_8_l[pos_hw + 0] = ( (DxTxx*metric_l[pos_m + 0] + DxTxy*metric_l[pos_m + 1] + DxTxz*metric_l[pos_m + 2])*rb1
				                         + (DyTxx*metric_l[pos_m + 3] + DyTxy*metric_l[pos_m + 4] + DyTxz*metric_l[pos_m + 5])*rb2
				                         + (DzTxx*metric_l[pos_m + 6] + DzTxy*metric_l[pos_m + 7] + DzTxz*metric_l[pos_m + 8])*rb3 )*rrho;

				        hwave_8_l[pos_hw + 1] = ( (DxTxy*metric_l[pos_m + 0] + DxTyy*metric_l[pos_m + 1] + DxTyz*metric_l[pos_m + 2])*rb1
				                         + (DyTxy*metric_l[pos_m + 3] + DyTyy*metric_l[pos_m + 4] + DyTyz*metric_l[pos_m + 5])*rb2
				                         + (DzTxy*metric_l[pos_m + 6] + DzTyy*metric_l[pos_m + 7] + DzTyz*metric_l[pos_m + 8])*rb3 )*rrho;

				        hwave_8_l[pos_hw + 2] = ( (DxTxz*metric_l[pos_m + 0] + DxTyz*metric_l[pos_m + 1] + DxTzz*metric_l[pos_m + 2])*rb1
				                         + (DyTxz*metric_l[pos_m + 3] + DyTyz*metric_l[pos_m + 4] + DyTzz*metric_l[pos_m + 5])*rb2
				                         + (DzTxz*metric_l[pos_m + 6] + DzTyz*metric_l[pos_m + 7] + DzTzz*metric_l[pos_m + 8])*rb3 )*rrho;

				        // Hooke's law
				        hwave_8_l[pos_hw + 3] = (lam2mu*DxVx*metric_l[pos_m + 0] + lam*DxVy*metric_l[pos_m + 1] + lam*DxVz*metric_l[pos_m + 2])*rb1
				                        + (lam2mu*DyVx*metric_l[pos_m + 3] + lam*DyVy*metric_l[pos_m + 4] + lam*DyVz*metric_l[pos_m + 5])*rb2
				                        + (lam2mu*DzVx*metric_l[pos_m + 6] + lam*DzVy*metric_l[pos_m + 7] + lam*DzVz*metric_l[pos_m + 8])*rb3;

				        hwave_8_l[pos_hw + 4] = (lam*DxVx*metric_l[pos_m + 0] + lam2mu*DxVy*metric_l[pos_m + 1] + lam*DxVz*metric_l[pos_m + 2])*rb1
				                        + (lam*DyVx*metric_l[pos_m + 3] + lam2mu*DyVy*metric_l[pos_m + 4] + lam*DyVz*metric_l[pos_m + 5])*rb2
				                        + (lam*DzVx*metric_l[pos_m + 6] + lam2mu*DzVy*metric_l[pos_m + 7] + lam*DzVz*metric_l[pos_m + 8])*rb3;

				        hwave_8_l[pos_hw + 5] = (lam*DxVx*metric_l[pos_m + 0] + lam*DxVy*metric_l[pos_m + 1] + lam2mu*DxVz*metric_l[pos_m + 2])*rb1
				                        + (lam*DyVx*metric_l[pos_m + 3] + lam*DyVy*metric_l[pos_m + 4] + lam2mu*DyVz*metric_l[pos_m + 5])*rb2
				                        + (lam*DzVx*metric_l[pos_m + 6] + lam*DzVy*metric_l[pos_m + 7] + lam2mu*DzVz*metric_l[pos_m + 8])*rb3;

				        hwave_8_l[pos_hw + 6] = mu * ( (DxVx*metric_l[pos_m + 1] + DxVy*metric_l[pos_m + 0])*rb1
				                               + (DyVx*metric_l[pos_m + 4] + DyVy*metric_l[pos_m + 3])*rb2
				                               + (DzVx*metric_l[pos_m + 7] + DzVy*metric_l[pos_m + 6])*rb3 );

				        hwave_8_l[pos_hw + 7] = mu * ( (DxVx*metric_l[pos_m + 2] + DxVz*metric_l[pos_m + 0])*rb1
				                               + (DyVx*metric_l[pos_m + 5] + DyVz*metric_l[pos_m + 3])*rb2
				                               + (DzVx*metric_l[pos_m + 8] + DzVz*metric_l[pos_m + 6])*rb3 );

				        hwave_1_l[tpos_hw] = mu * ( (DxVy*metric_l[pos_m + 2] + DxVz*metric_l[pos_m + 1])*rb1
				                               + (DyVy*metric_l[pos_m + 5] + DyVz*metric_l[pos_m + 4])*rb2
				                               + (DzVy*metric_l[pos_m + 8] + DzVz*metric_l[pos_m + 7])*rb3 );

#ifdef FreeSurface
				        /**************************************** Calculate step 2  ***********************************/
				        if (freenode && iz >= izn) {
                            n = izn_final - iz - 1;
							pos_m = plane * xstep_l_m5 + iy * ystep_l_m5 + iz * MSIZE;
							pos_mm = iy * ystep_l_m + iz * MSIZE;
							rrhojac = 1.0/metric_c5[pos_m + 12]/metric_c5[pos_m + 9];

							// hVx -----------------------------------------
							for (l = 0; l < 5; l++) {
                                ii = addn(plane_comp, l, wx + x0 + x1 + 1); //hd
                                pos = ii * xstep_l_w_8 + iy * ystep_l_w_8 + iz * WSIZE_V;
								pos_m = ii * xstep_l_m5 + iy * ystep_l_m5 + iz * MSIZE;
                                vecTx[l] = metric_c5[pos_m + 9] * (metric_c5[pos_m + 0] * wave_8_c[pos + 3] + metric_c5[pos_m + 1] * wave_8_c[pos + 6] + metric_c5[pos_m + 2] * wave_8_c[pos + 7]);

								jj = -y0 + l + iy;
                                pos = plane * xstep_l_w_8 + jj * ystep_l_w_8 + iz * WSIZE_V;
								pos_m = plane * xstep_l_m5 + jj * ystep_l_m5 + iz * MSIZE;
								vecTy[l] = metric_c5[pos_m + 9] * (metric_c5[pos_m + 3] * wave_8_c[pos + 3] + metric_c5[pos_m + 4] * wave_8_c[pos + 6] + metric_c5[pos_m + 5] * wave_8_c[pos + 7]);

								kk = -z0 + l + iz;
                                pos = plane * xstep_l_w_8 + iy * ystep_l_w_8 + kk * WSIZE_V;
								pos_m = plane * xstep_l_m5 + iy * ystep_l_m5 + kk * MSIZE;
                                vecTz[l] = metric_c5[pos_m + 9] * (metric_c5[pos_m + 6] * wave_8_c[pos + 3] + metric_c5[pos_m + 7] * wave_8_c[pos + 6] + metric_c5[pos_m + 8] * wave_8_c[pos + 7]);
							}

				        	if (FlagZ == 1) {                                   // forward direction
				        		if (n == 2) {                                   // save vec_tmp data for vecTz[-2] which will be used when n=0
									vecTz_tmp1[iy] = vecTz[0];
									vecTz_tmp2[iy] = vecTz[1];
				        		}
				        		vecTz[-2] = vecTz_tmp1[iy];
				        		vecTz[-1] = vecTz_tmp2[iy];
								for (l = 0; l < 3 - n; l++)                     // traction image in different condition point
									vecTz[4-l] = -vecTz[l+2*(n-1)];

				        		vecTz[n + 1] = 0.0;                             // free surface set to 0
				        	} else {                                            // back forward
				        		if (n == 0) vecTz[4] = -vecTz[2];               // traction image when n=0, otherwise no need
				        		if (n < 2) vecTz[n+3] = 0.0;
				        	}

					        DxTx = vec_LL(vecTx,0,FlagX);
					        DyTy = vec_LL(vecTy,0,FlagY);
					        DzTz = vec_LL(vecTz,0,FlagZ);

                            pos_hw = iz*WSIZE_V;
					        hwave_8_l[pos_hw + 0] = (DxTx*rb1 + DyTy*rb2 + DzTz*rb3) * rrhojac;

					        // hVy ----------------------------------------------
							for (l = 0; l < 5; l++) {
                                ii = addn(plane_comp, l, wx + x0 + x1 + 1);  //hd
                                tpos8 = ii * xstep_l_w_8 + iy * ystep_l_w_8 + iz * WSIZE_V;
                                tpos1 = ii * xstep_l_w_1 + iy * ystep_l_w_1 + iz;
								pos_m = ii * xstep_l_m5 + iy * ystep_l_m5 + iz * MSIZE;
                                vecTx[l] = metric_c5[pos_m + 9] * (metric_c5[pos_m + 0] * wave_8_c[tpos8 + 6] + metric_c5[pos_m + 1] * wave_8_c[tpos8 + 4] + metric_c5[pos_m + 2] * wave_1_c[tpos1]);

								jj = -y0 + l + iy;
                                tpos8 = plane * xstep_l_w_8 + jj * ystep_l_w_8 + iz * WSIZE_V;
                                tpos1 = plane * xstep_l_w_1 + jj * ystep_l_w_1 + iz;
								pos_m = plane * xstep_l_m5 + jj * ystep_l_m5 + iz * MSIZE;
                                vecTy[l] = metric_c5[pos_m + 9] * (metric_c5[pos_m + 3] * wave_8_c[tpos8 + 6] + metric_c5[pos_m + 4] * wave_8_c[tpos8 + 4] + metric_c5[pos_m + 5] * wave_1_c[tpos1]);

								kk = -z0 + l + iz;
								// pos = plane * xstep_l_w + iy * ystep_l_w + kk * WSIZE;
                                tpos8 = plane * xstep_l_w_8 + iy * ystep_l_w_8 + kk * WSIZE_V;
                                tpos1 = plane * xstep_l_w_1 + iy * ystep_l_w_1 + kk;
								pos_m = plane * xstep_l_m5 + iy * ystep_l_m5 + kk * MSIZE;
								vecTz[l] = metric_c5[pos_m + 9] * (metric_c5[pos_m + 6] * wave_8_c[tpos8 + 6] + metric_c5[pos_m + 7] * wave_8_c[tpos8 + 4] + metric_c5[pos_m + 8] * wave_1_c[tpos1]);
							}

				        	if (FlagZ == 1) {                                   // forward direction
				        		if (n == 2) {                                   // save vec_tmp data for vecTz[-2] which will be used when n=0
									vecTz_tmp3[iy] = vecTz[0];
									vecTz_tmp4[iy] = vecTz[1];
				        		}
				        		vecTz[-2] = vecTz_tmp3[iy];
				        		vecTz[-1] = vecTz_tmp4[iy];
								for (l = 0; l < 3 - n; l++)                     // traction image in different condition point
									vecTz[4-l] = -vecTz[l+2*(n-1)];

				        		vecTz[n + 1] = 0.0;                             // free surface set to 0
				        	} else {                                            // back forward
				        		if (n == 0) vecTz[4] = -vecTz[2];               // traction image when n=0, otherwise no need
				        		if (n < 2) vecTz[n+3] = 0.0;
				        	}

					        DxTx = vec_LL(vecTx,0,FlagX);
					        DyTy = vec_LL(vecTy,0,FlagY);
					        DzTz = vec_LL(vecTz,0,FlagZ);

                            pos_hw = iz*WSIZE_V;
					        hwave_8_l[pos_hw + 1] = (DxTx*rb1 + DyTy*rb2 + DzTz*rb3) * rrhojac;

					        // hVz ----------------------------------------------
							for (l = 0; l < 5; l++) {
                                ii = addn(plane_comp, l, wx + x0 + x1 + 1);  //hd
                                tpos8 = ii * xstep_l_w_8 + iy * ystep_l_w_8 + iz * WSIZE_V;
                                tpos1 = ii * xstep_l_w_1 + iy * ystep_l_w_1 + iz;
								pos_m = ii * xstep_l_m5 + iy * ystep_l_m5 + iz * MSIZE;
                                vecTx[l] = metric_c5[pos_m + 9] * (metric_c5[pos_m + 0] * wave_8_c[tpos8 + 7] + metric_c5[pos_m + 1] * wave_1_c[tpos1] + metric_c5[pos_m + 2] * wave_8_c[tpos8 + 5]);

								jj = -y0 + l + iy;
                                tpos8 = plane * xstep_l_w_8 + jj * ystep_l_w_8 + iz * WSIZE_V;
                                tpos1 = plane * xstep_l_w_1 + jj * ystep_l_w_1 + iz;
								pos_m = plane * xstep_l_m5 + jj * ystep_l_m5 + iz * MSIZE;
                                vecTy[l] = metric_c5[pos_m + 9] * (metric_c5[pos_m + 3] * wave_8_c[tpos8 + 7] + metric_c5[pos_m + 4] * wave_1_c[tpos1] + metric_c5[pos_m + 5] * wave_8_c[tpos8]);

								kk = -z0 + l + iz;
                                tpos8 = plane * xstep_l_w_8 + iy * ystep_l_w_8 + kk * WSIZE_V;
                                tpos1 = plane * xstep_l_w_1 + iy * ystep_l_w_1 + kk;
								pos_m = plane * xstep_l_m5 + iy * ystep_l_m5 + kk * MSIZE;
                                vecTz[l] = metric_c5[pos_m + 9] * (metric_c5[pos_m + 6] * wave_8_c[tpos8 + 7] + metric_c5[pos_m + 7] * wave_1_c[tpos1] + metric_c5[pos_m + 8] * wave_8_c[tpos8 + 5]);
							}
				        	if (FlagZ == 1) {                                   // forward direction
				        		if (n == 2) {                                   // save vec_tmp data for vecTz[-2] which will be used when n=0
									vecTz_tmp5[iy] = vecTz[0];
									vecTz_tmp6[iy] = vecTz[1];
				        		}
				        		vecTz[-2] = vecTz_tmp5[iy];
				        		vecTz[-1] = vecTz_tmp6[iy];
								for (l = 0; l < 3 - n; l++)                     // traction image in different condition point
									vecTz[4-l] = -vecTz[l+2*(n-1)];

				        		vecTz[n + 1] = 0.0;                             // free surface set to 0
				        	} else {                                            // back forward
				        		if (n == 0) vecTz[4] = -vecTz[2];               // traction image when n=0, otherwise no need
				        		if (n < 2) vecTz[n+3] = 0.0;
				        	}
					        DxTx = vec_LL(vecTx,0,FlagX);
					        DyTy = vec_LL(vecTy,0,FlagY);
					        DzTz = vec_LL(vecTz,0,FlagZ);

                            pos_hw = iz*WSIZE_V;
					        hwave_8_l[pos_hw + 2] = (DxTx*rb1 + DyTy*rb2 + DzTz*rb3) * rrhojac;

				           /*********************************** Calculate step 3  *************************************/
                            pos = plane * xstep_l_w_8 + iy * ystep_l_w_8 + iz * WSIZE_V;
							pos_m = plane * xstep_l_m5 + iy * ystep_l_m5 + iz * MSIZE;

					        if (n == 0) {
                                   DzVx = matVx2Vz[3*0 + 0]*DxVx
                                        + matVx2Vz[3*0 + 1]*DxVy
                                        + matVx2Vz[3*0 + 2]*DxVz
                                        + matVy2Vz[3*0 + 0]*DyVx
                                        + matVy2Vz[3*0 + 1]*DyVy
                                        + matVy2Vz[3*0 + 2]*DyVz;

                                   DzVy = matVx2Vz[3*1 + 0]*DxVx
                                        + matVx2Vz[3*1 + 1]*DxVy
                                        + matVx2Vz[3*1 + 2]*DxVz
                                        + matVy2Vz[3*1 + 0]*DyVx
                                        + matVy2Vz[3*1 + 1]*DyVy
                                        + matVy2Vz[3*1 + 2]*DyVz;

                                   DzVz = matVx2Vz[3*2 + 0]*DxVx
                                        + matVx2Vz[3*2 + 1]*DxVy
                                        + matVx2Vz[3*2 + 2]*DxVz
                                        + matVy2Vz[3*2 + 0]*DyVx
                                        + matVy2Vz[3*2 + 1]*DyVy
                                        + matVy2Vz[3*2 + 2]*DyVz;
					        } else if (n == 1) {
						        // Get Lz directly
                                DzVx = L22(wave_8_c,(pos + 0), 8, FlagZ);
						        DzVy = L22(wave_8_c,(pos + 1), 8, FlagZ);
						        DzVz = L22(wave_8_c,(pos + 2), 8, FlagZ);
					        } else {
						        // Get Lz directly
                                DzVx = L24(wave_8_c,(pos + 0), 8, FlagZ);
						        DzVy = L24(wave_8_c,(pos + 1), 8, FlagZ);
						        DzVz = L24(wave_8_c,(pos + 2), 8, FlagZ);
					        }

                            tpos_hw8 = iz*WSIZE_V;
                            tpos_hw1 = iz;

					        hwave_8_l[tpos_hw8 + 3] = (lam2mu*metric_c5[pos_m + 0]*DxVx + lam*metric_c5[pos_m + 1]*DxVy + lam*metric_c5[pos_m + 2]*DxVz)*rb1
						                        + (lam2mu*metric_c5[pos_m + 3]*DyVx + lam*metric_c5[pos_m + 4]*DyVy + lam*metric_c5[pos_m + 5]*DyVz)*rb2
						                        + (lam2mu*metric_c5[pos_m + 6]*DzVx + lam*metric_c5[pos_m + 7]*DzVy + lam*metric_c5[pos_m + 8]*DzVz)*rb3;

					        hwave_8_l[tpos_hw8 + 4] = (lam*metric_c5[pos_m + 0]*DxVx + lam2mu*metric_c5[pos_m + 1]*DxVy + lam*metric_c5[pos_m + 2]*DxVz)*rb1
						                        + (lam*metric_c5[pos_m + 3]*DyVx + lam2mu*metric_c5[pos_m + 4]*DyVy + lam*metric_c5[pos_m + 5]*DyVz)*rb2
						                        + (lam*metric_c5[pos_m + 6]*DzVx + lam2mu*metric_c5[pos_m + 7]*DzVy + lam*metric_c5[pos_m + 8]*DzVz)*rb3;

					        hwave_8_l[tpos_hw8 + 5] = (lam*metric_c5[pos_m + 0]*DxVx + lam*metric_c5[pos_m + 1]*DxVy + lam2mu*metric_c5[pos_m + 2]*DxVz)*rb1
						                        + (lam*metric_c5[pos_m + 3]*DyVx + lam*metric_c5[pos_m + 4]*DyVy + lam2mu*metric_c5[pos_m + 5]*DyVz)*rb2
						                        + (lam*metric_c5[pos_m + 6]*DzVx + lam*metric_c5[pos_m + 7]*DzVy + lam2mu*metric_c5[pos_m + 8]*DzVz)*rb3;

					        hwave_8_l[tpos_hw8 + 6] = mu*( (metric_c5[pos_m + 1]*DxVx + metric_c5[pos_m + 0]*DxVy)*rb1
						                             + (metric_c5[pos_m + 4]*DyVx + metric_c5[pos_m + 3]*DyVy)*rb2
						                             + (metric_c5[pos_m + 7]*DzVx + metric_c5[pos_m + 6]*DzVy)*rb3 );

					        hwave_8_l[tpos_hw8 + 7] = mu*( (metric_c5[pos_m + 2]*DxVx + metric_c5[pos_m + 0]*DxVz)*rb1
						                             + (metric_c5[pos_m + 5]*DyVx + metric_c5[pos_m + 3]*DyVz)*rb2
						                             + (metric_c5[pos_m + 8]*DzVx + metric_c5[pos_m + 6]*DzVz)*rb3 );

					        hwave_1_l[tpos_hw1] = mu*( (metric_c5[pos_m + 2]*DxVy + metric_c5[pos_m + 1]*DxVz)*rb1
						                             + (metric_c5[pos_m + 5]*DyVy + metric_c5[pos_m + 4]*DyVz)*rb2
						                             + (metric_c5[pos_m + 8]*DzVy + metric_c5[pos_m + 7]*DzVz)*rb3 );
						}
#endif

#ifdef usePML
                        if((isx1 && ix < PML_ND) || (isx2 && ix >= ni - PML_ND) || (isy1 && iybeg + iy < PML_ND) || (isy2 && iybeg + iy >= nj - PML_ND) || (isz1 && izbeg + iz < PML_ND) || (isz2 && izbeg + iz >= nk - PML_ND)) {
                            //pos = plane * xstep_l_w_8 + iy * ystep_l_w_8 + iz * WSIZE_V;
						    //pos_m = iz * MSIZE;
                            float hW_add[WSIZE] __attribute__((aligned(32)));
                            memset(hW_add, 0, sizeof(float) * WSIZE);
                            float *aux_w, *aux_hw, *aux_tw, *aux_mw;
                            float aux_wx[WSIZE] __attribute__((aligned(32))), aux_hwx[WSIZE] __attribute__((aligned(32))), aux_twx[WSIZE] __attribute__((aligned(32))), aux_mwx[WSIZE] __attribute__((aligned(32)));
                            float aux_wy[WSIZE] __attribute__((aligned(32))), aux_hwy[WSIZE] __attribute__((aligned(32))), aux_twy[WSIZE] __attribute__((aligned(32))), aux_mwy[WSIZE] __attribute__((aligned(32)));
                            float aux_wz[WSIZE] __attribute__((aligned(32))), aux_hwz[WSIZE] __attribute__((aligned(32))), aux_twz[WSIZE] __attribute__((aligned(32))), aux_mwz[WSIZE] __attribute__((aligned(32)));

                            pos_hw = iz * WSIZE_V;
                            tpos_hw = iz; 

                            if (isx1 && ix < PML_ND) {
                                int pni = PML_ND;
                                int pnj = ny;
                                int pnk = lnz;
                                idx0 = 0;
                                Aux_Wx = Aux->Wx1;
                                Aux_hWx = Aux->hWx1;

                                pos_a = (nk1 + izbeg + iz + (nj1 + iybeg + iy) * pnk + (ix - idx0) * pnj * pnk) * WSIZE;
                                iget_reply = 0;
                                athread_dma_iget(aux_wx, Aux_Wx + pos_a, sizeof(float) * WSIZE, &iget_reply);

                                aux_w = aux_wx;
                                aux_hw = aux_hwx;
                                aux_tw = aux_twx;
                                aux_mw = aux_mwx;
                                Aux_W = Aux->Wx1;
                                Aux_hW = Aux->hWx1;

                                iget_reply_pml = 0;

                                if (flag == 1) {
                                    Aux_tW = Aux->tWx1;
                                    Aux_mW = Aux->Wx1;
                                }
                                else if (flag == 2) {
                                    Aux_tW = Aux->tWx1;
                                    Aux_mW = Aux->mWx1;
                                    athread_dma_iget(aux_tw, Aux_tW + pos_a, sizeof(float) * WSIZE, &iget_reply_pml);
                                }
                                else if (flag == 3) {
                                    Aux_tW = Aux->mWx1;
                                    Aux_mW = Aux->tWx1;
                                }

                                athread_dma_iget(aux_mw, Aux_mW + pos_a, sizeof(float) * WSIZE, &iget_reply_pml);

                                b1 = bx[ix & DX];
                                d1 = dx[ix % DX];
                                ad1 = -(ax[ix % DX] + d1);

                                rb1 = 1.0f/b1;
                                
                                float Deriv[WSIZE];

                                Deriv[0] = (DxTxx*metric_l[pos_m + 0] + DxTxy*metric_l[pos_m + 1] + DxTxz*metric_l[pos_m + 2])*rrho;
                                Deriv[1] = (DxTxy*metric_l[pos_m + 0] + DxTyy*metric_l[pos_m + 1] + DxTyz*metric_l[pos_m + 2])*rrho;
                                Deriv[2] = (DxTxz*metric_l[pos_m + 0] + DxTyz*metric_l[pos_m + 1] + DxTzz*metric_l[pos_m + 2])*rrho;
                                Deriv[3] = (DxVx*metric_l[pos_m + 0]*lam2mu + DxVy*metric_l[pos_m + 1]*lam    + DxVz*metric_l[pos_m + 2]*lam   );
                                Deriv[4] = (DxVx*metric_l[pos_m + 0]*lam    + DxVy*metric_l[pos_m + 1]*lam2mu + DxVz*metric_l[pos_m + 2]*lam   );
                                Deriv[5] = (DxVx*metric_l[pos_m + 0]*lam    + DxVy*metric_l[pos_m + 1]*lam    + DxVz*metric_l[pos_m + 2]*lam2mu);
                                Deriv[6] = (DxVx*metric_l[pos_m + 1] + DxVy*metric_l[pos_m + 0])*mu;
                                Deriv[7] = (DxVx*metric_l[pos_m + 2] + DxVz*metric_l[pos_m + 0])*mu;
                                Deriv[8] = (DxVy*metric_l[pos_m + 2] + DxVz*metric_l[pos_m + 1])*mu;

                                // Auxiliary Equations
                                athread_dma_wait_value(&iget_reply, 1);

                                simd_load(va, &aux_wx[0]);
                                simd_load(vb, &hW_add[0]);
                                va = vb - va * rb1;
                                simd_store(va, &hW_add[0]);
                                hW_add[8] += -aux_wx[8] * rb1;
                                    
                                simd_load(va, &aux_wx[0]);
                                simd_load(vb, &Deriv[0]);
                                va = va * ad1 + vb * d1;
                                simd_store(va, &aux_hwx[0]);
                                aux_hwx[8] = ad1 * aux_wx[8] + d1 * Deriv[8];

#ifdef FreeSurface
                                if (freenode && izbeg + iz == nk-1) {
                                    // Get Lz by Lx and Ly
                                    DzVx1 = matVx2Vz[3*0+0] * DxVx
                                        + matVx2Vz[3*0+1] * DxVy
                                        + matVx2Vz[3*0+2] * DxVz;

                                    DzVy1 = matVx2Vz[3*1+0] * DxVx
                                        + matVx2Vz[3*1+1] * DxVy
                                        + matVx2Vz[3*1+2] * DxVz;

                                    DzVz1 = matVx2Vz[3*2+0] * DxVx
                                        + matVx2Vz[3*2+1] * DxVy
                                        + matVx2Vz[3*2+2] * DxVz;

                                    aux_hwx[3] += d1*b1*(lam2mu*metric_l[pos_m + 6]*DzVx1 + lam   *metric_l[pos_m + 7]*DzVy1 + lam   *metric_l[pos_m + 8]*DzVz1);
                                    aux_hwx[4] += d1*b1*(lam   *metric_l[pos_m + 6]*DzVx1 + lam2mu*metric_l[pos_m + 7]*DzVy1 + lam   *metric_l[pos_m + 8]*DzVz1);
                                    aux_hwx[5] += d1*b1*(lam   *metric_l[pos_m + 6]*DzVx1 + lam   *metric_l[pos_m + 7]*DzVy1 + lam2mu*metric_l[pos_m + 8]*DzVz1);
                                    aux_hwx[6] += d1*b1*mu*(metric_l[pos_m + 7]*DzVx1 + metric_l[pos_m + 6]*DzVy1);
                                    aux_hwx[7] += d1*b1*mu*(metric_l[pos_m + 8]*DzVx1 + metric_l[pos_m + 6]*DzVz1);
                                    aux_hwx[8] += d1*b1*mu*(metric_l[pos_m + 8]*DzVy1 + metric_l[pos_m + 7]*DzVz1);
                                }
#endif

                                if(flag == 1) {
                                    athread_dma_wait_value(&iget_reply_pml, 1);

                                    simd_load(va, &aux_hw[0]);
                                    simd_load(vb, &aux_mw[0]);
                                    va = va * rkb + vb;
                                    simd_store(va, &aux_tw[0]);

                                    aux_tw[8] = aux_mw[8] + rkb * aux_hw[8];

                                    athread_dma_iput(Aux_tW + pos_a, aux_tw, sizeof(float) * WSIZE, &iput_reply_pml);
                                    iput_count++;
                                }
                                else if(flag == 2) {
                                    athread_dma_wait_value(&iget_reply_pml, 2);

                                    simd_load(va, &aux_hw[0]);
                                    simd_load(vb, &aux_tw[0]);
                                    va = va * rkb + vb;
                                    simd_store(va, &aux_tw[0]);

                                    aux_tw[8] += rkb * aux_hw[8];

                                    athread_dma_iput(Aux_tW + pos_a, aux_tw, sizeof(float) * WSIZE, &iput_reply_pml);
                                    iput_count++;
                                }
                                else if (flag == 3) {
                                    athread_dma_wait_value(&iget_reply_pml, 1);
                                }

                                simd_load(va, &aux_hw[0]);
                                simd_load(vb, &aux_mw[0]);
                                va = va * rka + vb;
                                simd_store(va, &aux_hw[0]);

                                aux_hw[8] = aux_mw[8] + rka * aux_hw[8];

                                athread_dma_iput(Aux_hW + pos_a, aux_hw, sizeof(float)*WSIZE, &iput_reply_pml);
                                iput_count++;

                            }
                            if (isx2 && ix >= ni - PML_ND) {
                                int pni = PML_ND;
                                int pnj = ny;
                                int pnk = lnz;
                                idx0 = ni - PML_ND;
                                Aux_Wx = Aux->Wx2;
                                Aux_hWx = Aux->hWx2;

                                pos_a = (nk1 + izbeg + iz + (nj1 + iybeg + iy) * pnk + (ix - idx0) * pnj * pnk) * WSIZE;
                                iget_reply = 0;
                                athread_dma_iget(aux_wx, Aux_Wx + pos_a, sizeof(float) * WSIZE, &iget_reply);

                                aux_w = aux_wx;
                                aux_hw = aux_hwx;
                                aux_tw = aux_twx;
                                aux_mw = aux_mwx;
                                Aux_W = Aux->Wx2;
                                Aux_hW = Aux->hWx2;

                                iget_reply_pml = 0;

                                if(flag == 1) {
                                    Aux_tW = Aux->tWx2;
                                    Aux_mW = Aux->Wx2;
                                }
                                else if (flag == 2) {
                                    Aux_tW = Aux->tWx2;
                                    Aux_mW = Aux->mWx2;
                                    athread_dma_iget(aux_tw, Aux_tW + pos_a, sizeof(float) * WSIZE, &iget_reply_pml);
                                }
                                else if (flag == 3) {
                                    Aux_tW = Aux->mWx2;
                                    Aux_mW = Aux->tWx2;
                                }

                                athread_dma_iget(aux_mw, Aux_mW + pos_a, sizeof(float)*WSIZE, &iget_reply_pml);

                                b1 = bx[ix % DX];
                                d1 = dx[ix % DX];
                                ad1 = -(ax[ix % DX] + d1);

                                rb1 = 1.0f/b1;

                                float Deriv[WSIZE];

                                Deriv[0] = (DxTxx*metric_l[pos_m + 0] + DxTxy*metric_l[pos_m + 1] + DxTxz*metric_l[pos_m + 2])*rrho;
                                Deriv[1] = (DxTxy*metric_l[pos_m + 0] + DxTyy*metric_l[pos_m + 1] + DxTyz*metric_l[pos_m + 2])*rrho;
                                Deriv[2] = (DxTxz*metric_l[pos_m + 0] + DxTyz*metric_l[pos_m + 1] + DxTzz*metric_l[pos_m + 2])*rrho;
                                Deriv[3] = (DxVx*metric_l[pos_m + 0]*lam2mu + DxVy*metric_l[pos_m + 1]*lam    + DxVz*metric_l[pos_m + 2]*lam   );
                                Deriv[4] = (DxVx*metric_l[pos_m + 0]*lam    + DxVy*metric_l[pos_m + 1]*lam2mu + DxVz*metric_l[pos_m + 2]*lam   );
                                Deriv[5] = (DxVx*metric_l[pos_m + 0]*lam    + DxVy*metric_l[pos_m + 1]*lam    + DxVz*metric_l[pos_m + 2]*lam2mu);
                                Deriv[6] = (DxVx*metric_l[pos_m + 1] + DxVy*metric_l[pos_m + 0])*mu;
                                Deriv[7] = (DxVx*metric_l[pos_m + 2] + DxVz*metric_l[pos_m + 0])*mu;
                                Deriv[8] = (DxVy*metric_l[pos_m + 2] + DxVz*metric_l[pos_m + 1])*mu;

                                // Auxiliary Equations
                                athread_dma_wait_value(&iget_reply, 1);

                                simd_load(va, &aux_wx[0]);
                                simd_load(vb, &hW_add[0]);
                                va = vb - va * rb1;
                                simd_store(va, &hW_add[0]);
                                hW_add[8] += -aux_wx[8] * rb1;
                                    
                                simd_load(va, &aux_wx[0]);
                                simd_load(vb, &Deriv[0]);
                                va = va * ad1 + vb * d1;
                                simd_store(va, &aux_hwx[0]);
                                aux_hwx[8] = ad1 * aux_wx[8] + d1 * Deriv[8];

#ifdef FreeSurface
                                if (freenode && izbeg + iz == nk-1) {
                                    // Get Lz by Lx and Ly
                                    DzVx1 = matVx2Vz[3*0+0] * DxVx
                                        + matVx2Vz[3*0+1] * DxVy
                                        + matVx2Vz[3*0+2] * DxVz;

                                    DzVy1 = matVx2Vz[3*1+0] * DxVx
                                        + matVx2Vz[3*1+1] * DxVy
                                        + matVx2Vz[3*1+2] * DxVz;

                                    DzVz1 = matVx2Vz[3*2+0] * DxVx
                                        + matVx2Vz[3*2+1] * DxVy
                                        + matVx2Vz[3*2+2] * DxVz;

                                    aux_hwx[3] += d1*b1*(lam2mu*metric_l[pos_m + 6]*DzVx1 + lam   *metric_l[pos_m + 7]*DzVy1 + lam   *metric_l[pos_m + 8]*DzVz1);
                                    aux_hwx[4] += d1*b1*(lam   *metric_l[pos_m + 6]*DzVx1 + lam2mu*metric_l[pos_m + 7]*DzVy1 + lam   *metric_l[pos_m + 8]*DzVz1);
                                    aux_hwx[5] += d1*b1*(lam   *metric_l[pos_m + 6]*DzVx1 + lam   *metric_l[pos_m + 7]*DzVy1 + lam2mu*metric_l[pos_m + 8]*DzVz1);
                                    aux_hwx[6] += d1*b1*mu*(metric_l[pos_m + 7]*DzVx1 + metric_l[pos_m + 6]*DzVy1);
                                    aux_hwx[7] += d1*b1*mu*(metric_l[pos_m + 8]*DzVx1 + metric_l[pos_m + 6]*DzVz1);
                                    aux_hwx[8] += d1*b1*mu*(metric_l[pos_m + 8]*DzVy1 + metric_l[pos_m + 7]*DzVz1);
                                }
#endif 

                                if(flag == 1) {
                                    athread_dma_wait_value(&iget_reply_pml, 1);

                                    simd_load(va, &aux_hw[0]);
                                    simd_load(vb, &aux_mw[0]);
                                    va = va * rkb + vb;
                                    simd_store(va, &aux_tw[0]);

                                    aux_tw[8] = aux_mw[8] + rkb * aux_hw[8];

                                    athread_dma_iput(Aux_tW + pos_a, aux_tw, sizeof(float) * WSIZE, &iput_reply_pml);
                                    iput_count++;
                                }
                                else if(flag == 2) {
                                    athread_dma_wait_value(&iget_reply_pml, 2);

                                    simd_load(va, &aux_hw[0]);
                                    simd_load(vb, &aux_tw[0]);
                                    va = va * rkb + vb;
                                    simd_store(va, &aux_tw[0]);

                                    aux_tw[8] += rkb * aux_hw[8];

                                    athread_dma_iput(Aux_tW + pos_a, aux_tw, sizeof(float) * WSIZE, &iput_reply_pml);
                                    iput_count++;
                                }
                                else if (flag == 3) {
                                    athread_dma_wait_value(&iget_reply_pml, 1);
                                }

                                simd_load(va, &aux_hw[0]);
                                simd_load(vb, &aux_mw[0]);
                                va = va * rka + vb;
                                simd_store(va, &aux_hw[0]);

                                aux_hw[8] = aux_mw[8] + rka * aux_hw[8];

                                athread_dma_iput(Aux_hW + pos_a, aux_hw, sizeof(float)*WSIZE, &iput_reply_pml);
                                iput_count++;

                            }
                            if (isy1 && iybeg + iy < PML_ND) {
                                int pnj = PML_ND;
                                int pnk = lnz;
                                idx0 = 0;
                                Aux_Wy = Aux->Wy1;
                                Aux_hWy = Aux->hWy1;

                                pos_a = (nk1 + izbeg + iz + (iybeg + iy - idx0) * pnk + (ni1 + ix) * pnj * pnk) * WSIZE;
                                iget_reply = 0;
                                athread_dma_iget(aux_wy, Aux_Wy + pos_a, sizeof(float) * WSIZE, &iget_reply);

                                aux_w = aux_wy;
                                aux_hw = aux_hwy;
                                aux_tw = aux_twy;
                                aux_mw = aux_mwy;
                                Aux_W = Aux->Wy1;
                                Aux_hW = Aux->hWy1;

                                iget_reply_pml = 0;

                                if (flag == 1) {
                                    Aux_tW = Aux->tWy1;
                                    Aux_mW = Aux->Wy1;
                                }
                                else if (flag == 2) {
                                    Aux_tW = Aux->tWy1;
                                    Aux_mW = Aux->mWy1;
                                    athread_dma_iget(aux_tw, Aux_tW + pos_a, sizeof(float) * WSIZE, &iget_reply_pml);
                                }
                                else if (flag == 3) {
                                    Aux_tW = Aux->mWy1;
                                    Aux_mW = Aux->tWy1;
                                }

                                athread_dma_iget(aux_mw, Aux_mW + pos_a, sizeof(float)*WSIZE, &iget_reply_pml);

                                b1 = by[iy];
                                d1 = dy[iy];
                                ad1 = -(ay[iy] + d1);

                                rb1 = 1.0f/b1;

                                float Deriv[WSIZE];

                                Deriv[0] = (DyTxx*metric_l[pos_m + 3] + DyTxy*metric_l[pos_m + 4] + DyTxz*metric_l[pos_m + 5])*rrho;
                                Deriv[1] = (DyTxy*metric_l[pos_m + 3] + DyTyy*metric_l[pos_m + 4] + DyTyz*metric_l[pos_m + 5])*rrho;
                                Deriv[2] = (DyTxz*metric_l[pos_m + 3] + DyTyz*metric_l[pos_m + 4] + DyTzz*metric_l[pos_m + 5])*rrho;
                                Deriv[3] = (DyVx*metric_l[pos_m + 3]*lam2mu + DyVy*metric_l[pos_m + 4]*lam    + DyVz*metric_l[pos_m + 5]*lam   );
                                Deriv[4] = (DyVx*metric_l[pos_m + 3]*lam    + DyVy*metric_l[pos_m + 4]*lam2mu + DyVz*metric_l[pos_m + 5]*lam   );
                                Deriv[5] = (DyVx*metric_l[pos_m + 3]*lam    + DyVy*metric_l[pos_m + 4]*lam    + DyVz*metric_l[pos_m + 5]*lam2mu);
                                Deriv[6] = (DyVx*metric_l[pos_m + 4] + DyVy*metric_l[pos_m + 3])*mu;
                                Deriv[7] = (DyVx*metric_l[pos_m + 5] + DyVz*metric_l[pos_m + 3])*mu;
                                Deriv[8] = (DyVy*metric_l[pos_m + 5] + DyVz*metric_l[pos_m + 4])*mu;

                                // Auxiliary Equations
                                athread_dma_wait_value(&iget_reply, 1);

                                simd_load(va, &aux_wy[0]);
                                simd_load(vb, &hW_add[0]);
                                va = vb - va * rb1;
                                simd_store(va, &hW_add[0]);
                                hW_add[8] += -aux_wy[8] * rb1;
                                simd_load(va, &aux_wy[0]);
                                simd_load(vb, &Deriv[0]);
                                va = va * ad1 + vb * d1;
                                simd_store(va, &aux_hwy[0]);
                                aux_hwy[8] = ad1 * aux_wy[8] + d1 * Deriv[8];

#ifdef FreeSurface
                                if (freenode && izbeg + iz == nk-1) {
                                    // Get Lz by Lx and Ly
                                    DzVx1 = matVy2Vz[3*0+0] * DyVx
                                            + matVy2Vz[3*0+1] * DyVy
                                            + matVy2Vz[3*0+2] * DyVz;

                                    DzVy1 = matVy2Vz[3*1+0] * DyVx
                                            + matVy2Vz[3*1+1] * DyVy
                                            + matVy2Vz[3*1+2] * DyVz;

                                    DzVz1 = matVy2Vz[3*2+0] * DyVx
                                            + matVy2Vz[3*2+1] * DyVy
                                            + matVy2Vz[3*2+2] * DyVz;

                                    aux_hwy[3] += d1*b1*(lam2mu*metric_l[pos_m + 6]*DzVx1 + lam   *metric_l[pos_m + 7]*DzVy1 + lam   *metric_l[pos_m + 8]*DzVz1);
                                    aux_hwy[4] += d1*b1*(lam   *metric_l[pos_m + 6]*DzVx1 + lam2mu*metric_l[pos_m + 7]*DzVy1 + lam   *metric_l[pos_m + 8]*DzVz1);
                                    aux_hwy[5] += d1*b1*(lam   *metric_l[pos_m + 6]*DzVx1 + lam   *metric_l[pos_m + 7]*DzVy1 + lam2mu*metric_l[pos_m + 8]*DzVz1);
                                    aux_hwy[6] += d1*b1*mu*(metric_l[pos_m + 7]*DzVx1 + metric_l[pos_m + 6]*DzVy1);
                                    aux_hwy[7] += d1*b1*mu*(metric_l[pos_m + 8]*DzVx1 + metric_l[pos_m + 6]*DzVz1);
                                    aux_hwy[8] += d1*b1*mu*(metric_l[pos_m + 8]*DzVy1 + metric_l[pos_m + 7]*DzVz1);
                                }
#endif

                                if(flag == 1) {
                                    athread_dma_wait_value(&iget_reply_pml, 1);

                                    simd_load(va, &aux_hw[0]);
                                    simd_load(vb, &aux_mw[0]);
                                    va = va * rkb + vb;
                                    simd_store(va, &aux_tw[0]);

                                    aux_tw[8] = aux_mw[8] + rkb * aux_hw[8];

                                    athread_dma_iput(Aux_tW + pos_a, aux_tw, sizeof(float) * WSIZE, &iput_reply_pml);
                                    iput_count++;
                                }
                                else if(flag == 2) {
                                    athread_dma_wait_value(&iget_reply_pml, 2);

                                    simd_load(va, &aux_hw[0]);
                                    simd_load(vb, &aux_tw[0]);
                                    va = va * rkb + vb;
                                    simd_store(va, &aux_tw[0]);

                                    aux_tw[8] += rkb * aux_hw[8];

                                    athread_dma_iput(Aux_tW + pos_a, aux_tw, sizeof(float) * WSIZE, &iput_reply_pml);
                                    iput_count++;
                                }
                                else if (flag == 3) {
                                    athread_dma_wait_value(&iget_reply_pml, 1);
                                }

                                simd_load(va, &aux_hw[0]);
                                simd_load(vb, &aux_mw[0]);
                                va = va * rka + vb;
                                simd_store(va, &aux_hw[0]);

                                aux_hw[8] = aux_mw[8] + rka * aux_hw[8];

                                athread_dma_iput(Aux_hW + pos_a, aux_hw, sizeof(float)*WSIZE, &iput_reply_pml);
                                iput_count++;

                            }
                            if (isy2 && iybeg + iy >= nj - PML_ND) {
                                int pnj = PML_ND;
                                int pnk = lnz;
                                idx0 = nj - PML_ND;
                                Aux_Wy = Aux->Wy2;
                                Aux_hWy = Aux->hWy2;

                                pos_a = (nk1 + izbeg + iz + (iybeg + iy - idx0) * pnk + (ni1 + ix) * pnj * pnk) * WSIZE;
                                iget_reply = 0;
                                athread_dma_iget(aux_wy, Aux_Wy + pos_a, sizeof(float) * WSIZE, &iget_reply);

                                aux_w = aux_wy;
                                aux_hw = aux_hwy;
                                aux_tw = aux_twy;
                                aux_mw = aux_mwy;
                                Aux_W = Aux->Wy2;
                                Aux_hW = Aux->hWy2;

                                iget_reply_pml = 0;

                                if (flag == 1) {
                                    Aux_tW = Aux->tWy2;
                                    Aux_mW = Aux->Wy2;
                                }
                                else if (flag == 2) {
                                    Aux_tW = Aux->tWy2;
                                    Aux_mW = Aux->mWy2;
                                    athread_dma_iget(aux_tw, Aux_tW + pos_a, sizeof(float) * WSIZE, &iget_reply_pml);
                                }
                                else if (flag == 3) {
                                    Aux_tW = Aux->mWy2;
                                    Aux_mW = Aux->tWy2;
                                }

                                athread_dma_iget(aux_mw, Aux_mW + pos_a, sizeof(float)*WSIZE, &iget_reply_pml);

                                b1 = by[iy];
                                d1 = dy[iy];
                                ad1 = -(ay[iy] + d1);

                                rb1 = 1.0f/b1;

                                float Deriv[WSIZE];

                                Deriv[0] = (DyTxx*metric_l[pos_m + 3] + DyTxy*metric_l[pos_m + 4] + DyTxz*metric_l[pos_m + 5])*rrho;
                                Deriv[1] = (DyTxy*metric_l[pos_m + 3] + DyTyy*metric_l[pos_m + 4] + DyTyz*metric_l[pos_m + 5])*rrho;
                                Deriv[2] = (DyTxz*metric_l[pos_m + 3] + DyTyz*metric_l[pos_m + 4] + DyTzz*metric_l[pos_m + 5])*rrho;
                                Deriv[3] = (DyVx*metric_l[pos_m + 3]*lam2mu + DyVy*metric_l[pos_m + 4]*lam    + DyVz*metric_l[pos_m + 5]*lam   );
                                Deriv[4] = (DyVx*metric_l[pos_m + 3]*lam    + DyVy*metric_l[pos_m + 4]*lam2mu + DyVz*metric_l[pos_m + 5]*lam   );
                                Deriv[5] = (DyVx*metric_l[pos_m + 3]*lam    + DyVy*metric_l[pos_m + 4]*lam    + DyVz*metric_l[pos_m + 5]*lam2mu);
                                Deriv[6] = (DyVx*metric_l[pos_m + 4] + DyVy*metric_l[pos_m + 3])*mu;
                                Deriv[7] = (DyVx*metric_l[pos_m + 5] + DyVz*metric_l[pos_m + 3])*mu;
                                Deriv[8] = (DyVy*metric_l[pos_m + 5] + DyVz*metric_l[pos_m + 4])*mu;

                                // Auxiliary Equations
                                athread_dma_wait_value(&iget_reply, 1);

                                simd_load(va, &aux_wy[0]);
                                simd_load(vb, &hW_add[0]);
                                va = vb - va * rb1;
                                simd_store(va, &hW_add[0]);
                                hW_add[8] += -aux_wy[8] * rb1;
                                simd_load(va, &aux_wy[0]);
                                simd_load(vb, &Deriv[0]);
                                va = va * ad1 + vb * d1;
                                simd_store(va, &aux_hwy[0]);
                                aux_hwy[8] = ad1 * aux_wy[8] + d1 * Deriv[8];

#ifdef FreeSurface
                                if (freenode && izbeg + iz == nk-1) {
                                    // Get Lz by Lx and Ly
                                    DzVx1 = matVy2Vz[3*0+0] * DyVx
                                            + matVy2Vz[3*0+1] * DyVy
                                            + matVy2Vz[3*0+2] * DyVz;

                                    DzVy1 = matVy2Vz[3*1+0] * DyVx
                                            + matVy2Vz[3*1+1] * DyVy
                                            + matVy2Vz[3*1+2] * DyVz;

                                    DzVz1 = matVy2Vz[3*2+0] * DyVx
                                            + matVy2Vz[3*2+1] * DyVy
                                            + matVy2Vz[3*2+2] * DyVz;

                                    aux_hwy[3] += d1*b1*(lam2mu*metric_l[pos_m + 6]*DzVx1 + lam   *metric_l[pos_m + 7]*DzVy1 + lam   *metric_l[pos_m + 8]*DzVz1);
                                    aux_hwy[4] += d1*b1*(lam   *metric_l[pos_m + 6]*DzVx1 + lam2mu*metric_l[pos_m + 7]*DzVy1 + lam   *metric_l[pos_m + 8]*DzVz1);
                                    aux_hwy[5] += d1*b1*(lam   *metric_l[pos_m + 6]*DzVx1 + lam   *metric_l[pos_m + 7]*DzVy1 + lam2mu*metric_l[pos_m + 8]*DzVz1);
                                    aux_hwy[6] += d1*b1*mu*(metric_l[pos_m + 7]*DzVx1 + metric_l[pos_m + 6]*DzVy1);
                                    aux_hwy[7] += d1*b1*mu*(metric_l[pos_m + 8]*DzVx1 + metric_l[pos_m + 6]*DzVz1);
                                    aux_hwy[8] += d1*b1*mu*(metric_l[pos_m + 8]*DzVy1 + metric_l[pos_m + 7]*DzVz1);
                                }
#endif

                                if(flag == 1) {
                                    athread_dma_wait_value(&iget_reply_pml, 1);

                                    simd_load(va, &aux_hw[0]);
                                    simd_load(vb, &aux_mw[0]);
                                    va = va * rkb + vb;
                                    simd_store(va, &aux_tw[0]);

                                    aux_tw[8] = aux_mw[8] + rkb * aux_hw[8];

                                    athread_dma_iput(Aux_tW + pos_a, aux_tw, sizeof(float) * WSIZE, &iput_reply_pml);
                                    iput_count++;
                                }
                                else if(flag == 2) {
                                    athread_dma_wait_value(&iget_reply_pml, 2);

                                    simd_load(va, &aux_hw[0]);
                                    simd_load(vb, &aux_tw[0]);
                                    va = va * rkb + vb;
                                    simd_store(va, &aux_tw[0]);

                                    aux_tw[8] += rkb * aux_hw[8];

                                    athread_dma_iput(Aux_tW + pos_a, aux_tw, sizeof(float) * WSIZE, &iput_reply_pml);
                                    iput_count++;
                                }
                                else if (flag == 3) {
                                    athread_dma_wait_value(&iget_reply_pml, 1);
                                }

                                simd_load(va, &aux_hw[0]);
                                simd_load(vb, &aux_mw[0]);
                                va = va * rka + vb;
                                simd_store(va, &aux_hw[0]);

                                aux_hw[8] = aux_mw[8] + rka * aux_hw[8];

                                athread_dma_iput(Aux_hW + pos_a, aux_hw, sizeof(float)*WSIZE, &iput_reply_pml);
                                iput_count++;

                            }
                            if (isz1 && izbeg + iz < PML_ND) {
                                int pnj = ny;
                                int pnk = PML_ND;
                                idx0 = 0;
                                Aux_Wz = Aux->Wz1;
                                Aux_hWz = Aux->hWz1;

                                pos_a = (izbeg + iz - idx0 + (nj1 + iybeg + iy) * pnk + (ni1 + ix) * pnj * pnk) * WSIZE;
                                iget_reply = 0;
                                athread_dma_iget(aux_wz, Aux_Wz + pos_a, sizeof(float)*WSIZE, &iget_reply);

                                aux_w = aux_wz;
                                aux_hw = aux_hwz;
                                aux_tw = aux_twz;
                                aux_mw = aux_mwz;
                                Aux_W = Aux->Wz1;
                                Aux_hW = Aux->hWz1;

                                iget_reply_pml = 0;

                                if(flag == 1) {
                                    Aux_tW = Aux->tWz1;
                                    Aux_mW = Aux->Wz1;
                                }
                                else if (flag == 2) {
                                    Aux_tW = Aux->tWz1;
                                    Aux_mW = Aux->mWz1;
                                    athread_dma_iget(aux_tw, Aux_tW + pos_a, sizeof(float) * WSIZE, &iget_reply_pml);
                                }
                                else if (flag == 3) {
                                    Aux_tW = Aux->mWz1;
                                    Aux_mW = Aux->tWz1;
                                }

                                athread_dma_iget(aux_mw, Aux_mW + pos_a, sizeof(float)*WSIZE, &iget_reply_pml);

                                b1 = bz[iz];
                                d1 = dz[iz];
                                ad1 = -(az[iz] + d1);

                                rb1 = 1.0f/b1;

                                float Deriv[WSIZE];

                                Deriv[0] = (DzTxx*metric_l[pos_m + 6] + DzTxy*metric_l[pos_m + 7] + DzTxz*metric_l[pos_m + 8])*rrho;
                                Deriv[1] = (DzTxy*metric_l[pos_m + 6] + DzTyy*metric_l[pos_m + 7] + DzTyz*metric_l[pos_m + 8])*rrho;
                                Deriv[2] = (DzTxz*metric_l[pos_m + 6] + DzTyz*metric_l[pos_m + 7] + DzTzz*metric_l[pos_m + 8])*rrho;
                                Deriv[3] = (DzVx*metric_l[pos_m + 6]*lam2mu + DzVy*metric_l[pos_m + 7]*lam    + DzVz*metric_l[pos_m + 8]*lam   );
                                Deriv[4] = (DzVx*metric_l[pos_m + 6]*lam    + DzVy*metric_l[pos_m + 7]*lam2mu + DzVz*metric_l[pos_m + 8]*lam   );
                                Deriv[5] = (DzVx*metric_l[pos_m + 6]*lam    + DzVy*metric_l[pos_m + 7]*lam    + DzVz*metric_l[pos_m + 8]*lam2mu);
                                Deriv[6] = (DzVx*metric_l[pos_m + 7] + DzVy*metric_l[pos_m + 6])*mu;
                                Deriv[7] = (DzVx*metric_l[pos_m + 8] + DzVz*metric_l[pos_m + 6])*mu;
                                Deriv[8] = (DzVy*metric_l[pos_m + 8] + DzVz*metric_l[pos_m + 7])*mu;

                                // Auxiliary Equations
                                athread_dma_wait_value(&iget_reply, 1);
                                
                                simd_load(va, &aux_wz[0]);
                                simd_load(vb, &hW_add[0]);
                                va = vb - va * rb1;
                                simd_store(va, &hW_add[0]);
                                hW_add[8] += -aux_wz[8] * rb1;

                                simd_load(va, &aux_wz[0]);
                                simd_load(vb, &Deriv[0]);
                                va = va * ad1 + vb * d1;
                                simd_store(va, &aux_hwz[0]);
                                aux_hwz[8] = ad1 * aux_wz[8] + d1 * Deriv[8];

                                if(flag == 1) {
                                    athread_dma_wait_value(&iget_reply_pml, 1);

                                    simd_load(va, &aux_hw[0]);
                                    simd_load(vb, &aux_mw[0]);
                                    va = va * rkb + vb;
                                    simd_store(va, &aux_tw[0]);

                                    aux_tw[8] = aux_mw[8] + rkb * aux_hw[8];

                                    athread_dma_iput(Aux_tW + pos_a, aux_tw, sizeof(float) * WSIZE, &iput_reply_pml);
                                    iput_count++;
                                }
                                else if(flag == 2) {
                                    athread_dma_wait_value(&iget_reply_pml, 2);

                                    simd_load(va, &aux_hw[0]);
                                    simd_load(vb, &aux_tw[0]);
                                    va = va * rkb + vb;
                                    simd_store(va, &aux_tw[0]);

                                    aux_tw[8] += rkb * aux_hw[8];

                                    athread_dma_iput(Aux_tW + pos_a, aux_tw, sizeof(float) * WSIZE, &iput_reply_pml);
                                    iput_count++;
                                }
                                else if (flag == 3) {
                                    athread_dma_wait_value(&iget_reply_pml, 1);
                                }

                                simd_load(va, &aux_hw[0]);
                                simd_load(vb, &aux_mw[0]);
                                va = va * rka + vb;
                                simd_store(va, &aux_hw[0]);

                                aux_hw[8] = aux_mw[8] + rka * aux_hw[8];

                                athread_dma_iput(Aux_hW + pos_a, aux_hw, sizeof(float)*WSIZE, &iput_reply_pml);
                                iput_count++;
                                
                            }
                            if (isz2 && izbeg + iz >= nk - PML_ND) {
                                int pnj = ny;
                                int pnk = PML_ND;
                                idx0 = nk - PML_ND;
                                Aux_Wz = Aux->Wz2;
                                Aux_hWz = Aux->hWz2;

                                pos_a = (izbeg + iz - idx0 + (nj1 + iybeg + iy) * pnk + (ni1 + ix) * pnj * pnk) * WSIZE;
                                iget_reply = 0;
                                athread_dma_iget(aux_wz, Aux_Wz + pos_a, sizeof(float)*WSIZE, &iget_reply);

                                aux_w = aux_wz;
                                aux_hw = aux_hwz;
                                aux_tw = aux_twz;
                                aux_mw = aux_mwz;
                                Aux_W = Aux->Wz2;
                                Aux_hW = Aux->hWz2;

                                iget_reply_pml = 0;

                                if (flag == 1) {
                                    Aux_tW = Aux->tWz2;
                                    Aux_mW = Aux->Wz2;
                                }
                                else if (flag == 2) {
                                    Aux_tW = Aux->tWz2;
                                    Aux_mW = Aux->mWz2;
                                    athread_dma_iget(aux_tw, Aux_tW + pos_a, sizeof(float) * WSIZE, &iget_reply_pml);
                                }
                                else if (flag == 3) {
                                    Aux_tW = Aux->mWz2;
                                    Aux_mW = Aux->tWz2;
                                }

                                athread_dma_iget(aux_mw, Aux_mW + pos_a, sizeof(float)*WSIZE, &iget_reply_pml);

                                b1 = bz[iz];
                                d1 = dz[iz];
                                ad1 = -(az[iz] + d1);

                                rb1 = 1.0f/b1;
                                
                                float Deriv[WSIZE];

                                Deriv[0] = (DzTxx*metric_l[pos_m + 6] + DzTxy*metric_l[pos_m + 7] + DzTxz*metric_l[pos_m + 8])*rrho;
                                Deriv[1] = (DzTxy*metric_l[pos_m + 6] + DzTyy*metric_l[pos_m + 7] + DzTyz*metric_l[pos_m + 8])*rrho;
                                Deriv[2] = (DzTxz*metric_l[pos_m + 6] + DzTyz*metric_l[pos_m + 7] + DzTzz*metric_l[pos_m + 8])*rrho;
                                Deriv[3] = (DzVx*metric_l[pos_m + 6]*lam2mu + DzVy*metric_l[pos_m + 7]*lam    + DzVz*metric_l[pos_m + 8]*lam   );
                                Deriv[4] = (DzVx*metric_l[pos_m + 6]*lam    + DzVy*metric_l[pos_m + 7]*lam2mu + DzVz*metric_l[pos_m + 8]*lam   );
                                Deriv[5] = (DzVx*metric_l[pos_m + 6]*lam    + DzVy*metric_l[pos_m + 7]*lam    + DzVz*metric_l[pos_m + 8]*lam2mu);
                                Deriv[6] = (DzVx*metric_l[pos_m + 7] + DzVy*metric_l[pos_m + 6])*mu;
                                Deriv[7] = (DzVx*metric_l[pos_m + 8] + DzVz*metric_l[pos_m + 6])*mu;
                                Deriv[8] = (DzVy*metric_l[pos_m + 8] + DzVz*metric_l[pos_m + 7])*mu;

                                // Auxiliary Equations
                                athread_dma_wait_value(&iget_reply, 1);

                                simd_load(va, &aux_wz[0]);
                                simd_load(vb, &hW_add[0]);
                                va = vb - va * rb1;
                                simd_store(va, &hW_add[0]);
                                hW_add[8] += -aux_wz[8] * rb1;

                                simd_load(va, &aux_wz[0]);
                                simd_load(vb, &Deriv[0]);
                                va = va * ad1 + vb * d1;
                                simd_store(va, &aux_hwz[0]);
                                aux_hwz[8] = ad1 * aux_wz[8] + d1 * Deriv[8]; 

                                if(flag == 1) {
                                    athread_dma_wait_value(&iget_reply_pml, 1);

                                    simd_load(va, &aux_hw[0]);
                                    simd_load(vb, &aux_mw[0]);
                                    va = va * rkb + vb;
                                    simd_store(va, &aux_tw[0]);

                                    aux_tw[8] = aux_mw[8] + rkb * aux_hw[8];

                                    athread_dma_iput(Aux_tW + pos_a, aux_tw, sizeof(float) * WSIZE, &iput_reply_pml);
                                    iput_count++;
                                }
                                else if(flag == 2) {
                                    athread_dma_wait_value(&iget_reply_pml, 2);

                                    simd_load(va, &aux_hw[0]);
                                    simd_load(vb, &aux_tw[0]);
                                    va = va * rkb + vb;
                                    simd_store(va, &aux_tw[0]);

                                    aux_tw[8] += rkb * aux_hw[8];

                                    athread_dma_iput(Aux_tW + pos_a, aux_tw, sizeof(float) * WSIZE, &iput_reply_pml);
                                    iput_count++;
                                }
                                else if (flag == 3) {
                                    athread_dma_wait_value(&iget_reply_pml, 1);
                                }

                                simd_load(va, &aux_hw[0]);
                                simd_load(vb, &aux_mw[0]);
                                va = va * rka + vb;
                                simd_store(va, &aux_hw[0]);

                                aux_hw[8] = aux_mw[8] + rka * aux_hw[8];

                                athread_dma_iput(Aux_hW + pos_a, aux_hw, sizeof(float)*WSIZE, &iput_reply_pml);
                                iput_count++;

                            }

                            for(iii = 0 ; iii < WSIZE_V; iii ++) {
                                hwave_8_l[pos_hw + iii] += hW_add[iii];
                            } 
                            hwave_1_l[tpos_hw] += hW_add[8];

                            //wait iput value
                            //athread_dma_wait_value(&iput_reply_pml, iput_count);
                        }
#endif
                    }

                    athread_dma_wait_value(&reply_mw_8, 1);
                    athread_dma_wait_value(&reply_mw_1, 1);
                    athread_dma_wait_value(&reply_tw_8, 1);
                    athread_dma_wait_value(&reply_tw_1, 1);

#ifdef usePML
                    //wait iput value
                    athread_dma_wait_value(&iput_reply_pml, iput_count);
#endif
                    
				        /************************************ Calculate step 4  ***************************************/
                    for (iz = 0; iz < izn_final; iz++) {
                        tpos8 = iz*WSIZE_V;
                        tpos1 = iz;
						if (flag == 1) {  //RK_begin
                            simd_load(va, &hwave_8_l[tpos8 + 0]);
                            simd_load(vb, &mwave_8_l[tpos8 + 0]);
                            va = va * rkb + vb;
                            simd_store(va, &twave_8_l[tpos8 + 0]);

				        	twave_1_l[tpos1] = mwave_1_l[tpos1] + rkb * hwave_1_l[tpos1];
						}
						if (flag == 2) { // RK_inner
                            simd_load(va, &hwave_8_l[tpos8 + 0]);
                            simd_load(vb, &twave_8_l[tpos8 + 0]);
                            va = va * rkb + vb;
                            simd_store(va, &twave_8_l[tpos8 + 0]);
                         
					        twave_1_l[tpos1] += rkb * hwave_1_l[tpos1];
						}
                        simd_load(va, &hwave_8_l[tpos8 + 0]);
                        simd_load(vb, &mwave_8_l[tpos8 + 0]);
                        va = va * rka + vb;
                        simd_store(va, &hwave_8_l[tpos8 + 0]);

				        hwave_1_l[tpos1] = mwave_1_l[tpos1] + rka * hwave_1_l[tpos1];
					} // end iz

                    if (izn_final > 0) {
                        athread_dma_put(hwave_8_s_tmp, hwave_8_l, sizeof(float)*izn_final*WSIZE_V);
                        athread_dma_put(hwave_1_s_tmp, hwave_1_l, sizeof(float)*izn_final);

                        athread_dma_put(twave_8_s_tmp, twave_8_l, sizeof(float)*izn_final*WSIZE_V);
                        athread_dma_put(twave_1_s_tmp, twave_1_l, sizeof(float)*izn_final);
                    }

					hwave_8_s_tmp += ystep_w_8;
					hwave_1_s_tmp += ystep_w_1;

					twave_8_s_tmp += ystep_w_8;
					twave_1_s_tmp += ystep_w_1;

					mwave_8_s_tmp += ystep_w_8;
					mwave_1_s_tmp += ystep_w_1;

					metric_s_tmp += ystep_m;
				} // end iy
                
                athread_dma_wait_value(&reply_w_8, y0 + wy + y1);
                athread_dma_wait_value(&reply_w_1, y0 + wy + y1);
                if (freenode && (izend >= (nz - nk1 - 6))) {
                    athread_dma_wait_value(&reply_m, y0 + wy + y1);
                }
                
                plane_last = addn(plane_last, 1, wx + x0 + x1 + 1);  //hd
                plane_comp = addn(plane_comp, 1, wx + x0 + x1 + 1);  //hd             // plane ++
			} // end of loop of ix
		} // end of loop of iiy
	} // end of loop of iiz
}
#endif
