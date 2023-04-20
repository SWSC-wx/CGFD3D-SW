#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "params.h"
#include "common.h"
int construct_media2(float *C, float *M){

    float vp, vs, rho, lam, mu;
    int pos_c, pos_m;
    int i, j, k;
    for (i = ni1; i < ni2; i++){
        for (j = nj1; j < nj2; j++){
            for (k = nk1; k < nk2; k++){

                pos_c = (i*ny*lnz + j*lnz + k)*CSIZE;
                pos_m = (i*ny*lnz + j*lnz + k)*MSIZE;

                vp  = 6.000e3;
                vs  = 3.461e3;
                rho = 2.670e3;
                //vp  = 3.00e3;
                //vs  = 1.50e3;
                //rho = 1.80e3;
                //float tmpx = C[pos_c + 0];
                //float tmpy = C[pos_c + 1];
                float tmpz = C[pos_c + 2];
                /*
                   if( tmpz > -3e3){
                   vp  = vp / 1.2f;
                   vs  = vs / 1.2f;
                   rho = rho/ 1.2f;

                   }
                   */
                lam = rho * (vp*vp - 2.0f*vs*vs);
                mu  = rho * vs * vs;


                M[pos_m + 10] = lam;
                M[pos_m + 11] = mu ;
                M[pos_m + 12] = rho;

            }
        }
    }

    return 0;
}

int construct_media(float *C, float *M){

    char filename[100];
    int num = PX * PY * PZ;
    sprintf(filename, "./medium%d_%d/Vs_mpi%d_%d_%d", NX, num, thisid[0], thisid[1], thisid[2]);
    FILE *file_vs = fopen(filename, "rb");
    sprintf(filename, "./medium%d_%d/Vp_mpi%d_%d_%d", NX, num, thisid[0], thisid[1], thisid[2]);
    FILE *file_vp = fopen(filename, "rb");
    sprintf(filename, "./medium%d_%d/rho_mpi%d_%d_%d", NX, num, thisid[0], thisid[1], thisid[2]);
    FILE *file_rho = fopen(filename, "rb");

    float vs, vp, rho, qs, lam, mu; 

    for (int i = nk1; i < nk2; i++)
    {   
        for (int j = nj1; j < nj2; j++)
        {   
            int offset = (i * (ny * nx) + j * nx + 3) * sizeof(float);//debug
            fseek(file_vs, offset, SEEK_SET);
            fseek(file_vp, offset, SEEK_SET);
            fseek(file_rho, offset, SEEK_SET);
            for (int k = ni1; k < ni2; k++)
            {   
                //int loc = i * (NY * NX) + j * NX + k;//实际位置 
                int pos_m = (k * (ny * lnz) + j * lnz + i) * MSIZE; //坐标变换 
                fread(&vs, sizeof(float), 1, file_vs);
                fread(&vp, sizeof(float), 1, file_vp);
                fread(&rho, sizeof(float), 1, file_rho);
                qs = 1;
                lam = rho * (vp * vp - 2.0f * vs * vs);
                mu = rho * vs * vs; 
                M[pos_m + 10] = lam;
                M[pos_m + 11] = mu;
                M[pos_m + 12] = rho;
                M[pos_m + 13] = qs;

            }
        }
    }
    return 0;



}

int cal_range_media(float *M, float *range){

    float rho_min = 1.0e30 , rho_max = -1.0e30;
    float vp_min = 1.0e30 , vp_max = -1.0e30;
    float vs_min = 1.0e30 , vs_max = -1.0e30;
    float rho, vp, vs, lam, mu;
    int pos;//, pos1, pos2, pos3;
    //float L;
    int i, j, k;
    if ( masternode )
        printf( "==================M = %f\n", M[MSIZE * 10000 + 11] );

    for ( i = ni1; i < ni2; i++)
        for ( j = nj1; j < nj2; j++)
            for ( k = nk1; k < nk2; k++){

                pos = (i*ny*lnz+j*lnz+k)*MSIZE;
                //lam = M[pos + 10];
                //mu  = M[pos + 11];
                //rho = M[pos + 12];
                vs  = M[pos + 10];
                vp  = M[pos + 11];
                rho = M[pos + 12];
                lam = rho * (vp*vp - 2.0f*vs*vs);
                mu  = rho * vs * vs;
                M[pos + 10] = lam;
                M[pos + 11] = mu;

                //vp = sqrtf((lam + 2.0f * mu)/rho);
                //vs = sqrtf(mu/rho);

                rho_min = _min(rho_min, rho);
                rho_max = _max(rho_max, rho);
                vp_min = _min(vp_min, vp);
                vp_max = _max(vp_max, vp);
                vs_min = _min(vs_min, vs);
                vs_max = _max(vs_max, vs);

            }

    range[0] = vp_min;
    range[1] = vp_max;
    range[2] = vs_min;
    range[3] = vs_max;
    range[4] = rho_min;
    range[5] = rho_max;
    return 0;
}
