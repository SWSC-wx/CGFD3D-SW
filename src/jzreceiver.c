// made by wyw @ 2018.4.8
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "params.h"

void geo2cart(float *lat, float* lon, float lat0, float lon0, float alpha, float x0, float y0){
  // init
  int i;
  float a = alpha / 180 * PI;
  float D = 111319.5;
 
  if(masternode) for(i=0;i<nrec;i++) printf("lat = %f, lon = %f\n", lat[i], lon[i]);

  // cal
  for(i = 0; i < nrec; i++){
    lon[i] = x0 + D * ( cos(a)*(lat[i]-lat0) + sin(a)*(lon[i]-lon0) * cos(lat[i]/180*PI) );
    lat[i] = y0 + D * ( sin(a)*(lat[i]-lat0) - cos(a)*(lon[i]-lon0) * cos(lat[i]/180*PI) );
  }
  if(masternode) for(i=0;i<nrec;i++) printf("lat = %f, lon = %f\n", lat[i], lon[i]);
  return ;
}

void seek_index(int *rec_index, float *lat, float *lon, float x0, float y0){
  
  // init
  int i, tmp_ipx, tmp_ipy, tmp_rank;
  int *tmp_count = (int*)malloc(sizeof(int)*PX*PY*PZ);
  for(i = 0; i < PX*PY*PZ; i++) tmp_count[i] = 0;

  // cal
  for(i = 0; i < nrec; i++){
    tmp_ipx = round( (lon[i] - x0 ) / DH ) / ni;
    tmp_ipy = round( (lat[i] - y0 ) / DH ) / nj;
    tmp_rank = tmp_ipx * PY + tmp_ipy;
    rec_index[i*4+0] = tmp_rank;
    rec_index[i*4+1] = tmp_count[tmp_rank];
    tmp_count[tmp_rank] ++;//size:PX*PY*PZ
    rec_index[i*4+2] = round( (lon[i] - x0 ) / DH ) - tmp_ipx * ni; 
    rec_index[i*4+3] = round( (lat[i] - y0 ) / DH ) - tmp_ipy * nj;
  }
  return ;
}

void get_nrec(){
  
  FILE *fp;
  //printf("DEBUG: IN_REC = %s\n", IN_REC);
  fp = fopen(IN_REC, "r");
  if(!fp) { printf("can not open receiver file: %f\n", IN_REC); exit(-1);}

  //fread(&nrec, sizeof(int), 1, fp);
  fscanf(fp, "%d", &nrec);
  fclose(fp);
  return ;
}

int Load_Receiver(int *rec_index, int this_rank){

  // init 
//  int nrec, nrec_local;//global
  int i;
  float *laty = (float *)malloc(sizeof(float)*nrec);
  float *lonx = (float *)malloc(sizeof(float)*nrec);
  //float lat0 = 30;// for trans itself 
  //float lon0 = 103;// for trans itself
  float x0 = 0.0;
  float y0 = 0.0;
  //float alpha = 90;// for trans itself

  // read from file
  FILE *fp;
  fp = fopen(IN_REC, "r");
  if(!fp) { printf("can not open receiver file: %f\n", IN_REC); exit(-1);}

  //fread(&nrec, sizeof(int), 1, fp);
  //fread(laty, sizeof(float), nrec, fp);
  //fread(lonx, sizeof(float), nrec, fp);
  fscanf(fp, "%d", &nrec);
  for (i = 0; i < nrec; i++) fscanf(fp, "%f\t%f\n", &lonx[i], &laty[i]);


  // read in x0, y0 in topo  
  FILE *ftopo = NULL;
  ftopo = fopen(IN_TOPO, "rb");
  if(ftopo == NULL) {printf("can not open topo file in getting receiver!\n"); exit(-1);}

  fread(&x0, sizeof(float), 1, ftopo);
  fread(&y0, sizeof(float), 1, ftopo);
  fclose(ftopo);
 
  if(masternode) printf("receiver check: x0 y0 = %f %f \n", x0, y0);

  //geo2cart(lat, lon, lat0, lon0, alpha, 0.0, 0.0);//

  // seek for rec_index
  seek_index(rec_index, laty, lonx, x0, y0);

  return 0;
}

int keep_seismo(float *R, int *rec_index, float *W, int it){

  // init
  int i, pos, tmp;
  // keep
  // size of R: RNT*nrec*3
  for (i = 0; i < nrec; i++){
    if(rec_index[i*4+0] == this_rank){
      pos = ((rec_index[i*4+2] + ni1) * lnz * ny + (rec_index[i*4+3]+nj1) * lnz + nk2 - 1) * WSIZE;
      tmp = rec_index[i*4+1];
      R[nrec*RNT*0 + tmp*RNT + it] = W[pos + 0]; // vx
      R[nrec*RNT*1 + tmp*RNT + it] = W[pos + 1]; // vy
      R[nrec*RNT*2 + tmp*RNT + it] = W[pos + 2]; // vz
    }
  }

  return 0;
}

int export_seismo(float *R, int *rec_index, int tnow){

  // init 
  int i;
  FILE *fp = NULL;
  char filename[400];
  // export
  for(i = 0; i < nrec; i++){
    if(rec_index[i*4 + 0] == this_rank){
    // Vx
    sprintf(filename, "./%s/seismo_Vx_id%03d.bin", OUT_DIR, i);
    //sprintf(filename, "./output/seismo_Vx_id%03d.bin", i);
    fp = fopen(filename, "wb");
    fwrite( (R + nrec*RNT*0 + rec_index[i*4+1]*RNT ), sizeof(float), tnow, fp);
    fclose(fp);
    // Vy
    sprintf(filename, "./%s/seismo_Vy_id%03d.bin", OUT_DIR, i);
    //sprintf(filename, "./output/seismo_Vy_id%03d.bin", i);
    fp = fopen(filename, "wb");
    fwrite( (R + nrec*RNT*1 + rec_index[i*4+1]*RNT), sizeof(float), tnow, fp);
    fclose(fp);
    // Vz
    sprintf(filename, "./%s/seismo_Vz_id%03d.bin", OUT_DIR, i);
    //sprintf(filename, "./output/seismo_Vz_id%03d.bin", i);
    fp = fopen(filename, "wb");
    fwrite( (R + nrec*RNT*2 + rec_index[i*4+1]*RNT), sizeof(float), tnow, fp);
    fclose(fp);
    }
  }
  return 0;
}
