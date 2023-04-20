#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "params.h"
#include "common.h"
#include "sort_int.h"
#include <string.h>
#include <sys/stat.h>
#include <mpi.h>

#define POW2(x) ((x)*(x))
#define GAUSS_FUN(t, a, t0) (exp(-POW2( ( (t)-(t0) ) / (a) )) / (a*1.772453850905516 ))
#define NSRCEXT 3

int angle2moment(float strike1, float dip1, float rake1, float Moment[]){

  double s,d,r;
  s = strike1 / 180.0 * PI;
  d = dip1    / 180.0 * PI;
  r = rake1   / 180.0 * PI;

  double M11, M22, M33, M12, M13, M23;

  // in Aki and Richard's
  M11 = -(sin(d)*cos(r)*sin(2.0*s) + sin(2.0*d)*sin(r)*sin(s)*sin(s));
  M22 = sin(d)*cos(r)*sin(2.0*s) - sin(2.0*d)*sin(r)*cos(s)*cos(s);
  M33 = -(M11+M22);
  M12 = sin(d)*cos(r)*cos(2.0*s) + 0.5*sin(2.0*d)*sin(r)*sin(2.0*s);
  M13 = -(cos(d)*cos(r)*cos(s) + cos(2.0*d)*sin(r)*sin(s));
  M23 = -(cos(d)*cos(r)*sin(s) - cos(2.0*d)*sin(r)*cos(s));

  // transform North-East-Down to x-y-z
  // 1 <-> y
  // 2 <-> x
  // 3 <-> -z
  Moment[0] = (float) ( M22 );    // Mxx = M22
  Moment[1] = (float) ( M11 );    // Myy = M11
  Moment[2] = (float) ( M33 );    // Mzz = M33
  Moment[3] = (float) ( M12 );    // Mxy = M12
  Moment[4] = (float) (-M23 );   // Mxz = -M23
  Moment[5] = (float) (-M13 );   // Myz = -M13

  return 0;
}


float ricker_deriv(float rickerfc, float t, float tdelay){
  float f0 = sqrt(PI)*0.5f;
  float r = PI * rickerfc * (t-tdelay);
  float rr = r*r;
  //float s = (rr-0.5f)*exp(-rr) * f0;
  float s = r*(3.0f - 2.0f*rr) * exp(-rr) * f0 * PI * rickerfc;

  return s;
}

// int AddSourceRicker(float *W, float *M, int it) {

//   int i, j, k, gi, gj, gk, pos, pos_m;
//   float V;
//   float t = it*DT;

//   float rickerfc = 1.5f;
//   float tdelay = 1.2/rickerfc;

//   float s = 0.0f;

//   s = ricker_deriv(rickerfc, t, tdelay);

//   if(it == 0) s = 0.0;
//   //if(masternode) printf("SOURCE s = %e\n", s);

//   float Mxx, Myy, Mzz, Mxy, Mxz, Myz;

//   float M0 = 1e18;
//   Mxx = 1.0f; Myy = 1.0f; Mzz = 1.0f;
//   Mxy = 0.0f; Mxz = 0.0f; Myz = 0.0f;

//   //if( nsrclocal == 0) return 0;
//   //if( t > 3.0/rickerfc) return 0;

//   for (i = ni1; i < ni2; i++){
//     for (j = nj1; j < nj2; j++){
//       for (k = nk1; k < nk2; k++){

//         gi = thisid[0]*ni + i - ni1;
//         gj = thisid[1]*nj + j - nj1;
//         gk = thisid[2]*nk + k - nk1;

//         int isrc = NX/2 - 1;
//         int jsrc = NY/2 - 1;
//         // int ksrc = NZ - 10;
//         int ksrc = NZ/2 - 1;

//         if (
//               gi >= isrc - NSRCEXT && gi <= isrc + NSRCEXT &&
//               gj >= jsrc - NSRCEXT && gj <= jsrc + NSRCEXT &&
//               gk >= ksrc - NSRCEXT && gk <= ksrc + NSRCEXT ){

//           float ra = 0.5*NSRCEXT;
//           float D1 = GAUSS_FUN(gi-isrc, ra, 0.0);
//           float D2 = GAUSS_FUN(gj-jsrc, ra, 0.0);
//           float D3 = GAUSS_FUN(gk-ksrc, ra, 0.0);
//           float amp = D1*D2*D3;
          
//           amp /= 0.9951563131100551;

//           pos = (i * ny * lnz + j * lnz + k) * WSIZE;
//           pos_m = (i * ny * lnz + j * lnz + k) * MSIZE;
//           V = amp * DT /(M[pos_m + 9] * DH * DH * DH);

//           W[pos + 3] += -M0 * Mxx * s * V;
//           W[pos + 4] += -M0 * Myy * s * V;
//           W[pos + 5] += -M0 * Mzz * s * V;
//           W[pos + 6] += -M0 * Mxy * s * V;
//           W[pos + 7] += -M0 * Mxz * s * V;
//           W[pos + 8] += -M0 * Myz * s * V;
//         }
//       }
//     }
//   }

//   return 0;
// }

int AddSourceRicker(float *W_8, float *W_1, float *M, int it) {

  int i, j, k, gi, gj, gk, pos, pos_m;
  int pos_8, pos_1;
  float V;
  float t = it*DT;

  float rickerfc = 1.5f;
  float tdelay = 1.2/rickerfc;

  float s = 0.0f;

  s = ricker_deriv(rickerfc, t, tdelay);

  if(it == 0) s = 0.0;
  //if(masternode) printf("SOURCE s = %e\n", s);

  float Mxx, Myy, Mzz, Mxy, Mxz, Myz;

  float M0 = 1e18;
  Mxx = 1.0f; Myy = 1.0f; Mzz = 1.0f;
  Mxy = 0.0f; Mxz = 0.0f; Myz = 0.0f;

  //if( nsrclocal == 0) return 0;
  //if( t > 3.0/rickerfc) return 0;

  for (i = ni1; i < ni2; i++){
    for (j = nj1; j < nj2; j++){
      for (k = nk1; k < nk2; k++){

        gi = thisid[0]*ni + i - ni1;
        gj = thisid[1]*nj + j - nj1;
        gk = thisid[2]*nk + k - nk1;

        int isrc = NX/2 - 1;
        int jsrc = NY/2 - 1;
        // int ksrc = NZ - 10;
        int ksrc = NZ/2 - 1;

        if (
              gi >= isrc - NSRCEXT && gi <= isrc + NSRCEXT &&
              gj >= jsrc - NSRCEXT && gj <= jsrc + NSRCEXT &&
              gk >= ksrc - NSRCEXT && gk <= ksrc + NSRCEXT ){

          float ra = 0.5*NSRCEXT;
          float D1 = GAUSS_FUN(gi-isrc, ra, 0.0);
          float D2 = GAUSS_FUN(gj-jsrc, ra, 0.0);
          float D3 = GAUSS_FUN(gk-ksrc, ra, 0.0);
          float amp = D1*D2*D3;
          
          amp /= 0.9951563131100551;

          // pos = (i * ny * lnz + j * lnz + k) * WSIZE;
          pos_8 = (i * ny * lnz + j * lnz + k) * WSIZE_V;
          pos_1 = (i * ny * lnz + j * lnz + k);
          pos_m = (i * ny * lnz + j * lnz + k) * MSIZE;
          V = amp * DT /(M[pos_m + 9] * DH * DH * DH);

          W_8[pos_8 + 3] += -M0 * Mxx * s * V;
          W_8[pos_8 + 4] += -M0 * Myy * s * V;
          W_8[pos_8 + 5] += -M0 * Mzz * s * V;
          W_8[pos_8 + 6] += -M0 * Mxy * s * V;
          W_8[pos_8 + 7] += -M0 * Mxz * s * V;
          W_1[pos_1] += -M0 * Myz * s * V;
        }
      }
    }
  }

  return 0;
}

int read_source(char *topofile, char *srcfile){
//int read_source(char *topofile, char *srcfile, 
//int *srci, int *srcj, int *srck, float *area, float *strike, float *dip, float *rake, float *rate){
  //int src_npts, src_nt; float src_dt;
  int src_npts; // global
  //int nsrc; // local
  char partfile[1024];
  FILE *fp;

  int *gsrci = NULL, *gsrcj = NULL, *gsrck = NULL;
  float *garea = NULL, *gstrike = NULL, *gdip = NULL;
  float *grake = NULL, *grate = NULL;

  int *gsrci_sort = NULL, *gsrcj_sort = NULL, *gsrck_sort = NULL;
  float *garea_sort = NULL, *gstrike_sort = NULL, *gdip_sort = NULL;
  float *grake_sort = NULL, *grate_sort = NULL;

  //int *srci = NULL, *srcj = NULL, *srck = NULL;
  //float *area = NULL, *strike, *dip = NULL;
  //float *rake = NULL, *rate = NULL;

  int *src_npts_local = NULL;
  int *sendcounts = NULL;
  int *displs = NULL, *displs_nt = NULL;
  MPI_Comm comm = MPI_COMM_WORLD;

  int i, j, k, l, pos;

  if(masternode){
    // reading source
    fp = fopen(srcfile, "rb");
    struct stat statbuf;
    stat(srcfile, &statbuf);
    int filesize = statbuf.st_size;

    fread(&src_npts, sizeof(int), 1, fp);
    fread(&src_nt,   sizeof(int), 1, fp);
    fread(&src_dt, sizeof(float), 1, fp);

    printf("* SOURCE global alloc ...\n");
    float *sx = (float *) malloc(sizeof(float) * src_npts);
    float *sy = (float *) malloc(sizeof(float) * src_npts);
    float *sz = (float *) malloc(sizeof(float) * src_npts);
    garea   = (float *) malloc(sizeof(float) * src_npts);
    gstrike = (float *) malloc(sizeof(float) * src_npts);
    gdip    = (float *) malloc(sizeof(float) * src_npts);

    grake = (float *) malloc(sizeof(float) * src_npts * src_nt);
    grate = (float *) malloc(sizeof(float) * src_npts * src_nt);

    printf("* SOURCE reading ...\n");
    for (i = 0; i < src_npts; i++){
      fread(&sx[i], sizeof(float), 1, fp);
      fread(&sy[i], sizeof(float), 1, fp);
      fread(&sz[i], sizeof(float), 1, fp);
      fread(&garea[i],   sizeof(float), 1, fp);
      fread(&gstrike[i], sizeof(float), 1, fp);
      fread(&gdip[i],    sizeof(float), 1, fp);
      fread(&grake[i*src_nt], sizeof(float), src_nt, fp);
      fread(&grate[i*src_nt], sizeof(float), src_nt, fp);
      //printf("src coord [%06d]= %.2e %.2e %.2e\n", i, sx[i], sy[i], sz[i]);
    }

    //fseek(fp, 0, SEEK_END);
    //int filesize = ftell(fp);
    fclose(fp);

    int correct_size = 2*sizeof(int) + sizeof(float)*(1+src_npts*(6+2*src_nt));

    printf("* SOURCE file is '%s': filesize = %d (%.1e Mb)\n* source npts = %d, nt = %d, dt = %10.2e\n",
        srcfile, filesize, filesize/(1024.*1024.), src_npts, src_nt, src_dt);

    if(!(filesize == correct_size)){
      printf("! Error! SOURCE input file size is not right, please check it carefully!\n");
    }

    // reading topo
    //float *topo = (float *) malloc(sizeof(float) * NX * NY * 3);
    //float *coord = (float *) malloc(sizeof(float) * NX * NY * NZ * 3);
    // this way cost too much memory !!!

    float r, h, a;
    float x0, y0, z0;
    if(0 == icoord){
      // construct coord by the code
      //for (i = 0; i < NX; i++)
      //  for (j = 0; j < NY; j++){
      //    pos = (i*NY + j)*CSIZE;
      //    topo[pos + 0] = (i-NX/2) * DH;
      //    topo[pos + 1] = (j-NY/2) * DH;
      //    //float r = powf(C[pos + 0], 2) + powf(C[pos + 1], 2);
      //    r = powf(topo[pos + 0], 2);
      //    r += powf(topo[pos + 1], 2);
      //    h = gauss_height;
      //    a = gauss_width;;
      //    topo[pos + 2] = h * expf(-r/powf(a,2));
         
      x0 = (0-NX/2)*DH;
      y0 = (0-NY/2)*DH;
      z0 = 0.0;
          
      
    }else{
      // construct coord by reading topo data
      fp = fopen(topofile, "rb");
      //fread(topo, sizeof(float), NX*NY*3, fp);
      fread(&x0, sizeof(float), 1, fp);
      fread(&y0, sizeof(float), 1, fp);
      fread(&z0, sizeof(float), 1, fp);
      fclose(fp);

      printf("x0 y0 z0 = %f %f %f\n", x0, y0, z0);
    }

    // constructing coordinate
    //float xmin = 1e30, xmax = -1e30;
    //float ymin = 1e30, ymax = -1e30;
    //float zmin = 1e30, zmax = -1e30;
    //int pos_t, pos_c;
    //for (i = 0; i < NX; i++)
    //  for (j = 0; j < NY; j++)
    //    for (k = 0; k < NZ; k++){
    //      pos_c = (i*NY*NZ + j*NZ + k)*CSIZE;
    //      pos_t = (i*NY + j)*CSIZE;

    //      coord[pos_c + 0] = topo[pos_t + 0];
    //      coord[pos_c + 1] = topo[pos_t + 1];
    //      coord[pos_c + 2] = topo[pos_t + 2] + (k-NZ+1)*DH;

    //      xmin = min(xmin, coord[pos_c + 0]);
    //      ymin = min(ymin, coord[pos_c + 1]);
    //      zmin = min(zmin, coord[pos_c + 2]);
    //      xmax = max(xmax, coord[pos_c + 0]);
    //      ymax = max(ymax, coord[pos_c + 1]);
    //      zmax = max(zmax, coord[pos_c + 2]);
    //    }

    //printf("global x = %.2e ~ %.2e, y = %.2e ~ %.2e, z = %.2e ~ %.2e\n", 
    //    xmin, xmax, ymin, ymax, zmin, zmax);

    gsrci = (int *) malloc(sizeof(int) * src_npts);
    gsrcj = (int *) malloc(sizeof(int) * src_npts);
    gsrck = (int *) malloc(sizeof(int) * src_npts);

    // searching source index
//    float rmin, tmpx, tmpy, tmpz;
//    for (l = 0; l < src_npts; l++){
//      gsrci[l] = (int)((sx[l] - topo[0])/DH + 0.5f);
//      gsrcj[l] = (int)((sy[l] - topo[1])/DH + 0.5f);
//      rmin = 1e30;
//      for (k = 0; k < NZ; k++){
//        i = gsrci[l];
//        j = gsrcj[l];
//        pos = (i*NY*NZ + j*NZ + k)*CSIZE;
//
//        tmpz = coord[pos + 2];
//
//        if(fabs(tmpz - sz[l]) < rmin){
//          rmin = (float) fabs(tmpz - sz[l]);
//          gsrck[l] = k;
//        }
//      }

    long offset;
    //float *zline = (float *) malloc(sizeof(float) * NZ);
    float zline;

    float rmin = 1e30, topox, topoy, topoz;

    printf("* SOURCE index search begin\n");
    if(0 == icoord){
      for (l = 0; l < src_npts; l++){
        gsrci[l] = (int)((sx[l] - x0)/DH + 0.5f);
        gsrcj[l] = (int)((sy[l] - y0)/DH + 0.5f);

        i = gsrci[l];
        j = gsrcj[l];

        topox = (i-NX/2)*DH;
        topoy = (j-NY/2)*DH;
        r = powf(topox, 2) +  powf(topoy, 2);
        h = gauss_height;
        a = gauss_width;;
        topoz = h * expf(-r/powf(a,2));

        rmin = 1e30; // reset when searching next source
        for (k = 0; k < NZ; k++){
          zline = topoz + DH*(k-NZ+1);
          if(fabs(zline - sz[l]) < rmin){
            rmin = (float) fabs(zline - sz[l]);
            gsrck[l] = k;
          }
        }
      }
    }else{
      fp = fopen(topofile, "rb");
      for (l = 0; l < src_npts; l++){
        gsrci[l] = (int)((sx[l] - x0)/DH + 0.5f);
        gsrcj[l] = (int)((sy[l] - y0)/DH + 0.5f);

        i = gsrci[l];
        j = gsrcj[l];

        // find z
        offset = ( (i*NY+j) * CSIZE + 0 ) * sizeof(float);

        fseek(fp, offset, SEEK_SET);
        fread(&topox, sizeof(float), 1, fp);
        fread(&topoy, sizeof(float), 1, fp);
        fread(&topoz, sizeof(float), 1, fp);

        //printf("topox = %f, topoy = %f, topoz = %f\n", topox, topoy, topoz);

        rmin = 1e30; // reset when searching next source
        for (k = 0; k < NZ; k++){
          zline = topoz + DH*(k-NZ+1);
          if(fabs(zline - sz[l]) < rmin){
            rmin = (float) fabs(zline - sz[l]);
            gsrck[l] = k;
          }
        }
      }
      fclose(fp);
    }
    printf("* SOURCE index search finish\n");

    printf("* SOURCE index check write begin\n");
    fp = fopen("source_index_check.txt", "w");
    //fprintf(fp, "x         y         z         i     j     k     shiftx    shifty    shiftz\n");
    fprintf(fp, "x         y         z         i     j     k\n");
    //float shift[3];
    for (l = 0; l < src_npts; l++){
      i = gsrci[l];
      j = gsrcj[l];
      k = gsrck[l];
      pos = (i*NY*NZ + j*NZ + k)*CSIZE;
      //shift[0] = sx[l] - coord[pos + 0];
      //shift[1] = sy[l] - coord[pos + 1];
      //shift[2] = sz[l] - coord[pos + 2];
      //fprintf(fp, "%-9.1e %-9.1e %-9.1e %-5d %-5d %-5d %-9.1e %-9.1e %-9.1e\n",
      fprintf(fp, "%-9.1e %-9.1e %-9.1e %-5d %-5d %-5d\n",
          sx[l], sy[l], sz[l], 
          i,j,k);
      //    ,shift[0],shift[1],shift[2] );
    }
    fclose(fp);
    printf("* SOURCE index check write finish\n");
    //free(coord); coord = NULL;
    //free(topo); topo = NULL;

    src_npts_local = (int *) malloc(sizeof(int) * PX * PY * PZ);
    sendcounts     = (int *) malloc(sizeof(int) * PX * PY * PZ);
    displs         = (int *) malloc(sizeof(int) * PX * PY * PZ);
    displs_nt      = (int *) malloc(sizeof(int) * PX * PY * PZ);
    int *src_rank  = (int *) malloc(sizeof(int) * src_npts);
    int *src_sort  = (int *) malloc(sizeof(int) * src_npts);
    int src_npts_local_sum = 0;
    displs[0] = 0;
    displs_nt[0] = 0;
    int myid;
    int gi1, gj1, gk1;
    int gi2, gj2, gk2;
    for (i = 0; i < PX; i++)
      for (j = 0; j < PY; j++)
        for (k = 0; k < PZ; k++){

          myid = i * PY * PZ + j * PZ + k;

          gi1 = i * ni, gi2 = i*ni + ni - 1;
          gj1 = j * nj, gj2 = j*nj + nj - 1;
          gk1 = k * nk, gk2 = k*nk + nk - 1;

          src_npts_local[myid] = 0;
          for (l = 0; l < src_npts; l++){
            if( gsrci[l] >= gi1 && gsrci[l] <= gi2 &&
                gsrcj[l] >= gj1 && gsrcj[l] <= gj2 &&
                gsrck[l] >= gk1 && gsrck[l] <= gk2) {
              src_npts_local[myid]++;
              src_rank[l] = myid;
            }
          }

          if(myid > 0){
            displs[myid] = displs[myid-1] + src_npts_local[myid-1];
            displs_nt[myid] = displs_nt[myid-1] + src_npts_local[myid-1] * src_nt;
          }
          src_npts_local_sum += src_npts_local[myid];

        }

    for(i = 0; i < PX*PY*PZ; i++) sendcounts[i] = src_npts_local[i] * src_nt;


    printf("* SOURCE sort begin\n");
    int *src_rank_sort = (int *) malloc(sizeof(int) * src_npts);
    memcpy(src_rank_sort, src_rank, sizeof(int)*src_npts);
    int i;
    for(i = 0 ; i < src_npts ; i ++) 
      src_sort[i] = i;
    //for(i = 0 ; i < src_npts ; i ++) 
    //  printf("before src_rank_sort[%d] = %d, src_sort[%d] = %d\n", i, src_rank_sort[i], i, src_sort[i]);
    //bubble_sort_int(src_rank_sort, src_sort, src_npts);
    sort_int(src_rank_sort, src_sort, src_npts);
    printf("* SOURCE sort end\n");
    //for(i = 0 ; i < src_npts ; i ++) 
    //  printf("after  src_rank_sort[%d] = %d, src_sort[%d] = %d\n", i, src_rank_sort[i], i, src_sort[i]);

    if(src_npts_local_sum != src_npts){
      printf("! Warning! There are %d/%d source points outside of the compute domain.\n",
          src_npts-src_npts_local_sum, src_npts);
    }

    gsrci_sort = (int *) malloc(sizeof(int)*src_npts);
    gsrcj_sort = (int *) malloc(sizeof(int)*src_npts);
    gsrck_sort = (int *) malloc(sizeof(int)*src_npts);
    garea_sort   = (float *) malloc(sizeof(float)*src_npts);
    gstrike_sort = (float *) malloc(sizeof(float)*src_npts);
    gdip_sort    = (float *) malloc(sizeof(float)*src_npts);
    grake_sort = (float *) malloc(sizeof(float)*src_npts*src_nt);
    grate_sort = (float *) malloc(sizeof(float)*src_npts*src_nt);

    int ii;
    printf("* SOURCE reset begin\n");
    for (i = 0; i < src_npts; i++){
      ii = src_sort[i];
      gsrci_sort[i] = gsrci[ii];
      gsrcj_sort[i] = gsrcj[ii];
      gsrck_sort[i] = gsrck[ii];
      garea_sort[i] = garea[ii];
      gstrike_sort[i] = gstrike[ii];
      gdip_sort[i] = gdip[ii];

      memcpy(&grake_sort[i*src_nt], &grake[ii*src_nt], sizeof(float)*src_nt);
      memcpy(&grate_sort[i*src_nt], &grate[ii*src_nt], sizeof(float)*src_nt);
      //for (int j = 0; j < src_nt; j++){
      //  grake_sort[i*src_nt + j] = grake[ii*src_nt + j];
      //  grate_sort[i*src_nt + j] = grate[ii*src_nt + j];
      //}

    }
    printf("* SOURCE reset end\n");

    if(LOAD_CKPT != 1){ // if LOAD_CKPT == 1, do not need to write source again
      int l;
      int pos_s = 0;
      for (i = 0; i < PX; i++){
        printf("* SOURCE partition write %d/%d\n", i, PX);
        for (j = 0; j < PY; j++){
          myid = i*PY+j;

          nsrc = src_npts_local[myid];
          
          pos_s = displs[myid];

          //if(nsrc>0) printf("pos_s = %07d at rank %06d; nsrc = %d\n", pos_s, myid, nsrc);
          
          sprintf(partfile, "TOPO/tp-PX%d/src-PX%d-PY%d.bin", i, i, j);
          fp = fopen(partfile, "wb");
          fwrite(&nsrc, sizeof(int), 1, fp);
          fwrite(&gsrci_sort[pos_s], sizeof(int), nsrc, fp);
          fwrite(&gsrcj_sort[pos_s], sizeof(int), nsrc, fp);
          fwrite(&gsrck_sort[pos_s], sizeof(int), nsrc, fp);
          fwrite(&garea_sort[pos_s],   sizeof(float), nsrc, fp);
          fwrite(&gstrike_sort[pos_s], sizeof(float), nsrc, fp);
          fwrite(&gdip_sort[pos_s],    sizeof(float), nsrc, fp);
          fwrite(&grake_sort[pos_s*src_nt], sizeof(float), nsrc*src_nt, fp);
          fwrite(&grate_sort[pos_s*src_nt], sizeof(float), nsrc*src_nt, fp);
          fclose(fp);
        }
      }
    }

    free(sx);
    free(sy);
    free(sz);
    free(src_rank);
    free(src_sort);

  } // end masternode

  //MPI_Bcast(&src_npts, 1, MPI_INT, 0, comm);
  MPI_Bcast(&src_nt, 1, MPI_INT,   0, comm);
  MPI_Bcast(&src_dt, 1, MPI_FLOAT, 0, comm);
  if(masternode) printf("* SOURCE scatter begin\n");
  MPI_Scatter(src_npts_local, 1, MPI_INT, &nsrc, 1, MPI_INT, 0, comm);
  if(masternode) printf("* SOURCE scatter finish\n");

  //  for (int i = 0; i < PX*PY*PZ; i++) 
  //    printf("rank = %5d, sendcounts = %5d, displs = %5d, displs_nt = %5d\n",
  //        i, sendcounts[i], displs[i], displs_nt[i]);

  //  DO NOT USE MPI_SCATTERV IN SUNWAY !!!
  //MPI_Barrier(comm);
  //if(masternode) printf("* SOURCE scatterv begin\n");
  //MPI_Scatterv(gsrci_sort, src_npts_local, displs, MPI_INT, srci, nsrc, MPI_INT, 0, comm);
  //MPI_Scatterv(gsrcj_sort, src_npts_local, displs, MPI_INT, srcj, nsrc, MPI_INT, 0, comm);
  //MPI_Scatterv(gsrck_sort, src_npts_local, displs, MPI_INT, srck, nsrc, MPI_INT, 0, comm);
  //MPI_Scatterv(garea_sort,   src_npts_local, displs, MPI_FLOAT, area,   nsrc, MPI_FLOAT, 0, comm);
  //MPI_Scatterv(gstrike_sort, src_npts_local, displs, MPI_FLOAT, strike, nsrc, MPI_FLOAT, 0, comm);
  //MPI_Scatterv(gdip_sort,    src_npts_local, displs, MPI_FLOAT, dip,    nsrc, MPI_FLOAT, 0, comm);
  //MPI_Scatterv(grake_sort, sendcounts, displs_nt, MPI_FLOAT, rake, nsrc*src_nt, MPI_FLOAT, 0, comm);
  //MPI_Scatterv(grate_sort, sendcounts, displs_nt, MPI_FLOAT, rate, nsrc*src_nt, MPI_FLOAT, 0, comm);
  //if(masternode) printf("* SOURCE scatterv finish\n");

  MPI_Barrier(MPI_COMM_WORLD);
  if(masternode) printf("* SOURCE Barrier and partition read begin\n");

  sprintf(partfile, "TOPO/tp-PX%d/src-PX%d-PY%d.bin", thisid[0], thisid[0], thisid[1]);
  fp = fopen(partfile, "rb");
  if(!fp) printf("Error: open %s\n", partfile);
  fread(&nsrc, sizeof(int), 1, fp);
//  if(nsrc > 0){ // malloc 0 is OK
  srci   = (int   *) malloc(sizeof(int)   * nsrc);
  srcj   = (int   *) malloc(sizeof(int)   * nsrc);
  srck   = (int   *) malloc(sizeof(int)   * nsrc);
  area   = (float *) malloc(sizeof(float) * nsrc);
  strike = (float *) malloc(sizeof(float) * nsrc);
  dip    = (float *) malloc(sizeof(float) * nsrc);
  rake   = (float *) malloc(sizeof(float) * nsrc * src_nt);
  rate   = (float *) malloc(sizeof(float) * nsrc * src_nt);
// }
  fread(srci  , sizeof(int),   nsrc, fp);
  fread(srcj  , sizeof(int),   nsrc, fp);
  fread(srck  , sizeof(int),   nsrc, fp);
  fread(area  , sizeof(float), nsrc, fp);
  fread(strike, sizeof(float), nsrc, fp);
  fread(dip   , sizeof(float), nsrc, fp);
  fread(rake  , sizeof(float), nsrc*src_nt, fp);
  fread(rate  , sizeof(float), nsrc*src_nt, fp);
  fclose(fp);
  
  if(masternode) printf("* SOURCE partition read finish\n");

  // transform global src index to local
  for (i = 0; i < nsrc; i++){
    srci[i] = srci[i] % ni + ni1;
    srcj[i] = srcj[i] % nj + nj1;
    srck[i] = srck[i] % nk + nk1;
  }

  if(masternode) printf("* SOURCE global free\n");
  if(NULL != gsrci  ) free(gsrci  );
  if(NULL != gsrcj  ) free(gsrcj  );
  if(NULL != gsrck  ) free(gsrck  );
  if(NULL != garea  ) free(garea  );
  if(NULL != gstrike) free(gstrike);
  if(NULL != gdip   ) free(gdip   );
  if(NULL != grake  ) free(grake  );
  if(NULL != grate  ) free(grate  );
  if(NULL != gsrci_sort  ) free(gsrci_sort  );
  if(NULL != gsrcj_sort  ) free(gsrcj_sort  );
  if(NULL != gsrck_sort  ) free(gsrck_sort  );
  if(NULL != garea_sort  ) free(garea_sort  );
  if(NULL != gstrike_sort) free(gstrike_sort);
  if(NULL != gdip_sort   ) free(gdip_sort   );
  if(NULL != grake_sort  ) free(grake_sort  );
  if(NULL != grate_sort  ) free(grate_sort  );
  if(NULL != src_npts_local ) free(src_npts_local);
  if(NULL != sendcounts ) free(sendcounts);
  if(NULL != displs ) free(displs);
  if(NULL != displs_nt ) free(displs_nt);

  return 0;
}

//int cal_source(int it, float *area, float *strike, float *dip, float *rake, float *rate, float *S){
int cal_source(int it, float *S){

  if(nsrc<1) return 0;
  // S = rate * area * moment tensor ( 6 components)
  float t = it * DT;
  int it1 = (int)(it*DT/src_dt);
  int it2 = it1 + 1;
  float t1 = it1 * src_dt;
  float t2 = it2 * src_dt;
  int i;

  if(it2 > src_nt-1) {
    for (i = 0; i < 6*nsrc; i++) S[i] = 0.0f;
  }else{

    for (i = 0; i < nsrc; i++){
      float s = strike[i];
      float d = dip[i];

      // linear interpolation
      float y1 = rake[i*src_nt + it1];
      float y2 = rake[i*src_nt + it2];
      float r = y1 + (y2 - y1) * (t - t1) / (t2 - t1);
      float moment[6];

      angle2moment(s, d, r, moment);

      y1 = rate[i*src_nt + it1];
      y2 = rate[i*src_nt + it2];
      float rate0 = y1 + (y2 - y1) * (t - t1) / (t2 - t1);
      rate0 *= area[i];

      S[i*6 + 0] = moment[0] * rate0;
      S[i*6 + 1] = moment[1] * rate0;
      S[i*6 + 2] = moment[2] * rate0;
      S[i*6 + 3] = moment[3] * rate0;
      S[i*6 + 4] = moment[4] * rate0;
      S[i*6 + 5] = moment[5] * rate0;
    }
  }

  return 0;
}

int add_source(float *W, float *M, float *S){
//int add_source(float *W, float *M, float *D, float *S){

  int i, j, k, l;
  int pos, pos_m;//, pos_d;
  float jac, mu, V;
  for (l = 0; l < nsrc; l++){
    i = srci[l];
    j = srcj[l];
    k = srck[l];

//    printf("i,j,k = %d %d %d\n", i-3+thisid[0]*ni, j-3+thisid[1]*nj, k-3+thisid[2]*nk);

    pos   = (i * ny * lnz + j * lnz + k) * WSIZE;
    pos_m = (i * ny * lnz + j * lnz + k) * MSIZE;
//    pos_d = (i * ny * nz + j * nz + k) * DSIZE;
    jac = M[pos_m + 9];
//    mu  = D[pos_d + 1];
    mu  = M[pos_m + 11];

    V = mu * DT /(jac * DH * DH * DH);

    W[pos + 3] += -S[l*6 + 0] * V;    // Mxx
    W[pos + 4] += -S[l*6 + 1] * V;    // Myy
    W[pos + 5] += -S[l*6 + 2] * V;    // Mzz
    W[pos + 6] += -S[l*6 + 3] * V;    // Mxy
    W[pos + 7] += -S[l*6 + 4] * V;    // Mxz
    W[pos + 8] += -S[l*6 + 5] * V;    // Myz

  }

  return 0;
}
