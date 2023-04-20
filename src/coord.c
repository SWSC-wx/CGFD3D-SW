#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "params.h"
#include "common.h"
#include "snapshot.h"
#define root 0

//int gauss_coord(float *C, float gauss_height, float gauss_width);
//int read_coord(float *C);

int gauss_coord(float *C, float gauss_height, float gauss_width){

  int i, j, k;
  for ( i = ni1; i < ni2; i++){
    for ( j = nj1; j < nj2; j++){
      for ( k = nk1; k < nk2; k++){
        int gi = thisid[0]*ni+i-ni1;
        int gj = thisid[1]*nj+j-nj1;
        int gk = thisid[2]*nk+k-nk1;

        int pos = (i*ny*lnz+j*lnz+k)*CSIZE;

        C[pos + 0] = (gi-NX/2) * DH;
        C[pos + 1] = (gj-NY/2) * DH;

        //float r = powf(C[pos + 0], 2) + powf(C[pos + 1], 2);
        float r = powf(C[pos + 0], 2);
        r += powf(C[pos + 1], 2);
        float h = gauss_height;
        float a = gauss_width;
        C[pos + 2] = (gk-NZ+1) * DH + h * expf(-r/powf(a,2));
      }
    }
  }

  return 0;
}

int read_coord_part(float *C){

  if(masternode) printf("read coord APART begin\n");
  char fpname[2014];
  //mkpath("TOPO", 0700);
  //sprintf(fpname, "devideC/tp-PX%s/tp-PX%d-PY%d.bin", );
  sprintf(fpname, "TOPO/tp-PX%d/tp-PX%d-PY%d.bin", thisid[0] , thisid[0], thisid[1]);
  //printf("open file name%s\n", fpname);
  FILE *fp = NULL;
  fp = fopen(fpname, "rb");
  if(fp == NULL) printf("read coord: %s open failed\n", fpname);

  float *local = (float*)malloc(3*ni*nj*sizeof(float));
  fread(local, sizeof(float), 3*ni*nj, fp);
  fclose(fp);
  int i,j,k;

 // printf("thisid:, %d, %d, %d, %s, local[0], local[1], local[2] = %f, %f, %f\n", thisid[0], thisid[1], thisid[2], fpname, local[3*ni*nj-3], local[3*ni*nj-2], local[3*ni*nj-1]);
  for( i = 0; i < ni; i++){
    for ( j = 0; j < nj; j++){
      for ( k = 0; k < nk; k++){
        float tmpx = local[3*(i*nj+j)+0];
        float tmpy = local[3*(i*nj+j)+1];
        float tmpz = local[3*(i*nj+j)+2];
        //printf("tmpx = %f\n", tmpx);
  //      xmin = (xmin<tmpx)?xmin:tmpx;
  //      xmax = (xmax>tmpx)?xmax:tmpx;
  //      ymin = (ymin<tmpy)?ymin:tmpy;
  //      ymax = (ymax>tmpy)?ymax:tmpy;
  //      zmin = (zmin<tmpz)?zmin:tmpz;
  //      zmax = (zmax>tmpz)?zmax:tmpz;

        int gk = thisid[2]*nk+k;

        int pos = ((i+ni1) * ny * lnz + (j+nj1) * lnz + (k+nk1) )*CSIZE;

        C[pos + 0] = tmpx;
        C[pos + 1] = tmpy;
        C[pos + 2] = tmpz + (gk-NZ+1) * DH;
  //       printf("C[pos+2] = %f\n", C[pos+2]);
      }
    }
  }

  if(masternode) printf("read coord: read finished\n");
  return 0;
}

int read_coord2(float *C){
  if(masternode) printf("read coord begin\n");

  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm comm2d, commvert;
  MPI_Group grp, worldgrp;
  float *global, *local; // topo(x,y)

  int rank2d[PX*PY];
  int i, j, k;
  for ( i = 0; i < PX; i++)
    for ( j = 0; j < PY; j++)
      rank2d[i*PY+j] = i*PY*PZ + j*PZ + 0;

  MPI_Comm_group(comm, &worldgrp);
  MPI_Group_incl(worldgrp, PX*PY, rank2d, &grp);
  MPI_Comm_create(comm, grp, &comm2d);

  local = (float*)malloc(3*ni*nj*sizeof(float));

  FILE *fp = NULL;
  if (masternode){
    printf("alloc global %d Bytes begin\n", 3*NX*NY*sizeof(float));
    global = (float*)malloc(3*NX*NY*sizeof(float));
    printf("read coord: alloc global finished\n");
    fp = fopen(IN_TOPO, "rb");
    if(!fp) printf("read coord: IN_TOPO open failed\n");
    fread(global, sizeof(float), 3*NX*NY, fp);
    fclose(fp);
    printf("read coord: read global finished\n");
  }

  int sizes[2] = {NX, 3*NY};
  int subsizes[2] = {ni, 3*nj};
  int starts[2] = {0, 0};
  MPI_Aint lb = 0;
  MPI_Aint extend = 3*nj*sizeof(float);
  MPI_Datatype type, subarrtype;
  MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_FLOAT, &type);
  MPI_Type_create_resized(type, lb, extend, &subarrtype);
  MPI_Type_commit(&subarrtype);
  if(masternode) printf("read coord: subarrtype commit finished\n");

  int counts[PX*PY];
  int displs[PX*PY];

  for ( i = 0; i < PX*PY; i++) counts[i] = 1;
  int disp = 0;
  //if(masternode) printf("read coord: displs:\n");
  for ( i = 0; i < PX; i++){
   // if(masternode) printf("iPX = %d: ", i);
    for ( j = 0; j < PY; j++){
      displs[i*PY+j] = disp;
      disp ++;
     // if(masternode) printf("%d ", displs[i * PY + j]);
    }
    disp += (ni-1)*PY;
   /// if(masternode) printf("\n");
  }
 // if(masternode) printf("\n");

  if(masternode) printf("read coord: MPI_Scatter global begin\n");
  if(MPI_COMM_NULL != comm2d){
    //MPI_Scatterv(global, counts, displs, subarrtype, local,3*ni*nj, MPI_FLOAT, root, comm2d);
    MPI_Scatterv(global, counts, displs, subarrtype, local, 1 , subarrtype, root, comm2d);
    if(masternode) printf("read coord: MPI_Scatter global finished\n");
  }
  if(masternode) printf("read coord: MPI_Scatter global finished v2\n");


  MPI_Barrier(comm);
  return 0;
}
int read_coord(float *C){
  if(masternode) printf("read coord begin\n");

  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm comm2d, commvert;
  MPI_Group grp, worldgrp;
  //int rank, nprocs;
  float *global, *local; // topo(x,y)

  int rank2d[PX*PY];
  int i, j, k;
  for ( i = 0; i < PX; i++)
    for ( j = 0; j < PY; j++)
      rank2d[i*PY+j] = i*PY*PZ + j*PZ + 0;

  MPI_Comm_group(comm, &worldgrp);
  MPI_Group_incl(worldgrp, PX*PY, rank2d, &grp);
  MPI_Comm_create(comm, grp, &comm2d);

  local = (float*)malloc(3*ni*nj*sizeof(float));

  FILE *fp = NULL;
  if (masternode){
    global = (float*)malloc(3*NX*NY*sizeof(float));
    if(masternode) printf("read coord: alloc global finished\n");
    fp = fopen(IN_TOPO, "rb");
    if(!fp) printf("read coord: IN_TOPO open failed\n");
    fread(global, sizeof(float), 3*NX*NY, fp);
    fclose(fp);
    if(masternode) printf("read coord: read global finished\n");
  }

  int sizes[2] = {NX, 3*NY};
  int subsizes[2] = {ni, 3*nj};
  int starts[2] = {0, 0};
  MPI_Aint lb = 0;
  MPI_Aint extend = 3*nj*sizeof(float);
  MPI_Datatype type, subarrtype;
  MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_FLOAT, &type);
  MPI_Type_create_resized(type, lb, extend, &subarrtype);
  MPI_Type_commit(&subarrtype);
  if(masternode) printf("read coord: subarrtype commit finished\n");

  int counts[PX*PY];
  int displs[PX*PY];

  for ( i = 0; i < PX*PY; i++) counts[i] = 1;
  int disp = 0;
 // if(masternode) printf("read coord: displs:\n");
  for ( i = 0; i < PX; i++){
  //  if(masternode) printf("iPX = %d: ", i);
    for ( j = 0; j < PY; j++){
      displs[i*PY+j] = disp;
      disp ++;
  //    if(masternode) printf("%d ", displs[i * PY + j]);
    }
    disp += (ni-1)*PY;
  //  if(masternode) printf("\n");
  }
 // if(masternode) printf("\n");

  if(masternode) printf("read coord: MPI_Scatter global begin\n");
  if(MPI_COMM_NULL != comm2d){
    MPI_Scatterv(global, counts, displs, subarrtype, local,3*ni*nj, MPI_FLOAT, root, comm2d);
    if(masternode) printf("read coord: MPI_Scatter global finished\n");
  }
  if(masternode) printf("read coord: MPI_Scatter global finished v2\n");


  MPI_Barrier(comm);
  int color = thisid[0] * PY + thisid[1];
  int key = thisid[2];
  MPI_Comm_split(comm, color, key, &commvert);

  MPI_Bcast(local, 3*ni*nj, MPI_FLOAT, root, commvert);
  if(masternode) printf("read coord: MPI_Bcast local finished\n");

  //MPI_Barrier(comm);
  //float xmin = 1.0e30, xmax = -1.0e30;
  //float ymin = 1.0e30, ymax = -1.0e30;
  //float zmin = 1.0e30, zmax = -1.0e30;
  for( i = 0; i < ni; i++){
    for ( j = 0; j < nj; j++){
      for ( k = 0; k < nk; k++){
        float tmpx = local[3*(i*nj+j)+0];
        float tmpy = local[3*(i*nj+j)+1];
        float tmpz = local[3*(i*nj+j)+2];
  //      xmin = (xmin<tmpx)?xmin:tmpx;
  //      xmax = (xmax>tmpx)?xmax:tmpx;
  //      ymin = (ymin<tmpy)?ymin:tmpy;
  //      ymax = (ymax>tmpy)?ymax:tmpy;
  //      zmin = (zmin<tmpz)?zmin:tmpz;
  //      zmax = (zmax>tmpz)?zmax:tmpz;

        int gk = thisid[2]*nk+k;

        int pos = ((i+ni1) * ny * lnz + (j+nj1) * lnz + (k+nk1) )*CSIZE;

        C[pos + 0] = tmpx;
        C[pos + 1] = tmpy;
        C[pos + 2] = tmpz + (gk-NZ+1) * DH;

      }
    }
  }
  if(masternode) printf("read coord: construct C finished\n");

//  printf("thisid = %d %d %d, x = %10.2e ~ %10.2e, y = %10.2e ~ %10.2e, z = %10.2e ~ %10.2e\n", 
//      thisid[0], thisid[1], thisid[2], xmin, xmax, ymin, ymax, zmin, zmax);
//
  if(masternode){
    free(global);
  }
  free(local);
  if(masternode) printf("read coord: free finished\n");
  if(MPI_GROUP_NULL != grp) MPI_Group_free(&grp);
  if(MPI_GROUP_NULL != worldgrp) MPI_Group_free(&worldgrp);
  if(MPI_COMM_NULL != comm2d) MPI_Comm_free(&comm2d);
  if(MPI_COMM_NULL != commvert) MPI_Comm_free(&commvert);
  if(MPI_COMM_NULL != comm2d) MPI_Type_free(&subarrtype);
  if(masternode) printf("read coord: MPI free finished\n");


//  printf("thisid = %d %d %d, x = %10.2e ~ %10.2e, y = %10.2e ~ %10.2e, z = %10.2e ~ %10.2e\n", 
//      thisid[0], thisid[1], thisid[2], 
//      xmin, xmax, ymin, ymax, zmin, zmax);

  return 0;
}

int construct_coord(float *C){

  if(icoord==0){
    //float h = 1.0e3;
    //float a = 1.0e3;
    float h = gauss_height;
    float a = gauss_width;
    gauss_coord(C, h, a);
  }else{
    //read_coord(C);
    read_coord_part(C);
  }

  return 0;
}

int cal_range_steph(float *C, float *range){

  float hmin = 1.0e30 , hmax = -1.0e30;
  int pos, pos1, pos2, pos3;
  float L;

  int i, j, k, ii, jj, kk;
  for ( i = ni1; i < ni2; i++)
    for ( j = nj1; j < nj2; j++)
      for ( k = nk1; k < nk2; k++){

        pos = (i*ny*lnz+j*lnz+k)*CSIZE;

        float x0[3], x1[3], x2[3], x3[3];
        x0[0] = C[pos + 0];
        x0[1] = C[pos + 1];
        x0[2] = C[pos + 2];

        for ( ii = -1; ii <= 1; ii = ii+2)
          for ( jj = -1; jj <= 1; jj = jj+2)
            for ( kk = -1; kk <= 1; kk = kk+2){

              pos1 = ((i+ii)*ny*lnz+j*lnz+k)*CSIZE;
              pos2 = (i*ny*lnz+(j+jj)*lnz+k)*CSIZE;
              pos3 = (i*ny*lnz+j*lnz+(k+kk))*CSIZE;

              x1[0] = C[pos1 + 0]; x2[0] = C[pos2 + 0]; x3[0] = C[pos3 + 0];
              x1[1] = C[pos1 + 1]; x2[1] = C[pos2 + 1]; x3[1] = C[pos3 + 1];
              x1[2] = C[pos1 + 2]; x2[2] = C[pos2 + 2]; x3[2] = C[pos3 + 2];

              L = dist_point2plane(x0, x1, x2, x3);
              hmin = _min(hmin, L);
              hmax = _max(hmax, L);

            }
      }

  //printf("rank %4d %4d %4d distance of point2plane = %10.2e ~ %10.2e\n", thisid[0], thisid[1], thisid[2], hmin, hmax);
  range[0] = hmin;
  range[1] = hmax;
  return 0;
}












