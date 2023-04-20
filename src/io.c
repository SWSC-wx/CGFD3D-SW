#include <stdio.h>
#include <stdlib.h>
//#include <mpi.h>
#include "common.h"
#include "params.h"

/*

void paras_msg(void) {

  printf( "_________________________________________________________\n"
          "|   Parameters setting:                                  |\n"
          "|                                                        |\n"
          "| TMAX = %g\n"
          "| DT   = %g\n"
          "| DH   = %g\n"
          "| NX   = %6d;   NY = %6d;   NZ = %6d\n"
          "| PX   = %6d;   PY = %6d;   PZ = %6d\n"
          "| Topography = %s\n"
          "| Interfaces (%d) = %s\n"
          "| Velocity = %s\n"
          "| Source = %s\n"
          "| icoord = %d\n"
          "| isource = %d\n"
          "| TSKP0 = %d\n"
          "| TSKP = %d\n"
          "|________________________________________________________|\n\n"
          ,
          TMAX,DT,DH,NX,NY,NZ,PX,PY,PZ,
          IN_TOPO,
          NUM_INTER,
          IN_INTER,
          IN_VELO,
          IN_SOURCE,
          icoord,
          isource,
          TSKP0,
          TSKP
          );

  return;
}
*/
void paras_msg(void) {

  printf( "_________________________________________________________\n"
          "|   Parameters setting:                                  |\n"
          "|                                                        |\n"
          "| TMAX = %g\n"
          "| DT   = %g\n"
          "| DH   = %g\n"
          "| NX   = %6d;   NY = %6d;   NZ = %6d\n"
          "| PX   = %6d;   PY = %6d;   PZ = %6d\n"
          "| TSKP0 = %d\n"
          "| TSKP = %d\n"
          "|________________________________________________________|\n\n"
          ,
          TMAX,DT,DH,NX,NY,NZ,PX,PY,PZ,
          TSKP0,
          TSKP
          );

  return;
}
/*
int read_Coord(struct Coord *C) {

  Grid2D topox;
  Grid2D topoy;
  Grid2D topoz;

  topox = Alloc2D(ni, nj);
  topoy = Alloc2D(ni, nj);
  topoz = Alloc2D(ni, nj);

  int i, j, k;
  char filex[456];
  char filey[456];
  char filez[456];

  sprintf(filex, "%s/topox.bin", INGRD);
  sprintf(filey, "%s/topoy.bin", INGRD);
  sprintf(filez, "%s/topoz.bin", INGRD);

  Grid2D tpx;
  Grid2D tpy;
  Grid2D tpz;

  tpx = Alloc2D(NX, NY);
  tpy = Alloc2D(NX, NY);
  tpz = Alloc2D(NX, NY);

  FILE *fp;
  fp = fopen(filex, "rb");
  fread(&tpx[0][0], sizeof(float), NX*NY, fp); 
  fclose(fp);

  fp = fopen(filey, "rb");
  fread(&tpy[0][0], sizeof(float), NX*NY, fp); 
  fclose(fp);

  fp = fopen(filez, "rb");
  fread(&tpz[0][0], sizeof(float), NX*NY, fp); 
  fclose(fp);

  for (i = 0; i < ni; i++)
    for (j = 0; j < nj; j++){
      topox[i][j] = tpx[thisid[0]*ni+i][thisid[1]*nj+j];
      topoy[i][j] = tpy[thisid[0]*ni+i][thisid[1]*nj+j];
      topoz[i][j] = tpz[thisid[0]*ni+i][thisid[1]*nj+j];
    }

  //MPI_File_open(SWMPI_COMM, filez, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  //MPI_File_set_view(fh, 0, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
  //MPI_File_read_all(fh, tmpvar, ni*nj, MPI_FLOAT, MPI_STATUS_IGNORE);
  //MPI_File_close(&fh);
  //
  for (i = 0; i < ni; i++)
    for (j = 0; j < nj; j++)
      for (k = 0; k < nk; k++){
        C->x[i+ni1][j+nj1][k+nk1] = topox[i][j];
        C->y[i+ni1][j+nj1][k+nk1] = topoy[i][j];
        C->z[i+ni1][j+nj1][nk2-k-1] = topoz[i][j] - ((PZ-thisid[2]-1)*nk+k) * DH;
      }

  //MPI_Type_free(&filetype);
  //free(tmpvar);
  Delloc2D(tpx);
  Delloc2D(tpy);
  Delloc2D(tpz);
  Delloc2D(topox);
  Delloc2D(topoy);
  Delloc2D(topoz);

  return 0;
}
*/
int write_xplane(char *prefix, Grid3D U, int gi){

  // gi is global index
  char filename[400];
  sprintf(filename, "%s_xplane_%03d.bin", prefix, gi);

  int n_i, i, ii, jj, kk, j, k;
  n_i = gi/ni;
  i = gi%ni; // transform global indx to local indx

  MPI_Group grp, worldgrp;
  MPI_Comm comm2d;

  int *ranks = malloc(sizeof(int) * PY * PZ);

  for (jj = 0; jj < PY; jj++)
    for (kk = 0; kk < PZ; kk++)
      ranks[jj*PZ+kk] = n_i * PY * PZ + jj * PZ + kk;

  MPI_Comm_group(SWMPI_COMM, &worldgrp);
  MPI_Group_incl(worldgrp, PY*PZ, ranks, &grp);
  MPI_Comm_create(SWMPI_COMM, grp, &comm2d);

  if (MPI_COMM_NULL != comm2d) {
    MPI_Datatype filetype;

    int gsize[2] = { NY, NZ};
    int lsize[2] = { nj, nk};
    int offset[2] = {thisid[1]*nj, thisid[2]*nk};

    MPI_Type_create_subarray(2, gsize, lsize, offset, MPI_ORDER_C, MPI_FLOAT, &filetype);
    MPI_Type_commit(&filetype);

    MPI_File fh;

    Grid2D tmpvar;
    tmpvar = Alloc2D(nj, nk);
    for (j = 0; j < nj; j++)
      for (k = 0; k < nk; k++){
        tmpvar[j][k] = U[i+ni1][j+nj1][k+nk1];
      }

    MPI_File_open(comm2d, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY,
        MPI_INFO_NULL, &fh);
    MPI_File_set_view(fh, 0, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
    MPI_File_write_all(fh, &tmpvar[0][0], nj*nk, MPI_FLOAT, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);
    //
    MPI_Type_free(&filetype);
    Delloc2D(tmpvar);
  }
  free(ranks);

  if(MPI_GROUP_NULL != grp) MPI_Group_free(&grp);
  if(MPI_GROUP_NULL != worldgrp) MPI_Group_free(&worldgrp);
  if(MPI_COMM_NULL != comm2d) MPI_Comm_free(&comm2d);

  return 0;
}

int write_yplane(char *prefix, Grid3D U, int gj){
  char filename[400];
  sprintf(filename, "%s_yplane_%03d.bin", prefix, gj);

  int n_j, j, ii, kk, i, k;
  n_j = gj/nj;
  j = gj%nj; // transform global indx to local indx

  Grid2D tmpvar;
  tmpvar = Alloc2D(ni, nk);

  MPI_Group grp, worldgrp;
  MPI_Comm comm2d;

  int *ranks = malloc(sizeof(int) * PX * PZ);

  for (ii = 0; ii < PX; ii++)
    for (kk = 0; kk < PZ; kk++)
      ranks[ii*PZ+kk] = ii * PY * PZ + n_j * PZ + kk;

  MPI_Comm_group(MPI_COMM_WORLD, &worldgrp);
  MPI_Group_incl(worldgrp, PX*PZ, ranks, &grp);
  MPI_Comm_create(MPI_COMM_WORLD, grp, &comm2d);

  if ( MPI_COMM_NULL != comm2d){
    for (i = 0; i < ni; i++)
      for (k = 0; k < nk; k++){
        tmpvar[i][k] = U[i+ni1][j+nj1][k+nk1];
      }

    MPI_Datatype filetype;

    int gsize[2] = { NX, NZ};
    int lsize[2] = { ni, nk};
    int offset[2] = {thisid[0]*ni, thisid[2]*nk};

    MPI_Type_create_subarray(2, gsize, lsize, offset, MPI_ORDER_C, MPI_FLOAT, &filetype);
    MPI_Type_commit(&filetype);

    MPI_File fh;

    MPI_File_open(comm2d, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY,
        MPI_INFO_NULL, &fh);
    MPI_File_set_view(fh, 0, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
    MPI_File_write_all(fh, &tmpvar[0][0], ni*nk, MPI_FLOAT, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);
    //
    Delloc2D(tmpvar);
    free(ranks);
    MPI_Type_free(&filetype);
  }

  if(MPI_GROUP_NULL != grp) MPI_Group_free(&grp);
  if(MPI_GROUP_NULL != worldgrp) MPI_Group_free(&worldgrp);
  if(MPI_COMM_NULL != comm2d) MPI_Comm_free(&comm2d);

  return 0;
}

int write_zplane(char *prefix, Grid3D U, int gk){
  char filename[400];
  sprintf(filename, "%s_zplane_%03d.bin", prefix, gk);

  int n_k, k, ii, jj, i, j;
  n_k = gk/nk;
  k = gk%nk; // transform global indx to local indx

  Grid2D tmpvar;
  tmpvar = Alloc2D(ni, nj);

  MPI_Group grp, worldgrp;
  MPI_Comm comm2d;

  int *ranks = malloc(sizeof(int) * PX * PY);

  for (ii = 0; ii < PX; ii++)
    for (jj = 0; jj < PY; jj++)
      ranks[ii*PY+jj] = ii * PY * PZ + jj * PZ + n_k;

  MPI_Comm_group(MPI_COMM_WORLD, &worldgrp);
  MPI_Group_incl(worldgrp, PX*PY, ranks, &grp);
  MPI_Comm_create(MPI_COMM_WORLD, grp, &comm2d);

  if ( MPI_COMM_NULL != comm2d) {
    for (i = 0; i < ni; i++)
      for (j = 0; j < nj; j++){
        tmpvar[i][j] = U[i+ni1][j+nj1][k+nk1];
      }

    MPI_Datatype filetype;

    int gsize[2] = { NX, NY};
    int lsize[2] = { ni, nj};
    int offset[2] = {thisid[0]*ni, thisid[1]*nj};

    MPI_Type_create_subarray(2, gsize, lsize, offset, MPI_ORDER_C, MPI_FLOAT, &filetype);
    MPI_Type_commit(&filetype);

    MPI_File fh;

    MPI_File_open(comm2d, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY,
        MPI_INFO_NULL, &fh);
    MPI_File_set_view(fh, 0, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
    MPI_File_write_all(fh, &tmpvar[0][0], ni*nj, MPI_FLOAT, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);
    //
    MPI_Type_free(&filetype);
    Delloc2D(tmpvar);
    free(ranks);
  }

  if(MPI_GROUP_NULL != worldgrp) MPI_Group_free(&worldgrp);
  if(MPI_GROUP_NULL != grp) MPI_Group_free(&grp);
  if(MPI_COMM_NULL != comm2d) MPI_Comm_free(&comm2d);

  return 0;
}


int write_xplane_tmp(char *prefix, float *W, int shift, int Unit, int gi){

  // gi is global index
  char filename[400];
  sprintf(filename, "%s_xplane_%03d.bin", prefix, gi);

  int n_i, i, jj, kk, j, k;
  n_i = gi/ni;
  i = gi%ni; // transform global indx to local indx

  MPI_Group grp, worldgrp;
  MPI_Comm comm2d;

  int *ranks = malloc(sizeof(int) * PY * PZ);

  for (jj = 0; jj < PY; jj++)
    for (kk = 0; kk < PZ; kk++)
      ranks[jj*PZ+kk] = n_i * PY * PZ + jj * PZ + kk;

  MPI_Comm_group(SWMPI_COMM, &worldgrp);
  MPI_Group_incl(worldgrp, PY*PZ, ranks, &grp);
  MPI_Comm_create(SWMPI_COMM, grp, &comm2d);

  if (MPI_COMM_NULL != comm2d) {
    MPI_Datatype filetype;

    int gsize[2] = { NY, NZ};
    int lsize[2] = { nj, nk};
    int offset[2] = {thisid[1]*nj, thisid[2]*nk};

    MPI_Type_create_subarray(2, gsize, lsize, offset, MPI_ORDER_C, MPI_FLOAT, &filetype);
    MPI_Type_commit(&filetype);

    MPI_File fh;

    Grid2D tmpvar;
    tmpvar = Alloc2D(nj, nk);
    for (j = 0; j < nj; j++)
      for (k = 0; k < nk; k++){
        //tmpvar[j][k] = U[i+ni1][j+nj1][k+nk1];
    	  int pos = ((i + ni1) * ny * lnz + (j + nj1) * lnz + (k + nk1)) * Unit;
    	  tmpvar[j][k] = W[pos + shift];

      }

    MPI_File_open(comm2d, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY,
        MPI_INFO_NULL, &fh);
    MPI_File_set_view(fh, 0, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
    MPI_File_write_all(fh, &tmpvar[0][0], nj*nk, MPI_FLOAT, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);
    //
    MPI_Type_free(&filetype);
    Delloc2D(tmpvar);
  }
  free(ranks);

  if(MPI_GROUP_NULL != grp) MPI_Group_free(&grp);
  if(MPI_GROUP_NULL != worldgrp) MPI_Group_free(&worldgrp);
  if(MPI_COMM_NULL != comm2d) MPI_Comm_free(&comm2d);

  return 0;
}


int write_yplane_tmp(char *prefix, float *W, int shift, int Unit, int gj){
  char filename[400];
  sprintf(filename, "%s_yplane_%03d.bin", prefix, gj);

  int n_j, j, ii, kk, i, k;
  n_j = gj/nj;
  j = gj%nj; // transform global indx to local indx

  Grid2D tmpvar;
  tmpvar = Alloc2D(ni, nk);

  MPI_Group grp, worldgrp;
  MPI_Comm comm2d;

  int *ranks = malloc(sizeof(int) * PX * PZ);

  for (ii = 0; ii < PX; ii++)
    for (kk = 0; kk < PZ; kk++)
      ranks[ii*PZ+kk] = ii * PY * PZ + n_j * PZ + kk;

  MPI_Comm_group(MPI_COMM_WORLD, &worldgrp);
  MPI_Group_incl(worldgrp, PX*PZ, ranks, &grp);
  MPI_Comm_create(MPI_COMM_WORLD, grp, &comm2d);

  if ( MPI_COMM_NULL != comm2d){
    for (i = 0; i < ni; i++)
      for (k = 0; k < nk; k++){
          //tmpvar[i][k] = U[i+ni1][j+nj1][k+nk1];
      	  int pos = ((i + ni1) * ny * lnz + (j + nj1) * lnz + (k + nk1)) * Unit;
      	  tmpvar[i][k] = W[pos + shift];
      }

    MPI_Datatype filetype;

    int gsize[2] = { NX, NZ};
    int lsize[2] = { ni, nk};
    int offset[2] = {thisid[0]*ni, thisid[2]*nk};

    MPI_Type_create_subarray(2, gsize, lsize, offset, MPI_ORDER_C, MPI_FLOAT, &filetype);
    MPI_Type_commit(&filetype);

    MPI_File fh;

    MPI_File_open(comm2d, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY,
        MPI_INFO_NULL, &fh);
    MPI_File_set_view(fh, 0, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
    MPI_File_write_all(fh, &tmpvar[0][0], ni*nk, MPI_FLOAT, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);
    //
    Delloc2D(tmpvar);
    free(ranks);
    MPI_Type_free(&filetype);
  }

  if(MPI_GROUP_NULL != grp) MPI_Group_free(&grp);
  if(MPI_GROUP_NULL != worldgrp) MPI_Group_free(&worldgrp);
  if(MPI_COMM_NULL != comm2d) MPI_Comm_free(&comm2d);

  return 0;
}



int write_zplane_tmp(char *prefix, float *W, int shift, int Unit, int gk){
  char filename[400];
  sprintf(filename, "%s_zplane_%03d.bin", prefix, gk);

  int n_k, k, ii, jj, i, j;
  n_k = gk/nk;
  k = gk%nk; // transform global indx to local indx

  Grid2D tmpvar;
  tmpvar = Alloc2D(ni, nj);

  MPI_Group grp, worldgrp;
  MPI_Comm comm2d;

  int *ranks = malloc(sizeof(int) * PX * PY);

  for (ii = 0; ii < PX; ii++)
    for (jj = 0; jj < PY; jj++)
      ranks[ii*PY+jj] = ii * PY * PZ + jj * PZ + n_k;

  MPI_Comm_group(MPI_COMM_WORLD, &worldgrp);
  MPI_Group_incl(worldgrp, PX*PY, ranks, &grp);
  MPI_Comm_create(MPI_COMM_WORLD, grp, &comm2d);

  if ( MPI_COMM_NULL != comm2d) {
    for (i = 0; i < ni; i++)
      for (j = 0; j < nj; j++){
          //tmpvar[i][j] = U[i+ni1][j+nj1][k+nk1];
    	  int pos = ((i + ni1) * ny * lnz + (j + nj1) * lnz + (k + nk1)) * Unit;
    	  tmpvar[i][j] = W[pos + shift];
      }

    MPI_Datatype filetype;

    int gsize[2] = { NX, NY};
    int lsize[2] = { ni, nj};
    int offset[2] = {thisid[0]*ni, thisid[1]*nj};

    MPI_Type_create_subarray(2, gsize, lsize, offset, MPI_ORDER_C, MPI_FLOAT, &filetype);
    MPI_Type_commit(&filetype);

    MPI_File fh;

    MPI_File_open(comm2d, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY,
        MPI_INFO_NULL, &fh);
    MPI_File_set_view(fh, 0, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
    MPI_File_write_all(fh, &tmpvar[0][0], ni*nj, MPI_FLOAT, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);
    //
    MPI_Type_free(&filetype);
    Delloc2D(tmpvar);
    free(ranks);
  }

  if(MPI_GROUP_NULL != worldgrp) MPI_Group_free(&worldgrp);
  if(MPI_GROUP_NULL != grp) MPI_Group_free(&grp);
  if(MPI_COMM_NULL != comm2d) MPI_Comm_free(&comm2d);

  return 0;
}



/*
int write_coord_xplane(struct Coord *C, int gi){
  write_xplane("output/x", C->x, gi);
  write_xplane("output/y", C->y, gi);
  write_xplane("output/z", C->z, gi);
  return 0;
}

int write_coord_yplane(struct Coord *C, int gi){
  write_yplane("output/x", C->x, gi);
  write_yplane("output/y", C->y, gi);
  write_yplane("output/z", C->z, gi);
  return 0;
}

int write_coord_zplane(struct Coord *C, int gi){
  write_zplane("output/x", C->x, gi);
  write_zplane("output/y", C->y, gi);
  write_zplane("output/z", C->z, gi);
  return 0;
}
*/
/*
int write_metric_xplane(struct Metric *M, int gi){
  write_xplane("output/jac",  M->jac,  gi);
  write_xplane("output/xi_x", M->xi_x, gi);
  write_xplane("output/xi_y", M->xi_y, gi);
  write_xplane("output/xi_z", M->xi_z, gi);
  write_xplane("output/et_x", M->et_x, gi);
  write_xplane("output/et_y", M->et_y, gi);
  write_xplane("output/et_z", M->et_z, gi);
  write_xplane("output/zt_x", M->zt_x, gi);
  write_xplane("output/zt_y", M->zt_y, gi);
  write_xplane("output/zt_z", M->zt_z, gi);
  return 0;
}
int write_metric_yplane(struct Metric *M, int gi){
  write_yplane("output/jac",  M->jac,  gi);
  write_yplane("output/xi_x", M->xi_x, gi);
  write_yplane("output/xi_y", M->xi_y, gi);
  write_yplane("output/xi_z", M->xi_z, gi);
  write_yplane("output/et_x", M->et_x, gi);
  write_yplane("output/et_y", M->et_y, gi);
  write_yplane("output/et_z", M->et_z, gi);
  write_yplane("output/zt_x", M->zt_x, gi);
  write_yplane("output/zt_y", M->zt_y, gi);
  write_yplane("output/zt_z", M->zt_z, gi);
  return 0;
}
int write_metric_zplane(struct Metric *M, int gi){
  write_zplane("output/jac",  M->jac,  gi);
  write_zplane("output/xi_x", M->xi_x, gi);
  write_zplane("output/xi_y", M->xi_y, gi);
  write_zplane("output/xi_z", M->xi_z, gi);
  write_zplane("output/et_x", M->et_x, gi);
  write_zplane("output/et_y", M->et_y, gi);
  write_zplane("output/et_z", M->et_z, gi);
  write_zplane("output/zt_x", M->zt_x, gi);
  write_zplane("output/zt_y", M->zt_y, gi);
  write_zplane("output/zt_z", M->zt_z, gi);
  return 0;
}
*/
