#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include "recorder.h"

// Calculates recording points for each core
// rec_nbgxyz rec_nedxyz...
// WARNING: Assumes NPZ = 1! Only surface outputs are needed!
void calcRecordingPoints(int *rec_nbgx, int *rec_nedx, int *rec_nbgy, int *rec_nedy, int *rec_nbgz, int *rec_nedz, int *rec_nxt, int *rec_nyt,
                         int *rec_nzt, MPI_Offset *displacement, long int nxt, long int nyt, long int nzt, int rec_NX, int rec_NY, int rec_NZ,
                         int NBGX, int NEDX, int NSKPX, int NBGY, int NEDY, int NSKPY, int NBGZ, int NEDZ, int NSKPZ, int *coord) {
  *displacement = 0;

  if (NBGX > nxt * (coord[0] + 1))
    *rec_nxt = 0;
  else if (NEDX < nxt * coord[0] + 1)
    *rec_nxt = 0;
  else {
    if (nxt * coord[0] >= NBGX) {
      int mod = (nxt * coord[0] + NBGX - 1) % NSKPX;
      *rec_nbgx = mod == 0 ? 0 : NSKPX - mod;
      *displacement += (nxt * coord[0] - NBGX) / NSKPX + 1;
    } else
      *rec_nbgx = NBGX - nxt * coord[0] - 1;  // since rec_nbgx is 0-based
    if (nxt * (coord[0] + 1) <= NEDX) {
      /* *rec_nedx = (nxt * (coord[0] + 1) + NBGX - 1) % NSKPX - NSKPX + nxt; */
      *rec_nedx = *rec_nbgx;
      while (*rec_nedx + NSKPX < nxt) {
        *rec_nedx += NSKPX;
      }

    } else {
      *rec_nedx = NEDX - nxt * coord[0] - 1;
    }
    *rec_nxt = (*rec_nedx - *rec_nbgx) / NSKPX + 1;
  }

  if (NBGY > nyt * (coord[1] + 1))
    *rec_nyt = 0;
  else if (NEDY < nyt * coord[1] + 1)
    *rec_nyt = 0;
  else {
    if (nyt * coord[1] >= NBGY) {
      int mod = (nyt * coord[1] + NBGY - 1) % NSKPY;
      *rec_nbgy = mod == 0 ? 0 : NSKPY - mod;
      *displacement += ((nyt * coord[1] - NBGY) / NSKPY + 1) * rec_NX;
    } else
      *rec_nbgy = NBGY - nyt * coord[1] - 1;  // since rec_nbgy is 0-based

    if (nyt * (coord[1] + 1) <= NEDY) {
      /* *rec_nedy = (nyt * (coord[1] + 1) + NBGY - 1) % NSKPY - NSKPY + nyt; */
      *rec_nedy = *rec_nbgy;
      while (*rec_nedy + NSKPY < nyt) {
        *rec_nedy += NSKPY;
      }
    } else {
      *rec_nedy = NEDY - nyt * coord[1] - 1;
    }
    *rec_nyt = (*rec_nedy - *rec_nbgy) / NSKPY + 1;
  }

  if (NBGZ > nzt)
    *rec_nzt = 0;
  else {
    *rec_nbgz = NBGZ - 1;  // since rec_nbgz is 0-based
    *rec_nedz = NEDZ - 1;
    *rec_nzt = (*rec_nedz - *rec_nbgz) / NSKPZ + 1;
  }

  if (*rec_nxt == 0 || *rec_nyt == 0 || *rec_nzt == 0) {
    *rec_nxt = 0;
    *rec_nyt = 0;
    *rec_nzt = 0;
  }

  // displacement assumes NPZ=1!
  *displacement *= sizeof(float);

  return;
}

void write_metadata(const char *OUT, const char *fn, unsigned int *disp_xyzw, int max_send_size, int RNX, int RNY, int RNZ, int WRITE_STEP, int size) {
  char file_name[1024];
  sprintf(file_name, "%s/%s", OUT, fn);

  FILE *fp = fopen(file_name, "w");
  if (fp == NULL) {
    fprintf(stderr, "Error: cannot open file %s, errno: %s\n", file_name, strerror(errno));
    exit(0);
  }

  unsigned int buf[5] = {max_send_size, RNX, RNY, RNZ, WRITE_STEP};
  fwrite(buf, sizeof(unsigned int), 5, fp);
  fwrite(disp_xyzw, sizeof(unsigned int), size, fp);
  fclose(fp);
}

void write_model(int this_rank, long int cur_step, char *filenamebase, float *buf, int max_send_size, float *snapshots, int snapshots_size, MPI_Comm MCW) {
  char file_name[1024];
  MPI_Gather(buf, max_send_size, MPI_FLOAT, snapshots, max_send_size, MPI_FLOAT, 0, MCW);

  if (this_rank == 0) {
    //// write SX
    sprintf(file_name, "%s%07ld.dat", filenamebase, cur_step);
    FILE *fp = fopen(file_name, "w");
    if (fp == NULL) {
      fprintf(stderr, "Error: cannot open file %s, errno: %s\n", file_name, strerror(errno));
      exit(0);
    }
    fwrite(snapshots, sizeof(int), snapshots_size, fp);
    fclose(fp);
  }
}
