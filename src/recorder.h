#ifndef __RECORD_H_JFKSKAJDFKEJZKCV
#define __RECORD_H_JFKSKAJDFKEJZKCV
#include <mpi.h>

void calcRecordingPoints(int *rec_nbgx, int *rec_nedx, int *rec_nbgy, int *rec_nedy, int *rec_nbgz, int *rec_nedz, int *rec_nxt, int *rec_nyt,
                         int *rec_nzt, MPI_Offset *displacement, long int nxt, long int nyt, long int nzt, int rec_NX, int rec_NY, int rec_NZ,
                         int NBGX, int NEDX, int NSKPX, int NBGY, int NEDY, int NSKPY, int NBGZ, int NEDZ, int NSKPZ, int *coord);

void write_metadata(const char *OUT, const char *fn, unsigned int *disp_xyzw, int max_send_size, int RNX, int RNY, int RNZ, int WRITE_STEP, int size);
void write_model(int this_rank, long int cur_step, char *filenamebase, float *buf, int max_send_size, float *snapshots, int snapshots_size, MPI_Comm MCW);

#endif
