#ifndef ___SNAPSHOTS_H_JFX93KX8AK20KI
#define ___SNAPSHOTS_H_JFX93KX8AK20KI
#include <sys/stat.h>
#include <sys/types.h>

void write_full_snapshot(unsigned long cur_step, int NX, int NY, int PX, int PY, const char *OUT, int this_rank, MPI_Comm MCW, const int *coord, unsigned int *disp_xyzw, int size, int TARGET_OUTPUT_AT_DEPTH, const float *W);
int mkpath(const char *path, mode_t mode);

#endif
