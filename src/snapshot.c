#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
//#include "recorder.h"
#include "params.h"
#include "snapshot.h"
#include "logger.h"
#include "common.h"

static int do_mkdir(const char *path, mode_t mode)
{
    struct stat     st;
    int             status = 0;

    if (stat(path, &st) != 0)
    {
        /* Directory does not exist. EEXIST for race condition */
        if (mkdir(path, mode) != 0 && errno != EEXIST)
            status = -1;
    }
    else if (!S_ISDIR(st.st_mode))
    {
        errno = ENOTDIR;
        status = -1;
    }

    return(status);
}

/**
** mkpath - ensure all directories in path exist
** Algorithm takes the pessimistic view and works top-down to ensure
** each directory in path exists, rather than optimistically creating
** the last element and working backwards.
*/
int mkpath(const char *path, mode_t mode)
{
    char           *pp;
    char           *sp;
    int             status;
    char           *copypath = strdup(path);

    status = 0;
    pp = copypath;
    while (status == 0 && (sp = strchr(pp, '/')) != 0)
    {
        if (sp != pp)
        {
            /* Neither root nor double slash in path */
            *sp = '\0';
            status = do_mkdir(copypath, mode);
            *sp = '/';
        }
        pp = sp + 1;
    }
    if (status == 0)
        status = do_mkdir(path, mode);
    free(copypath);
    return (status);
}

/* if (cur_step == (long int)(FULL_IMG_LIST[*next_full_output_timestep_idx] / DT)) { */
/* (*next_full_output_timestep_idx)++; */

// void write_full_snapshot(unsigned long cur_step, int NX, int NY, int PX, int PY, const char *OUT, int this_rank, MPI_Comm MCW, const int *coord, unsigned int *disp_xyzw, int size, int TARGET_OUTPUT_AT_DEPTH, const float *W) {
//   int NUM_DIRS = ceilf(1.0 * PX * PY * PZ / NUM_RANK_PER_DIR);
//   char ofvdir[1024];
//   sprintf(ofvdir, "%s/px%dpy%dts%ldndir%d-%d", OUT, PX, PY, cur_step, NUM_DIRS, this_rank % NUM_DIRS);
//   if (this_rank < NUM_DIRS) mkpath(ofvdir, 0700);
//   MPI_Barrier(MCW);


//   if (this_rank == 0) LOG_INFO("writing full output volume\n");
//   static int is_init_metadata = 0;
//   if (!is_init_metadata) {
//     is_init_metadata = 1;
//     /// calculate the displacement of each this_rank
//     unsigned int disp = 0;
//     if (ni * coord[0] >= 1) {
//       disp += (ni * coord[0] - 1) + 1;
//     }
//     if (nj * coord[1] >= 1) {
//       disp += ((nj * coord[1] - 1) + 1) * NX;
//     }
//     disp *= sizeof(float);

//     unsigned int dxyzw_one[5] = {disp, ni, nj, 1, 1};
//     MPI_Allgather(dxyzw_one, 5, MPI_UNSIGNED, disp_xyzw, 5, MPI_UNSIGNED, MCW);
//     if (this_rank == 0) write_metadata(OUT, "full-output-metadata", disp_xyzw, ni * nj, NX, NY, 1, 1, size * 5);
//   }

//   /// every process write its own data
//   char ofvfn[1024];
//   sprintf(ofvfn, "%s/xyz-rank%d.dat", ofvdir, this_rank);

//   int n2 = ny;
//   int n1 = lnz;
//   int k = lnz - nk1 - 1 - TARGET_OUTPUT_AT_DEPTH;
//   int iround;
//   for (iround = 0; iround < NUM_ROUND_OUTPUT; iround++) {
//     if (this_rank % NUM_ROUND_OUTPUT == iround) {
//       if (this_rank == iround) LOG_INFO("writing full volume of round %d/%d", iround+1, NUM_ROUND_OUTPUT);

//       FILE *ofvfp = fopen(ofvfn, "w"); /// output full volume file pointer
//       if (ofvfp == NULL) {
//         LOG_ERROR("[ERROR]: rank %d cannot open file %s (%s)", this_rank, ofvfn, strerror(errno));
//         exit(1);
//       }

//       int i, j;
//       for (j = nj1 ; j <= nj2 - 1 ; j = j + 1) {
//         for (i = ni1 ; i <= ni2 - 1 ; i = i + 1) {
//           unsigned int idx = (i * n1 * n2 + j * n1 + k) * WSIZE;
//           fwrite(&W[idx], sizeof(float), 3, ofvfp);
//         }
//       }

//       fclose(ofvfp);
//     } /// end of round
//     MPI_Barrier(MCW);
//   }

//   if (this_rank == 0) LOG_INFO("writing full output volume finished!\n");
// }

void write_full_snapshot(unsigned long cur_step, int NX, int NY, int PX, int PY, const char *OUT, int this_rank, MPI_Comm MCW, const int *coord, unsigned int *disp_xyzw, int size, int TARGET_OUTPUT_AT_DEPTH, const float *W) {
  int NUM_DIRS = ceilf(1.0 * PX * PY * PZ / NUM_RANK_PER_DIR);
  char ofvdir[1024];
  sprintf(ofvdir, "%s/px%dpy%dts%ldndir%d-%d", OUT, PX, PY, cur_step, NUM_DIRS, this_rank % NUM_DIRS);
  if (this_rank < NUM_DIRS) mkpath(ofvdir, 0700);
  MPI_Barrier(MCW);


  if (this_rank == 0) LOG_INFO("writing full output volume\n");
  static int is_init_metadata = 0;
  if (!is_init_metadata) {
    is_init_metadata = 1;
    /// calculate the displacement of each this_rank
    unsigned int disp = 0;
    if (ni * coord[0] >= 1) {
      disp += (ni * coord[0] - 1) + 1;
    }
    if (nj * coord[1] >= 1) {
      disp += ((nj * coord[1] - 1) + 1) * NX;
    }
    disp *= sizeof(float);

    unsigned int dxyzw_one[5] = {disp, ni, nj, 1, 1};
    MPI_Allgather(dxyzw_one, 5, MPI_UNSIGNED, disp_xyzw, 5, MPI_UNSIGNED, MCW);
    if (this_rank == 0) write_metadata(OUT, "full-output-metadata", disp_xyzw, ni * nj, NX, NY, 1, 1, size * 5);
  }

  /// every process write its own data
  char ofvfn[1024];
  sprintf(ofvfn, "%s/xyz-rank%d.dat", ofvdir, this_rank);

  int n2 = ny;
  int n1 = lnz;
  int k = lnz - nk1 - 1 - TARGET_OUTPUT_AT_DEPTH;
  int iround;
  for (iround = 0; iround < NUM_ROUND_OUTPUT; iround++) {
    if (this_rank % NUM_ROUND_OUTPUT == iround) {
      if (this_rank == iround) LOG_INFO("writing full volume of round %d/%d", iround+1, NUM_ROUND_OUTPUT);

      FILE *ofvfp = fopen(ofvfn, "w"); /// output full volume file pointer
      if (ofvfp == NULL) {
        LOG_ERROR("[ERROR]: rank %d cannot open file %s (%s)", this_rank, ofvfn, strerror(errno));
        exit(1);
      }

      int i, j;
      for (j = nj1 ; j <= nj2 - 1 ; j = j + 1) {
        for (i = ni1 ; i <= ni2 - 1 ; i = i + 1) {
          // unsigned int idx = (i * n1 * n2 + j * n1 + k) * WSIZE;
          unsigned int idx = (i * n1 * n2 + j * n1 + k) * WSIZE_V;
          fwrite(&W[idx], sizeof(float), 3, ofvfp);
        }
      }

      fclose(ofvfp);
    } /// end of round
    MPI_Barrier(MCW);
  }

  if (this_rank == 0) LOG_INFO("writing full output volume finished!\n");
}
