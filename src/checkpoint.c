#include <string.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "params.h"
#include "logger.h"
#include "snapshot.h"
#include "fastlz.h"

#define BLOCK_SIZE (2*64*1024)
#define MIN(a, b) ((a) < (b) ? (a) : (b))

const static int METHOD = 1;
const static int COMPRESS_LEVEL = 1;

typedef struct {
  void *dat;
  size_t pos;   /// current pos
  size_t size;  /// total size in bytes
} Array;

static Array *new_Array(void *dat, size_t size) {
  Array *ret = (Array *)malloc(sizeof *ret);
  ret->dat = dat;
  ret->pos = 0;
  ret->size = size;

  return ret;
}

static void delete_Array(Array **arr) {
  free(*arr);
  *arr = NULL;
}

static size_t aread(void *ptr, size_t size, size_t nmemb, Array *arr) {
  size_t request_bytes = size * nmemb;
  size_t remain_bytes = arr->size - arr->pos;
  size_t actual_get_bytes = MIN(request_bytes, remain_bytes);

  memcpy(ptr, arr->dat + arr->pos, actual_get_bytes);

  /// update pos
  arr->pos += actual_get_bytes;

  return actual_get_bytes;
}

static size_t awrite(void *ptr, size_t size, size_t nmemb, Array *arr) {
  size_t request_bytes = size * nmemb;
  size_t remain_bytes = arr->size - arr->pos;
  size_t actual_get_bytes = MIN(request_bytes, remain_bytes);

  memcpy(arr->dat + arr->pos, ptr, actual_get_bytes);

  /// update pos
  arr->pos += actual_get_bytes;

  return actual_get_bytes;
}

/* for Adler-32 checksum algorithm, see RFC 1950 Section 8.2 */
#define ADLER32_BASE 65521
static inline unsigned long update_adler32(unsigned long checksum, const void *buf, int len)
{
  const unsigned char* ptr = (const unsigned char*)buf;
  unsigned long s1 = checksum & 0xffff;
  unsigned long s2 = (checksum >> 16) & 0xffff;

  while(len>0)
  {
    unsigned k = len < 5552 ? len : 5552;
    len -= k;

    while(k >= 8)
    {
      s1 += *ptr++; s2 += s1;
      s1 += *ptr++; s2 += s1;
      s1 += *ptr++; s2 += s1;
      s1 += *ptr++; s2 += s1;
      s1 += *ptr++; s2 += s1;
      s1 += *ptr++; s2 += s1;
      s1 += *ptr++; s2 += s1;
      s1 += *ptr++; s2 += s1;
      k -= 8;
    }

    while(k-- > 0)
    {
      s1 += *ptr++; s2 += s1;
    }
    s1 = s1 % ADLER32_BASE;
    s2 = s2 % ADLER32_BASE;
  }
  return (s2 << 16) + s1;
}

static void write_chunk_header(FILE* f, int id, int options, unsigned long size, unsigned long checksum, unsigned long extra) {
  unsigned char buffer[16];

  buffer[0] = id & 255;
  buffer[1] = id >> 8;
  buffer[2] = options & 255;
  buffer[3] = options >> 8;
  buffer[4] = size & 255;
  buffer[5] = (size >> 8) & 255;
  buffer[6] = (size >> 16) & 255;
  buffer[7] = (size >> 24) & 255;
  buffer[8] = checksum & 255;
  buffer[9] = (checksum >> 8) & 255;
  buffer[10] = (checksum >> 16) & 255;
  buffer[11] = (checksum >> 24) & 255;
  buffer[12] = extra & 255;
  buffer[13] = (extra >> 8) & 255;
  buffer[14] = (extra >> 16) & 255;
  buffer[15] = (extra >> 24) & 255;

  fwrite(buffer, 16, 1, f);
}

static inline unsigned long readU16( const unsigned char* ptr )
{
  return ptr[0]+(ptr[1]<<8);
}

static inline unsigned long readU32( const unsigned char* ptr )
{
  return ptr[0]+(ptr[1]<<8)+(ptr[2]<<16)+(ptr[3]<<24);
}

static void read_chunk_header(FILE* f, int* id, int* options, unsigned long* size, unsigned long* checksum, unsigned long* extra) {
  unsigned char buffer[16];
  fread(buffer, 1, 16, f);

  *id = readU16(buffer) & 0xffff;
  *options = readU16(buffer+2) & 0xffff;
  *size = readU32(buffer+4) & 0xffffffff;
  *checksum = readU32(buffer+8) & 0xffffffff;
  *extra = readU32(buffer+12) & 0xffffffff;
}

static int pack_file(const void *ptr, size_t size, int method, int level, FILE *f) {
  unsigned long checksum;
  unsigned char buffer[BLOCK_SIZE];
  unsigned char result[BLOCK_SIZE*2]; /* FIXME twice is too large */
  int chunk_size;

  Array *arr = new_Array((void *)ptr, size);

  /* read file and place in archive */
  for(;;)
  {
    int compress_method = method;
    size_t bytes_read = aread(buffer, 1, BLOCK_SIZE, arr);
    if(bytes_read == 0)
      break;

    /* too small, don't bother to compress */
    if(bytes_read < 32)
      compress_method = 0;

    /* write to output */
    switch(compress_method)
    {
      /* FastLZ */
      case 1:
        chunk_size = fastlz_compress_level(level, buffer, bytes_read, result);
        checksum = update_adler32(1L, result, chunk_size);
        write_chunk_header(f, 17, 1, chunk_size, checksum, bytes_read);
        fwrite(result, 1, chunk_size, f);
        break;

        /* uncompressed, also fallback method */
      case 0:
      default:
        checksum = 1L;
        checksum = update_adler32(checksum, buffer, bytes_read);
        write_chunk_header(f, 17, 0, bytes_read, checksum, bytes_read);
        fwrite(buffer, 1, bytes_read, f);
        break;
    }
  }

  delete_Array(&arr);

  return 0;
}

static int unpack_file(void *ptr, size_t size, size_t fsize, FILE *in) {
  int chunk_id;
  int chunk_options;
  unsigned long chunk_size;
  unsigned long chunk_checksum;
  unsigned long chunk_extra;
  unsigned char buffer[BLOCK_SIZE];
  unsigned long checksum;

  unsigned long decompressed_size;
  unsigned long total_extracted;

  unsigned char* compressed_buffer;
  unsigned char* decompressed_buffer;
  unsigned long compressed_bufsize;
  unsigned long decompressed_bufsize;

  Array *arr = new_Array(ptr, size);

  /* initialize */
  total_extracted = 0;
  decompressed_size = 0;
  compressed_buffer = 0;
  decompressed_buffer = 0;
  compressed_bufsize = 0;
  decompressed_bufsize = 0;

  /* printf("======size: %lu\n", size); */

  /* main loop */
  for(;;)
  {
    /* end of file? */
    size_t pos = ftell(in);
    /* printf("pos %lu\n", pos); */
    if(pos >= fsize)
      break;

    read_chunk_header(in, &chunk_id, &chunk_options,
        &chunk_size, &chunk_checksum, &chunk_extra);

    if((chunk_id == 17))
    {
      unsigned long remaining;

      /* uncompressed */
      switch(chunk_options)
      {
        /* stored, simply copy to output */
        case 0:
          /* read one block at at time, write and update checksum */
          total_extracted += chunk_size;
          remaining = chunk_size;
          checksum = 1L;
          for(;;)
          {
            unsigned long r = (BLOCK_SIZE < remaining) ? BLOCK_SIZE: remaining;
            size_t bytes_read = fread(buffer, 1, r, in);
            if(bytes_read == 0)
              break;
            awrite(buffer, 1, bytes_read, arr);

            checksum = update_adler32(checksum, buffer, bytes_read);
            remaining -= bytes_read;
          }

          /* verify everything is written correctly */
          if(checksum != chunk_checksum)
          {
            printf("\nError: checksum mismatch. Aborted.\n");
            printf("Got %08lX Expecting %08lX\n", checksum, chunk_checksum);
          }
          break;

        /* compressed using FastLZ */
        case 1:
          /* enlarge input buffer if necessary */
          if(chunk_size > compressed_bufsize)
          {
            compressed_bufsize = chunk_size;
            free(compressed_buffer);
            compressed_buffer = (unsigned char*)malloc(compressed_bufsize);
          }

          /* enlarge output buffer if necessary */
          if(chunk_extra > decompressed_bufsize)
          {
            decompressed_bufsize = chunk_extra;
            free(decompressed_buffer);
            decompressed_buffer = (unsigned char*)malloc(decompressed_bufsize);
          }

          /* read and check checksum */
          fread(compressed_buffer, 1, chunk_size, in);
          checksum = update_adler32(1L, compressed_buffer, chunk_size);
          total_extracted += chunk_extra;

          /* verify that the chunk data is correct */
          if(checksum != chunk_checksum)
          {
            printf("\nError: checksum mismatch. Skipped.\n");
            printf("Got %08lX Expecting %08lX\n", checksum, chunk_checksum);
          }
          else
          {
            /* decompress and verify */
            remaining = fastlz_decompress(compressed_buffer, chunk_size, decompressed_buffer, chunk_extra);
            if(remaining != chunk_extra)
            {
              printf("\nError: decompression failed. Skipped.\n");
            }
            else {
              awrite(decompressed_buffer, 1, chunk_extra, arr);
            }
          }
          break;

        default:
          printf("\nError: unknown compression method (%d)\n", chunk_options);
          break;
      }

    }

    /* position of next chunk */
    fseek(in, pos + 16 + chunk_size, SEEK_SET);

    if (arr->pos >= arr->size) {
      break;
    }
  }

  /* free allocated stuff */
  free(compressed_buffer);
  free(decompressed_buffer);

  /* close working files */
  /* fclose(in); */

  delete_Array(&arr);
  /* so far so good */
  return 0;
}

static int safe_fwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream, int this_rank) {
  int num_save = fwrite(ptr, size, nmemb, stream);

  int curr_save_has_error = (num_save != nmemb);
  if (curr_save_has_error) {
    LOG_ERROR("[WARNING]: this_rank %d save %d elemement (expected %lu)", this_rank, num_save, nmemb);
  }
  return curr_save_has_error;
}

static int safe_fread(void *ptr, size_t size, size_t nmemb, FILE *stream, int this_rank) {
  int num_load = fread(ptr, size, nmemb, stream);

  int curr_load_has_error = (num_load != nmemb);
  if (curr_load_has_error) {
    LOG_ERROR("[WARNING]: this_rank %d load %d elemement (expected %lu)", this_rank, num_load, nmemb);
  }
  return curr_load_has_error;
}

//void load_ckpt(int *global_step, int PX, int PY, int PZ, int this_rank, int n1, int n2, int n3, int WSIZE, int MSIZE, int CUR_CKPT_IDX, int NUM_CKPT_DIR, const float *W, const float *hW, const float *mW, const float *tW, const float *M, const float *D, MPI_Comm MCW) {
//void load_ckpt(int *global_step, int PX, int PY, int PZ, int this_rank, int n1, int n2, int n3, int WSIZE, int MSIZE, int CUR_CKPT_IDX, int NUM_CKPT_DIR, const float *W, const float *hW, const float *mW, const float *tW, const float *M, MPI_Comm MCW) {
void load_ckpt(int *global_step, int PX, int PY, int PZ, int this_rank, int n1, int n2, int n3, int WSIZE, int MSIZE, int CUR_CKPT_IDX, int NUM_CKPT_DIR, const float *W_8, const float *W_1, MPI_Comm MCW) {

  //// create the necesary directory
  int NUM_DIRS = ceilf(1.0 * PX * PY * PZ / NUM_RANK_PER_DIR);

  char ofvdir[1024];
  sprintf(ofvdir, "restart/ts-%d-ndir%d-%d", CUR_CKPT_IDX % NUM_CKPT_DIR, NUM_DIRS, this_rank % NUM_DIRS);

  /// every process write its own data
  char ofvfn[1024];
  sprintf(ofvfn, "%s/W-hW-mW-tW-M-D-%d.dat", ofvdir, this_rank);

  /// save the content
  int curr_load_has_error = 0;
  int iround = 0;
  for (iround = 0; iround < NUM_ROUND_OUTPUT; iround++) {
    if (this_rank % NUM_ROUND_OUTPUT == iround) {
      FILE *ofvfp = fopen(ofvfn, "r"); /// output full volume file pointer
      if (ofvfp == NULL) {
        LOG_ERROR("[ERROR]: this_rank %d cannot open file %s (%s)\n", this_rank, ofvfn, strerror(errno));
        exit(1);
      }

      fseek(ofvfp, 0, SEEK_END);
      size_t fsize = ftell(ofvfp);
      fseek(ofvfp, 0, SEEK_SET);

      curr_load_has_error += unpack_file(global_step, sizeof(int) * 1 , fsize, ofvfp);
      // curr_load_has_error += unpack_file(W, sizeof(float) * n1 * n2 * n3 * WSIZE , fsize, ofvfp);
      curr_load_has_error += unpack_file(W_8, sizeof(float) * n1 * n2 * n3 * (WSIZE - 1) , fsize, ofvfp);
      curr_load_has_error += unpack_file(W_1, sizeof(float) * n1 * n2 * n3 , fsize, ofvfp);

//      curr_load_has_error += unpack_file(hW, sizeof(float) * n1 * n2 * n3 * WSIZE , fsize, ofvfp);
//      curr_load_has_error += unpack_file(mW, sizeof(float) * n1 * n2 * n3 * WSIZE , fsize, ofvfp);
//      curr_load_has_error += unpack_file(tW, sizeof(float) * n1 * n2 * n3 * WSIZE , fsize, ofvfp);
//      curr_load_has_error += unpack_file(M, sizeof(float) * n1 * n2 * n3 * MSIZE , fsize, ofvfp);
//      curr_load_has_error += unpack_file(D, sizeof(float) * n1 * n2 * n3, fsize, ofvfp);

      if (this_rank == iround) LOG_INFO("reading checkpoint from %s of round %d/%d @ timestep %d", ofvdir, iround+1, NUM_ROUND_OUTPUT, global_step);

      fclose(ofvfp);
    } /// end of round
    MPI_Barrier(MCW);
  }

  /// reduce the curr_save_has_error
  int sum_error;
  MPI_Reduce(&curr_load_has_error,&sum_error,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
  if (this_rank == 0) {
    if (sum_error) {
      LOG_ERROR("failed to load checkpoint");
      exit(1);
    } else {
      LOG_DEBUG("successfully load checkpoint.");
    }
  }
  MPI_Barrier(MCW);
}

//void save_ckpt(int global_step, int PX, int PY, int PZ, int this_rank, int n1, int n2, int n3, int WSIZE, int MSIZE, int CUR_CKPT_IDX, int NUM_CKPT_DIR, const float *W, const float *hW, const float *mW, const float *tW, const float *M, const float *D, MPI_Comm MCW) {
//void save_ckpt(int global_step, int PX, int PY, int PZ, int this_rank, int n1, int n2, int n3, int WSIZE, int MSIZE, int CUR_CKPT_IDX, int NUM_CKPT_DIR, const float *W, const float *hW, const float *mW, const float *tW, const float *M, MPI_Comm MCW) {
void save_ckpt(int global_step, int PX, int PY, int PZ, int this_rank, int n1, int n2, int n3, int WSIZE, int MSIZE, int CUR_CKPT_IDX, int NUM_CKPT_DIR, const float *W_8, const float *W_1, MPI_Comm MCW) {

  //// create the necesary directory
  int NUM_DIRS = ceilf(1.0 * PX * PY * PZ / NUM_RANK_PER_DIR);
  char ofvdir[1024];
  sprintf(ofvdir, "restart/ts-%d-ndir%d-%d", CUR_CKPT_IDX % NUM_CKPT_DIR, NUM_DIRS, this_rank % NUM_DIRS);
  if (this_rank < NUM_DIRS) mkpath(ofvdir, 0700);
  MPI_Barrier(MCW);

  /// every process write its own data
  char ofvfn[1024];
  sprintf(ofvfn, "%s/W-hW-mW-tW-M-D-%d.dat", ofvdir, this_rank);

  int curr_save_has_error = 0;
  int iround = 0;
  for (iround = 0; iround < NUM_ROUND_OUTPUT; iround++) {
    if (this_rank % NUM_ROUND_OUTPUT == iround) {
      if (this_rank == iround) LOG_INFO("writing checkpoint to %s of round %d/%d @ timestep %ld", ofvdir, iround+1, NUM_ROUND_OUTPUT, global_step);
      FILE *ofvfp = fopen(ofvfn, "w"); /// output full volume file pointer

      if (ofvfp == NULL) {
        LOG_ERROR("[ERROR]: this_rank %d cannot open file %s (%s). This directory is invalid. skipping\n", this_rank, ofvfn, strerror(errno));
        curr_save_has_error++;
      } else { /// open file sucessfully
        curr_save_has_error += pack_file(&global_step, sizeof(int) * 1, METHOD, COMPRESS_LEVEL, ofvfp);
        // curr_save_has_error += pack_file(W, sizeof(float) * n1 * n2 * n3 * WSIZE, METHOD, COMPRESS_LEVEL, ofvfp);
        curr_save_has_error += pack_file(W_8, sizeof(float) * n1 * n2 * n3 * (WSIZE - 1), METHOD, COMPRESS_LEVEL, ofvfp);
        curr_save_has_error += pack_file(W_1, sizeof(float) * n1 * n2 * n3, METHOD, COMPRESS_LEVEL, ofvfp);

//        curr_save_has_error += pack_file(hW, sizeof(float) * n1 * n2 * n3 * WSIZE, METHOD, COMPRESS_LEVEL, ofvfp);
//        curr_save_has_error += pack_file(mW, sizeof(float) * n1 * n2 * n3 * WSIZE, METHOD, COMPRESS_LEVEL, ofvfp);
//        curr_save_has_error += pack_file(tW, sizeof(float) * n1 * n2 * n3 * WSIZE, METHOD, COMPRESS_LEVEL, ofvfp);
//        curr_save_has_error += pack_file(M, sizeof(float) * n1 * n2 * n3 * MSIZE, METHOD, COMPRESS_LEVEL, ofvfp);
//      curr_save_has_error += pack_file(D, sizeof(float) * n1 * n2 * n3, METHOD, COMPRESS_LEVEL, ofvfp);

        fclose(ofvfp);
      } /// end of write
    } /// end of round

    MPI_Barrier(MCW); /// MPI_Barrier for each round

  } /// end of iround

  /// reduce the curr_save_has_error
  int sum_error;
  MPI_Reduce(&curr_save_has_error,&sum_error,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);

  if (this_rank == 0) {
    time_t timer;
    char tbuf[26];
    struct tm* tm_info;

    time(&timer);
    tm_info = localtime(&timer);
    strftime(tbuf, 26, "%Y-%m-%d %H:%M:%S", tm_info);

    const char *logfn = "restart/ckpt.log";
    FILE *fp = fopen(logfn, "a");
    if (fp == NULL) {
      LOG_ERROR("Cannot open file %s (%s)", logfn, strerror(errno));
    }

    if (sum_error == 0) {
      LOG_INFO("writing checkpoint to %s @ timestep %ld successfully.", ofvdir, global_step);
      if (fp != NULL) {
        fprintf(fp, "[%s]: saving checkpoint of %ld timestep to %s [OK]!\n", tbuf, global_step, ofvdir);
      }
    } else {
      if (fp != NULL) {
        fprintf(fp, "[%s]: saving checkpoint of %ld timestep to %s [FAILED]!\n", tbuf, global_step, ofvdir);
      }
      LOG_ERROR("[WARNING] writing checkpoint to %s @ timestep %ld failed! Please don't use this checkpoint!", ofvdir, global_step);
    }

    if (fp != NULL) {
      fclose(fp);
    }
  }
  MPI_Barrier(MCW);
}
