#ifndef PMCL3D_H_JDI39F0EKFJVH3USH
#define PMCL3D_H_JDI39F0EKFJVH3USH
#include <time.h>
typedef int int_pt;
typedef float real;

#define BLOCK_SIZE_X 2
#define BLOCK_SIZE_Y 2
#define BLOCK_SIZE_Z 32
#define align 32
#define loop 1
#define ngsl 8   /* number of ghost cells x loop */
#define ngsl2 16 /* ngsl * 2 */

#define Both 0
#define Left 1
#define Right 2
#define Front 3
#define Back 4

#define NEDZ_EP 160 /*max k to save final plastic strain*/

#define d_c1 (9.0 / 8.0)
#define d_c2 (-1.0 / 24.0)

#define NUM_RANK_PER_DIR (5000)
#define NUM_ROUND_OUTPUT (30) /// 16 is suggested by Wang, Xiyang

#endif /* end of include guard */
