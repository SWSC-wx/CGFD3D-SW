#ifndef __PGV__
#define __PGV__
#include<stdio.h>
#include<stdlib.h>
#include "struct.h"
#include "readJson.h"
#include "cJSON.h"
#include "params.h"
#include "macro.h"
#include "common.h"

#define CHECK( call ) call

typedef struct PGV {
	float * pgvh;
	float * pgv;
}PGV;


void allocatePGV(GRID grid, PGV * pgv);
void compare_pgv(PGV pgv, float *W_8, int nx, int ny, int nz);	
void outputPgvData(PARAMS params, MPI_COORD thisMPICoord, PGV Pgv, int nx, int ny);
void comparePGV(GRID grid, MPI_COORD thisMPICoord, float *W_8, PGV pgv);
void outputPGV(PARAMS params, GRID grid, MPI_COORD thisMPICoord, PGV pgv);
void freePGV(PGV pgv);

#endif //__PGV__