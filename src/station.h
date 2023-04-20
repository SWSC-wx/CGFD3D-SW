/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:struct.h
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-09-04
*   Discription:
*
================================================================*/
#ifndef __STATION__
#define __STATION__
#include<stdio.h>
#include<stdlib.h>
#include "struct.h"
#include "readJson.h"
#include "cJSON.h"
#include "params.h"
#include "macro.h"
#include "common.h"

typedef struct STATION 
{
	int *X;
	int *Y;
	int *Z;
	float *wave;

}STATION;


void allocStation( STATION * station, int stationNum, int NT );
void freeStation( STATION station );
int readStationIndex( GRID grid );
void initStationIndex( GRID grid, STATION station );
void storageStation( GRID grid, int NT, int stationNum, STATION station, float *W_8, float*W_1, int it );
void stationWrite(PARAMS params, GRID grid, MPI_COORD thisMPICoord, STATION station, int NT, int stationNum);

#endif //__STATION__