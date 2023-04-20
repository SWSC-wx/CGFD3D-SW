#include "struct.h"
#include "readJson.h"
#include <mpi.h>

typedef struct SLICE
{
	int X;
	int Y;
	int Z;
}SLICE;

typedef struct SLICE_DATA
{
	float * x;
	float * y;
	float * z;

}SLICE_DATA;

#define HALO 3
#define WSIZE 9
#define WSIZE_V 8
#define MSIZE 14
#define CSIZE 3

#define Memset memset
#define Free free
#define Memcpy memcpy
#define CHECK( call ) call



#define INDEX( i, j, k ) ( k + ( j ) * ( _nz_ ) + ( i ) * ( _nz_ ) * (  _ny_ ) )
#define INDEX3D( i, j, k, nx, ny, nz ) ( k + ( j ) * ( nz ) + ( i ) * ( nz ) * (  ny ) )
#define INDEX2D( i, j, nx, ny, nz ) ( j + ( i ) * ( ny )) //3D array with only one slice
#define Index2D( i, j, nx, ny ) ( j + ( i ) * ( ny ) )  //2D array


#define CALCULATE2D( i, j, startI, endI, startJ, endJ ) \
for ( i = startI; i < endI; i ++ ) {    \
for ( j = startJ; j < endJ; j ++ ) {    \


#define END_CALCULATE2D( ) }}

#define CALCULATE3D( i, j, k, startI, endI, startJ, endJ, startK, endK ) \
for ( i = startI; i < endI; i ++ ) {     \
for ( j = startJ; j < endJ; j ++ ) {     \
for ( k = startK; k < endK; k ++ ) {     \

#define END_CALCULATE3D( ) }}}

void isMaster(int rank);
//int Malloc( void ** mem, long long  size);
void locateSlice(PARAMS params, GRID grid, SLICE *slice);
void locateFreeSurfSlice(GRID grid, SLICE *slice);
void allocDataout( GRID grid, char XYZ, float ** dataout);
void freeDataout(float *dataout);
void allocSliceData(GRID grid, SLICE slice, SLICE_DATA * sliceData);
//void pack_iodata_x(float *datain, float *dataout, int nx, int ny, int nz, int I);
//void pack_iodata_y(float *datain, float *dataout, int nx, int ny, int nz, int J);
//void pack_iodata_z(float * datain, float * dataout, int nx, int ny, int nz, int K);
//void data2D_output_bin(GRID grid, SLICE slice,  MPI_COORD thisMPICoord, float * datain, SLICE_DATA sliceData, const char *name);
// void data2D_XYZ_out(MPI_COORD thisMPICoord, PARAMS params, GRID grid, float *W, SLICE slice, SLICE_DATA sliceData, char var, int it);
void data2D_XYZ_out( MPI_COORD thisMPICoord, PARAMS params, GRID grid, float * my_W, float * my_W_Tyz, SLICE slice, SLICE_DATA sliceData, char var, int it, int SAMP, int SAMPLE_SIZE ,int sub_rank_x , MPI_Comm sub_comm_x,int sub_rank_y , MPI_Comm sub_comm_y, int Gather );


void data2D_Model_out( MPI_COORD thisMPICoord, PARAMS params, GRID grid, float * coord, float * M, SLICE slice, SLICE_DATA sliceData, int SAMP, int SAMPLE_SIZE ,int sub_rank_x , MPI_Comm sub_comm_x,int sub_rank_y , MPI_Comm sub_comm_y, int Gather );

