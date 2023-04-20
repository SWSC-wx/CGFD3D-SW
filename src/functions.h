/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:functions.h
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-09-05
*   Discription:
*
================================================================*/
#ifndef __FUNCTIONS__
#define __FUNCTIONS__

/*extern*/ void getParams( PARAMS * params);
double bilinear( double x, double y, double x1, double x2, double y1, double y2, double f11, double f12, double f21, double f22 ); 
double bilinearInterp( double x, double y, double x1, double y1, double x2, double y2, double A, double B, double C, double D );
double interp2d(double x[2], double y[2], double z[4], double x_, double y_ );
void cart2LonLat( GRID grid, PJ * P,  PJ_DIRECTION PD, COORD coord, LONLAT LonLat );
void projTrans( double lon_0, double lat_0, GRID grid, float * coord, LONLAT LonLat  );
/*extern*/ void preprocessTerrain( PARAMS params, MPI_Comm comm_cart, MPI_COORD thisMPICoord, GRID grid, float * C );
/*extern*/ void init_grid( PARAMS params, GRID * grid, MPI_COORD thisMPICoord );

/*extern*/ void allocCoord( GRID grid, COORD * coord );
/*extern*/ void freeCoord( COORD coord );
/*extern*/ void constructCoord( GRID grid, MPI_COORD thisMPICoord, float * coord );

/*extern*/ void finalize_MPI( MPI_Comm * comm_cart );
/*extern*/ void init_MPI( int *argc, char *** argv, PARAMS params, MPI_Comm * comm_cart, MPI_COORD * thisMPICoord, MPI_NEIGHBOR * mpiNeigbor );

/*extern*/ void data2D_output_bin( PARAMS params, GRID grid, MPI_COORD thisMPICoord, float * data, const char * name );
/*extern*/ void data2D_output_nc( PARAMS params, GRID grid, MPI_COORD thisMPICoord, COORD coord );
void dealSource( PARAMS params, GRID grid, float * coord, MPI_COORD thisMPICoord, long long ** srcIndex, long long ** srcSlicePos, MOMENT_RATE *momentRate, SOURCE_FILE_INPUT * src_in_);
/*extern*/ void freeSourceIndex( SOURCE_INDEX srcIndex );
/*extern*/ void readSource( PARAMS params, GRID grid, MPI_COORD thisMPICoord );
int constructMedium(GRID grid, float *M);
void dealWeisenShenMedium( PARAMS params, GRID grid, float * medium, float * coord, MPI_COORD thisMPICoord);
void allocReadMomentRate( SOURCE_FILE_INPUT src_in, MOMENT_RATE momentRate );
void freeReadMomentRate( MOMENT_RATE momentRate );
MOMENT_RATE allocMomentRateSlice( SOURCE_FILE_INPUT src_in);
MOMENT_RATE allocMomentRate(long long pointNum, int nt );
void freeMomentRateSlice( MOMENT_RATE momentRateSlice  );
void interpMomentRate( SOURCE_FILE_INPUT src_in, long long * srcSlicePos, MOMENT_RATE momentRate, MOMENT_RATE momentRateSlice, float t_weight, int srcIt );
// void addMomentRate( float * W, GRID grid, MOMENT_RATE momentRateSlice, long long * srcIndex, long long * srcSlicePos, int npts, float DH, float DT, float * M, int it );
void addMomentRate( float * my_W, float * my_W_Tyz, GRID grid,  MOMENT_RATE momentRateSlice, long long * srcIndex, long long * srcSlicePos, int npts, float DH, float DT, float * M, int it );
void verifySrcIndex( MPI_COORD thisMPICoord, long long * srcIndex, int npts  );
#endif //__FUNCTIONS__
