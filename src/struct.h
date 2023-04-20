/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:struct.h
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-09-04
*   Discription:
*
================================================================*/
#ifndef __STRUCT__
#define __STRUCT__

typedef struct COORD
{	
	float * x;
	float * y;
	float * z;
}COORD;


typedef struct LONLAT
{
	double * lon;
	double * lat;
	double * depth;
}LONLAT;

typedef struct GRID
{

	int PX;
	int PY;
	int PZ;

	int _NX_;
	int _NY_;
	int _NZ_;

	int _NX;
	int _NY;
	int _NZ;

	int NX;
	int NY;
	int NZ;

	int _nx_;
	int _ny_;
	int _nz_;

	int _nx;
	int _ny;
	int _nz;

	int nx;
	int ny;
	int nz;

	int frontNX;
	int frontNY;
	int frontNZ;

	int _frontNX;
	int _frontNY;
	int _frontNZ;

	int originalX;
	int originalY;

	int _originalX;
	int _originalY;
	//int originalZ;

	int halo;
	
	int DH;

}GRID;

typedef struct MPI_COORD
{
	int X;
	int Y;
	int Z;

}MPI_COORD;

typedef struct MPI_NEIGHBOR
{
	int X1; //left
	int X2; //right

	int Y1; //front
	int Y2; //back

	int Z1; //down
	int Z2; //up


	
}MPI_NEIGHBOR;

typedef struct NCFILE
{
	int ncID;

	int ntDimID;	
	int nzDimID;	
	int nyDimID;	
	int nxDimID;
	
	int VxVarID;
	int VyVarID;
	int VzVarID;
	
	int coordXVarID;
	int coordYVarID;
	int coordZVarID;
	
	
	int lonVarID;
	int latVarID;


}NCFILE;

typedef struct SOURCE_FILE_INPUT
{
	
	long long npts; //source point number
	int nt;   //number of source time sequences of every point
	float dt; //time sample interval

	float * lon;
	float * lat;
	float * coordZ;

	float * area;
	float * strike;
	float * dip;

	float * rake;
	float * rate;


}SOURCE_FILE_INPUT;

typedef struct SOURCE_INDEX
{
	int * X;
	int * Y;
	int * Z;

}SOURCE_INDEX;

typedef struct POINT_INDEX
{
	int X;
	int Y;
	int Z;

}POINT_INDEX;

/*
typedef struct MOMENT_RATE
{
	float * Mxx;
	float * Myy;
	float * Mzz;
	float * Mxy;
	float * Mxz;
	float * Myz;
}MOMENT_RATE;
*/
typedef float* MOMENT_RATE;


typedef struct WAVE
{
	float * Vx;
	float * Vy;
	float * Vz;

	float * Txx;
	float * Tyy;
	float * Tzz;
	float * Txy;
	float * Txz;
	float * Tyz;
}WAVE;

typedef struct MEDIUM
{
	float * Vp;
	float * Vs;
	float * rho;
}MEDIUM;


#endif //__STRUCT__
