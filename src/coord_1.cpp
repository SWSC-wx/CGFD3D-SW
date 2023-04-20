/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:coord.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-09-06
*   Discription:
*
================================================================*/
#include "headers.h"


void constructCoord( GRID grid, MPI_COORD thisMPICoord, float * C )
{
	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;

	int _nx = grid._nx;
	int _ny = grid._ny;
	int _nz = grid._nz;

	int frontNX = grid.frontNX;
	int frontNY = grid.frontNY;
	int frontNZ = grid.frontNZ;

	int NZ = grid.NZ;

	int originalX = grid.originalX;
	int originalY = grid.originalY;
	

	int halo = grid.halo;

	float DH = grid.DH;


	int i = 0, j = 0, k = 0;
	long long index = 0, pos;
	int I = 0, J = 0, K = 0;
	
	int nx = grid.nx;
	int ny = grid.ny;
	int nz = grid.nz;

	FOR_LOOP3D( i, j, k, 0, _nx_, 0, _ny_, 0, _nz_ )
        //index = INDEX( i, j, k );
        index = INDEX( i, j, k ) * CSIZE;
		I = frontNX + i;
		J = frontNY + j;
		K = frontNZ + k;
		C[index + 0] = ( I - halo ) * DH - originalX * DH;
		C[index + 1] = ( J - halo ) * DH - originalY * DH;
		C[index + 2] = ( K - halo ) * DH - NZ * DH ; 
    END_LOOP3D( )

}
