/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:init_grid.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-09-06
*   Discription:
*
================================================================*/
#include "headers.h"

void init_grid( PARAMS params, GRID * grid, MPI_COORD thisMPICoord )
{
	int resX = 0;	
	int resY = 0;
	int resZ = 0;	
		
	grid->PX = params.PX;
	grid->PY = params.PY;
	grid->PZ = params.PZ;

	grid->_NX_ = params.NX + 2 * HALO;
	grid->_NY_ = params.NY + 2 * HALO;
	grid->_NZ_ = params.NZ + 2 * HALO;
	
	grid->_NX = params.NX + HALO;
	grid->_NY = params.NY + HALO;
	grid->_NZ = params.NZ + HALO;

	grid->NX = params.NX;
	grid->NY = params.NY;
	grid->NZ = params.NZ;
	//partition strategy
	int marginX = params.NX / params.PX;
	int marginY = params.NY / params.PY;
	int marginZ = params.NZ / params.PZ;
	// int marginX = 165;
	// int marginY = 165;
	// int marginZ = 165;
	// int marginX = (params.PX == 1) ? params.NX : (3*params.NX)/(4*params.PX - 2);
	// int marginY = (params.PY == 1) ? params.NY : (3*params.NY)/(4*params.PY - 2);
	// int marginZ = (params.PZ == 1) ? params.NZ : (3*params.NZ)/(4*params.PZ - 2);
	int partitionX = (params.PX > 2) ? (params.NX - 2*marginX) / (params.PX - 2) : 0;
	int partitionY = (params.PY > 2) ? (params.NY - 2*marginY) / (params.PY - 2) : 0;
	int partitionZ = (params.PZ > 2) ? (params.NZ - 2*marginZ) / (params.PZ - 2) : 0;
	if(thisMPICoord.X == 0 && thisMPICoord.Y == 0 && thisMPICoord.Z == 0) {
		printf("marginX = %d, marginY = %d, marginZ = %d, partitionX = %d, partitionY = %d, partitionZ = %d\n", marginX, marginY, marginZ, partitionX, partitionY, partitionZ);
	}
	//end of partition
	if(thisMPICoord.X == 0 || thisMPICoord.X == params.PX - 1) {
		grid->nx = marginX;
	}
	else {
		grid->nx = partitionX;
	}

	if(thisMPICoord.Y == 0 || thisMPICoord.Y == params.PY - 1) {
		grid->ny = marginY;
	}
	else {
		grid->ny = partitionY;
	}

	if(thisMPICoord.Z == 0 || thisMPICoord.Z == params.PZ - 1) {
		grid->nz = marginZ;
	}
	else {
		grid->nz = partitionZ;
	}

	// grid->nx = params.NX / params.PX;
	// grid->ny = params.NY / params.PY;
	// grid->nz = params.NZ / params.PZ;
	
	resX = params.PX > 2 ? params.NX - 2*marginX - (params.PX - 2)*partitionX : params.NX % marginX;
	resY = params.PY > 2 ? params.NY - 2*marginY - (params.PY - 2)*partitionY : params.NY % marginY;
	resZ = params.PZ > 2 ? params.NZ - 2*marginZ - (params.PZ - 2)*partitionZ : params.NZ % marginZ;

	if ( thisMPICoord.X < resX )
	{
		grid->nx ++;
		grid->frontNX = (thisMPICoord.X == 0) ? 0 : marginX + 1 + (thisMPICoord.X - 1)*(partitionX + 1);
	}	
	else
	{	
		grid->frontNX = (thisMPICoord.X == 0) ? 0 : marginX + 1 + (resX - 1) * ( partitionX + 1 ) + ( thisMPICoord.X - resX ) * partitionX;
	}

	if ( thisMPICoord.Y < resY )
	{
		grid->ny ++;
		grid->frontNY = (thisMPICoord.Y == 0) ? 0 : marginY + 1 + (thisMPICoord.Y - 1)*(partitionY + 1);
	}	
	else
	{	
		grid->frontNY = (thisMPICoord.Y == 0) ? 0 : marginY + 1 + (resY - 1) * ( partitionY + 1 ) + ( thisMPICoord.Y - resY ) * partitionY;
	}

	if ( thisMPICoord.Z < resZ )
	{
		grid->nz ++;
		grid->frontNZ = (thisMPICoord.Z == 0) ? 0 : marginZ + 1 + (thisMPICoord.Z - 1)*(partitionZ + 1);
	}	
	else
	{	
		grid->frontNZ = (thisMPICoord.Z == 0) ? 0 : marginZ + 1 + (resZ - 1) * ( partitionZ + 1 ) + ( thisMPICoord.Z - resZ ) * partitionZ;
	}

	MPI_Barrier( MPI_COMM_WORLD );
	//cout << "X: " << thisMPICoord.X << ", Y: " << thisMPICoord.Y << ", Z: " << thisMPICoord.Z << " ====> frontNX = " << grid->frontNX << ", frontNY = " << grid->frontNY << ", frontNZ = " << grid->frontNZ << endl;
	//cout << "X: " << thisMPICoord.X << ", Y: " << thisMPICoord.Y << ", Z: " << thisMPICoord.Z << " ====> nx = " << grid->nx << ", ny = " << grid->ny << ", nz = " << grid->nz << endl;
	MPI_Barrier( MPI_COMM_WORLD );

	grid->_frontNX = grid->frontNX + HALO;
	grid->_frontNY = grid->frontNY + HALO;
	grid->_frontNZ = grid->frontNZ + HALO;


	grid->_nx = grid->nx + HALO;
	grid->_ny = grid->ny + HALO;
	grid->_nz = grid->nz + HALO;

	grid->_nx_ = grid->nx + 2 * HALO;
	grid->_ny_ = grid->ny + 2 * HALO;
	grid->_nz_ = grid->nz + 2 * HALO;
	
	grid->originalX = params.centerX;
	grid->originalY = params.centerY;

	grid->_originalX = grid->originalX + HALO;
	grid->_originalY = grid->originalY + HALO;

	grid->halo = HALO;

	grid->DH = params.DH;
}
