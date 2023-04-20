/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:propagate.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-11-20
*   Discription:
*
================================================================*/
#include "pgv.h"

void allocatePGV(GRID grid, PGV * pgv)
{
	int nx = grid.nx;
	int ny = grid.ny;


	int len = sizeof(float)*nx*ny*2;
	float * pPgv = NULL;

	CHECK( Malloc( ( void ** )&pPgv, len) );
	CHECK( memset( pPgv, 0, len ));

	int num = nx * ny;

	pgv->pgvh = pPgv;
	pgv->pgv  = pgv->pgvh + num;
}

void freePGV(PGV pgv)
{
	free( pgv.pgvh );
}

void compare_pgv(PGV pgv, float *W_8, int nx, int ny, int nz)			
{
	int _nx_ = nx + 2 * HALO;
	int _ny_ = ny + 2 * HALO;
	int _nz_ = nz + 2 * HALO;

	int i0 = 0;
	int j0 = 0;

	int i = i0 + HALO;
	int j = j0 + HALO;
	int k = _nz_ - HALO - 1;

	long long index_v, index, pos;

	float Vx = 0.0f, Vy = 0.0f, Vz = 0.0f, Vh = 0.0f, V = 0.0f;

	double c = 1.0 / Cv;
	CALCULATE2D( i0, j0, 0, nx, 0, ny )
		i = i0 + HALO;
		j = j0 + HALO;
		index_v = INDEX( i, j, k )*WSIZE_V;	
		index = INDEX( i, j, k );	
		pos = Index2D( i0, j0, nx, ny );
			
		//if ( i0 == nx - 1 && j0 == ny - 1 )
		//	printf( "nx = %d, ny = %d, pos = %d, pgvh = %p, pgv = %p\n", nx, ny, pos, pgv.pgvh, pgv.pgv );

		Vx = W_8[index_v + 0] * c;
		Vy = W_8[index_v + 1] * c;
		Vz = W_8[index_v + 2] * c;
		
		Vh = sqrtf( Vx * Vx + Vy * Vy );
		V  = sqrtf( Vx * Vx + Vy * Vy + Vz * Vz);
		
		if ( pgv.pgvh[pos] < Vh )
		{
			pgv.pgvh[pos] = Vh;
		}
		if ( pgv.pgv[pos] < V )
		{
			pgv.pgv[pos] = V;
		}
	END_CALCULATE2D( )
}

void outputPgvData(PARAMS params, MPI_COORD thisMPICoord, PGV Pgv, int nx, int ny)
{
	char fileName1[1024] = { 0 };
	char fileName2[1024] = { 0 };
	FILE * filePgvh = NULL;
	FILE * filePgv  = NULL;

	int i, j;

	sprintf( fileName1, "%s/PGVh_Z_mpi_%d_%d_%d.bin", params.OUT, thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );
	sprintf( fileName2, "%s/PGV_Z_mpi_%d_%d_%d.bin",  params.OUT, thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );
	filePgvh = fopen( fileName1, "wb" );
	filePgv  = fopen( fileName2, "wb" );

	float vm = 0.0;

	fwrite(Pgv.pgvh, sizeof( float ), nx * ny, filePgvh);
	fwrite(Pgv.pgv , sizeof( float ), nx * ny, filePgv);

	fclose(filePgvh);
	fclose(filePgv);
}



void comparePGV(GRID grid, MPI_COORD thisMPICoord, float *W_8, PGV pgv)
{
	int nx = grid.nx;
	int ny = grid.ny;
	int nz = grid.nz;
	compare_pgv(pgv, W_8, nx, ny, nz);
}

void outputPGV(PARAMS params, GRID grid, MPI_COORD thisMPICoord, PGV pgv)
{

	int nx = grid.nx;
	int ny = grid.ny;

	int size;
	size = sizeof( float ) * nx * ny * 2;
	outputPgvData( params, thisMPICoord, pgv, nx, ny );
}	