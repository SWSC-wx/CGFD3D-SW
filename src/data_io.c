/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:data_io.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-11-06
*   Discription:
*
================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "data_io.h"
static int masternode;
static int myrank;
void isMaster(int rank) 
{
	myrank = rank;
	if (rank == 0) {
	    masternode == 1;
		printf("this node is master\n");
	} else {
	    masternode == 0;
	}
}

int Malloc( void ** mem, long long  size)
{
	*mem = malloc( size );
	if ( *mem == NULL )
	{
		printf( "can not malloc, Error: %s:%d\n", __FILE__, __LINE__ );
	}
	return 0;
}

void locateSlice( PARAMS params, GRID grid, SLICE * slice )
{
	int sliceX = params.sliceX;
	int sliceY = params.sliceY;
	int sliceZ = params.sliceZ;
	
	int frontNX = grid.frontNX;
	int frontNY = grid.frontNY;
	int frontNZ = grid.frontNZ;
	
 	slice->X = sliceX - frontNX + HALO;
 	slice->Y = sliceY - frontNY + HALO;
 	slice->Z = sliceZ - frontNZ + HALO;
	
	//printf( "slice.X = %d, slice.Y = %d, slice.Z = %d\n", slice->X, slice->Y, slice->Z );
	
}

void locateFreeSurfSlice( GRID grid, SLICE * slice )
{
	int sliceX = -1;
	int sliceY = -1;
	int sliceZ = grid.NZ - 1;
	
	int frontNX = grid.frontNX;
	int frontNY = grid.frontNY;
	int frontNZ = grid.frontNZ;
	
 	slice->X = sliceX - frontNX + HALO;
 	slice->Y = sliceY - frontNY + HALO;
 	slice->Z = sliceZ - frontNZ + HALO;
	
	//printf( "slice.X = %d, slice.Y = %d, slice.Z = %d\n", slice->X, slice->Y, slice->Z );
	
}


void allocDataout( GRID grid, char XYZ, float ** dataout )
{
	int nx = grid.nx;
	int ny = grid.ny;
	int nz = grid.nz;

	long long num = 0;
	
	switch( XYZ )
	{
		case 'X':
			num = ny * nz;
			break;
		case 'Y':
			num = nx * nz;
			break;
		case 'Z':
			num = nx * ny;
			break;
	}
		
	float * pData = NULL;
	long long size = sizeof( float ) * num;

	CHECK( Malloc( ( void ** )&pData, size ) );
	CHECK( Memset(  pData, 0, size ) ); 
	
	*dataout = pData;
}

void freeDataout( float * dataout  )
{

	Free( dataout  );
}


void pack_iodata_x( float * datain, float * dataout, int nx, int ny, int nz, int I, int SIZE  )
{
	int _nx_ = nx + 2 * HALO;
	int _ny_ = ny + 2 * HALO;
	int _nz_ = nz + 2 * HALO;

	//printf( "datain = %p\n", datain  );

	int j0 = 0;
	int k0 = 0;
	
	int i = I;
	int j = j0 + HALO;
	int k = k0 + HALO;

	long long index, pos;
	
	CALCULATE2D( j0, k0, 0, ny, 0, nz )
		j = j0 + HALO;
		k = k0 + HALO;
		index = INDEX( i, j, k ) * SIZE;	
		pos = Index2D( j0, k0, ny, nz );
		dataout[pos] = datain[index];
		//printf( "1:datain = %e\n", datain[pos]  );
	END_CALCULATE2D( )

}


void pack_iodata_y( float * datain, float * dataout, int nx, int ny, int nz, int J, int SIZE  )
{
	int _nx_ = nx + 2 * HALO;
	int _ny_ = ny + 2 * HALO;
	int _nz_ = nz + 2 * HALO;

	int i0 = 0;
	int k0 = 0;
	
	int i = i0 + HALO;
	int j = J;
	int k = k0 + HALO;

	long long index, pos;
	
	CALCULATE2D( i0, k0, 0, nx, 0, nz )
		i = i0 + HALO;
		k = k0 + HALO;
		index = INDEX( i, j, k ) * SIZE;	
		pos = Index2D( i0, k0, nx, nz );
		dataout[pos] = datain[index];
		//printf( "2:datain = %e\n", datain[pos]  );
	END_CALCULATE2D( )

}


void pack_iodata_z( float * datain, float * dataout, int nx, int ny, int nz, int K, int SIZE  )
{
	int _nx_ = nx + 2 * HALO;
	int _ny_ = ny + 2 * HALO;
	int _nz_ = nz + 2 * HALO;

	int i0 = 0;
	int j0 = 0;
	
	int i = i0 + HALO;
	int j = j0 + HALO;
	int k = K;

	long long index, pos;
	
	CALCULATE2D( i0, j0, 0, nx, 0, ny )
		i = i0 + HALO;
		j = j0 + HALO;
		index = INDEX( i, j, k ) * SIZE;	
		pos = Index2D( i0, j0, nx, ny );
		dataout[pos] = datain[index];
		//printf( "3:datain = %e\n", datain[index]  );
	END_CALCULATE2D( )


}

void allocSliceData( GRID grid, SLICE slice, SLICE_DATA * sliceData )
{
	int _nx = grid._nx;
	int _ny = grid._ny;
	int _nz = grid._nz;

	if ( slice.X >= HALO && slice.X < _nx )
		allocDataout( grid, 'X', &( sliceData->x ) );
	if ( slice.Y >= HALO && slice.Y < _ny )
		allocDataout( grid, 'Y', &( sliceData->y ) );
	if ( slice.Z >= HALO && slice.Z < _nz )
		allocDataout( grid, 'Z', &( sliceData->z ) );
}



void freeSliceData( GRID grid, SLICE slice, SLICE_DATA sliceData )
{
	int _nx = grid._nx;
	int _ny = grid._ny;
	int _nz = grid._nz;

	if ( slice.X >= HALO && slice.X < _nx )
		freeDataout( sliceData.x );
	if ( slice.Y >= HALO && slice.Y < _ny )
		freeDataout( sliceData.y );
	if ( slice.Z >= HALO && slice.Z < _nz )
		freeDataout( sliceData.z );
}

/*
void data2D_output_bin( GRID grid, SLICE slice, 
						MPI_COORD thisMPICoord, 
						float * datain, int SIZE, SLICE_DATA sliceData, 
						const char * name )
{

	int _nx = grid._nx;
	int _ny = grid._ny;
	int _nz = grid._nz;

	int nx = grid.nx;
	int ny = grid.ny;
	int nz = grid.nz;


	if ( slice.X >= HALO && slice.X < _nx )
	{

		pack_iodata_x( datain, sliceData.x, nx, ny, nz, slice.X, SIZE  );


		FILE * fp;
		char fileName[256];
		sprintf( fileName, "%s_X_mpi_%d_%d_%d.bin", name, thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );

		fp = fopen( fileName, "wb" ); 
		
		fwrite( sliceData.x, sizeof( float ), ny * nz, fp );
		
		fclose( fp );
	}

	if ( slice.Y >= HALO && slice.Y < _ny )
	{

		pack_iodata_y( datain, sliceData.y, nx, ny, nz, slice.Y, SIZE  );

		FILE * fp;
		char fileName[256];
		sprintf( fileName, "%s_Y_mpi_%d_%d_%d.bin", name, thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );

		fp = fopen( fileName, "wb" ); 

		fwrite( sliceData.y, sizeof( float ), nx * nz, fp );
		
		fclose( fp );
	}

	if ( slice.Z >= HALO && slice.Z < _nz )
	{

		pack_iodata_z( datain, sliceData.z, nx, ny, nz, slice.Z, SIZE  );

		FILE * fp;
		char fileName[256];
		sprintf( fileName, "%s_Z_mpi_%d_%d_%d.bin", name, thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );

		fp = fopen( fileName, "wb" ); 
		
		fwrite( sliceData.z, sizeof( float ), nx * ny, fp );

		fclose( fp );
	}

}
*/
void data2D_output_bin( GRID grid, SLICE slice, 
						MPI_COORD thisMPICoord, 
						float * datain, int SIZE, SLICE_DATA sliceData, 
						const char * name, int SAMP, int SAMPLE_SIZE , int sub_rank_x , MPI_Comm sub_comm_x,int sub_rank_y , MPI_Comm sub_comm_y, int Gather)
{

	int _nx = grid._nx;
	int _ny = grid._ny;
	int _nz = grid._nz;

	int nx = grid.nx;
	int ny = grid.ny;
	int nz = grid.nz;
	int i;
	int j;
	float *sub_recv_x;
	float *sub_recv_y;
	float *sub_recv_z;
	int sub_size_x;
	int sub_size_y;
    MPI_Comm_size(sub_comm_x, &sub_size_x);
	MPI_Comm_size(sub_comm_y, &sub_size_y);

	if ( slice.X >= HALO && slice.X < _nx )
	{

		pack_iodata_x( datain, sliceData.x, nx, ny, nz, slice.X, SIZE  );
		if(Gather){
			FILE * fp;
			char fileName[256];
			if(sub_rank_x == 0){
				sprintf( fileName, "%s_X_mpi_%d_%d_%d.bin", name, thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );
				fp = fopen( fileName, "wb" ); 
				sub_recv_x = (float *)malloc(sub_size_x*ny*nz*sizeof(float));
			}

			MPI_Gather(sliceData.x,ny*nz,MPI_FLOAT,sub_recv_x,ny*nz,MPI_FLOAT,0,sub_comm_x);


			if(sub_rank_x == 0){
				if(SAMP){
					for(j = 0;j<ny*sub_size_x;j+=SAMPLE_SIZE){
						for(i = 0;i<nz;i+=SAMPLE_SIZE){
							fwrite( sub_recv_x + i + j*nz, sizeof( float ), 1, fp );
						}
					}
				}
				else{
					fwrite( sub_recv_x, sizeof( float ), ny*nz*sub_size_x , fp );
				}
				fclose( fp );
				free(sub_recv_x);
			}
		}else{	
			FILE * fp;
			char fileName[256];
			sprintf( fileName, "%s_X_mpi_%d_%d_%d.bin", name, thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );
			fp = fopen( fileName, "wb" ); 		
			if(SAMP){
				for(j = 0;j<ny;j+=SAMPLE_SIZE){
					for(i = 0;i<nz;i+=SAMPLE_SIZE){
						fwrite( sliceData.x + i + j*nz, sizeof( float ), 1, fp );
					}
				}
			}
			else{
				fwrite( sliceData.x, sizeof( float ), ny * nz, fp );	
			}
			fclose( fp );
		}
		
	}

	if ( slice.Y >= HALO && slice.Y < _ny )
	{

		pack_iodata_y( datain, sliceData.y, nx, ny, nz, slice.Y, SIZE  );
		if(Gather){
			FILE * fp;
			char fileName[256];
			if(sub_rank_y == 0){
				sprintf( fileName, "%s_Y_mpi_%d_%d_%d.bin", name, thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );
				fp = fopen( fileName, "wb" ); 
				sub_recv_y = (float *)malloc(sub_size_y*nx*nz*sizeof(float));
			}

			MPI_Gather(sliceData.y,nx*nz,MPI_FLOAT,sub_recv_y,nx*nz,MPI_FLOAT,0,sub_comm_y);

			if(sub_rank_y == 0){
				if(SAMP){
					for(j = 0;j<nx*sub_size_y;j+=SAMPLE_SIZE){
						for(i = 0;i<nz;i+=SAMPLE_SIZE){
							fwrite( sub_recv_y + i + j*nz, sizeof( float ), 1, fp );
						}
					}
				}
				else{
					fwrite( sub_recv_y, sizeof( float ), nx*nz*sub_size_y , fp );
				}
				fclose( fp );
				free(sub_recv_y);
			}
		}else{
			FILE * fp;
			char fileName[256];
			sprintf( fileName, "%s_Y_mpi_%d_%d_%d.bin", name, thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );
			fp = fopen( fileName, "wb" ); 
			if(SAMP){
				for(j = 0;j<nx;j+=SAMPLE_SIZE){
					for(i = 0;i<nz;i+=SAMPLE_SIZE){
						fwrite( sliceData.y + i + j*nz, sizeof( float ), 1, fp );
					}
				}
			}
			else{
				fwrite( sliceData.y, sizeof( float ), nx * nz, fp );
			}
			fclose( fp );
		}
	}

	if ( slice.Z >= HALO && slice.Z < _nz )
	{

		pack_iodata_z( datain, sliceData.z, nx, ny, nz, slice.Z, SIZE  );
		if(Gather){
			FILE * fp;
			char fileName[256];
			if(sub_rank_y == 0){
				sprintf( fileName, "%s_Z_mpi_%d_%d_%d.bin", name, thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );
				fp = fopen( fileName, "wb" ); 
				sub_recv_z = (float *)malloc(sub_size_y*nx*ny*sizeof(float));
			}

			MPI_Gather(sliceData.z,nx*ny,MPI_FLOAT,sub_recv_z,nx*ny,MPI_FLOAT,0,sub_comm_y);

			if(sub_rank_y == 0){
				if(SAMP){
					for(j = 0;j<nx*sub_size_y;j+=SAMPLE_SIZE){
						for(i = 0;i<ny;i+=SAMPLE_SIZE){
							fwrite( sub_recv_z + i + j*ny, sizeof( float ), 1, fp );
						}
					}
				}
				else{
					fwrite( sub_recv_z, sizeof( float ), nx*ny*sub_size_y , fp );
				}
				fclose( fp );
				free(sub_recv_z);
			}
		}else{
			FILE * fp;
			char fileName[256];
			sprintf( fileName, "%s_Z_mpi_%d_%d_%d.bin", name, thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );

			fp = fopen( fileName, "wb" ); 
			if(SAMP){
				for(j = 0;j<nx;j+=SAMPLE_SIZE){
					for(i = 0;i<ny;i+=SAMPLE_SIZE){
						fwrite( sliceData.z + i + j*ny, sizeof( float ), 1, fp );
					}
				}
			}
			else{
				fwrite( sliceData.z, sizeof( float ), nx * ny, fp );
			}
			fclose( fp );
		}
	}

}


// void data2D_XYZ_out( MPI_COORD thisMPICoord, PARAMS params, GRID grid, float * W, SLICE slice, SLICE_DATA sliceData, char var, int it )
// {

// 	switch ( var )
// 	{
// 		case 'V':
// 			{

// 				char VxFileName[128], VyFileName[128], VzFileName[128];

// 				sprintf( VxFileName, "%s/Vx_%d", params.OUT, it );
// 				sprintf( VyFileName, "%s/Vy_%d", params.OUT, it );
// 				sprintf( VzFileName, "%s/Vz_%d", params.OUT, it );

// 				data2D_output_bin( grid, slice, thisMPICoord, &W[0], WSIZE, sliceData, VxFileName );
// 				data2D_output_bin( grid, slice, thisMPICoord, &W[1], WSIZE, sliceData, VyFileName );
// 				data2D_output_bin( grid, slice, thisMPICoord, &W[2], WSIZE, sliceData, VzFileName );

				
// 			}	
// 			break;	
// 		case 'T':
// 			{

// 				char TxxFileName[128], TyyFileName[128], TzzFileName[128];
// 				char TxyFileName[128], TxzFileName[128], TyzFileName[128];

// 				sprintf( TxxFileName, "%s/Txx_%d", params.OUT, it );
// 				sprintf( TyyFileName, "%s/Tyy_%d", params.OUT, it );
// 				sprintf( TzzFileName, "%s/Tzz_%d", params.OUT, it );
                                             
// 				sprintf( TxyFileName, "%s/Txy_%d", params.OUT, it );
// 				sprintf( TxzFileName, "%s/Txz_%d", params.OUT, it );
// 				sprintf( TyzFileName, "%s/Tyz_%d", params.OUT, it );

// 				data2D_output_bin( grid, slice, thisMPICoord, &W[3], WSIZE, sliceData, TxxFileName );
// 				data2D_output_bin( grid, slice, thisMPICoord, &W[4], WSIZE, sliceData, TyyFileName );
// 				data2D_output_bin( grid, slice, thisMPICoord, &W[5], WSIZE, sliceData, TzzFileName );
// 				data2D_output_bin( grid, slice, thisMPICoord, &W[6], WSIZE, sliceData, TxyFileName );
// 				data2D_output_bin( grid, slice, thisMPICoord, &W[7], WSIZE, sliceData, TxzFileName );
// 				data2D_output_bin( grid, slice, thisMPICoord, &W[8], WSIZE, sliceData, TyzFileName );
				
// 			}	
// 			break;	

// 		case 'F':
// 			{
// 				char VxFileName[128], VyFileName[128], VzFileName[128];
// 				sprintf( VxFileName, "%s/FreeSurfVx_%d", params.OUT, it );
// 				sprintf( VyFileName, "%s/FreeSurfVy_%d", params.OUT, it );
// 				sprintf( VzFileName, "%s/FreeSurfVz_%d", params.OUT, it );

// 				data2D_output_bin( grid, slice, thisMPICoord, &W[0], WSIZE, sliceData, VxFileName );
// 				data2D_output_bin( grid, slice, thisMPICoord, &W[1], WSIZE, sliceData, VyFileName );
// 				data2D_output_bin( grid, slice, thisMPICoord, &W[2], WSIZE, sliceData, VzFileName );
// 			}

// 	}

// }

void data2D_XYZ_out( MPI_COORD thisMPICoord, PARAMS params, GRID grid, float * my_W, float * my_W_Tyz, SLICE slice, SLICE_DATA sliceData, char var, int it, int SAMP, int SAMPLE_SIZE,int sub_rank_x , MPI_Comm sub_comm_x,int sub_rank_y , MPI_Comm sub_comm_y, int Gather )
// void data2D_XYZ_out( MPI_COORD thisMPICoord, PARAMS params, GRID grid, float * W, SLICE slice, SLICE_DATA sliceData, char var, int it )
{

	switch ( var )
	{
		case 'V':
			{

				char VxFileName[128], VyFileName[128], VzFileName[128];

				sprintf( VxFileName, "%s/Vx_%d", params.OUT, it );
				sprintf( VyFileName, "%s/Vy_%d", params.OUT, it );
				sprintf( VzFileName, "%s/Vz_%d", params.OUT, it );

				// data2D_output_bin( grid, slice, thisMPICoord, &W[0], WSIZE, sliceData, VxFileName );
				// data2D_output_bin( grid, slice, thisMPICoord, &W[1], WSIZE, sliceData, VyFileName );
				// data2D_output_bin( grid, slice, thisMPICoord, &W[2], WSIZE, sliceData, VzFileName );
				data2D_output_bin( grid, slice, thisMPICoord, &my_W[0], WSIZE_V, sliceData, VxFileName, SAMP, SAMPLE_SIZE, sub_rank_x, sub_comm_x, sub_rank_y, sub_comm_y, Gather  );
				data2D_output_bin( grid, slice, thisMPICoord, &my_W[1], WSIZE_V, sliceData, VyFileName, SAMP, SAMPLE_SIZE, sub_rank_x, sub_comm_x, sub_rank_y, sub_comm_y, Gather  );
				data2D_output_bin( grid, slice, thisMPICoord, &my_W[2], WSIZE_V, sliceData, VzFileName, SAMP, SAMPLE_SIZE, sub_rank_x, sub_comm_x, sub_rank_y, sub_comm_y, Gather  );

				
			}	
			break;	
		case 'T':
			{

				char TxxFileName[128], TyyFileName[128], TzzFileName[128];
				char TxyFileName[128], TxzFileName[128], TyzFileName[128];

				sprintf( TxxFileName, "%s/Txx_%d", params.OUT, it );
				sprintf( TyyFileName, "%s/Tyy_%d", params.OUT, it );
				sprintf( TzzFileName, "%s/Tzz_%d", params.OUT, it );
                                             
				sprintf( TxyFileName, "%s/Txy_%d", params.OUT, it );
				sprintf( TxzFileName, "%s/Txz_%d", params.OUT, it );
				sprintf( TyzFileName, "%s/Tyz_%d", params.OUT, it );

				// data2D_output_bin( grid, slice, thisMPICoord, &W[3], WSIZE, sliceData, TxxFileName );
				// data2D_output_bin( grid, slice, thisMPICoord, &W[4], WSIZE, sliceData, TyyFileName );
				// data2D_output_bin( grid, slice, thisMPICoord, &W[5], WSIZE, sliceData, TzzFileName );
				// data2D_output_bin( grid, slice, thisMPICoord, &W[6], WSIZE, sliceData, TxyFileName );
				// data2D_output_bin( grid, slice, thisMPICoord, &W[7], WSIZE, sliceData, TxzFileName );
				// data2D_output_bin( grid, slice, thisMPICoord, &W[8], WSIZE, sliceData, TyzFileName );
				data2D_output_bin( grid, slice, thisMPICoord, &my_W[3], WSIZE - 1, sliceData, TxxFileName, SAMP, SAMPLE_SIZE, sub_rank_x, sub_comm_x, sub_rank_y, sub_comm_y, Gather  );
				data2D_output_bin( grid, slice, thisMPICoord, &my_W[4], WSIZE - 1, sliceData, TyyFileName, SAMP, SAMPLE_SIZE, sub_rank_x, sub_comm_x, sub_rank_y, sub_comm_y, Gather );
				data2D_output_bin( grid, slice, thisMPICoord, &my_W[5], WSIZE - 1, sliceData, TzzFileName, SAMP, SAMPLE_SIZE, sub_rank_x, sub_comm_x, sub_rank_y, sub_comm_y, Gather  );
				data2D_output_bin( grid, slice, thisMPICoord, &my_W[6], WSIZE - 1, sliceData, TxyFileName, SAMP, SAMPLE_SIZE, sub_rank_x, sub_comm_x, sub_rank_y, sub_comm_y, Gather  );
				data2D_output_bin( grid, slice, thisMPICoord, &my_W[7], WSIZE - 1, sliceData, TxzFileName, SAMP, SAMPLE_SIZE, sub_rank_x, sub_comm_x, sub_rank_y, sub_comm_y, Gather  );
				data2D_output_bin( grid, slice, thisMPICoord, &my_W_Tyz[0], 1, sliceData, TyzFileName, SAMP, SAMPLE_SIZE, sub_rank_x, sub_comm_x, sub_rank_y, sub_comm_y, Gather  );
				
			}	
			break;	

		case 'F':
			{
				char VxFileName[128], VyFileName[128], VzFileName[128];
				sprintf( VxFileName, "%s/FreeSurfVx_%d", params.OUT, it );
				sprintf( VyFileName, "%s/FreeSurfVy_%d", params.OUT, it );
				sprintf( VzFileName, "%s/FreeSurfVz_%d", params.OUT, it );

				// data2D_output_bin( grid, slice, thisMPICoord, &W[0], WSIZE, sliceData, VxFileName );
				// data2D_output_bin( grid, slice, thisMPICoord, &W[1], WSIZE, sliceData, VyFileName );
				// data2D_output_bin( grid, slice, thisMPICoord, &W[2], WSIZE, sliceData, VzFileName );
				data2D_output_bin( grid, slice, thisMPICoord, &my_W[0], WSIZE - 1, sliceData, VxFileName, SAMP, SAMPLE_SIZE, sub_rank_x, sub_comm_x, sub_rank_y, sub_comm_y, Gather  );
				data2D_output_bin( grid, slice, thisMPICoord, &my_W[1], WSIZE - 1, sliceData, VyFileName, SAMP, SAMPLE_SIZE, sub_rank_x, sub_comm_x, sub_rank_y, sub_comm_y, Gather  );
				data2D_output_bin( grid, slice, thisMPICoord, &my_W[2], WSIZE - 1, sliceData, VzFileName, SAMP, SAMPLE_SIZE, sub_rank_x, sub_comm_x, sub_rank_y, sub_comm_y, Gather  );
			}

	}

}

void data2D_Model_out( MPI_COORD thisMPICoord, PARAMS params, GRID grid, float * C, float * M, SLICE slice, SLICE_DATA sliceData, int SAMP, int SAMPLE_SIZE,int sub_rank_x , MPI_Comm sub_comm_x,int sub_rank_y , MPI_Comm sub_comm_y, int Gather )
{
	char XName[256], YName[256], ZName[256];

	{
		sprintf( XName, "%s/coordX", params.OUT );
		sprintf( YName, "%s/coordY", params.OUT );
		sprintf( ZName, "%s/coordZ", params.OUT );

		data2D_output_bin( grid, slice, thisMPICoord, &C[0], CSIZE, sliceData, XName, SAMP, SAMPLE_SIZE, sub_rank_x, sub_comm_x, sub_rank_y, sub_comm_y, Gather  );
		data2D_output_bin( grid, slice, thisMPICoord, &C[1], CSIZE, sliceData, YName, SAMP, SAMPLE_SIZE, sub_rank_x, sub_comm_x, sub_rank_y, sub_comm_y, Gather  );
		data2D_output_bin( grid, slice, thisMPICoord, &C[2], CSIZE, sliceData, ZName, SAMP, SAMPLE_SIZE, sub_rank_x, sub_comm_x, sub_rank_y, sub_comm_y, Gather  );

		memset( XName, 0, 256 );
		memset( YName, 0, 256 );
		memset( ZName, 0, 256 );
	}
	{
		sprintf( XName, "%s/Vs", params.OUT );
		sprintf( YName, "%s/Vp", params.OUT );
		sprintf( ZName, "%s/rho", params.OUT );

		data2D_output_bin( grid, slice, thisMPICoord, &M[10], MSIZE, sliceData, XName, SAMP, SAMPLE_SIZE, sub_rank_x, sub_comm_x, sub_rank_y, sub_comm_y, Gather  );
		data2D_output_bin( grid, slice, thisMPICoord, &M[11], MSIZE, sliceData, YName, SAMP, SAMPLE_SIZE, sub_rank_x, sub_comm_x, sub_rank_y, sub_comm_y, Gather  );
		data2D_output_bin( grid, slice, thisMPICoord, &M[12], MSIZE, sliceData, ZName, SAMP, SAMPLE_SIZE, sub_rank_x, sub_comm_x, sub_rank_y, sub_comm_y, Gather  );
	}

}
