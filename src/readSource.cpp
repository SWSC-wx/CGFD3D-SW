/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:readSource.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-09-17
*   Discription:
*
================================================================*/

#include "headers.h"


int readSourceInfo( MPI_COORD thisMPICoord, SOURCE_FILE_INPUT * src_in )
{
	char fileName[256];
	sprintf( fileName, "output/source_mpi_%d_%d_%d.bin", thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );
	
	FILE * file = fopen( fileName, "rb" );

	if ( file == NULL )
	{
		return 0;
	}
	
	fread( &( src_in->npts ), sizeof( int ), 1, file );
	fread( &( src_in->nt ), sizeof( int ), 1, file );
	fread( &( src_in->dt ), sizeof( float ), 1, file );
	
	fclose( file );

	//printf("npts = %d\n", src_in->npts );
	//printf("nt = %d\n", src_in->nt );
	//printf("dt = %f\n", src_in->dt );
	return 1;


}

void allocReadMomentRate( SOURCE_FILE_INPUT src_in, MOMENT_RATE momentRate )
{

	int size = src_in.npts * src_in.nt;
	momentRate = (float *) malloc(6*size*sizeof(float));
	// momentRate->Mxx = ( float * ) malloc( sizeof( float ) * size );
	// momentRate->Myy = ( float * ) malloc( sizeof( float ) * size );
	// momentRate->Mzz = ( float * ) malloc( sizeof( float ) * size );
	// momentRate->Mxy = ( float * ) malloc( sizeof( float ) * size );
	// momentRate->Mxz = ( float * ) malloc( sizeof( float ) * size );
	// momentRate->Myz = ( float * ) malloc( sizeof( float ) * size );

}

void freeReadMomentRate( MOMENT_RATE momentRate )
{
    free(momentRate);
	// free( momentRate.Mxx );
	// free( momentRate.Myy );
	// free( momentRate.Mzz );
	// free( momentRate.Mxy );
	// free( momentRate.Mxz );
	// free( momentRate.Myz );

}

/*
void allocMomentRateDiff(SOURCE_FILE_INPUT src_in, MOMENT_RATE momentRateDiff) 
{
	int size = src_in.npts*src_in.nt;
	momentRateDiff = (float *) malloc(6*size*sizeof(float));
}

void freeMomentRateDiff(MOMENT_RATE momentRateDiff) 
{
    free(momentRateDiff);
}
*/

void allocReadIndex( SOURCE_FILE_INPUT src_in, long long ** srcIndex )
{
	int npts = src_in.npts;

	*srcIndex = ( long long * ) malloc( sizeof( long long ) * npts ); 

}

void freeReadIndex( long long * srcIndex )
{
	free( srcIndex );
}

/*
void readIndexMomentRate( MPI_COORD thisMPICoord,  SOURCE_FILE_INPUT src_in, long long * srcIndex, MOMENT_RATE momentRate )
{
	char fileName[256];
	sprintf( fileName, "output/source_mpi_%d_%d_%d.bin", thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );
	
	// float * Mxx = momentRate.Mxx;
	// float * Myy = momentRate.Myy;
	// float * Mzz = momentRate.Mzz;
	// float * Mxy = momentRate.Mxy;
	// float * Mxz = momentRate.Mxz;
	// float * Myz = momentRate.Myz;
	
	FILE * file = fopen( fileName, "rb" );

	fseek( file, sizeof( int ) + sizeof( int ) + sizeof( float ), SEEK_SET );

	int npts = src_in.npts;
	int nt = src_in.nt;

	int p = 0;
	for ( p = 0; p < npts; p ++ )
	{
		fread( srcIndex + p, sizeof( long long ), 1, file );
		fread( Mxx + p, sizeof( float ), nt, file );
		fread( Myy + p, sizeof( float ), nt, file );
		fread( Mzz + p, sizeof( float ), nt, file );
		fread( Mxy + p, sizeof( float ), nt, file );
		fread( Mxz + p, sizeof( float ), nt, file );
		fread( Myz + p, sizeof( float ), nt, file );
	}

	fclose( file );
}
*/



void verifySrcIndex( MPI_COORD thisMPICoord, long long * srcIndex, int npts  )
{

	char fileName[256];
	sprintf( fileName, "output/srcIndex_%d_%d_%d.bin",  thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );
	FILE * file;

	file = fopen( fileName, "wb" ); 
	if(file == NULL){
		mkdir("output", 0777);
		file = fopen( fileName, "wb" );
	}
		
	fwrite( srcIndex, sizeof( long long ), npts, file );


	fclose( file );

}


MOMENT_RATE allocMomentRateSlice( SOURCE_FILE_INPUT src_in)
{
	
	int size = src_in.npts;
	MOMENT_RATE momentRateSlice = (float *) malloc(6*size*sizeof(float));
	return momentRateSlice;
	// momentRateSlice->Mxx = ( float * ) malloc( sizeof( float ) * size );
	// momentRateSlice->Myy = ( float * ) malloc( sizeof( float ) * size );
	// momentRateSlice->Mzz = ( float * ) malloc( sizeof( float ) * size );
	// momentRateSlice->Mxy = ( float * ) malloc( sizeof( float ) * size );
	// momentRateSlice->Mxz = ( float * ) malloc( sizeof( float ) * size );
	// momentRateSlice->Myz = ( float * ) malloc( sizeof( float ) * size );

}

void freeMomentRateSlice( MOMENT_RATE momentRateSlice  )
{
	free(momentRateSlice);
	// free( momentRateSlice.Mxx );
	// free( momentRateSlice.Myy );
	// free( momentRateSlice.Mzz );
	// free( momentRateSlice.Mxy );
	// free( momentRateSlice.Mxz );
	// free( momentRateSlice.Myz );
	
}

void interpMomentRate( SOURCE_FILE_INPUT src_in, long long  * srcSlicePos, MOMENT_RATE momentRate,  MOMENT_RATE momentRateSlice, float t_weight, int srcIt )
{
	int npts = src_in.npts;
	int nt = src_in.nt;
	
	int i = 0;

	long long ii = 0;
	long long pos = 0;

	for ( i = 0; i < npts; i ++ )
	{
        ii = srcSlicePos[i];
	    pos = ii * nt + srcIt;
	    // printf("debug for interpMomentRate npts:%d\n", npts);	
	    long long pos_n0 = pos*6 + 0;
	    long long pos_n1 = pos*6 + 1;
	    long long pos_n2 = pos*6 + 2;
	    long long pos_n3 = pos*6 + 3;
	    long long pos_n4 = pos*6 + 4;
	    long long pos_n5 = pos*6 + 5;
	    long long posP1_n0 = (pos + 1)*6 + 0;
	    long long posP1_n1 = (pos + 1)*6 + 1;
	    long long posP1_n2 = (pos + 1)*6 + 2;
	    long long posP1_n3 = (pos + 1)*6 + 3;
	    long long posP1_n4 = (pos + 1)*6 + 4;
	    long long posP1_n5 = (pos + 1)*6 + 5;
		long long pos_s0 = i*6 + 0;
		long long pos_s1 = i*6 + 1;
		long long pos_s2 = i*6 + 2;
		long long pos_s3 = i*6 + 3;
		long long pos_s4 = i*6 + 4;
		long long pos_s5 = i*6 + 5;
		momentRateSlice[pos_s0] = (momentRate[posP1_n0] - momentRate[pos_n0])*t_weight + momentRate[pos_n0];
		momentRateSlice[pos_s1] = (momentRate[posP1_n1] - momentRate[pos_n1])*t_weight + momentRate[pos_n1];
		momentRateSlice[pos_s2] = (momentRate[posP1_n2] - momentRate[pos_n2])*t_weight + momentRate[pos_n2];
		momentRateSlice[pos_s3] = (momentRate[posP1_n3] - momentRate[pos_n3])*t_weight + momentRate[pos_n3];
		momentRateSlice[pos_s4] = (momentRate[posP1_n4] - momentRate[pos_n4])*t_weight + momentRate[pos_n4];
		momentRateSlice[pos_s5] = (momentRate[posP1_n5] - momentRate[pos_n5])*t_weight + momentRate[pos_n5];

		// momentRateSlice.Mxx[i] = ( momentRate.Mxx[pos + 1] - momentRate.Mxx[pos] ) * t_weight + momentRate.Mxx[pos];
		// momentRateSlice.Myy[i] = ( momentRate.Myy[pos + 1] - momentRate.Myy[pos] ) * t_weight + momentRate.Myy[pos];
		// momentRateSlice.Mzz[i] = ( momentRate.Mzz[pos + 1] - momentRate.Mzz[pos] ) * t_weight + momentRate.Mzz[pos];
		// momentRateSlice.Mxy[i] = ( momentRate.Mxy[pos + 1] - momentRate.Mxy[pos] ) * t_weight + momentRate.Mxy[pos];
		// momentRateSlice.Mxz[i] = ( momentRate.Mxz[pos + 1] - momentRate.Mxz[pos] ) * t_weight + momentRate.Mxz[pos];
		// momentRateSlice.Myz[i] = ( momentRate.Myz[pos + 1] - momentRate.Myz[pos] ) * t_weight + momentRate.Myz[pos];
		
	}
}

double rickerFun( float t, float t0, float a )
{
    return exp( -( t - t0  ) * ( t - t0 ) / ( a * a  )  ) / ( sqrt( 3.141592657 ) * a  );   
}
	


// void addMomentRate( float * W, GRID grid,  MOMENT_RATE momentRateSlice, long long * srcIndex, long long * srcSlicePos, int npts, float DH, float DT, float * M, int it )
// {
// 	long long i = 0;
//     int x = 0, y = 0, z = 0;
//     long long idx = 0, pos;
	
//     float V = 0.0;
//     float Jac, Mu;
    
//     int	_nx_ = grid._nx_;
//     int	_ny_ = grid._ny_;
//     int	_nz_ = grid._nz_;
    
//     int	frontNX = grid.frontNX;
//     int	frontNY = grid.frontNY;
//     int	frontNZ = grid.frontNZ;


//     float t = 0.0, t0 = 3.0, a = 1.5;
//     t = it * DT;
    
// 	for ( i = 0; i < npts; i ++ )
// 	{
// 		idx = srcIndex[i];
// 		pos = srcSlicePos[i];
// 		//z = idx / ( _nz_ * _ny_ ) + frontNX;
// 		//y = idx % ( _nz_ * _ny_ ) / _nz_ + frontNY;
// 		//x = idx % ( _nz_ * _ny_ ) % _nz_ + frontNZ;

//         //idx = INDEX( x, y, z );

// 		Jac = M[idx * MSIZE + 9];
// 		Mu  = M[idx * MSIZE + 11];
// 		V = - DT * Mu / ( Jac * DH * DH * DH );
//         //if ( i == 200 )
//         //    printf( "%f\n", momentRateSlice.Mxx[pos]  );
// 		long long pos_s0 = i*6 + 0;
// 		long long pos_s1 = i*6 + 1;
// 		long long pos_s2 = i*6 + 2;
// 		long long pos_s3 = i*6 + 3;
// 		long long pos_s4 = i*6 + 4;
// 		long long pos_s5 = i*6 + 5;

// 		W[idx * WSIZE + 3] += momentRateSlice[pos_s0] * V;
// 		W[idx * WSIZE + 4] += momentRateSlice[pos_s1] * V;
// 		W[idx * WSIZE + 5] += momentRateSlice[pos_s2] * V;
// 		W[idx * WSIZE + 6] += momentRateSlice[pos_s3] * V;
// 		W[idx * WSIZE + 7] += momentRateSlice[pos_s4] * V;
// 		W[idx * WSIZE + 8] += momentRateSlice[pos_s5] * V;

// 		// W[idx * WSIZE + 3] += momentRateSlice.Mxx[i] * V;
// 		// W[idx * WSIZE + 4] += momentRateSlice.Myy[i] * V;
// 		// W[idx * WSIZE + 5] += momentRateSlice.Mzz[i] * V;
// 		// W[idx * WSIZE + 6] += momentRateSlice.Mxy[i] * V;
// 		// W[idx * WSIZE + 7] += momentRateSlice.Mxz[i] * V;
// 		// W[idx * WSIZE + 8] += momentRateSlice.Myz[i] * V;
// 	}
// }

void addMomentRate( float * my_W, float * my_W_Tyz, GRID grid,  MOMENT_RATE momentRateSlice, long long * srcIndex, long long * srcSlicePos, int npts, float DH, float DT, float * M, int it )
{
	long long i = 0;
    int x = 0, y = 0, z = 0;
    long long idx = 0, pos;
	
    float V = 0.0;
    float Jac, Mu;
    
    int	_nx_ = grid._nx_;
    int	_ny_ = grid._ny_;
    int	_nz_ = grid._nz_;
    
    int	frontNX = grid.frontNX;
    int	frontNY = grid.frontNY;
    int	frontNZ = grid.frontNZ;


    float t = 0.0, t0 = 3.0, a = 1.5;
    t = it * DT;
    
	for ( i = 0; i < npts; i ++ )
	{
		idx = srcIndex[i];
		pos = srcSlicePos[i];
		//z = idx / ( _nz_ * _ny_ ) + frontNX;
		//y = idx % ( _nz_ * _ny_ ) / _nz_ + frontNY;
		//x = idx % ( _nz_ * _ny_ ) % _nz_ + frontNZ;

        //idx = INDEX( x, y, z );

		Jac = M[idx * MSIZE + 9];
		Mu  = M[idx * MSIZE + 11];
		V = - DT * Mu / ( Jac * DH * DH * DH );
        //if ( i == 200 )
        //    printf( "%f\n", momentRateSlice.Mxx[pos]  );
		long long pos_s0 = i*6 + 0;
		long long pos_s1 = i*6 + 1;
		long long pos_s2 = i*6 + 2;
		long long pos_s3 = i*6 + 3;
		long long pos_s4 = i*6 + 4;
		long long pos_s5 = i*6 + 5;

		my_W[idx * WSIZE_V + 3] += momentRateSlice[pos_s0] * V;
		my_W[idx * WSIZE_V + 4] += momentRateSlice[pos_s1] * V;
		my_W[idx * WSIZE_V + 5] += momentRateSlice[pos_s2] * V;
		my_W[idx * WSIZE_V + 6] += momentRateSlice[pos_s3] * V;
		my_W[idx * WSIZE_V + 7] += momentRateSlice[pos_s4] * V;
		my_W_Tyz[idx] += momentRateSlice[pos_s5] * V;

		// W[idx * WSIZE + 3] += momentRateSlice.Mxx[i] * V;
		// W[idx * WSIZE + 4] += momentRateSlice.Myy[i] * V;
		// W[idx * WSIZE + 5] += momentRateSlice.Mzz[i] * V;
		// W[idx * WSIZE + 6] += momentRateSlice.Mxy[i] * V;
		// W[idx * WSIZE + 7] += momentRateSlice.Mxz[i] * V;
		// W[idx * WSIZE + 8] += momentRateSlice.Myz[i] * V;
	}
}







