/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:source.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-09-14
*   Discription:
*
================================================================*/
#include "headers.h"

void allocSourceLocation( SOURCE_FILE_INPUT * src_in  )
{
	long long npts = src_in->npts;
	int nt = src_in->nt;

	src_in->lon = ( float * )malloc( npts * sizeof( float ) );
	src_in->lat = ( float * )malloc( npts * sizeof( float ) );
	src_in->coordZ = ( float * )malloc( npts * sizeof( float ) );
	
	memset( src_in->lon,   0, npts * sizeof( float ) );
	memset( src_in->lat,   0, npts * sizeof( float ) );
	memset( src_in->coordZ,0, npts * sizeof( float ) );

}


void freeSourceLocation( SOURCE_FILE_INPUT src_in  )
{
	free( src_in.lon );
	free( src_in.lat );
	free( src_in.coordZ );
}


void readSourceInfo  ( SOURCE_FILE_INPUT * src_in,  char * sourceFileName )
{
	FILE * sourceFile = fopen( sourceFileName, "rb" );
	if(sourceFile == NULL){
		printf( "There is not %s file!\n", sourceFileName );
		exit( 1 );
	}
	
#ifndef SOURCE_NPTS_LONG_LONG
	int nnnnnn;
	fread( &( nnnnnn ), sizeof( int ), 1, sourceFile ); 
	src_in->npts = nnnnnn;
#else
	fread( &( src_in->npts ), sizeof( long long  ), 1, sourceFile ); 
#endif
	fread( &( src_in->nt ),   sizeof( int ), 1, sourceFile ); 
	fread( &( src_in->dt ),   sizeof( float ), 1, sourceFile ); 
	
	fclose( sourceFile );

}

void readSourceLocation( SOURCE_FILE_INPUT src_in,  char * sourceFileName )
{
	
	FILE * sourceFile = fopen( sourceFileName, "rb" );

	if ( NULL == sourceFile ){
        cout << sourceFileName << ", is source file" << endl;
		exit( 1 );
	}
		
	fseek( sourceFile, sizeof( int ) + sizeof( int ) + sizeof( float ), SEEK_CUR );

	long long npts = src_in.npts;
	int nt = src_in.nt;
	
	long long i = 0;
	for ( i = 0; i < npts; i ++ )
	{
		fread( &( src_in.lon[i]    ), 1, sizeof( float ), sourceFile ); 
		fread( &( src_in.lat[i]    ), 1, sizeof( float ), sourceFile ); 
		fread( &( src_in.coordZ[i] ), 1, sizeof( float ), sourceFile ); 
		
		fseek( sourceFile, 3 * sizeof( float ) + 2 * nt * sizeof( float ), SEEK_CUR );
	}

	//int nnnn = 1000;

	//cout << "lon[" << nnnn  <<"] = " << src_in.lon[nnnn] << endl;
	//cout << "npts = " << npts << ", nt = " << nt << ", dt = " << src_in.dt << endl;
	fclose( sourceFile );
}


void LonLat2cart( PJ * P,  PJ_DIRECTION PD, SOURCE_FILE_INPUT src_in )
{
	long long npts = src_in.npts;

	PJ_COORD * pj_coord;


	pj_coord = ( PJ_COORD * )malloc( sizeof( PJ_COORD ) * src_in.npts );

	long long i = 0;
	for ( i = 0; i < npts; i ++ )
	{
		pj_coord[i].lp.lam = src_in.lon[i] * DEGREE2RADIAN;
		pj_coord[i].lp.phi = src_in.lat[i] * DEGREE2RADIAN;
	}
	

	proj_trans_array( P, PD, npts, pj_coord );


	for ( i = 0; i < npts; i ++ )
	{
		src_in.lon[i] = pj_coord[i].xy.x;
		src_in.lat[i] = pj_coord[i].xy.y;
	}
	

	free( pj_coord );
}

void projTrans( PARAMS params, SOURCE_FILE_INPUT src_in )
{
	float lon_0 = params.centerLongitude;
	float lat_0 = params.centerLatitude;

	long long i = 0;
	float maxLon = - 180.0, minLon = 180.0, maxLat = -90.0, minLat = 90.0;
	for ( i = 0; i < src_in.npts; i ++ )
	{	
		if ( maxLon < src_in.lon[i] )
		{
			maxLon = src_in.lon[i];
		}
		if ( minLon > src_in.lon[i] )
		{
			minLon = src_in.lon[i];
		}

		if ( maxLat < src_in.lat[i] )
		{
			maxLat = src_in.lat[i];
		}
		if ( minLat > src_in.lat[i] )
		{
			minLat = src_in.lat[i];
		}
	}

	
	//cout << "maxLon = " << maxLon << ", minLon = " << minLon << ", maxLat = " << maxLat << ", minLat = " << minLat << endl;
	
	PJ_CONTEXT * C;
	PJ * P;
	
	C = proj_context_create( );
	
	char projStr[256];//""
	sprintf( projStr, "+proj=aeqd +lon_0=%lf +lat_0=%lf +x_0=0.0 +y_0=0.0 +ellps=WGS84", lon_0, lat_0 );

	//printf( projStr  );
	//printf( "\n"  );
	P = proj_create( C, projStr );
	if ( NULL == P )
	{
		printf( "Failed to create projection\n"  );
	}


	LonLat2cart( P, PJ_FWD, src_in );
	proj_destroy( P );
	proj_context_destroy( C );
}

void allocateSourceIndex( SOURCE_INDEX * srcIndex, long long npts )
{
	srcIndex->X = ( int * )malloc( npts * sizeof( int ) );
	srcIndex->Y = ( int * )malloc( npts * sizeof( int ) );
	srcIndex->Z = ( int * )malloc( npts * sizeof( int ) );
}


void freeSourceIndex( SOURCE_INDEX srcIndex )
{
	free( srcIndex.X );
	free( srcIndex.Y );
	free( srcIndex.Z );
}

long long srcCoord2Index( GRID grid, float Depth, float * coord, SOURCE_FILE_INPUT src_in, map< long long, long long > & point_index )//map<int, POINT_INDEX> pos_pointIndex )
{
	long long i = 0;
	long long npts = src_in.npts;
	int srcX, srcY, srcZ;

	float DH = grid.DH;
	double DZ = 0.0;

	float * x = src_in.lon;
	float * y = src_in.lat;
	float * z = src_in.coordZ;

	int originalX = grid.originalX;
	int originalY = grid.originalY;

	int frontNX = grid.frontNX;
	int frontNY = grid.frontNY;
	int frontNZ = grid.frontNZ;

	int halo = grid.halo;
	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;

	int _nx = grid._nx;
	int _ny = grid._ny;
	int _nz = grid._nz;

	long long  pointNum = 0;

	//cout << "x = " << x[npts-1] << ", y = " << y[npts-1] << ", z = " << z[npts-1] << endl;
	float maxCoordZ = -Depth;
	float minCoordZ = Depth;
	for ( int jj = 0; jj < npts; jj ++  )
	{
		if( maxCoordZ < z[jj] )
			maxCoordZ = z[jj];

		if( minCoordZ > z[jj] )
			minCoordZ = z[jj];
	}

	for ( int jj = 0; jj < npts; jj ++  )
	{
		z[jj] = z[jj] - maxCoordZ;
	}
	//cout << "maxCoordZ = " << maxCoordZ << ", minCoordZ = " << minCoordZ  << endl;
	
	int extendNZ = 0;
	
	//cout << npts << endl;

	int k = 0;

	long long index1 = 0, index2 = 0, index = 0;
	double dz1 = 0.0, dz2 = 0.0;
	for ( i = 0; i < npts; i ++ )
	{
		srcX = int( x[i] / DH + 0.5 ) + originalX - frontNX + halo; 
		srcY = int( y[i] / DH + 0.5 ) + originalY - frontNY + halo; 
		//if ( srcX >= 0 && srcX < _nx_ && srcY >= 0 && srcY < _ny_ )

#define NO_SOURCE_SMOOTH

#ifdef NO_SOURCE_SMOOTH
		if ( srcX >= halo && srcX < _nx && srcY >= halo && srcY < _ny )
		{
			//srcIndex.X[i] = srcX;
			//srcIndex.Y[i] = srcY;
			for ( k = halo; k < _nz; k ++ )
			{
				index1 = INDEX( srcX, srcY, k - 1 ) * CSIZE + 2;
				index2 = INDEX( srcX, srcY, k + 1 ) * CSIZE + 2;
				index  = INDEX( srcX, srcY, k )     * CSIZE + 2;
			
				dz1 = ( coord[index ] - coord[index1] ) * 0.5;
				dz2 = ( coord[index2] - coord[index ] ) * 0.5;
				
				if ( coord[index] - dz1 <= z[i] && z[i] <  dz2 + coord[index] )
				{	
					srcZ = k;
					point_index[i] = INDEX( srcX, srcY, srcZ );
					pointNum ++;					

				}
					

			}

		}
#else
		
		if ( srcX >= 0 && srcX < _nx_ && srcY >= 0 && srcY < _ny_ )
		{
			for ( k = 0; k < _nz_; k ++ )
			{
				if ( k - 1 == -1 )
				{	
					index  = INDEX( srcX, srcY, 0 ) * CSIZE + 2;
					index2 = INDEX( srcX, srcY, 1 ) * CSIZE + 2;
					dz2 = ( coord[index2] - coord[index] ) * 0.5;
					dz1 = dz2;
				}
				if ( k + 1 == _nz_ )
				{
					index1 = INDEX( srcX, srcY, _nz_ - 2 ) * CSIZE + 2;
					index  = INDEX( srcX, srcY, _nz_ - 1 ) * CSIZE + 2;
					dz1 = ( coord[index] - coord[index1] ) * 0.5;
					dz2 = dz1;
					
				}
				if ( k - 1 != -1 && k + 1 != _nz_ )
				{
					index1 = INDEX( srcX, srcY, k - 1 ) * CSIZE + 2;
					index  = INDEX( srcX, srcY, k )     * CSIZE + 2;
					index2 = INDEX( srcX, srcY, k + 1 ) * CSIZE + 2;
					dz1 = ( coord[index] - coord[index1] ) * 0.5;
					dz2 = ( coord[index2] - coord[index] ) * 0.5;
					
				}
				
				
				if ( coord[index] - dz1 <= z[i] && z[i] <  dz2 + coord[index] )
				{	
					srcZ = k;
					point_index[i] = INDEX( srcX, srcY, srcZ );
					pointNum ++;					

				}
			}
		}
#endif // NO_SOURCE_SMOOTH



	}

	return pointNum;
}






void allocSourceParams( SOURCE_FILE_INPUT * src_in, long long pointNum  )
{
	if ( 0 == pointNum )
		return;
	int nt = src_in->nt;
	
	src_in->area   = ( float * )malloc( pointNum * sizeof( float ) );
	src_in->strike = ( float * )malloc( pointNum * sizeof( float ) );
	src_in->dip    = ( float * )malloc( pointNum * sizeof( float ) );

	src_in->rake = ( float * )malloc( pointNum * nt * sizeof( float ) );
	src_in->rate = ( float * )malloc( pointNum * nt * sizeof( float ) );


}
void freeSourceParams( SOURCE_FILE_INPUT src_in, long long pointNum  )
{
	if ( 0 == pointNum )
		return;


	free( src_in.area   );
	free( src_in.strike );
	free( src_in.dip    );

	free( src_in.rake );
	free( src_in.rate );
}


void readSourceParams( SOURCE_FILE_INPUT src_in,  char * sourceFileName, map<long long, long long> & point_index  )
{

	int size = point_index.size( );
	if ( 0 == size )
	{
		return;

	}
	FILE * sourceFile = fopen( sourceFileName, "rb" );

	if ( NULL == sourceFile ){
        printf( "%s is source file!\n", sourceFileName ); 
		exit( 1 );
		//cout << sourceFileName << ", is source file" << endl;
	}
	
	int npts = src_in.npts;
	int nt = src_in.nt;
	long long i = 0;
	long long nptsIndex; 




	int headSize = sizeof( int ) + sizeof( int ) + sizeof( int ); // npts nt dt
	long long byteSize = 6 * sizeof( float ) + //lon lat coordZ area strike dip
						 2 * sizeof( float ) * nt;

	for ( map< long long, long long>::iterator it = point_index.begin( ); it != point_index.end( ); it ++ )
	{
		nptsIndex = it->first;
		//cout << nptsIndex << endl;
		fseek( sourceFile, headSize + byteSize * nptsIndex + 3 * sizeof( float ), SEEK_SET )  ;  
		//	
		fread( src_in.area + i,   sizeof( float ), 1, sourceFile );
		fread( src_in.strike + i, sizeof( float ), 1, sourceFile );
		fread( src_in.dip + i,    sizeof( float ), 1, sourceFile );

		fread( src_in.rake + i * nt, sizeof( float ), nt, sourceFile );
		fread( src_in.rate + i * nt, sizeof( float ), nt, sourceFile );
		i ++;
	}

	fclose( sourceFile );
}


MOMENT_RATE allocMomentRate(long long pointNum, int nt )
{
	if ( 0 == pointNum )
		return NULL;

	int size = pointNum * nt;
	MOMENT_RATE momentRate = (float *) malloc(6*size*sizeof(float));
	return momentRate;
}

void freeMomentRate( MOMENT_RATE momentRate, long long pointNum )
{
	if ( 0 == pointNum )
		return;
	free(momentRate);
}

void solveMomentRate( SOURCE_FILE_INPUT src_in, MOMENT_RATE momentRate, long long pointNum )
{
	if ( 0 == pointNum )
		return;
	
	float s, d, r, a;
	float * strike = src_in.strike, * dip = src_in.dip, * area = src_in.area;
	float * rake = src_in.rake, * rate = src_in.rate;
	long long p = 0, index = 0, index_m = 0;
	int it = 0, nt = src_in.nt;

	for ( p = 0; p < pointNum; p ++ )
	{
		s = strike[p]*DEGREE2RADIAN;
		d = dip[p]*DEGREE2RADIAN;
		a = area[p];
		for ( it = 0; it <  nt; it ++ )
		{
			index = it + p*nt;
			index_m = 6*index;
            r = rake[index] * DEGREE2RADIAN; 
			momentRate[index_m + 0] = -(sin(d)*cos(r)*sin(s*2.0) + sin(d*2.0)*sin(r) * sin(s) * sin(s)) * a * rate[index];   
			momentRate[index_m + 1] = (sin(d) * cos(r) * sin(s * 2.0) - sin(d * 2.0) * sin(r) * cos(s) * cos(s)) * a * rate[index];    
			momentRate[index_m + 2] = -(momentRate[index_m + 0] + momentRate[index_m + 1]);                                                                                 
			momentRate[index_m + 3] = (sin(d) * cos(r) * cos(s * 2.0) + 0.5 * sin(d * 2.0) * sin(r) * sin(s * 2.0)) * a * rate[index]; 
			momentRate[index_m + 4] = (cos(d) * cos(r) * cos(s) + cos(d * 2.0) * sin(r) * sin(s)) * a * rate[index];                   
			momentRate[index_m + 5] = (cos(d) * cos(r) * sin(s) - cos(d * 2.0) * sin(r) * cos(s)) * a * rate[index];                   
		}
	}
}

void dealDuplicateIndex( MOMENT_RATE momentRate, 
						 map< long long, long long > & point_index, 
						 map< long long, long long > & index_point, 
						 SOURCE_FILE_INPUT src_in )
{
	int size = point_index.size( );
	if ( 0 == size )
	{
		return;

	}
	long long pnt = 0, idx = 0;

	//map < long long, long long > index_point;

	int i = 0;
	for ( map < long long, long long >::iterator it = point_index.begin( ); it != point_index.end( ); it ++ )
	{	
		idx = it->second;
		index_point[idx] = i;//It means that the hashmap index_point stores the max pnt since pnt is in a high order.
		i ++;	
	}

	int t = 0;
	int nt = src_in.nt;
	int ii, indexII, indexI;
	i = 0;

	// float * Mxx = momentRate.Mxx;
	// float * Myy = momentRate.Myy;
	// float * Mzz = momentRate.Mzz;
	// float * Mxy = momentRate.Mxy;
	// float * Mxz = momentRate.Mxz;
	// float * Myz = momentRate.Myz;


	for ( map < long long, long long >::iterator it = point_index.begin( ); it != point_index.end( ); it ++ )
	{	
		idx = it->second;
		ii = index_point[idx];
		if ( ii > i )
		{
			for ( t = 0; t < nt; t ++ )
			{
				indexI =  6*(t + i*nt);
				indexII = 6*(t + ii*nt);

				momentRate[indexII + 0] += momentRate[indexI + 0]; 
				momentRate[indexII + 1] += momentRate[indexI + 1]; 
				momentRate[indexII + 2] += momentRate[indexI + 2]; 
				momentRate[indexII + 3] += momentRate[indexI + 3]; 
				momentRate[indexII + 4] += momentRate[indexI + 4]; 
				momentRate[indexII + 5] += momentRate[indexI + 5]; 

				// Mxx[indexII] += Mxx[indexI]; 
				// Myy[indexII] += Myy[indexI]; 
				// Mzz[indexII] += Mzz[indexI]; 
				// Mxy[indexII] += Mxy[indexI]; 
				// Mxz[indexII] += Mxz[indexI]; 
				// Myz[indexII] += Myz[indexI]; 
			}
		}
		i ++;
	}
}
/*
void outputSourceData( PARAMS params,  SOURCE_FILE_INPUT src_in,   map< long long, long long> & index_point, MOMENT_RATE momentRate, MPI_COORD thisMPICoord )
{

	
	int size = index_point.size( );
	if ( 0 == size )
	{
		return;

	}
	char fileName[256];
	sprintf( fileName, "%s/source_mpi_%d_%d_%d.bin", params.OUT,  thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );
	
	FILE * file = fopen( fileName, "wb" );
	
	if ( NULL == file )
		printf( "The file %s can not be opened\n", fileName );
	
	int npts = src_in.npts;
	int nt = src_in.nt;
	float dt = src_in.dt;

	fwrite( &size, sizeof( int ), 1, file );
	fwrite( &nt, sizeof( int ), 1, file );
	fwrite( &dt, sizeof( float ), 1, file );

	long long index = 0;
	long long ii = 0;
	long long pos = 0;

	// float * Mxx = momentRate.Mxx;
	// float * Myy = momentRate.Myy;
	// float * Mzz = momentRate.Mzz;
	// float * Mxy = momentRate.Mxy;
	// float * Mxz = momentRate.Mxz;
	// float * Myz = momentRate.Myz;

	for( map<long long, long long>::iterator it = index_point.begin( ); it != index_point.end( ); it ++ )
	{
		index = it->first;
		ii = it->second;
		pos = ii * nt;
		fwrite( &index, sizeof( long long ), 1, file);
		fwrite( Mxx + pos, sizeof( float ), nt, file);
		fwrite( Myy + pos, sizeof( float ), nt, file);
		fwrite( Mzz + pos, sizeof( float ), nt, file);
		fwrite( Mxy + pos, sizeof( float ), nt, file);
		fwrite( Mxz + pos, sizeof( float ), nt, file);
		fwrite( Myz + pos, sizeof( float ), nt, file);
	}
	fclose( file );
}
*/

void verifyLocation( SOURCE_FILE_INPUT src_in )
{
	FILE * file[3];
	file[0] = fopen( "output/source_coord_x.bin", "wb" );
	file[1] = fopen( "output/source_coord_y.bin", "wb" );
	file[2] = fopen( "output/source_coord_z.bin", "wb" );
	if(file[0] == NULL){
		mkdir("output", 0777);
		file[0] = fopen( "output/source_coord_x.bin", "wb" );
	    file[1] = fopen( "output/source_coord_y.bin", "wb" );
	    file[2] = fopen( "output/source_coord_z.bin", "wb" );
	}

	long long npts = src_in.npts;

	float * x = src_in.lon;
	float * y = src_in.lat;
	float * z = src_in.coordZ;



	fwrite( x, sizeof( float ), npts, file[0] );
	fwrite( y, sizeof( float ), npts, file[1] );
	fwrite( z, sizeof( float ), npts, file[2] );

	fclose( file[0] );
	fclose( file[1] );
	fclose( file[2] );

}

void dealSource( PARAMS params, GRID grid, float * coord, MPI_COORD thisMPICoord, long long ** srcIndex, long long ** srcSlicePos, MOMENT_RATE *momentRate, SOURCE_FILE_INPUT * src_in_ )
{
	if(!params.useMultiSource){
		return;
	}
	SOURCE_FILE_INPUT src_in;
	char sourceFileName[256];
	sprintf( sourceFileName, "%s", params.sourceFile );
	//cout << sourceFileName << ", is source file" << endl;
	readSourceInfo( &src_in, sourceFileName );
	allocSourceLocation(&src_in);
	readSourceLocation(src_in, sourceFileName); 	
//	cout << "npts = " << src_in.npts << ", nt = " << src_in.nt << ", dt = " << src_in.dt << endl;
	projTrans( params, src_in );	
	//allocateSourceIndex( srcIndex, src_in.npts );
	float Depth = params.Depth * 1000;
	map< long long, long long >  point_index;
	map<long long, long long >  index_point;
	long long pointNum = srcCoord2Index( grid, Depth, coord, src_in, point_index );
	/*
	if( 0 == thisMPICoord.X && 0 == thisMPICoord.Y && 0 == thisMPICoord.Z )
	{
		verifyLocation( src_in );
	}
	*/
	freeSourceLocation( src_in );
	MPI_Barrier( MPI_COMM_WORLD );

	allocSourceParams( &src_in, pointNum );
	readSourceParams( src_in, sourceFileName, point_index );
	*momentRate = allocMomentRate(pointNum, src_in.nt );
	solveMomentRate( src_in, *momentRate, pointNum );
	freeSourceParams( src_in, pointNum );
	dealDuplicateIndex(*momentRate, point_index,	index_point, src_in);
	point_index.clear( );
	// outputSourceData( params, src_in, index_point, *momentRate, thisMPICoord );
	//freeMomentRate( momentRate, pointNum );
    int npts = index_point.size( );
    *srcIndex = ( long long * )malloc( npts * sizeof( long long ) );
    *srcSlicePos = ( long long * )malloc( npts * sizeof( long long ) );
    long long ii = 0;
	for( map<long long, long long>::iterator it = index_point.begin( ); it != index_point.end( ); it ++ )
	{
		//index = it->first;
		//ii = it->second;
        //(*srcIndex)[ii] = index;
        (*srcIndex)[ii] = it->first;
        (*srcSlicePos)[ii] = it->second;
        ii ++;
    }
    src_in_->nt = src_in.nt;
    src_in_->dt = src_in.dt;
    src_in_->npts = npts;
	index_point.clear( );
}



