#include"station.h"

#define CHECK( call ) call

int Malloc( void ** mem, long long  size  )
{
	*mem = malloc( size );
	if ( *mem == NULL )
	{
		printf( "can not malloc, Error: %s:%d\n", __FILE__, __LINE__ );
	}
	return 0;

}

void allocStation( STATION * station, int stationNum, int NT )
{
	int sizeIdx = sizeof(int)*stationNum*3;

	int *pIndex = NULL;

	CHECK( Malloc( ( void ** )&pIndex, sizeIdx ) );
	CHECK( memset( pIndex, 0, sizeIdx ) );

	station->X = pIndex;
	station->Y = pIndex + stationNum;
	station->Z = pIndex + stationNum * 2;

	
	long long sizeWave = sizeof(float)*NT*stationNum*WSIZE;
	float *pWave = NULL;
	
	CHECK( Malloc( ( void ** )&pWave, sizeWave ) );
	CHECK( memset( pWave, 0, sizeWave ) );

	station->wave = pWave;

}

void freeStation( STATION station )
{
	free( station.X );
	free( station.wave );
}

int readStationIndex( GRID grid ) 
{
	char jsonFile[1024] = { 0 };
	strcpy( jsonFile, "station.json" );
	FILE * fp;
	fp = fopen( jsonFile, "r" );

	if ( NULL == fp )
	{
		printf( "There is not %s file!\n", jsonFile );
		return 0;
	}
	
	fseek( fp, 0, SEEK_END );
	int len = ftell( fp );
	
	fseek( fp, 0, SEEK_SET );

	char * jsonStr = ( char * ) malloc( len * sizeof( char ) );

	if ( NULL == jsonStr )
	{
		printf( "Can't allocate json string memory\n" );
		return 0;
	}

	fread( jsonStr, sizeof( char ), len, fp );
	
	//printf( "%s\n", jsonStr );
	cJSON * object;
	cJSON * objArray;

	object = cJSON_Parse( jsonStr );
	if ( NULL == object )
	{
		printf( "Can't parse json file!\n");
		return 0;
	}

	fclose( fp );

	int stationCnt= 0;
	if (objArray = cJSON_GetObjectItem(object, "station(point)"))
	{
		stationCnt = cJSON_GetArraySize( objArray );
	}
	
	cJSON *stationObj, *stationItem;

	int i, j, stationNum;
	int X, Y, Z, thisX, thisY, thisZ;
	int frontNX = grid.frontNX;
	int frontNY = grid.frontNY;
	int frontNZ = grid.frontNZ;

	int _nx = grid._nx;
	int _ny = grid._ny;
	int _nz = grid._nz;

	stationNum = 0;
	for ( i = 0; i < stationCnt; i ++  )
	{
		stationObj = cJSON_GetArrayItem( objArray, i );
		int a = cJSON_GetArraySize( stationObj );
		if ( a != 3 )
		{
			printf( "In file %s, the coodinate index don't equal to 3. However, it equals to %d\n", jsonFile, a );
			return 0;
		}
	
		stationItem = cJSON_GetArrayItem( stationObj, 0 );
		X = stationItem->valueint;
		thisX = X - frontNX + HALO;

		stationItem = cJSON_GetArrayItem( stationObj, 1 );
		Y = stationItem->valueint;
		thisY = Y - frontNY + HALO;

		stationItem = cJSON_GetArrayItem( stationObj, 2 );
		Z = stationItem->valueint;
		thisZ = Z - frontNZ + HALO;
			
		if ( thisX >= HALO && thisX < _nx &&
			 thisY >= HALO && thisY < _ny &&
			 thisZ >= HALO && thisZ < _nz )
		{
			stationNum ++;
		}

	}
	
	return stationNum;
}


void initStationIndex( GRID grid, STATION station ) 
{
	char jsonFile[1024] = { 0 };
	strcpy( jsonFile, "station.json" );
	FILE * fp;
	fp = fopen( jsonFile, "r" );

	if ( NULL == fp )
	{
		printf( "There is not %s file!\n", jsonFile );
		return;
	}
	
	fseek( fp, 0, SEEK_END );
	int len = ftell( fp );
	
	fseek( fp, 0, SEEK_SET );

	char * jsonStr = ( char * ) malloc( len * sizeof( char ) );

	if ( NULL == jsonStr )
	{
		printf( "Can't allocate json string memory\n" );
		return;
	}

	fread( jsonStr, sizeof( char ), len, fp );
	
	//printf( "%s\n", jsonStr );
	cJSON * object;
	cJSON * objArray;

	object = cJSON_Parse( jsonStr );
	if ( NULL == object )
	{
		printf( "Can't parse json file!\n");
		return;
	}

	fclose( fp );

	int stationCnt= 0;

	if (objArray = cJSON_GetObjectItem(object, "station(point)"))
	{
		stationCnt = cJSON_GetArraySize( objArray );
		//printf( "stationCnt = %d\n", stationCnt );
	}

	cJSON *stationObj, *stationItem;
	int i, j;
	int thisX, thisY, thisZ;
	
	int frontNX = grid.frontNX;
	int frontNY = grid.frontNY;
	int frontNZ = grid.frontNZ;

	int _nx = grid._nx;
	int _ny = grid._ny;
	int _nz = grid._nz;
	
	int X, Y, Z;

	int stationIdx = 0;
	for ( i = 0; i < stationCnt; i ++  )
	{
		stationObj = cJSON_GetArrayItem( objArray, i );

		int a = cJSON_GetArraySize( stationObj );
		if ( a != 3 )
		{
			printf( "In file %s, the coodinate index don't equal to 3. However, it equals to %d\n", jsonFile, a );
			return;
		}
	
		stationItem = cJSON_GetArrayItem( stationObj, 0 );
		X = stationItem->valueint;
		thisX = X - frontNX + HALO;

		stationItem = cJSON_GetArrayItem( stationObj, 1 );
		Y = stationItem->valueint;
		thisY = Y - frontNY + HALO;

		stationItem = cJSON_GetArrayItem( stationObj, 2 );
		Z = stationItem->valueint;
		thisZ = Z - frontNZ + HALO;
			
		if ( thisX >= HALO && thisX < _nx &&
			 thisY >= HALO && thisY < _ny &&
			 thisZ >= HALO && thisZ < _nz )
		{
			station.X[stationIdx] = thisX;
			station.Y[stationIdx] = thisY;
			station.Z[stationIdx] = thisZ;
			stationIdx ++;
		}
	}
}

void storageStation( GRID grid, int NT, int stationNum, STATION station, float *W_8, float*W_1, int it )
{
	long long num = stationNum;

	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;

	long long index_v = 0, index = 0, pos = 0;
	int X, Y, Z;

	float C1 = 1.0 / Cv;
	float C2 = 1.0 / Cs;

    for (int i = 0; i < stationNum; i++) {
		X = station.X[i];
		Y = station.Y[i];
		Z = station.Z[i];

		index_v = INDEX( X, Y, Z ) * WSIZE_V;
		index = INDEX( X, Y, Z );
		pos = ( it + i * NT ) * WSIZE;
		station.wave[pos + 0] = W_8[index + 0] * C1;
		station.wave[pos + 1] = W_8[index + 1] * C1;
		station.wave[pos + 2] = W_8[index + 2] * C1;
		station.wave[pos + 3] = W_8[index + 3] * C2;
		station.wave[pos + 4] = W_8[index + 4] * C2;
		station.wave[pos + 5] = W_8[index + 5] * C2;
		station.wave[pos + 6] = W_8[index + 6] * C2;
		station.wave[pos + 7] = W_8[index + 7] * C2;
		station.wave[pos + 8] = W_1[index] * C2;
    }
}


void stationWrite(PARAMS params, GRID grid, MPI_COORD thisMPICoord, STATION station, int NT, int stationNum)
{	
	FILE * fp;
	char fileName[256];
	sprintf( fileName, "%s/station_mpi_%d_%d_%d.bin", params.OUT, thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );

	for (int i = 0; i < stationNum; i ++){
		station.X[i] = grid.frontNX + station.X[i] - HALO;
		station.Y[i] = grid.frontNY + station.Y[i] - HALO;
		station.Z[i] = grid.frontNZ + station.Z[i] - HALO;
	}

	fp = fopen( fileName, "wb" ); 

	fwrite( station.X, sizeof( int ), stationNum * 3, fp );
	fwrite( station.wave, sizeof( float ), NT * stationNum * WSIZE, fp );
	
	fclose( fp );
}