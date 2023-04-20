#ifndef __READJSON__
#define __READJSON__
typedef struct PARAMS
{
	
    double TMAX;
    
	double DT;
    double DH;
    
	int NX;
    int NY;
    int NZ;

    int PX;
    int PY;
    int PZ;

	int centerX;
	int centerY;

	double centerLatitude ; 
	double centerLongitude;
    
	int SourceX; 
	int SourceY;
	int SourceZ;

    int icoord;	
	int imedium;
    int useMedium;
    int useTerrain;
	int useMultiSource;
	int sliceX; 
	int sliceY;
	int sliceZ;
	int sliceFreeSurf;


	int nPML;
	
	int SliceX; 
	int SliceY;
	int SliceZ;

	int itSlice			;
	int itStep			;
	char waveOutput[64]	;
	char sliceName[64]	;
	int itStart			;
	int itEnd			;
	int igpu			;
    char OUT[256];
    
	
	char TerrainDir[256];
	
	
	int SRTM90;
	
	int lonStart;
	int latStart;
	int blockX;
	int blockY;
		
	float Depth;


	float MLonStart;
	float MLatStart;
	float MLonEnd;
	float MLatEnd;

	float MLonStep;
	float MLatStep;

	float MVeticalStep;
	
	char MediumDir[256];

	char sourceFile[256];


}PARAMS;

#endif
