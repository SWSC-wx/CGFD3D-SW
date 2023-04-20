/*
  _____________________________________________________________________________  
 | Name   | Type      | Command | Description                                  |   
 |________|___________|_________|______________________________________________| 
 | TMAX   | <FLOAT>   |   -T    | propagation time                             |   
 | DT     | <FLOAT>   |   -t    | time step (seconds)                          |   
 | TSKP   | <INTEGER> |   -s    | Time skip steps for snapshot output          |   
 |________|___________|_________|______________________________________________| 
 | DH     | <FLOAT>   |   -H    | spatial step for x, z (meters)               |  
 | NX     | <INTEGER> |   -X    | x dimension                                  |  
 | NY     | <INTEGER> |   -Y    | y dimension                                  |  
 | NZ     | <INTEGER> |   -Z    | z dimension                                  |  
 | PX     | <INTEGER> |   -x    | partition of x dimension                     |  
 | PY     | <INTEGER> |   -y    | partition of y dimension                     |  
 | PZ     | <INTEGER> |   -z    | partition of z dimension                     |  
 |________|___________|_________|______________________________________________| 
 | INGRD  | <STRING>  |         | grid input file                              |   
 | INVEL  | <STRING>  |         | velocity model input file                    |   
 | INSRC  | <STRING>  |         | source input file                            |   
 | INREC  | <STRING>  |         | receiver input file                          |   
 |________|___________|_________|______________________________________________| 
                                                                            
*/

//#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "params.h"
#include "macdrp.h"
#include "common.h"
#include "cJSON.h"
#include "snapshot.h"

// Default IN3D Values
const float def_TMAX  = 20.00;
const float def_DT    = 0.01;
const float def_TSKP  = 50;
const float def_FM  = 1.2;

const float def_DH = 100.0;
const int   def_NX = 32;
const int   def_NY = 32;
const int   def_NZ = 2;

const int def_NBGX = 1;
const int def_NBGY = 1;
const int def_NBGZ = 1;

const int def_NEDZ = 1;

const int def_NSKPX = 1;
const int def_NSKPY = 1;
const int def_NSKPZ = 1;

const int def_WRITE_STEP = 1;

const int def_SAVE_FULL_IMG = 0;
const int def_SAVE_FULL_IMG_START_STEP = 0;
const int def_SAVE_FULL_IMG_PER_STEP = 100;

const int def_SAVE_CKPT = 0;
const int def_SAVE_CKPT_START_STEP = 100;
const int def_SAVE_CKPT_PER_STEP = 100;
const int def_CUR_CKPT_IDX = 0;
const int def_NUM_CKPT_DIR = 3;
const int def_LOAD_CKPT = 0;

const int def_SAMP = 0;
const int def_SAMPLE_SIZE = 1;
const int def_Gather = 0;
const int def_Gather_SIZE_X = 1;
const int def_Gather_SIZE_Y = 1;

const char def_OUT[100] = "output_sfc";

const char def_INGRD[128] = "./input/coord";
const char def_INVEL[128] = "./input/media.bin";
const char def_INSRC[128] = "./input/source";
const char def_INREC[128] = "./input/receiver.bin";

void get_params( int argc, char **argv) {
  // Fill in default values
  TMAX = def_TMAX;
  DT   = def_DT;
  TSKP = def_TSKP;
  rickerfc = def_FM;

  DH = def_DH;
  NX = def_NX;
  NY = def_NY;
  NZ = def_NZ;

  strcpy(INGRD, def_INGRD);
  strcpy(INVEL, def_INVEL);
  strcpy(INSRC, def_INSRC);
  strcpy(INREC, def_INREC);

  NBGX = def_NBGX;
  NBGY = def_NBGY;
  NBGZ = def_NBGZ;
       
  NEDZ = def_NEDZ;
       
  NSKPX = def_NSKPX;
  NSKPY = def_NSKPY;
  NSKPZ = def_NSKPZ;

  WRITE_STEP = def_WRITE_STEP;

  SAVE_FULL_IMG = def_SAVE_FULL_IMG;
  SAVE_FULL_IMG_START_STEP = def_SAVE_FULL_IMG_START_STEP;
  SAVE_FULL_IMG_PER_STEP = def_SAVE_FULL_IMG_PER_STEP;

  SAVE_CKPT = def_SAVE_CKPT;
  SAVE_CKPT_START_STEP = def_SAVE_CKPT_START_STEP;
  SAVE_CKPT_PER_STEP = def_SAVE_CKPT_PER_STEP;
  CUR_CKPT_IDX = def_CUR_CKPT_IDX;
  NUM_CKPT_DIR = def_NUM_CKPT_DIR;
  LOAD_CKPT = def_LOAD_CKPT;

  SAMPLE_SIZE = def_SAMPLE_SIZE;
  SAMPLE_SIZE = def_SAMPLE_SIZE;

  Gather = def_Gather;
  Gather_SIZE_X = def_Gather_SIZE_X;
  Gather_SIZE_Y = def_Gather_SIZE_Y;

  FILE *fp = fopen("params.json","r");
  fseek(fp, 0, SEEK_END);
  long len = ftell(fp);
  fseek(fp, 0, SEEK_SET);
  char *str = (char*)malloc(len+1);
  fread(str, 1, len, fp);
  fclose(fp);

  //printf("configure json file:\n%s\n", str);

  cJSON *root = cJSON_Parse(str);
  if (NULL == root) printf("Error!\n");

  cJSON *item;
  item = cJSON_GetObjectItem(root, "IN_SOURCE");
  memcpy(IN_SOURCE, item->valuestring, strlen(item->valuestring));
  item = cJSON_GetObjectItem(root, "IN_TOPO");
  memcpy(IN_TOPO, item->valuestring, strlen(item->valuestring));
  item = cJSON_GetObjectItem(root, "IN_VELO");
  memcpy(IN_VELO, item->valuestring, strlen(item->valuestring));
  item = cJSON_GetObjectItem(root, "IN_INTER");
  memcpy(IN_INTER, item->valuestring, strlen(item->valuestring));
  item = cJSON_GetObjectItem(root, "IN_REC");
  memcpy(IN_REC, item->valuestring, strlen(item->valuestring));
	item = cJSON_GetObjectItem(root, "OUT_DIR");
  //cJSON *name_json;
  memcpy(OUT_DIR, item->valuestring, strlen(item->valuestring));
  //name_json = cJSON_GetObjectItem(root, "INGRD");
  //if (NULL != name_json)
  //{
  //    INGRD = cJSON_Print(name_json);
  //    printf("INGRD : %s\n", INGRD);
  //    free(name_json);
  //}

  //name_json = cJSON_GetObjectItem(root, "INVEL");
  //if (NULL != name_json)
  //{
  //    INVEL = cJSON_Print(name_json);
  //    printf("INVEL : %s\n", INVEL);
  //    free(name_json);
  //}

  //name_json = cJSON_GetObjectItem(root, "INSRC");
  //if (NULL != name_json)
  //{
  //    INSRC = cJSON_Print(name_json);
  //    printf("INSRC : %s\n", INSRC);
  //    free(name_json);
  //}

  //name_json = cJSON_GetObjectItem(root, "INREC");
  //if (NULL != name_json)
  //{
  //    INREC = cJSON_Print(name_json);
  //    printf("INREC : %s\n", INREC);
  //    free(name_json);
  //}

  NUM_INTER  = cJSON_GetObjectItem(root, "NUM_INTER")->valueint;
  TSKP0  = cJSON_GetObjectItem(root, "TSKP0")->valueint;
  TSKP  = cJSON_GetObjectItem(root, "TSKP")->valueint;
  NX  = cJSON_GetObjectItem(root, "NX")->valueint;
  NY  = cJSON_GetObjectItem(root, "NY")->valueint;
  NZ  = cJSON_GetObjectItem(root, "NZ")->valueint;
  PX  = cJSON_GetObjectItem(root, "PX")->valueint;
  PY  = cJSON_GetObjectItem(root, "PY")->valueint;
  PZ  = cJSON_GetObjectItem(root, "PZ")->valueint;
  icoord = cJSON_GetObjectItem(root, "icoord"  )->valueint;
  isource = cJSON_GetObjectItem(root, "isource"  )->valueint;
  DH  = cJSON_GetObjectItem(root, "DH")->valuedouble;
  TMAX  = cJSON_GetObjectItem(root, "TMAX")->valuedouble;
  DT  = cJSON_GetObjectItem(root, "DT")->valuedouble;
  rickerfc = (float)cJSON_GetObjectItem(root, "rickerfc" )->valuedouble;
  Qsf0 = (float)cJSON_GetObjectItem(root, "Qsf0" )->valuedouble;
  gauss_height = (float)cJSON_GetObjectItem(root, "gauss_height" )->valuedouble;
  gauss_width  = (float)cJSON_GetObjectItem(root, "gauss_width" )->valuedouble;
  cJSON *p;
  if(p = cJSON_GetObjectItem(root, "SAVE_CKPT" ))
    SAVE_CKPT = p->valueint;
  if(p = cJSON_GetObjectItem(root, "SAVE_CKPT_START_STEP" ))
    SAVE_CKPT_START_STEP = p->valueint;
  if(p = cJSON_GetObjectItem(root, "SAVE_CKPT_PER_STEP" ))
    SAVE_CKPT_PER_STEP = p->valueint;
  if(p = cJSON_GetObjectItem(root, "CUR_CKPT_IDX" ))
    CUR_CKPT_IDX = p->valueint;
  if(p = cJSON_GetObjectItem(root, "NUM_CKPT_DIR" ))
    NUM_CKPT_DIR = p->valueint;
  if(p = cJSON_GetObjectItem(root, "LOAD_CKPT" ))
    LOAD_CKPT = p->valueint;

  if(p = cJSON_GetObjectItem(root, "SAMP" ))
    SAMP = p->valueint;
  if(p = cJSON_GetObjectItem(root, "SAMPLE_SIZE" ))
    SAMPLE_SIZE = p->valueint;

  if(p = cJSON_GetObjectItem(root, "Gather" ))
    Gather = p->valueint;
  if(p = cJSON_GetObjectItem(root, "Gather_SIZE_X" ))
    Gather_SIZE_X = p->valueint;
  if(p = cJSON_GetObjectItem(root, "Gather_SIZE_Y" ))
    Gather_SIZE_Y = p->valueint;

  if(p = cJSON_GetObjectItem(root, "NBGX" ))
    NBGX = p->valueint;
  if(p = cJSON_GetObjectItem(root, "NBGY" ))
    NBGY = p->valueint;
  if(p = cJSON_GetObjectItem(root, "NBGZ" ))
    NBGZ = p->valueint;
  if(p = cJSON_GetObjectItem(root, "NEDX" ))
    NEDX = p->valueint;
  if(p = cJSON_GetObjectItem(root, "NEDY" ))
    NEDY = p->valueint;
  if(p = cJSON_GetObjectItem(root, "NEDZ" ))
    NEDZ = p->valueint;
  if(p = cJSON_GetObjectItem(root, "NSKPX" ))
    NSKPX = p->valueint;
  if(p = cJSON_GetObjectItem(root, "NSKPY" ))
    NSKPY = p->valueint;
  if(p = cJSON_GetObjectItem(root, "NSKPZ" ))
    NSKPZ = p->valueint;
  if(p = cJSON_GetObjectItem(root, "WRITE_STEP" ))
    WRITE_STEP= p->valueint;
  if(p = cJSON_GetObjectItem(root, "SAVE_FULL_IMG" ))
    SAVE_FULL_IMG = p->valueint;
  if(p = cJSON_GetObjectItem(root, "SAVE_FULL_IMG_START_STEP" ))
    SAVE_FULL_IMG_START_STEP = p->valueint;
  if(p = cJSON_GetObjectItem(root, "SAVE_FULL_IMG_PER_STEP" ))
    SAVE_FULL_IMG_PER_STEP = p->valueint;
  if(p = cJSON_GetObjectItem(root, "OUT" )) {
    strcpy(OUT, p->valuestring);
  }

  cJSON_Delete(root);

  mkpath(OUT, 0700);
  mkpath(OUT_DIR, 0700);

  //printf("---------------------------------------------\n");
  //printf("TSKP  : %d\n", TSKP);
  //printf("TMAX  : %f\n", TMAX);
  //printf("DT  : %f\n", DT);
  //printf("DH  : %f\n", DH);
  //printf("NX  : %d\n", NX);
  //printf("NY  : %d\n", NY);
  //printf("NZ  : %d\n", NZ);
  //printf("PX  : %d\n", PX);
  //printf("PY  : %d\n", PY);
  //printf("PZ  : %d\n", PZ);
  //printf("---------------------------------------------\n");

  rDH = 1.0 / DH;
  int i;
  for (i = 0; i < 5; i++){
    coeF_DH[i] = coeF[i] * rDH;
    coeB_DH[i] = coeB[i] * rDH;
  }
  for (i = 0; i < 3; i++){
    coe24F_DH[i] = coe24F[i] * rDH;
    coe24B_DH[i] = coe24B[i] * rDH;
  }

  return;
}
