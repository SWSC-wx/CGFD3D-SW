/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:headers.h
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-09-04
*   Discription:
*
================================================================*/
#ifndef __HEADERS__
#define __HEADERS__

#include <sys/types.h>
#include <sys/stat.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// #include <fstream>
#include <iostream>
#include <iomanip>
using namespace std;
#include <list>
#include <map>
#include <vector>

#include "macro.h"
#include "struct.h"
#include "readJson.h"
#include "cJSON.h"

#include <proj.h>
#include <netcdf.h>
#include <mpi.h>



#include "functions.h"
//#include "global_variable.h"
#include "params.h"

#ifdef __cplusplus
extern "C"
{
#endif

#include "pml.h"
#include "common.h"
#include "snapshot.h"
#include "timer.h"
#include "params.h"
#include "macdrp.h"
#include "logger.h"
#include "data_io.h"
#include "station.h"
#include "pgv.h"

#ifdef __cplusplus
}
#endif

#endif //__HEADERS__
