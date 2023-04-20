#!/bin/bash
###
 # @Descripttion: your project
 # @version: 1.0
 # @Author: liuxiaohui
 # @Date: 2023-04-18 15:50:31
 # @LastEditors: liuxiaohui
 # @LastEditTime: 2023-04-18 17:41:17
### 

#MPIHOME=
PROJHOME=
SQLITE3=
NETCDFHOME=

export LD_LIBRARY_PATH=${PROJHOME}/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${SQLITE3}/lib:${LD_LIBRARY_PATH}
#export LD_LIBRARY_PATH=${MPIHOME}/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${NETCDFHOME}/lib:${LD_LIBRARY_PATH}
export PROJ_LIB=${PROJHOME}/share/proj

#RUN=${MPIHOME}/bin/mpirun
PX=`cat params.json | grep "\"PX\"" | tr -cd "[0-9]"`
PY=`cat params.json | grep "\"PY\"" | tr -cd "[0-9]"`
PZ=`cat params.json | grep "\"PZ\"" | tr -cd "[0-9]"`

#${RUN} -np $1 ./bin/main
date=`date +%m%d-%H-%M-%S`
bsub -J cgfd3 -o ./log/${date}.log -debug  ./cgfd3d
