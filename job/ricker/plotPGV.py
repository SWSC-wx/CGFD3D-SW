'''
Author: Wenqiang Wang @ SUSTech on Sep 11, 2021
15:40
'''
import json
import numpy as np
from pyscripts.GRID import GRID
import matplotlib.pyplot as plt
import  sys
import math


FreeSurf = 1

it = 1000

var = 'Vy' #Vy Vz



if len( sys.argv ) > 1:
    it = int( sys.argv[1] )




jsonsFile = open( "params.json" )
params = json.load( jsonsFile )
grid = GRID( params )
FAST_AXIS = params["FAST_AXIS"]

sample = 5

outputPath = params["out"]
fileNameX = params["out"] + "/coordX/coordX"
fileNameY = params["out"] + "/coordY/coordY"

Gather = params["Gather"]
Gather_SIZE_X = params["Gather_SIZE_X"]
Gather_SIZE_Y = params["Gather_SIZE_Y"]

if Gather == 0:
    Gather_SIZE_X = 1
    Gather_SIZE_Y = 1


if FreeSurf == 1:
    #fileName = params["out"] + "/FreeSurf%s_%d" % ( var, it )
    fileName = params["out"] + "/PGV"
else:
    #fileName = params["out"] + "/%s_%d" % ( var, it )
    fileName = outputPath + "/it_%d/%s/%s_%d" % (it,var,var,it)



'''
for Z in range( grid.PZ ):
    for Y in range( grid.PY ):
        for X in range( grid.PX ):
            print( "nx = %d, ny = %d, nz = %d\n" % ( grid.nx[X], grid.ny[Y], grid.nz[Z] ) )
'''

sliceX = params["sliceX"] - grid.frontNX
sliceY = params["sliceY"] - grid.frontNY
sliceZ = params["sliceZ"] - grid.frontNZ

for mpiSliceX in range( grid.PX ):
    if sliceX[mpiSliceX] >= 0 and sliceX[mpiSliceX] < grid.nx[mpiSliceX]:
        break

for mpiSliceY in range( grid.PY ):
    if sliceY[mpiSliceY] >= 0 and sliceY[mpiSliceY] < grid.ny[mpiSliceY]:
        break

for mpiSliceZ in range( grid.PZ ):
    if sliceZ[mpiSliceZ] >= 0 and sliceZ[mpiSliceZ] < grid.nz[mpiSliceZ]:
        break



if FAST_AXIS == 'Z':
    dataX = np.zeros( [grid.NX, grid.NY] )
    dataY = np.zeros( [grid.NX, grid.NY] )
    data  = np.zeros( [grid.NX, grid.NY] )
    Intensity = np.zeros( [grid.NX, grid.NY] )
    print(grid.NX, grid.NY)
else:
    dataX = np.zeros( [grid.NY, grid.NX] )
    dataY = np.zeros( [grid.NY, grid.NX] )
    data  = np.zeros( [grid.NY, grid.NX] )
    Intensity = np.zeros( [grid.NY, grid.NX] )

offset_x = np.zeros( grid.PX, dtype = 'int32' )

ii = 0
jj = 0
kk = 0
for jj in range( grid.PX// Gather_SIZE_X ):
    begin = jj*Gather_SIZE_X
    for kk in range(begin,begin +Gather_SIZE_X):
        offset_x[jj] += grid.nx[kk]


NUM_RANK_PER_DIR = 5000

for mpiX in range( grid.PX ):
    for mpiY in range( grid.PY):
        mpiZ = mpiSliceZ
        ii = mpiX // Gather_SIZE_X
        if (mpiX % Gather_SIZE_X == 0):
            XFile = open( "%s_Z_mpi_%d_%d_%d.bin" % ( fileNameX, mpiX, mpiY, mpiZ ), "rb" )
            YFile = open( "%s_Z_mpi_%d_%d_%d.bin" % ( fileNameY, mpiX, mpiY, mpiZ ), "rb" )
        if FreeSurf:
            mpiZ = grid.PZ - 1
        else:
            mpiZ = mpiSliceZ
        NUM_DIRS = math.ceil(1.0 * grid.PX * grid.PY / NUM_RANK_PER_DIR)

        if (mpiX % Gather_SIZE_X == 0):
            nx = offset_x[ii]
            ny = grid.ny[mpiY]
            J  = grid.frontNY[mpiY]
            J_ = grid.frontNY[mpiY] + ny
            I  = grid.frontNX[mpiX]
            I_ = grid.frontNX[mpiX] + nx

            datax = np.fromfile( XFile, dtype='float32', count = ny * nx )
            datay = np.fromfile( YFile, dtype='float32', count = ny * nx )
            if FAST_AXIS == 'Z':
                dataX[I:I_, J:J_] = np.reshape( datax, (nx, ny) )
                dataY[I:I_, J:J_] = np.reshape( datay, (nx, ny) )
            else:
                dataX[J:J_, I:I_] = np.reshape( datax, (ny, nx) )
                dataY[J:J_, I:I_] = np.reshape( datay, (ny, nx) )

        nx = grid.nx[mpiX]
        ny = grid.ny[mpiY]
        J  = grid.frontNY[mpiY]
        J_ = grid.frontNY[mpiY] + ny
        I  = grid.frontNX[mpiX]
        I_ = grid.frontNX[mpiX] + nx
        FileNamePGV = "%s/PGV/ndir%d-%d/%s_Z_mpi_%d_%d_%d.bin" % (params['out'], NUM_DIRS, (mpiX*grid.PY + mpiY)%NUM_DIRS, "PGVh", mpiX, mpiY, mpiZ)
        print(FileNamePGV)
        File  = open( "%s" % (FileNamePGV), "rb" )
        data_ = np.fromfile( FileNamePGV, dtype='float32', count = ny * nx )
        print(max(data_))

        if FAST_AXIS == 'Z':
            data [I:I_, J:J_] = np.reshape( data_, (nx, ny) )
        else:
            data [J:J_, I:I_] = np.reshape( data_, (ny, nx) )

PGV = data/10.0

print(PGV.shape)
'''
for j in range( grid.NY ):
    for i in range( grid.NX ):
        if  PGV[j,i] >= 1.41:
            Intensity[j,i] = 11
        if  PGV[j,i] >= 0.72 and  PGV[j,i] < 1.41:
            Intensity[j,i] = 10
        if  PGV[j,i] >= 0.36 and  PGV[j,i] < 0.72:
            Intensity[j,i] = 9
        if  PGV[j,i] >= 0.19 and  PGV[j,i] < 0.36:
            Intensity[j,i] = 8
        if  PGV[j,i] >= 0.10 and  PGV[j,i] < 0.19:
            Intensity[j,i] = 7
        if  PGV[j,i] >= 0.05 and  PGV[j,i] < 0.1:
            Intensity[j,i] = 6
        if  PGV[j,i] >= 0.02 and  PGV[j,i] < 0.05:
            Intensity[j,i] = 5
        if  PGV[j,i] >= 0.01 and  PGV[j,i] < 0.02:
            Intensity[j,i] = 4
        if  PGV[j,i] >= 0.005 and  PGV[j,i] < 0.01:
            Intensity[j,i] = 3
        if  PGV[j,i] >= 0.001 and  PGV[j,i] < 0.005:
            Intensity[j,i] = 2
        if  PGV[j,i] < 0.001:
            Intensity[j,i] = 1
'''
for i in range( grid.NX ):
    print("----> %d/%d\t\r" %(i, grid.NX))
    for j in range( grid.NY ):
        if  PGV[i,j] >= 1.41:
            Intensity[i,j] = 11
        if  PGV[i,j] >= 0.72 and  PGV[i,j] < 1.41:
            Intensity[i,j] = 10
        if  PGV[i,j] >= 0.36 and  PGV[i,j] < 0.72:
            Intensity[i,j] = 9
        if  PGV[i,j] >= 0.19 and  PGV[i,j] < 0.36:
            Intensity[i,j] = 8
        if  PGV[i,j] >= 0.10 and  PGV[i,j] < 0.19:
            Intensity[i,j] = 7
        if  PGV[i,j] >= 0.05 and  PGV[i,j] < 0.1:
            Intensity[i,j] = 6
        if  PGV[i,j] >= 0.02 and  PGV[i,j] < 0.05:
            Intensity[i,j] = 5
        if  PGV[i,j] >= 0.01 and  PGV[i,j] < 0.02:
            Intensity[i,j] = 4
        if  PGV[i,j] >= 0.005 and PGV[i,j] < 0.01:
            Intensity[i,j] = 3
        if  PGV[i,j] >= 0.001 and PGV[i,j] < 0.005:
            Intensity[i,j] = 2
        if  PGV[i,j] < 0.001:
            Intensity[i,j] = 1


print(np.max(Intensity))
nPML = params["nPML"]
NX = grid.NX
NY = grid.NY

'''
np.save( "./npy/PGV.npy", PGV)
logPGV = np.log( PGV )
np.save( "./npy/logPGV.npy", logPGV)

#np.save( "Intensity.npy", Intensity[nPML:NY - nPML, nPML:NX-nPML] )
np.save( "PGV.npy", PGV [::sample, ::sample]  )
logPGV = np.log( PGV )
np.save( "logPGV.npy", logPGV [::sample, ::sample] )
'''

dpi = 300
#fig = plt.figure( dpi = dpi, figsize = ( 1920 // dpi, 1080 // dpi ) )
unit = 1 # 1km = 1000m
vm = np.max( np.abs( PGV ) )
print("start ploting\n")
plt.contourf( dataX[::sample, ::sample], dataY[::sample, ::sample], Intensity[::sample, ::sample], [5, 6, 7, 8, 9, 10, 11], cmap = "seismic" )
#plt.pcolormesh( data[::sample, ::sample] // unit, cmap = "seismic" )
#plt.pcolormesh(np.log( PGV ), cmap = "seismic" ) #origin= "lower" )
#plt.pcolormesh(PGV, cmap = "seismic" ) #origin= "lower" )
#plt.pcolormesh( dataX[::sample, ::sample], dataY[::sample, ::sample], logPGV[::sample, ::sample], cmap = "jet" ) #origin= "lower" )
#plt.pcolormesh( dataX[::sample, ::sample], dataY[::sample, ::sample], logPGV[::sample, ::sample], cmap = "jet", vmin=0, vmax=vm ) #origin= "lower" )
plt.colorbar( )
plt.axis( "image" )
plt.title( fileName )

plt.show( )
plt.savefig( fileName + "2.pdf" )
#plt.plot( dataY[grid.NZ - 1, :] )

#print( grid )


