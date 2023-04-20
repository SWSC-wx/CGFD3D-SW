'''
Author: Wenqiang Wang @ SUSTech on Sep 11, 2021
15:40
'''
import json
import numpy as np
from pyscripts.GRID import GRID
import matplotlib.pyplot as plt
import sys


it = 600
if len( sys.argv ) > 1:
	it = int( sys.argv[1] )

var = sys.argv[2] #Vy Vz


jsonsFile = open( "params.json" )
params = json.load( jsonsFile )
grid = GRID( params )
FAST_AXIS = params["FAST_AXIS"]

sample = 1

outputPath = params["out"]
fileNameY = outputPath + "/coordY/coordY"
fileNameZ = outputPath + "/coordZ/coordZ"

varname = "./png" + "/YZ%s_%d"%( var, it )

fileName = outputPath + "/it_%d/%s/%s_%d" % (it,var,var,it)
#fileName = params["out"] + "/Vs"

SAMP = params["SAMP"]
SAMPLE_SIZE = params["SAMPLE_SIZE"]
Gather = params["Gather"]
Gather_SIZE_X = params["Gather_SIZE_X"]
Gather_SIZE_Y = params["Gather_SIZE_Y"]

SUB_SIZE = Gather_SIZE_Y

if SAMP == 0:
	SAMPLE_SIZE = 1
if Gather == 0:
	Gather_SIZE_X = 1
	Gather_SIZE_Y = 1


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
	dataY = np.zeros( [grid.NY//SAMPLE_SIZE, grid.NZ//SAMPLE_SIZE] )
	dataZ = np.zeros( [grid.NY//SAMPLE_SIZE, grid.NZ//SAMPLE_SIZE] )
	data  = np.zeros( [grid.NY//SAMPLE_SIZE, grid.NZ//SAMPLE_SIZE] )
else:
	dataY = np.zeros( [grid.NZ//SAMPLE_SIZE, grid.NY//SAMPLE_SIZE] )
	dataZ = np.zeros( [grid.NZ//SAMPLE_SIZE, grid.NY//SAMPLE_SIZE] )
	data  = np.zeros( [grid.NZ//SAMPLE_SIZE, grid.NY//SAMPLE_SIZE] )

offset_x = np.zeros( grid.PX, dtype = 'int32' )
offset_y = np.zeros( grid.PY, dtype = 'int32' )
offset_z = np.zeros( grid.PZ, dtype = 'int32' )

jj = 0
for jj in range( grid.PY// Gather_SIZE_Y ):
	begin = jj*Gather_SIZE_Y
	for begin in range(begin,begin +Gather_SIZE_Y):
		offset_y[jj] += grid.ny[begin]

ii = 0
mpiX = mpiSliceX
for mpiZ in range( grid.PZ ):
	for ii in range( grid.PY // Gather_SIZE_Y ):
		mpiY = ii * Gather_SIZE_Y

		'''
		if Gather == 0 :
			fileNameZ_1 = fileNameZ + 'ndir%d/'%(mpiY%SUB_SIZE)
			fileNameY_1 = fileNameY + 'ndir%d/'%(mpiY%SUB_SIZE)
			fileName_1 = fileName + 'ndir%d/'%(mpiY%SUB_SIZE)
		else:
			fileNameZ_1 = fileNameZ 
			fileNameY_1 = fileNameY 
			fileName_1 = fileName 
		'''
		fileY = open( "%s_X_mpi_%d_%d_%d.bin" % ( fileNameY, mpiX, mpiY, mpiZ ), "rb" )
		fileZ = open( "%s_X_mpi_%d_%d_%d.bin" % ( fileNameZ, mpiX, mpiY, mpiZ ), "rb" )
		file  = open( "%s_X_mpi_%d_%d_%d.bin" % ( fileName , mpiX, mpiY, mpiZ ), "rb" )
		if Gather:
			nz = grid.nz[mpiZ]//SAMPLE_SIZE
			ny = offset_y[ii]
		else:
			ny = grid.ny[mpiY]//SAMPLE_SIZE
			nz = grid.nz[mpiZ]//SAMPLE_SIZE
		
		print( "ny = %d, nz = %d" % ( ny, nz ) )
		datay = np.fromfile( fileY, dtype='float32', count = int(ny * nz) )
		dataz = np.fromfile( fileZ, dtype='float32', count = int(ny * nz) )
		data_ = np.fromfile( file , dtype='float32', count = int(ny * nz) )

		J  = grid.frontNY[mpiY]//SAMPLE_SIZE
		J_ = grid.frontNY[mpiY]//SAMPLE_SIZE + ny
		K  = grid.frontNZ[mpiZ]//SAMPLE_SIZE
		K_ = grid.frontNZ[mpiZ]//SAMPLE_SIZE + nz

		if FAST_AXIS == 'Z':
			dataY[J:J_, K:K_] = np.reshape( datay, ( ny, nz ) )
			dataZ[J:J_, K:K_] = np.reshape( dataz, ( ny, nz ) )
			data [J:J_, K:K_] = np.reshape( data_, ( ny, nz ) )
		else:
			dataY[K:K_, J:J_] = np.reshape( datay, ( nz, ny ) )
			dataZ[K:K_, J:J_] = np.reshape( dataz, ( nz, ny ) )
			data [K:K_, J:J_] = np.reshape( data_, ( nz, ny ) )

dpi = 300
#fig = plt.figure( dpi = dpi, figsize = ( 1920 // dpi, 1080 // dpi ) )
unit = 1 # 1km = 1000m
vm = np.max( np.abs( data ) ) / 20
#plt.pcolormesh( dataY, dataZ, data, vmax = vm, vmin = -vm, cmap = "jet" )
#plt.imshow( data, vmax = vm, vmin = -vm, cmap = "jet", origin= "lower" )
plt.pcolormesh( dataY, dataZ, data, vmax = vm, vmin = -vm, cmap = "seismic" ) #origin= "lower" )
#plt.pcolormesh( dataY, cmap = "jet" )
#plt.pcolormesh( data[::sample, ::sample] // unit, cmap = "seismic" )
plt.colorbar( )
plt.axis( "image" )
plt.title( fileName )
plt.savefig( varname + ".png" )
#plt.plot( dataZ[grid.NZ - 1, :] )

#print( grid )


