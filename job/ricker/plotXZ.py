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

var = sys.argv[2] #Vx Vz


jsonsFile = open( "params.json" )
params = json.load( jsonsFile )
grid = GRID( params )

sample = 1

outputPath = params["out"]
#outputPath = 'output_mpe'
fileNameX = outputPath + "/coordX/coordX"
fileNameY = outputPath + "/coordY/coordY"
fileNameZ = outputPath + "/coordZ/coordZ"

varname = "./png" + "/XZ%s_%d"%( var, it )

fileName = outputPath + "/it_%d/%s/%s_%d" % (it,var,var,it)
#fileName = params["out"] + "/Vs"


FAST_AXIS = params["FAST_AXIS"]

SAMP = params["SAMP"]
SAMPLE_SIZE = params["SAMPLE_SIZE"]
Gather = params["Gather"]
Gather_SIZE_X = params["Gather_SIZE_X"]
Gather_SIZE_Y = params["Gather_SIZE_Y"]

SUB_SIZE = Gather_SIZE_X

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
	dataX = np.zeros( [grid.NX//SAMPLE_SIZE, grid.NZ//SAMPLE_SIZE] )
	dataZ = np.zeros( [grid.NX//SAMPLE_SIZE, grid.NZ//SAMPLE_SIZE] )
	data  = np.zeros( [grid.NX//SAMPLE_SIZE, grid.NZ//SAMPLE_SIZE] )
else:
	dataX = np.zeros( [grid.NZ//SAMPLE_SIZE, grid.NX//SAMPLE_SIZE] )
	dataZ = np.zeros( [grid.NZ//SAMPLE_SIZE, grid.NX//SAMPLE_SIZE] )
	data  = np.zeros( [grid.NZ//SAMPLE_SIZE, grid.NX//SAMPLE_SIZE] )


offset_x = np.zeros( grid.PX// Gather_SIZE_X, dtype = 'int32' )
offset_y = np.zeros( grid.PY, dtype = 'int32' )
offset_z = np.zeros( grid.PZ, dtype = 'int32' )
jj = 0
kk = 0
for jj in range( grid.PX// Gather_SIZE_X ):
	begin = jj*Gather_SIZE_X
	for kk in range(begin,begin + Gather_SIZE_X):
		offset_x[jj] += grid.nx[kk]


ii = 0
mpiY = mpiSliceY
for mpiZ in range( grid.PZ ):
	for ii in range( grid.PX // Gather_SIZE_X):
		mpiX = ii * Gather_SIZE_X
		'''
		if Gather == 0 :
			fileNameX_1 = fileNameX + 'ndir%d/'%(mpiX%SUB_SIZE)
			fileNameZ_1 = fileNameZ + 'ndir%d/'%(mpiX%SUB_SIZE)
			fileName_1 = fileName + 'ndir%d/'%(mpiX%SUB_SIZE)
		else:
			fileNameX_1 = fileNameX 
			fileNameZ_1 = fileNameZ 
			fileName_1 = fileName 
		'''
		fileX = open( "%s_Y_mpi_%d_%d_%d.bin" % ( fileNameX, mpiX, mpiY, mpiZ ), "rb" )
		fileZ = open( "%s_Y_mpi_%d_%d_%d.bin" % ( fileNameZ, mpiX, mpiY, mpiZ ), "rb" )
		file  = open( "%s_Y_mpi_%d_%d_%d.bin" % ( fileName , mpiX, mpiY, mpiZ ), "rb" )
		if Gather:
			nz = grid.nz[mpiZ]//SAMPLE_SIZE
			nx = offset_x[ii]
		else:
			nx = grid.nx[mpiX]//SAMPLE_SIZE
			nz = grid.nz[mpiZ]//SAMPLE_SIZE
		
		print( "nx = %d, nz = %d" % ( nx, nz ) )
		datax = np.fromfile( fileX, dtype='float32', count = int(nx* nz) )
		#print( np.shape( datax ) )
		dataz = np.fromfile( fileZ, dtype='float32', count = int(nx* nz) )
		data_  = np.fromfile( file, dtype='float32', count = int(nx* nz) )

		I  = grid.frontNX[mpiX]//SAMPLE_SIZE
		I_ = grid.frontNX[mpiX]//SAMPLE_SIZE + nx
		K  = grid.frontNZ[mpiZ]//SAMPLE_SIZE
		K_ = grid.frontNZ[mpiZ]//SAMPLE_SIZE + nz

		if FAST_AXIS == 'Z':
			#print("I = %d,I_ = %d,K = %d,K_ = %d,nx = %d,nz =%d\n"%(I,I_,K,K_,nx,nz))
			dataX[I:I_, K:K_] = np.reshape( datax, ( nx, nz ) )
			dataZ[I:I_, K:K_] = np.reshape( dataz, ( nx, nz ) )
			data [I:I_, K:K_] = np.reshape( data_, ( nx, nz ) )
		else:
			dataX[K:K_, I:I_] = np.reshape( datax, ( nz, nx ) )
			dataZ[K:K_, I:I_] = np.reshape( dataz, ( nz, nx ) )
			data [K:K_, I:I_] = np.reshape( data_, ( nz, nx ) )
		if it == 20000:
			print("datax = %f\n"%(datax[2000]))
			print("dataz = %f\n"%(dataz[2000]))
			print("data_ = %f\n"%(data_[2000]))
		

dpi = 300
#fig = plt.figure( dpi = dpi, figsize = ( 1920 // dpi, 1080 // dpi ) )
unit = 1 # 1km = 1000m
vm = np.max( np.abs( data ) ) / 20
#plt.pcolormesh( dataX, dataZ, dataX, cmap = "jet" )
plt.pcolormesh( dataX, dataZ, data, vmax = vm, vmin = -vm, cmap = "seismic" )
#plt.pcolormesh( data[::sample, ::sample] // unit, vmin = -vm / 2, vmax = vm, cmap = "jet" )
#plt.pcolormesh( data[5:grid.NZ -5:sample, 5:grid.NX - 5:sample] // unit, cmap = "seismic" )
#plt.pcolormesh( data, vmax = vm, vmin = -vm, cmap = "jet" )
#plt.imshow( data, cmap = "jet", origin= "lower" )
plt.colorbar( )
#plt.colorbar( orientation = "horizontal" )
plt.axis( "image" )
plt.title( fileName )
plt.savefig( varname + ".png" )
#plt.plot( dataZ[grid.NZ - 1, :] )

#print( grid )


