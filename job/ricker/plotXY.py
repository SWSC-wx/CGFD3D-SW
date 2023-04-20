'''
Author: Wenqiang Wang @ SUSTech on Sep 11, 2021
15:40
'''
import json
import numpy as np
from pyscripts.GRID import GRID
import matplotlib.pyplot as plt
import  sys


FreeSurf = 0

it =600


if len( sys.argv ) > 1:
	it = int( sys.argv[1] )

var = sys.argv[2] #Vx Vy



jsonsFile = open( "params.json" )
params = json.load( jsonsFile )
grid = GRID( params )
FAST_AXIS = params["FAST_AXIS"]

sample = 5

outputPath = params["out"]
fileNameX = outputPath + "/coordX/coordX"
fileNameY = outputPath + "/coordY/coordY"

SAMP = params["SAMP"]
SAMPLE_SIZE = params["SAMPLE_SIZE"]
Gather = params["Gather"]
Gather_SIZE_X = params["Gather_SIZE_X"]
Gather_SIZE_Y = params["Gather_SIZE_Y"]


varname = "./png" + "/XY%s_%d"%( var, it )

if FreeSurf == 1:
	fileName = outputPath + "/FreeSurf%s_%d" % ( var, it )
else:
	fileName = outputPath + "/it_%d/%s/%s_%d" % (it,var,var,it)

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
	dataX = np.zeros( [grid.NX//SAMPLE_SIZE, grid.NY//SAMPLE_SIZE] )
	dataY = np.zeros( [grid.NX//SAMPLE_SIZE, grid.NY//SAMPLE_SIZE] )
	data  = np.zeros( [grid.NX//SAMPLE_SIZE, grid.NY//SAMPLE_SIZE] )
else:
	dataX = np.zeros( [grid.NY//SAMPLE_SIZE, grid.NX//SAMPLE_SIZE] )
	dataY = np.zeros( [grid.NY//SAMPLE_SIZE, grid.NX//SAMPLE_SIZE] )
	data  = np.zeros( [grid.NY//SAMPLE_SIZE, grid.NX//SAMPLE_SIZE] )

offset_x = np.zeros( grid.PX, dtype = 'int32' )
offset_y = np.zeros( grid.PY, dtype = 'int32' )
offset_z = np.zeros( grid.PZ, dtype = 'int32' )

jj = 0
for jj in range( grid.PX// Gather_SIZE_X ):
	begin = jj*Gather_SIZE_X
	for begin in range(begin,begin +Gather_SIZE_X):
		offset_x[jj] += grid.nx[begin]
		#print("jj : %d,nx : %d"%(jj,offset_x[jj]))


ii = 0

for mpiY in range( grid.PY ):
	for ii in range( grid.PX // Gather_SIZE_X):
		mpiX = ii * Gather_SIZE_X
		mpiZ = mpiSliceZ
		'''
		if Gather == 0 :
			fileNameX_1 = fileNameX + 'ndir%d/'%(mpiX%SUB_SIZE)
			fileNameY_1 = fileNameY + 'ndir%d/'%(mpiX%SUB_SIZE)
			fileName_1 = fileName + 'ndir%d/'%(mpiX%SUB_SIZE)
		else:
			fileNameX_1 = fileNameX 
			fileNameY_1 = fileNameY 
			fileName_1 = fileName 
		'''
		XFile = open( "%s_Z_mpi_%d_%d_%d.bin" % ( fileNameX, mpiX, mpiY, mpiZ ), "rb" )
		YFile = open( "%s_Z_mpi_%d_%d_%d.bin" % ( fileNameY, mpiX, mpiY, mpiZ ), "rb" )
		if FreeSurf:
			mpiZ = grid.PZ - 1
		else:
			mpiZ = mpiSliceZ
		File  = open( "%s_Z_mpi_%d_%d_%d.bin" % ( fileName , mpiX, mpiY, mpiZ ), "rb" )
			
		#ny = grid.ny[mpiY]
		#nx = grid.nx[mpiX]
		if Gather:
			ny = grid.ny[mpiY]//SAMPLE_SIZE
			nx = offset_x[ii]
		else:
			nx = grid.nx[mpiX]//SAMPLE_SIZE
			ny = grid.ny[mpiY]//SAMPLE_SIZE
		print( "nx = %d, ny = %d" % ( nx, ny ) )
		datax = np.fromfile( XFile, dtype='float32', count = int(nx * ny) )
		datay = np.fromfile( YFile, dtype='float32', count = int(nx * ny) )
		data_ = np.fromfile(  File, dtype='float32', count = int(nx * ny) )
		

		J  = grid.frontNY[mpiY]//SAMPLE_SIZE
		J_ = grid.frontNY[mpiY]//SAMPLE_SIZE + ny
		I  = grid.frontNX[mpiX]//SAMPLE_SIZE
		I_ = grid.frontNX[mpiX]//SAMPLE_SIZE + nx

		if FAST_AXIS == 'Z':
			dataX[I:I_, J:J_] = np.reshape( datax, (nx, ny) )
			dataY[I:I_, J:J_] = np.reshape( datay, (nx, ny) )
			data [I:I_, J:J_] = np.reshape( data_, (nx, ny) )
		else:
			dataX[J:J_, I:I_] = np.reshape( datax, (ny, nx) )
			dataY[J:J_, I:I_] = np.reshape( datay, (ny, nx) )
			data [J:J_, I:I_] = np.reshape( data_, (ny, nx) )
		
		if it == 20000:
			print("data = %f\n"%(data_[2000]))

np.save( "./data_npy/%s-%d.npy"%(var, it), data[::sample, ::sample]  )
dpi = 300
#fig = plt.figure( dpi = dpi, figsize = ( 1920 // dpi, 1080 // dpi ) )
unit = 1 # 1km = 1000m
vm = np.max( np.abs( data ) ) / 2
print("vm = ", vm)
vm = 0.1
#plt.pcolormesh( data[::sample, ::sample] // unit, cmap = "seismic" )
#plt.pcolormesh( dataX, dataY, data, vmax = vm, vmin = -vm, shading='auto',  cmap = "seismic" )
#plt.pcolormesh( dataX, dataY, data, shading='auto',  cmap = "seismic" )
plt.pcolormesh( data.T, vmax = vm, vmin = -vm, cmap = "seismic" ) #origin= "lower" )
plt.colorbar( )
plt.axis( "image" )
plt.title( fileName )
#plt.show( )
plt.savefig( varname + ".png" )
#plt.plot( dataY[grid.NZ - 1, :] )

#print( grid )


