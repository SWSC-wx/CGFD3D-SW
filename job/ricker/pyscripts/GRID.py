#/usr/bin/env python
'''
Author: Wenqiang Wang @ SUSTech on Sep 11, 2021
15:40
'''
import json
import numpy as np

class GRID:
	def __init__( self, params, HALO = 3 ):
		#resX = 0 	
		#resY = 0 
		#resZ = 0 	
			
		self.PX = params["PX"] 
		self.PY = params["PY"] 
		self.PZ = params["PZ"] 

		self._NX_ = params["NX"] + 2 * HALO 
		self._NY_ = params["NY"] + 2 * HALO 
		self._NZ_ = params["NZ"] + 2 * HALO 

		self._NX = params["NX"] + HALO 
		self._NY = params["NY"] + HALO 
		self._NZ = params["NZ"] + HALO 

		self.NX = params["NX"] 
		self.NY = params["NY"] 
		self.NZ = params["NZ"] 

		#resX = self.NX % self.PX 
		#resY = self.NY % self.PY 
		#resZ = self.NZ % self.PZ

		self.nx = np.zeros( self.PX, dtype = 'int32' )
		self.ny = np.zeros( self.PY, dtype = 'int32' )
		self.nz = np.zeros( self.PZ, dtype = 'int32' )
		
		self.frontNX = np.zeros( self.PX, dtype = 'int32' )
		self.frontNY = np.zeros( self.PY, dtype = 'int32' )
		self.frontNZ = np.zeros( self.PZ, dtype = 'int32' )

		Gather = params["Gather"]
		
		marginX = self.NX // self.PX - params["PMLresX"]
		marginY = self.NY // self.PY - params["PMLresY"]
		marginZ = self.NZ // self.PZ - params["PMLresZ"]

		if self.PX > 2:
			partitionX = (self.NX - 2*marginX) // (self.PX - 2)
		else:
			partitionX = 0
		if self.PY > 2:
			partitionY = (self.NY - 2*marginY) // (self.PY - 2)
		else:
			partitionY = 0
		if self.PZ > 2:
			partitionZ = (self.NZ - 2*marginZ) // (self.PZ - 2)
		else:
			partitionZ = 0

		for mpiX in range( self.PX ):
			if (mpiX == 0 or mpiX == self.PX - 1) :
				self.nx[mpiX] = marginX
			else :
				self.nx[mpiX] = partitionX
			

		for mpiY in range( self.PY ):
			if (mpiY == 0 or mpiY == self.PY - 1) :
				self.ny[mpiY] = marginY
			else :
				self.ny[mpiY] = partitionY
			

		for mpiZ in range( self.PZ ):
			if (mpiZ == 0 or mpiZ == self.PZ - 1) :
				self.nz[mpiZ] = marginZ
			else :
				self.nz[mpiZ] = partitionZ
			

		if self.PX > 2:
			self.resX = self.NX - 2*marginX - (self.PX - 2)*partitionX
		else:
			self.resX = self.NX % marginX
		if self.PY > 2:
			self.resY = self.NY - 2*marginY - (self.PY - 2)*partitionY
		else:
			self.resY = self.NY % marginY
		if self.PZ > 2:
			self.resZ = self.NZ - 2*marginZ - (self.PZ - 2)*partitionZ
		else:
			self.resZ = self.NZ % marginZ

		for mpiX in range( self.PX ):
			if ( mpiX < self.resX ):
				self.nx[mpiX] += 1
				if mpiX != 0 :
					self.frontNX[mpiX] = marginX + 1 + (mpiX - 1)*(partitionX + 1) 
				else:
					self.frontNX[mpiX] = 0
			else:
				if mpiX != 0 :
					self.frontNX[mpiX] = marginX + 1 + (self.resX - 1) * ( partitionX + 1 ) + ( mpiX - self.resX ) * partitionX 
				else:
					self.frontNX[mpiX] = 0


		for mpiY in range( self.PY ):
			if ( mpiY < self.resY ):
				self.ny[mpiY] += 1
				if mpiY != 0 :
					self.frontNY[mpiY] = marginY + 1 + (mpiY - 1)*(partitionY + 1) 
				else:
					self.frontNY[mpiY] = 0
			else:
				if mpiY != 0 :
					self.frontNY[mpiY] = marginY + 1 + (self.resY - 1) * ( partitionY + 1 ) + ( mpiY - self.resY ) * partitionY 
				else:
					self.frontNY[mpiY] = 0 

		for mpiZ in range( self.PZ ):
			if ( mpiZ < self.resZ ):
				self.nz[mpiZ] += 1
				if mpiZ != 0 :
					self.frontNZ[mpiZ] = marginZ + 1 + (mpiZ - 1)*(partitionZ + 1) 
				else:
					self.frontNZ[mpiZ] = 0
			else:
				if mpiZ != 0 :
					self.frontNZ[mpiZ] = marginZ + 1 + (self.resZ - 1) * ( partitionZ + 1 ) + ( mpiZ - self.resZ ) * partitionZ 
				else:
					self.frontNZ[mpiZ] = 0
		'''
		if Gather:
			for mpiX in range( self.PX ):
				if mpiX != 0 :
					self.frontNX[mpiX] = marginX + (mpiX - 1)*(partitionX) 
				else:
					self.frontNX[mpiX] = 0
			
			for mpiY in range( self.PY ):
				if mpiY != 0 :
					self.frontNY[mpiY] = marginY + (mpiY - 1)*(partitionY) 
				else:
					self.frontNY[mpiY] = 0
			
			for mpiZ in range( self.PZ ):
				if mpiZ != 0 :
					self.frontNZ[mpiZ] = marginZ + (mpiZ - 1)*(partitionZ) 
				else:
					self.frontNZ[mpiZ] = 0
		'''
		self._frontNX = self.frontNX + HALO 
		self._frontNY = self.frontNY + HALO 
		self._frontNZ = self.frontNZ + HALO 


		self._nx = self.nx + HALO 
		self._ny = self.ny + HALO 
		self._nz = self.nz + HALO 

		self._nx_ = self.nx + 2 * HALO 
		self._ny_ = self.ny + 2 * HALO 
		self._nz_ = self.nz + 2 * HALO 

		self.originalX = params["centerX"] 
		self.originalY = params["centerY"] 

		self._originalX = self.originalX + HALO 
		self._originalY = self.originalY + HALO 

		self.halo = HALO 

		self.DH = params["DH"] 

def main( ):
	jsonsFile = open( "params.json" )
	params = json.load( jsonsFile )
	grid = GRID( params )
	print( grid )


if __name__ == '__main__':
	main()




