### 安装mpi
 yum install -y hwloc-libs libevent-devel
 tar jxvf openmpi-4.1.0.tar.bz2
 cd openmpi-4.1.0
 mkdir build && cd build
 ../configure --prefix=/usr/local/openmp --enable-shared
 make -j
 make install

tar -zxvf zlib-1.2.8.tar.gz
cd zlib-1.2.8
./configure --prefix=/$home/netcdf_install
make
make check
make install

yum install -y zlib-devel# wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.7/src/hdf5-1.10.7.tar.gz
tar zxvf hdf5-1.10.7.tar.gz
cd hdf5-1.10.7
mkdir build && cd build
CC=mpicc FC=mpifort F77=mpifort \
  ../configure --prefix=/usr/local/hdf5 --enable-fortran --enable-parallel --enable-shared \
  --enable-hl --with-szlib=/usr/local/szip
make -j
make install



yum install -y zlib-devel# wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.7/src/hdf5-1.10.7.tar.gz
tar zxvf hdf5-1.10.7.tar.gz
cd hdf5-1.10.7
mkdir build && cd build
CC=mpicc FC=mpifort F77=mpifort \
  ../configure --prefix=/usr/local/hdf5 --enable-fortran --enable-parallel --enable-shared \
  --enable-hl --with-szlib=/usr/local/szip
make -j
make install


wget https://github.com/Unidata/netcdf-c/archive/refs/tags/v4.7.4.tar.gz -O netcdf-c-4.7.4.tar.gz
tar zxvf netcdf-c-4.7.4.tar.gz
cd netcdf-c-4.7.4
mkdir build && cd build
CFLAGS="-I/usr/local/hdf5/include -I/usr/local/pnetcdf/include" \
  CPPFLAGS="-I/usr/local/hdf5/include -I/usr/local/pnetcdf/include" \
  LDFLAGS="-L/usr/local/hdf5/lib -L/usr/local/pnetcdf/lib" \
  CC=mpicc ../configure --prefix=/usr/local/netcdf \
  --enable-shared --enable-static --enable-pnetcdf --enable-netcdf-4 \
  --enable-largefile --enable-large-file-tests
make -j
make install