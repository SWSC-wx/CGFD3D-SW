/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:macro.h
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-09-04
*   Discription:
*
================================================================*/
#ifndef __MACRO__
#define __MACRO__

#define HALO 3
#define Cv 1.0e-1
#define Cs 1.0e-7
#define PI 3.1415926535898
#define RADIAN2DEGREE ( 180.0 / PI )
#define DEGREE2RADIAN ( PI / 180.0 )



//#define X_FAST
#ifdef X_FAST

#define INDEX( i, j, k ) ( i + ( j ) * ( _nx_ ) + ( k ) * ( _nx_ ) * (  _ny_ ) )
#define INDEX3D( i, j, k, nx, ny, nz ) ( i + ( j ) * ( nx ) + ( k ) * ( nx ) * (  ny ) )
#define INDEX2D( i, j, nx, ny, nz ) ( i + ( j ) * ( nx )) //3D array with only one slice
#define Index2D( i, j, nx, ny ) ( i + ( j ) * ( nx ) )  //2D array


#define FOR_LOOP2D( i, j, startI, endI, startJ, endJ ) \
for ( j = startJ; j < endJ; j ++ ) {    \
for ( i = startI; i < endI; i ++ ) {    \


#define END_LOOP2D( ) }}

#define FOR_LOOP3D( i, j, k, startI, endI, startJ, endJ, startK, endK ) \
for ( k = startK; k < endK; k ++ ) {     \
for ( j = startJ; j < endJ; j ++ ) {     \
for ( i = startI; i < endI; i ++ ) {     \

#define END_LOOP3D( ) }}}
#define LON_FAST

#else

#define CALCULATE2D( i, j, startI, endI, startJ, endJ ) \
for ( j = startJ; j < endJ; j ++ ) {    \
for ( i = startI; i < endI; i ++ ) {    

#define END_CALCULATE2D( ) }}

#define CALCULATE3D( i, j, k, startI, endI, startJ, endJ, startK, endK ) \
for ( k = startK; k < endK; k ++ ) {     \
for ( j = startJ; j < endJ; j ++ ) {     \
for ( i = startI; i < endI; i ++ ) {     

#define END_CALCULATE3D( ) }}}

#define CALCULATE1D( i, startI, endI )         \
    for ( i = ( startI ); i < endI; i ++ ) {
#define END_CALCULATE1D( ) } 

#define INDEX( i, j, k ) ( k + ( j ) * ( _nz_ ) + ( i ) * ( _nz_ ) * (  _ny_ ) )
#define INDEX3D( i, j, k, nx, ny, nz ) ( k + ( j ) * ( nz ) + ( i ) * ( nz ) * (  ny ) )
#define INDEX2D( i, j, nx, ny, nz ) ( j + ( i ) * ( ny )) //3D array with only one slice
#define Index2D( i, j, nx, ny ) ( j + ( i ) * ( ny ) )  //2D array


#define FOR_LOOP2D( i, j, startI, endI, startJ, endJ ) \
for ( i = startI; i < endI; i ++ ) {    \
for ( j = startJ; j < endJ; j ++ ) {    \


#define END_LOOP2D( ) }}

#define FOR_LOOP3D( i, j, k, startI, endI, startJ, endJ, startK, endK ) \
for ( i = startI; i < endI; i ++ ) {     \
for ( j = startJ; j < endJ; j ++ ) {     \
for ( k = startK; k < endK; k ++ ) {     \

#define END_LOOP3D( ) }}}






#endif

#endif //__MACRO__
