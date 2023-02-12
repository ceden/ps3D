
OSX = "" 

#########################################################################
# platform dependent options
#########################################################################

ifdef OSX

 FFTW = ${HOME}/src/0_Third_party
 P3DFFT = ${HOME}/src/0_Third_party/p3dfft-2.7.9
 
 F90 =  mpif90
 CDFFLAGS = -I/usr/local/include -lnetcdff -lnetcdf # -Dwith_netcdf_parallel 
 MPIFLAGS =
 F90FLAGS =  -Wall -O3 -fbounds-check -finit-real=nan -fallow-argument-mismatch
 P3DFFT_FLAGS =  -I${P3DFFT}/include  -L${P3DFFT}/build -lp3dfft -L${FFTW}/lib -lfftw3
endif






