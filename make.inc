
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


ifdef LEVANTE
 P3DFFT = /work/uo0780/p3dfft
 FFTW = /sw/spack-levante/fftw-3.3.10-iynhcf
 cdf = nf-config
 F90 = mpifort
 CDFFLAGS = -I`$(cdf) --includedir`  `$(cdf) --flibs` 
 CDFFLAGS = ${CDFFLAGS} -Wl,-rpath,/sw/spack-levante/netcdf-fortran-4.5.3-k6xq5g/lib  
 CDFFLAGS = ${CDFFLAGS} -Dwith_netcdf_parallel
 MPIFLAGS =
 F90FLAGS =
 P3DFFT_FLAGS = -I${P3DFFT}/include  -L${P3DFFT}/lib -lp3dfft 
 P3DFFT_FLAGS = ${P3DFFT_FLAGS} -Wl,-rpath,${FFTW}/lib -L$(FFTW)/lib -lfftw3
endif




