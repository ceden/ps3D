
include ../make.inc

all_files = config.o main_module.o timing_module.o advection.o restart.o diagnose.o\
            integrate.o mixing.o my_fft.o parallel_mpi.o pressure_solver.o  \
            diag_balance.o diag_balance_io.o eigenvectors.o \
            diag_opt_balance.o diag_opt_balance_io.o


all: main.o  $(all_files)
	$(F90) main.o $(all_files) $(F90FLAGS) $(CDFFLAGS) $(MPIFLAGS) $(P3DFFT_FLAGS) -o model.x

clean :
	 rm -f *.o *.mod 

timing_module.o    : timing_module.f90
	 	$(F90) $(F90FLAGS) -c timing_module.f90
        
main_module.o    : main_module.f90
	 	$(F90) $(F90FLAGS) $(MPIFLAGS) -c main_module.f90
main.o         : main.f90 main_module.o timing_module.o
	 	$(F90) $(F90FLAGS) $(MPIFLAGS) -c main.f90
integrate.o: integrate.f90 main_module.o  timing_module.o
		$(F90) $(F90FLAGS) $(MPIFLAGS) -c integrate.f90
mixing.o: mixing.f90 main_module.o 
		$(F90) $(F90FLAGS) $(MPIFLAGS) -c mixing.f90
advection.o: advection.f90 main_module.o 
		$(F90) $(F90FLAGS) $(MPIFLAGS) -c advection.f90				
parallel_mpi.o: parallel_mpi.f90 main_module.o 
		$(F90) $(F90FLAGS) $(MPIFLAGS) -c parallel_mpi.f90		
pressure_solver.o: pressure_solver.f90 main_module.o 
		$(F90) $(F90FLAGS) $(MPIFLAGS) -c pressure_solver.f90				
restart.o: restart.f90 main_module.o 
		$(F90) $(F90FLAGS) $(MPIFLAGS) -c restart.f90		
        
my_fft.o: my_fft.f90 main_module.o
		$(F90) $(F90FLAGS) $(MPIFLAGS) $(P3DFFT_FLAGS)  -c my_fft.f90       
       
diagnose.o: diagnose.F90 main_module.o
		$(F90) $(F90FLAGS) $(CDFFLAGS) $(MPIFLAGS)  -c diagnose.F90 
diag_balance.o: diag_balance.f90 main_module.o
		$(F90) $(F90FLAGS)  $(MPIFLAGS)  -c diag_balance.f90 
diag_balance_io.o: diag_balance_io.F90 main_module.o
		$(F90) $(F90FLAGS) $(CDFFLAGS) $(MPIFLAGS) -c diag_balance_io.F90 
eigenvectors.o: eigenvectors.f90 
		$(F90) $(F90FLAGS)  $(MPIFLAGS) -c eigenvectors.f90 
diag_opt_balance.o: diag_opt_balance.f90 main_module.o
		$(F90) $(F90FLAGS)  $(MPIFLAGS) -c diag_opt_balance.f90 
diag_opt_balance_io.o: diag_opt_balance_io.f90 main_module.o
		$(F90) $(F90FLAGS) $(CDFFLAGS) $(MPIFLAGS) -c diag_opt_balance_io.f90 

config.o: config.f90 main_module.o diag_balance.o diag_opt_balance.o
		$(F90) $(F90FLAGS) $(MPIFLAGS)  $(CDFFLAGS) -c config.f90

