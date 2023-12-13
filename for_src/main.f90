
program main
 use main_module
 use timing_module   
 implicit none
 include "mpif.h"
 integer :: ierr,iargc,n
 character (len=80) :: arg
 
 !----------------------------------------------------------------------
 ! initialize parallelisation with MPI
 !----------------------------------------------------------------------
 call mpi_init(ierr)
 call MPI_Comm_rank(MPI_COMM_WORLD, my_pe, ierr)
 call MPI_Comm_size(MPI_COMM_WORLD, n_pes, ierr)
 if (my_pe==0)  print *,' '
 if (my_pe==0)  print *,' '
 if (my_pe==0)  print *,' '
 
 if (n_pes>1) then
   if (iargc() < 2)  then
      if (my_pe==0) print*,'ERROR: not enough command line input'
      call MPI_ABORT(mpi_comm_world, 99, IERR)
   endif
   call getarg(1,arg); read(arg,*) n_pes_x
   call getarg(2,arg); read(arg,*) n_pes_y
 else 
  n_pes_x = 1; n_pes_y =1
 endif
 
 if (n_pes_x*n_pes_y /= n_pes) then
   if (my_pe==0) print*,' ERROR: n_pes_x times n_pes_y not equal number of PEs'
   call MPI_ABORT(MPI_COMM_WORLD, 99, ierr)
 endif
 if (my_pe==0)  print '(a,i4,a,i4)','Using processor grid ',n_pes_x,' x ',n_pes_y

 !----------------------------------------------------------------------
 ! setup the model
 !----------------------------------------------------------------------
 call tic('setup')
 if (my_pe==0)  print *,' setting parameters '
 call flush(6)
 call set_parameter()
 
 if (enable_vertical_boundaries) then
  if ( (mod(nx,2)/=0.and.nx>1).or.(mod(ny,2)/=0.and.ny>1) ) then
   if (my_pe==0) print*,'ERROR: nx,ny need to be even or one'
   call MPI_ABORT(MPI_COMM_WORLD, 99, ierr)
  endif
 else
  if ( (mod(nx,2)/=0.and.nx>1).or.(mod(ny,2)/=0.and.ny>1).or.(mod(nz,2)/=0.and.nz>1) ) then
   if (my_pe==0) print*,'ERROR: nx,ny,nz need to be even or one'
   call MPI_ABORT(MPI_COMM_WORLD, 99, ierr)
  endif 
 endif
 
 !----------------------------------------------------------------------
 ! P3DFFT initialization 
 !----------------------------------------------------------------------
 if (my_pe==0)  print *,' P3DFFT initializations '
 call flush(6)
 call p3dfft_initialization 


 if (enable_vertical_boundaries) then
   if (fstart(3) /= 1 .or. fend(3) /=nz )  then
      if (my_pe==0) print*,' ERROR: n_pes_k in spectral space >1'
      call MPI_ABORT(mpi_comm_world, 99, IERR)
   endif   
 endif

 if (my_pe==0)  print *,' allocate '
 call flush(6)
 call allocate_main_module()
 
 call wavenumber_initialisation

 if (my_pe==0)  print *,' initialization of diagnostic output  '
 call flush(6)

 if (enable_diag_snap) call init_snap_cdf
 if (enable_diag_balance.or.enable_diag_balance_chunks) call init_diag_balance
 if (enable_diag_opt_balance) call init_diag_opt_balance
   
 if (my_pe==0) then
     print*,' '
     print'(a,i4,a,i4,a,i4)',' nx x ny x nz = ',nx,' x ',ny,' x ',nz
     print'(a,f8.2,a,f8.2,a,f8.2,a)',' Lx x Ly x Lz = ',Lx,' x ',Ly,' x ',Lz,' m'
     print'(a,f8.2,a)',' Ro      = ',Ro
     print'(a,f8.2,a)',' delta   = ',sqrt(dsqr)
     print'(a,f8.2,a)',' Delta t = ',dt,' s'
     print'(a,f8.2,a)',' Delta x = ',dx,' m'
     print'(a,f8.2,a)',' Delta y = ',dy,' m'
     print'(a,f8.2,a)',' Delta z = ',dz,' m'
     print'(a,e12.6,a)',' Ah      = ',Ah,' m**2/s'
     print'(a,e12.6,a)',' Kh      = ',Kh,' m**2/s'
     print'(a,e12.6,a)',' Av      = ',Av,' m**2/s'
     print'(a,e12.6,a)',' Kv      = ',Kv,' m**2/s'
     print'(a,e12.6,a)',' Ahbi    = ',Ahbi,' m**4/s'
     print'(a,e12.6,a)',' Khbi    = ',Khbi,' m**4/s'
     print'(a,e12.6,a)',' Avbi    = ',Avbi,' m**4/s'
     print'(a,e12.6,a)',' Kvbi    = ',Kvbi,' m**4/s'
     print'(a,e12.6,a)',' f0      = ',f0,' 1/s'
     print'(a,e12.6,a)',' N0      = ',N0,' 1/s'
     print'(a,f8.2,a)',' runlen  = ',runlen,' s'
     print'(a,e12.6)' ,' eps_ab  = ',eps_ab
     print'(a,L)'     ,' enable_vertical_boundaries  = ',enable_vertical_boundaries
     print'(a,L)'     ,' enable_AB_3_order           = ',enable_AB_3_order
     print*,' '
 endif
 
 time = 0.; itt = 0
 taum2 = 1; taum1 = 2; tau = 3

 if (my_pe==0)  print *,' setting initial conditions '
 call flush(6)
 call set_initial_conditions
 
 if (enable_vertical_boundaries) then
   if (my_blk_k == 1      ) w(:,:,1-onx:0) = 0d0
   if (my_blk_k == n_pes_k) w(:,:,nz:nz+onx) = 0d0
 endif
 if (my_pe==0)  print *,' done setting initial conditions '
 
 call toc('setup')

 if (my_pe==0)  print *,' starting main loop '
 call flush(6)
  
 !----------------------------------------------------------------------
 ! reading restart
 !----------------------------------------------------------------------
 call read_restart 
 !----------------------------------------------------------------------
 ! start main integration loop
 !----------------------------------------------------------------------

 call tic('loop')   
 do while (time < runlen) 

   call integrate
      
   call tic('diagnostics')     
   if ( mod(itt,max(1,int(tsmonint/dt)))  == 0 .or. itt == 0)  call diagnose  
   if ( mod(itt,max(1,int( snapint/dt)))  == 0 .or. itt == 0)  then      
      if (enable_diag_snap)     call diag_snap()      
      if (enable_diag_balance.or.enable_diag_balance_chunks)  call diag_balance()
      if (enable_diag_opt_balance) then
         call diag_opt_balance()       
         call write_diag_opt_balance()
      endif 
      
   endif   
   call toc('diagnostics')
   
   tau   = mod(tau,3)+1
   taum1 = mod(taum1,3)+1
   taum2 = mod(taum2,3)+1
   time = time + dt       
   itt = itt + 1 
   
 enddo
 call toc('loop') 


 call write_restart
 
!--------------------------------------------------------------
!     show timing results here
!--------------------------------------------------------------
 !do n = 0,n_pes
  n=0
     call mpi_barrier(MPI_COMM_WORLD, ierr)
     if (my_pe == n) then
        print'(/,a,i4)','Timing summary for PE #',my_pe 
        print'(a,f15.2,a)',' costs for measuring      = ',timing_secs('tictoc'),' s'
        print'(a,f15.2,a)',' setup time summary       = ',timing_secs('setup'),' s'
        print'(a,f15.2,a)',' loop                     = ',timing_secs('loop'),' s'
        print'(a,f15.2,a)',' integrate                = ',timing_secs('integrate'),' s'
        print'(a,f15.2,a)',' pressure                 = ',timing_secs('pressure'),' s'
        print'(a,f15.2,a)',' boundary exchange        = ',timing_secs('boundary'),' s'
        print'(a,f15.2,a)',' diagnostics              = ',timing_secs('diagnostics'),' s'
       
       endif
  !enddo
 
 
 if (my_pe==0) then
  call get_free_iounit(n)
  open(n,file='ritt',form='formatted',status='unknown')
  write(n,*) itt
  close(n)
 endif
 
 !----------------------------------------------------------------------
 ! cancel parallelisation and quit
 !----------------------------------------------------------------------

 if (my_pe==0) print'(/a/)','cancelling MPI service'
 call mpi_finalize(ierr)


end program main





