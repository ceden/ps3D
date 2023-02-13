
program main
 use main_module
 use timing_module   
 implicit none
 include "mpif.h"
 integer :: ierr,i,j,k,k2,iargc,n
 character (len=80) :: arg
 real*8,external :: hat_ksqr
 
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
 
 !----------------------------------------------------------------------
 ! wavenumber grid
 !----------------------------------------------------------------------

 allocate( kx(nx/2+1), ky(ny), kz(nz) )
 
 do i=1,nx/2+1
   kx(i) =  2*(i-1)*pi/Lx
 enddo
 do j=1,(ny+1)/2
   ky(j) =  2*(j-1)*pi/Ly 
 enddo 
 do j=(ny+1)/2+1,ny
   ky(j) =  -2*(ny-j+1)*pi/Ly
 enddo
 
 if (enable_vertical_boundaries) then
  do k=1,nz
   kz(k) =  (k-1.)*pi/Lz  
  enddo 
   
  phi(:,1) = sqrt(1./Lz)
  psi(:,1) = 0d0
  psib(:,1) = 0d0
  do k2=2,nz
   do k=1,nz
     psi(k,k2) = sqrt(2./Lz)*sin((k2-1)*k*pi/nz)
   enddo  
   k=1
   phi(k,k2) = (psi(k,k2))/(2*sin((k2-1)*pi/(2*nz)) )
   psib(k,k2)= (psi(k,k2))/(2*cos((k2-1)*pi/(2*nz)) )
   do k=2,nz
    phi(k,k2) = (psi(k,k2)-psi(k-1,k2))/(2*sin((k2-1)*pi/(2*nz)) )
    psib(k,k2)= (psi(k,k2)+psi(k-1,k2))/(2*cos((k2-1)*pi/(2*nz)) )
   enddo
    
  enddo
  
 else
  do k=1,(nz+1)/2
   kz(k) =  2*(k-1)*pi/Lz 
  enddo 
  do k=(nz+1)/2+1,nz
   kz(k) =  -2*(nz-k+1)*pi/Lz
  enddo
 endif 
 

 !----------------------------------------------------------------------
 ! wavenumber related stuff
 !----------------------------------------------------------------------
 if (enable_vertical_boundaries) then
   do j=fstart(2),fend(2)
    do i=fstart(1),fend(1)
       ksqr(i,j,1) =  hat_ksqr(kx(i),dx) + hat_ksqr(ky(j),dy) 
    enddo
   enddo
 else
   do k=fstart(3),fend(3)
    do j=fstart(2),fend(2)
     do i=fstart(1),fend(1)
      ksqr(i,j,k) =  hat_ksqr(kx(i),dx) + hat_ksqr(ky(j),dy) + hat_ksqr(kz(k),dz)/dsqr
     enddo
    enddo   
   enddo 
 endif



  if (my_pe==0)  print *,' initialization of diagnostic output  '
  call flush(6)

  if (enable_diag_snap) call init_snap_cdf
  !if (enable_diag_spec) call init_diag_spec()
  !if (enable_diag_snap_chunks) call init_snap_cdf_chunks
  if (enable_diag_balance.or.enable_diag_balance_chunks) call init_diag_balance
  !if (enable_diag_opt_balance) call init_diag_opt_balance
   
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
 
 time = 0.
 itt = 0

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
      !if (enable_diag_snap_par) call diag_snap_par()
      !if (enable_diag_spec)     call diag_spec()
      if (enable_diag_balance.or.enable_diag_balance_chunks)  call diag_balance()
      !if (enable_diag_snap_chunks) call diag_snap_chunks()
      !if (enable_diag_opt_balance) then
      !   call diag_opt_balance()       
      !   call write_diag_opt_balance()
      !endif 
      
   endif   
   call toc('diagnostics')
        
   n = taum2; taum2 = taum1; taum1 = tau; tau = n
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






subroutine p3dfft_initialization
 use main_module
 implicit none
 include "mpif.h"
 integer :: n,ierr
 
 ! Set up work structures for P3DFFT
 if (my_pe==0)  print *,'Setting up P3DFFT '
 call p3dfft_setup ((/n_pes_x,n_pes_y/),nx,ny,nz,MPI_COMM_WORLD,overwrite=.false.)
 if (my_pe==0)  print *,'Done setting up P3DFFT '
 
 
 ! Get dimensions for the original array of real numbers, X-pencils
 call p3dfft_get_dims(istart,iend,isize,1)

! Get dimensions for the R2C-forward-transformed array of complex numbers
!   Z-pencils (depending on how the library was compiled, the first
!   dimension could be either X or Z)
 call p3dfft_get_dims(fstart,fend,fsize,2)


 is_pe = istart(1); ie_pe = iend(1)
 js_pe = istart(2); je_pe = iend(2)
 ks_pe = istart(3); ke_pe = iend(3)
 

 if (my_pe==0)  print *,' '
 if (my_pe==0)  print *,'Layout of physical domain decomposition '
 do n=0,n_pes-1
  call mpi_barrier(MPI_COMM_WORLD, ierr)
  if (my_pe==n) print*,'PE#',my_pe,': i=',istart(1),':',iend(1), &
        ' j=',istart(2),':',iend(2),' k=',istart(3),':',iend(3) 
 enddo
 
 
 call mpi_barrier(MPI_COMM_WORLD, ierr)
 call flush(6) 
 call mpi_barrier(MPI_COMM_WORLD, ierr)
 
 if (my_pe==0)  print *,' '
 if (my_pe==0)  print *,'Layout of spectral domain decomposition '
 do n=0,n_pes-1
  call mpi_barrier(MPI_COMM_WORLD, ierr)
  if (my_pe==n) print*,'PE#',my_pe,': i=',fstart(1),':',fend(1), &
        ' j=',fstart(2),':',fend(2),' k=',fstart(3),':',fend(3) 
 enddo

 
 call mpi_barrier(MPI_COMM_WORLD, ierr)
 call flush(6)
 call mpi_barrier(MPI_COMM_WORLD, ierr)
 
 n_pes_i = 1
 n_pes_j = n_pes_x
 n_pes_k = n_pes_y

 if (my_pe==0)  print *,' n_pes_i =',n_pes_i
 if (my_pe==0)  print *,' n_pes_j =',n_pes_j
 if (my_pe==0)  print *,' n_pes_k =',n_pes_k

 
 call mpi_barrier(MPI_COMM_WORLD, ierr)
 call flush(6)
 call mpi_barrier(MPI_COMM_WORLD, ierr)
 
  my_blk_i = 1
  my_blk_j = mod(my_pe,n_pes_j)+1
  my_blk_k = mod(my_pe/n_pes_j,n_pes_k)+1
   
 
 if (my_pe==0)  print *,' '
 do n=0,n_pes-1
  call mpi_barrier(MPI_COMM_WORLD, ierr)
  if (my_pe==n) print*,'PE#',my_pe,' my_blk_i=',my_blk_i,' my_blk_j=',my_blk_j,' my_blk_k=',my_blk_k
 enddo

end subroutine p3dfft_initialization


