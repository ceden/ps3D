



subroutine p3dfft_initialization
 use main_module
 use p3Dfft
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




subroutine wavenumber_initialisation
 use main_module
 implicit none
 integer :: i,j,k,k2
 real*8,external :: hat_ksqr
 
 !----------------------------------------------------------------------
 ! wavenumber grid
 !----------------------------------------------------------------------
 
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
end subroutine wavenumber_initialisation


subroutine my_ifft_2D(QS,QP)
 use main_module
 use p3Dfft
 implicit none
 real(real_type)    :: QP(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3))
 complex(real_type) :: QS(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3))
 call p3dfft_btran_c2r(QS,QP,'nff')
 QP=QP /(Lx*Ly) 
end subroutine my_ifft_2D


subroutine my_fft_2D(QP,QS)
 use main_module
 use p3Dfft
 implicit none
 real(real_type)    :: QP(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3))
 complex(real_type) :: QS(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3))
 call p3dfft_ftran_r2c(QP,QS,'ffn')
 QS = QS*dx*dy
end subroutine my_fft_2D





subroutine my_ifft_3D(QS,QP)
 use main_module
 use p3Dfft
 implicit none
 real(real_type)    :: QP(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3))
 complex(real_type) :: QS(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3))
 call p3dfft_btran_c2r(QS,QP,'fft')
 QP=QP /(Lx*Ly*Lz) 
end subroutine my_ifft_3D




subroutine my_fft_3D(QP,QS)
 use main_module
 use p3Dfft
 implicit none
 real(real_type)    :: QP(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3))
 complex(real_type) :: QS(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3))
 call p3dfft_ftran_r2c(QP,QS,'fft')
 QS = QS*dx*dy*dz 
end subroutine my_fft_3D





subroutine my_ifft_vert_mode(QS,QP,mode)
 use main_module
 use p3Dfft
 implicit none
 real(real_type)    :: QP(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3))
 complex(real_type) :: QS(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3))
 real(real_type)    :: mode(nz,nz)
 complex(real_type) :: Qloc(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) 
 integer :: k1,k2
   
 Qloc = 0
 do k1=1,nz ! mode
    do k2=1,nz ! z 
       Qloc(:,:,k2) = Qloc(:,:,k2) + QS(:,:,k1)*mode(k2,k1)
    enddo     
 enddo    
 call p3dfft_btran_c2r(Qloc,QP,'nff')
 QP=QP /(Lx*Ly) 
end subroutine my_ifft_vert_mode




subroutine my_fft_vert_mode(QP,QS,mode)
 use main_module
 use p3Dfft
 implicit none
 real(real_type)    :: QP(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3))
 complex(real_type) :: QS(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3))
 real(real_type)    :: mode(nz,nz)
 complex(real_type) :: Qloc(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) 
 integer :: k1,k2
 
 call p3dfft_ftran_r2c(QP,Qloc,'ffn')
 Qloc = Qloc*dx*dy
 QS = 0
 do k1=1,nz ! mode
    do k2=1,nz ! z 
       QS(:,:,k1) = QS(:,:,k1) + Qloc(:,:,k2)*mode(k2,k1)*dz
    enddo     
 enddo  
end subroutine my_fft_vert_mode

