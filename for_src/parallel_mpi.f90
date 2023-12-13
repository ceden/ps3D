



subroutine halt_stop(string)
!--------------------------------------------------------------
!     controlled stop, should not be called from python
!--------------------------------------------------------------
      implicit none
      character*(*) :: string
      integer :: ierr,code,my_pe
      include "mpif.h"
      call mpi_comm_rank(MPI_COMM_WORLD,my_pe,ierr)
      print*,' global pe #',my_pe,' : ',string
      print*,' global pe #',my_pe,' aborting '
      code=99
      call MPI_ABORT(mpi_comm_world, code, IERR)
end subroutine halt_stop


subroutine pe0_bcast(a,len)
!--------------------------------------------------------------
!     Broadcast a vector from pe0 to all other pe
!--------------------------------------------------------------
      use main_module   
      implicit none
      integer, intent(in) :: len
      real*8, intent(inout) :: a(len)
      integer :: ierr
      include "mpif.h"
      call mpi_bcast(a,len,mpi_real8,0,MPI_COMM_WORLD,ierr)
end subroutine pe0_bcast




subroutine global_max(x)
!--------------------------------------------------------------
!     Get the max of real x over all PEs in sub domain
!--------------------------------------------------------------
      use main_module   
      implicit none
      real(real_type),intent(inout)    :: x
      real(real_type)    :: x_sym,x_sym2
      integer :: ierr
      include "mpif.h"
      x_sym = x
      call mpi_allreduce(x_sym,x_sym2,1,mpi_real8,MPI_MAX,mpi_comm_world,ierr)
      x = x_sym2
 end subroutine global_max



subroutine global_sum(x)
!--------------------------------------------------------------
!     Do a sum of real x over all PEs in sub domain
!--------------------------------------------------------------
      use main_module   
      implicit none
      real(real_type),intent(inout)    :: x
      real(real_type)    :: x_sym,x_sym2
      integer :: ierr
      include "mpif.h"
      x_sym = x
      call mpi_allreduce(x_sym,x_sym2,1,mpi_real8,MPI_SUM,mpi_comm_world,ierr)
      x = x_sym2
end subroutine global_sum


subroutine border_exchg_3D(a)
!--------------------------------------------------------------
! Exchange overlapping areas of array a in all PEs 
! and also for periodic boundaries
!--------------------------------------------------------------
  use main_module   
  implicit none
  real(real_type), intent(inout)  :: a(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx)
  call border_exchg_in_z(a)
  call border_exchg_in_y(a)
  call border_exchg_in_x(a)
end subroutine border_exchg_3D


subroutine border_exchg_in_x(a)
!--------------------------------------------------------------
! Exchange overlapping areas of array a in all PEs 
! and also for periodic boundaries, only in x direction
!--------------------------------------------------------------
  use main_module   
  implicit none
  real(real_type), intent(inout)  :: a(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx)
  integer  ::  tag=0, ierr,i,len,east,west
  include "mpif.h"
  integer,dimension(MPI_STATUS_SIZE)  :: Status
  integer :: my_comm
  
  my_comm = MPI_COMM_WORLD
  west = my_pe-1
  if (my_blk_i==1      ) west = my_pe+ n_pes_i-1
  east = my_pe+1
  if (my_blk_i==n_pes_i) east = my_pe- (n_pes_i-1)

  if ( n_pes_i > 1) then
     len=(je_pe-js_pe+1+2*onx)*(ke_pe-ks_pe+1+2*onx)
     do i=1,onx
       if (my_blk_i /=1 )       call mpi_send(a(is_pe+i-1,:,:),len,mpi_real8,west,tag,my_comm,ierr)
       if (my_blk_i /= n_pes_i) call mpi_recv(a(ie_pe+i  ,:,:),len,mpi_real8,east,tag,my_comm,status,ierr)
       if (my_blk_i ==1 )       call mpi_send(a(is_pe+i-1,:,:),len,mpi_real8,west,tag,my_comm,ierr)
       if (my_blk_i == n_pes_i) call mpi_recv(a(ie_pe+i  ,:,:),len,mpi_real8,east,tag,my_comm,status,ierr)
     enddo
     do i=1,onx
       if (my_blk_i /= n_pes_i) call mpi_send(a(ie_pe-i+1,:,:),len,mpi_real8,east,tag,my_comm,ierr)
       if (my_blk_i /=1 )       call mpi_recv(a(is_pe-i  ,:,:),len,mpi_real8,west,tag,my_comm,status,ierr)
       if (my_blk_i == n_pes_i) call mpi_send(a(ie_pe-i+1,:,:),len,mpi_real8,east,tag,my_comm,ierr)
       if (my_blk_i ==1 )       call mpi_recv(a(is_pe-i  ,:,:),len,mpi_real8,west,tag,my_comm,status,ierr)
     enddo
  else
    do i=1,onx
     a(nx+i,:,:) = a(i     ,:,:)
     a(1-i ,:,:) = a(nx-i+1,:,:) 
    enddo
  endif
end subroutine border_exchg_in_x



subroutine border_exchg_in_y(a)
!--------------------------------------------------------------
! Exchange overlapping areas of array a in all PEs 
! and also for periodic boundaries, only in y direction
!--------------------------------------------------------------
  use main_module   
  implicit none
  real(real_type), intent(inout)  :: a(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx)
  integer  ::  tag=0, ierr,j,len,north,south
  include "mpif.h"
  integer,dimension(MPI_STATUS_SIZE)  :: Status
  integer :: my_comm
  
  my_comm = MPI_COMM_WORLD
  south = my_pe-n_pes_i 
  north = my_pe+n_pes_i
  if (my_blk_j==1)       south = my_pe+ n_pes_i*(n_pes_j-1)
  if (my_blk_j==n_pes_j) north = my_pe- n_pes_i*(n_pes_j-1)

  if ( n_pes_j > 1) then
     len=(ie_pe-is_pe+1+2*onx)*(ke_pe-ks_pe+1+2*onx)
     do j=1,onx
       if (my_blk_j /=1 )       call mpi_send(a(:,js_pe+j-1,:),len,mpi_real8,south,tag,my_comm,ierr)
       if (my_blk_j /= n_pes_j) call mpi_recv(a(:,je_pe+j  ,:),len,mpi_real8,north,tag,my_comm,status,ierr)
       if (my_blk_j ==1 )       call mpi_send(a(:,js_pe+j-1,:),len,mpi_real8,south,tag,my_comm,ierr)
       if (my_blk_j == n_pes_j) call mpi_recv(a(:,je_pe+j  ,:),len,mpi_real8,north,tag,my_comm,status,ierr)
     enddo
     do j=1,onx
       if (my_blk_j /= n_pes_j) call mpi_send(a(:,je_pe-j+1,:),len,mpi_real8,north,tag,my_comm,ierr)
       if (my_blk_j /=1 )       call mpi_recv(a(:,js_pe-j  ,:),len,mpi_real8,south,tag,my_comm,status,ierr)
       if (my_blk_j == n_pes_j) call mpi_send(a(:,je_pe-j+1,:),len,mpi_real8,north,tag,my_comm,ierr)
       if (my_blk_j ==1 )       call mpi_recv(a(:,js_pe-j  ,:),len,mpi_real8,south,tag,my_comm,status,ierr)
     enddo
  else
     do j=1,onx
      a(:,ny+j,:) = a(:,j  ,:)
      a(:,1-j ,:) = a(:,ny-j+1,:) 
     enddo
  endif
end subroutine border_exchg_in_y





subroutine border_exchg_in_z(a)
!--------------------------------------------------------------
! Exchange overlapping areas of array a in all PEs 
! and also for periodic boundaries, only in vertical direction
!--------------------------------------------------------------
  use main_module   
  implicit none
  real(real_type), intent(inout)  :: a(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx)
  integer  ::  tag=0, ierr,k,len,top,bottom
  include "mpif.h"
  integer,dimension(MPI_STATUS_SIZE)  :: Status
  integer :: my_comm
  
  my_comm = MPI_COMM_WORLD
  bottom = my_pe - n_pes_i*n_pes_j
  top    = my_pe + n_pes_i*n_pes_j
  if (my_blk_k==1)       bottom = my_pe+ n_pes_j*n_pes_i*(n_pes_k-1)
  if (my_blk_k==n_pes_k) top    = my_pe- n_pes_j*n_pes_i*(n_pes_k-1)

  if (n_pes_k > 1) then
     len=(ie_pe-is_pe+1+2*onx)*(je_pe-js_pe+1+2*onx)
     do k=1,onx
       if (my_blk_k /=1 )       call mpi_send(a(:,:,ks_pe+k-1),len,mpi_real8,bottom,tag,my_comm,ierr)
       if (my_blk_k /= n_pes_k) call mpi_recv(a(:,:,ke_pe+k  ),len,mpi_real8,top,tag,my_comm,status,ierr)
       if (my_blk_k ==1 )       call mpi_send(a(:,:,ks_pe+k-1),len,mpi_real8,bottom,tag,my_comm,ierr)
       if (my_blk_k == n_pes_k) call mpi_recv(a(:,:,ke_pe+k  ),len,mpi_real8,top,tag,my_comm,status,ierr)
     enddo
     do k=1,onx
       if (my_blk_k /= n_pes_k) call mpi_send(a(:,:,ke_pe-k+1),len,mpi_real8,top,tag,my_comm,ierr)
       if (my_blk_k /=1 )       call mpi_recv(a(:,:,ks_pe-k  ),len,mpi_real8,bottom,tag,my_comm,status,ierr)
       if (my_blk_k == n_pes_k) call mpi_send(a(:,:,ke_pe-k+1),len,mpi_real8,top,tag,my_comm,ierr)
       if (my_blk_k ==1 )       call mpi_recv(a(:,:,ks_pe-k  ),len,mpi_real8,bottom,tag,my_comm,status,ierr)
     enddo
  else 
   do k=1,onx 
     a(:,:,nz+k) = a(:,:,k)
     a(:,:,1-k)  = a(:,:,nz-k+1) 
   enddo
  endif
end subroutine border_exchg_in_z




subroutine cumsum_in_z(a)
!--------------------------------------------------------------
!  cumulative sum of 3D array a over vertical direction 
!--------------------------------------------------------------
  use main_module   
  implicit none
  real(real_type),intent(inout) :: a(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx)
  integer :: n,k,k1
  
  call border_exchg_in_z(a)
  do n=1,n_pes_k
    if (my_blk_k==n) then
     if (my_blk_k==1) then
       k1 = 2
     else
       k1=ks_pe
     endif
     do k=k1,ke_pe
       a(:,:,k) = a(:,:,k-1) + a(:,:,k)
     enddo
    endif
    call border_exchg_in_z(a)
  enddo
end subroutine cumsum_in_z




subroutine cumsum_in_y(a)
!--------------------------------------------------------------
!  cumulative sum of 3D array a over y direction 
!--------------------------------------------------------------
  use main_module   
  implicit none
  real(real_type),intent(inout) :: a(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx)
  integer :: n,j,j1
  
  call border_exchg_in_y(a)
  do n=1,n_pes_j
    if (my_blk_j==n) then
     if (my_blk_j==1) then
       j1 = 2
     else
       j1=js_pe
     endif
     do j=j1,je_pe
       a(:,j,:) = a(:,j-1,:) + a(:,j,:)
     enddo
    endif
    call border_exchg_in_y(a)
  enddo  
end subroutine cumsum_in_y






