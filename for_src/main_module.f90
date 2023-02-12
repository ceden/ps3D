

module main_module
 use p3DFFT
 implicit none

!---------------------------------------------------------------------------------
! some fixed parameter
!---------------------------------------------------------------------------------
  real*8, parameter :: version = 0.1
  real*8 :: pi = 3.14159265358979323846264338327950588 ! will be set below
!---------------------------------------------------------------------------------
! switches for configuration
!---------------------------------------------------------------------------------
  logical :: enable_nonlinear = .true.    
  logical :: enable_vertical_boundaries = .false.
  logical :: enable_AB_3_order = .false.
!---------------------------------------------------------------------------------
! model parameter
!---------------------------------------------------------------------------------
  integer:: nx,ny,nz    ! number of grid points   
  real*8 :: Lx,Ly,Lz    ! extent of domain in m
  real*8 :: dx,dy,dz    ! extent of grid cell in m
  real*8 :: dt          ! time step
  real*8 :: f0          ! Coriolis freq.
  real*8 :: N0          ! Stability freq.
 
  real*8 :: Kh = 0d0  ! harmonic diffusivity
  real*8 :: Kv = 0d0  ! harmonic diffusivity
  real*8 :: Ah = 0d0  ! lateral harmonic viscosity
  real*8 :: Av = 0d0  ! vertical harmonic viscosity
  
  real*8 :: Khbi = 0d0  ! biharmonic diffusivity
  real*8 :: Kvbi = 0d0  ! biharmonic diffusivity
  real*8 :: Ahbi = 0d0  ! lateral biharmonic viscosity
  real*8 :: Avbi = 0d0  ! vertical biharmonic viscosity
  real*8 :: Ro=1., dsqr=1.
!---------------------------------------------------------------------------------
! time stepping parameter
!---------------------------------------------------------------------------------
  integer ::  taum2 = 1, taum1 = 2, tau = 3 !  time step index
  real*8  :: eps_ab                         ! Adam-Bashforth weighting
  real*8, parameter :: AB3_a=  23d0/12d0 , AB3_b = -16d0/12d0, AB3_c = 5d0/12d0
  integer :: itt = 0                ! current time step number
  real*8 :: time,runlen = 0.        ! current time in s, length of integration in s

!---------------------------------------------------------------------------------
! model fields
!---------------------------------------------------------------------------------
  real*8, allocatable :: u(:,:,:),du(:,:,:,:)    ! velocity and tendency
  real*8, allocatable :: v(:,:,:),dv(:,:,:,:)    ! velocity and tendency
  real*8, allocatable :: w(:,:,:),dw(:,:,:,:)    ! velocity and tendency
  real*8, allocatable :: b(:,:,:),db(:,:,:,:)    ! buoyancy and tendency
  real*8, allocatable :: p(:,:,:)                ! pressure
  
  real*8, allocatable :: flux_east(:,:,:) ! auxilliary fields
  real*8, allocatable :: flux_north(:,:,:)
  real*8, allocatable :: flux_top(:,:,:)

  real*8, allocatable :: kx(:),ky(:),kz(:)  ! wavenumber grid
!---------------------------------------------------------------------------------
!     Parallel domain setup
!---------------------------------------------------------------------------------
  integer :: n_pes_x,n_pes_y
  integer :: n_pes     ! total number of processors
  integer :: my_pe     ! index of this processor from 0 to n_pes-1
  integer :: n_pes_i   ! total number of processors in x direction
  integer :: n_pes_j   ! total number of processors in y direction
  integer :: n_pes_k
  integer :: my_blk_i  ! index of this processor in x direction from 1 to n_pes_i
  integer :: my_blk_j  ! index of this processor in y direction from 1 to n_pes_j
  integer :: my_blk_k
  integer :: is_pe     ! start index of grid points in x direction of this processor
  integer :: ie_pe     ! end index of grid points in x direction of this processor
  integer :: js_pe     ! start index of grid points in y direction of this processor
  integer :: je_pe     ! end index of grid points in y direction of this processor
  integer :: ks_pe
  integer :: ke_pe
  integer :: onx=2     ! number of overlapping points in x and y direction
 
 
  integer :: istart(3),iend(3),isize(3) ! physical grid point index of this processor 
  integer :: fstart(3),fend(3),fsize(3) ! spectral grid point index of this processor
  
  complex(p3dfft_type), allocatable :: FS(:,:,:)
  real*8, allocatable :: ksqr(:,:,:)
  real*8, allocatable :: psi(:,:),phi(:,:),psib(:,:) ! vertical modes

!---------------------------------------------------------------------------------
! diagnostic setup
!---------------------------------------------------------------------------------
 logical :: enable_diag_snap = .false.
 real*8 :: tsmonint = 0. , snapint = 0.
end module main_module



subroutine  allocate_main_module
 use main_module
 implicit none
 
 pi = acos(0d0)*2d0
 
 allocate(u(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) )
 allocate(v(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) )
 allocate(w(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) )
 allocate(b(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) )
 allocate(p(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) )
 u=0;v=0;w=0;b=0;p=0
 allocate(du(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx,3) )
 allocate(dv(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx,3) )
 allocate(dw(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx,3) )
 allocate(db(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx,3) )
 db=0;du=0;dv=0;dw=0
 allocate( flux_east (is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) )
 allocate( flux_north(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) )
 allocate( flux_top  (is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) )
 flux_east=0;flux_north=0;flux_top=0
 
 if (enable_vertical_boundaries) then
  allocate(psi(nz,nz),phi(nz,nz),psib(nz,nz) )
  psi=0;phi=0;psib=0
  allocate( ksqr(fstart(1):fend(1),fstart(2):fend(2),1) )
  allocate( FS(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) )
  ksqr = 0.; FS = 0.
 else
  allocate( ksqr(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) )
  allocate( FS(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) )
  ksqr = 0.; FS = 0.
 endif

end subroutine  allocate_main_module
