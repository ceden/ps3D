

module main_module
 implicit none

 integer, parameter :: real_type = KIND(1.d0)

!---------------------------------------------------------------------------------
! some fixed parameter
!---------------------------------------------------------------------------------
  real(real_type), parameter :: version = 0.1
  real(real_type) :: pi = 3.14159265358979323846264338327950588 ! will be set below
!---------------------------------------------------------------------------------
! switches for configuration
!---------------------------------------------------------------------------------
  logical :: enable_nonlinear = .true.    
  logical :: enable_vertical_boundaries = .false.
  logical :: enable_AB_3_order = .false.
  
  logical :: enable_diag_balance        = .false.
  logical :: enable_diag_balance_chunks = .false.
  logical :: enable_diag_balance_filter = .false.
  integer :: diag_balance_filter_width  = 3 
  logical :: enable_diag_opt_balance    = .false.
  logical :: enable_diag_opt_time_ave   = .false.
  real(real_type)  :: opt_balance_period         = 1.
  real(real_type)  :: opt_balance_average        = 1.
  integer :: opt_balance_max_Itts       = 100
  integer :: opt_balance_average_times  = 1
  real(real_type)  :: opt_balance_tol            = 1d-9
 
!---------------------------------------------------------------------------------
! model parameter
!---------------------------------------------------------------------------------
  integer:: nx,ny,nz    ! number of grid points   
  real(real_type) :: Lx,Ly,Lz    ! extent of domain in m
  real(real_type) :: dx,dy,dz    ! extent of grid cell in m
  real(real_type) :: dt          ! time step
  real(real_type) :: f0          ! Coriolis freq.
  real(real_type) :: N0          ! Stability freq.
 
  real(real_type) :: Kh = 0d0  ! harmonic diffusivity
  real(real_type) :: Kv = 0d0  ! harmonic diffusivity
  real(real_type) :: Ah = 0d0  ! lateral harmonic viscosity
  real(real_type) :: Av = 0d0  ! vertical harmonic viscosity
  
  real(real_type) :: Khbi = 0d0  ! biharmonic diffusivity
  real(real_type) :: Kvbi = 0d0  ! biharmonic diffusivity
  real(real_type) :: Ahbi = 0d0  ! lateral biharmonic viscosity
  real(real_type) :: Avbi = 0d0  ! vertical biharmonic viscosity
  real(real_type) :: Ro=1., dsqr=1.
!---------------------------------------------------------------------------------
! time stepping parameter
!---------------------------------------------------------------------------------
  integer ::  taum2 = 1, taum1 = 2, tau = 3 !  time step index
  real(real_type)  :: eps_ab                         ! Adam-Bashforth weighting
  real(real_type), parameter :: AB3_a=  23d0/12d0 , AB3_b = -16d0/12d0, AB3_c = 5d0/12d0
  integer :: itt = 0                ! current time step number
  real(real_type) :: time,runlen = 0.        ! current time in s, length of integration in s

!---------------------------------------------------------------------------------
! model fields
!---------------------------------------------------------------------------------
  real(real_type), allocatable :: u(:,:,:),du(:,:,:,:)    ! velocity and tendency
  real(real_type), allocatable :: v(:,:,:),dv(:,:,:,:)    ! velocity and tendency
  real(real_type), allocatable :: w(:,:,:),dw(:,:,:,:)    ! velocity and tendency
  real(real_type), allocatable :: b(:,:,:),db(:,:,:,:)    ! buoyancy and tendency
  real(real_type), allocatable :: p(:,:,:)                ! pressure
  
  real(real_type), allocatable :: flux_east(:,:,:) ! auxilliary fields
  real(real_type), allocatable :: flux_north(:,:,:)
  real(real_type), allocatable :: flux_top(:,:,:)

  real(real_type), allocatable :: kx(:),ky(:),kz(:)  ! wavenumber grid
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
  
  complex(real_type), allocatable :: FS(:,:,:)
  real(real_type), allocatable :: ksqr(:,:,:)
  real(real_type), allocatable :: psi(:,:),phi(:,:),psib(:,:) ! vertical modes

!---------------------------------------------------------------------------------
! diagnostic setup
!---------------------------------------------------------------------------------
 logical :: enable_diag_snap = .false.
 real(real_type) :: tsmonint = 0. , snapint = 0.
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
  allocate( ksqr(fstart(1):fend(1),fstart(2):fend(2),1) ) ;ksqr = 0.
 else
  allocate( ksqr(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) );ksqr = 0.
 endif
 allocate( FS(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); FS = 0.
 allocate( kx(nx/2+1), ky(ny), kz(nz) )
 
end subroutine  allocate_main_module
