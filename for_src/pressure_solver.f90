
subroutine pressure_solver
 use main_module
 use timing_module   
 implicit none
 integer :: i,j,k
 real(real_type) :: aloc(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)
 real(real_type) :: a_tri(nz),b_tri(nz),c_tri(nz),fxa
 complex(real_type) :: bloc(fstart(3):fend(3))
 
 call tic('pressure')
 
 !---------------------------------------------------------------------------------
 ! non-hydrostatic pressure forcing =  (u_x + v_y + w_z)
 !---------------------------------------------------------------------------------
 do k=ks_pe,ke_pe
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
    aloc(i,j,k) = (u(i,j,k)-u(i-1,j,k))/dx  +(v(i,j,k)-v(i,j-1,k))/dy &
                 +(w(i,j,k)-w(i,j,k-1))/dz
   enddo
  enddo
 enddo
 aloc = aloc/dt

 !---------------------------------------------------------------------------------
 ! solve for non-hydrostatic pressure  
 !---------------------------------------------------------------------------------

  if (enable_vertical_boundaries) then
  
   call my_fft_2D(aloc,FS)
   fxa = 1d0/dz**2/dsqr
   a_tri(2:nz)   = fxa
   c_tri(1:nz-1) = fxa
   do j=fstart(2),fend(2)
    do i=fstart(1),fend(1)
       
       if (ksqr(i,j,1)>0) then      
        b_tri(1)      = -ksqr(i,j,1) - fxa
        b_tri(2:nz-1) = -ksqr(i,j,1) - 2*fxa 
        b_tri(nz)     = -ksqr(i,j,1) - fxa
        call solve_tridiag(a_tri,b_tri,c_tri,FS(i,j,1:nz),FS(i,j,1:nz),nz)
       else
        bloc = FS(i,j,:)
        FS(i,j,1) = 0.
        FS(i,j,2) = bloc(1)*dz**2*dsqr
        do k=2,nz-1
          FS(i,j,k+1) = 2*FS(i,j,k) - FS(i,j,k-1) + bloc(k)*dz**2*dsqr
        enddo
       endif 
    enddo
   enddo
   call my_ifft_2D(FS,aloc)
   
  else 
  
   call my_fft_3D(aloc,FS)
   do k=fstart(3),fend(3)
    do j=fstart(2),fend(2)
     do i=fstart(1),fend(1)
      if ( ksqr(i,j,k) > 0. ) then
       FS(i,j,k)= -FS(i,j,k)/ksqr(i,j,k)
      else
       FS(i,j,k) = 0.
      endif  
     enddo
    enddo   
   enddo 
   call my_ifft_3D(FS,aloc)
   
  endif ! enable_vertical_boundaries
  
  
  p(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) = aloc
  call border_exchg_3D(p) 
  
  ! with odd number of grid points it is necessary to iterate here in addition
  !call congrad(fc,p)
  !call border_exchg_3D(p) 

  if (enable_vertical_boundaries) then
   if (my_blk_k == 1 ) then
     do k=1-onx,0
       p(:,:,k) =  p(:,:,1)
     enddo
   endif
   if (my_blk_k == n_pes_k ) then
     do k=nz+1,nz+onx
       p(:,:,k) =  p(:,:,nz)
     enddo
   endif
  endif

 !---------------------------------------------------------------------------------
 ! remove non-hydrostatic pressure gradient for final result
 !---------------------------------------------------------------------------------
 do k=ks_pe,ke_pe
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
    u(i,j,k) = u(i,j,k) -      dt*( p(i+1,j,k)-p(i,j,k))/dx  
    v(i,j,k) = v(i,j,k) -      dt*( p(i,j+1,k)-p(i,j,k))/dy 
    w(i,j,k) = w(i,j,k) - dt/dsqr*( p(i,j,k+1)-p(i,j,k))/dz 
   enddo
  enddo
 enddo
 if (enable_vertical_boundaries) then
   if (my_blk_k == 1      ) w(:,:,1-onx:0) = 0d0
   if (my_blk_k == n_pes_k) w(:,:,nz:nz+onx) = 0d0
 endif

 call toc('pressure')

 call tic('boundary')
 call border_exchg_3D(u) 
 call border_exchg_3D(v) 
 call border_exchg_3D(w) 
 call toc('boundary')
 
end subroutine pressure_solver





subroutine solve_tridiag(a,b,c,d,x,n)
      implicit none
!---------------------------------------------------------------------------------
!        a - sub-diagonal (means it is the diagonal below the main diagonal)
!        b - the main diagonal
!        c - sup-diagonal (means it is the diagonal above the main diagonal)
!        d - right part
!        x - the answer
!        n - number of equations
!---------------------------------------------------------------------------------
        integer, parameter :: real_type = KIND(1.d0)
        integer,intent(in) :: n
        real(real_type),dimension(n),intent(in) :: a,b,c
        complex(real_type),dimension(n),intent(in ) :: d
        complex(real_type),dimension(n),intent(out) :: x
        complex(real_type),dimension(n) :: dp
        real(real_type),dimension(n) :: cp
        real(real_type) :: m,fxa
        integer i
 
! initialize c-prime and d-prime
        cp(1) = c(1)/b(1)
        dp(1) = d(1)/b(1)
! solve for vectors c-prime and d-prime
         do i = 2,n
           m = b(i)-cp(i-1)*a(i)
           fxa = 1D0/m
           cp(i) = c(i)*fxa
           dp(i) = (d(i)-dp(i-1)*a(i))*fxa
         enddo
! initialize x
         x(n) = dp(n) 
! solve for x from the vectors c-prime and d-prime
        do i = n-1, 1, -1
          x(i) = dp(i)-cp(i)*x(i+1)  
        end do
end subroutine solve_tridiag






