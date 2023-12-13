
! two balanced jets with small perturbation get unstable
! testing optimal balance with time averaging

module config_module
 ! use this module only locally in this file
 implicit none
 integer,parameter :: fac=1
 real(real_type),parameter :: Ri  = 2000    ! Richardson number 
 real(real_type),parameter :: delta = 0.02 ! aspect ratio
 real(real_type),parameter :: H0   = 1000.0       ! total depth
 real(real_type) :: U0_loc,M0,kmax,jet_scale
end module config_module



subroutine set_parameter
 use main_module 
 use config_module
 implicit none
 
 nx=30*fac; ny=30*2*fac; nz = 10*fac
 
 Ro = 1.; dsqr = 1.
 f0 = 1e-4
 N0 = f0/delta
 U0_loc = sqrt(1./Ri)*N0*H0   
 M0 = sqrt(f0*U0_loc/H0)  
 kmax   = 1./ ( sqrt(2./5.)*sqrt(1+Ri)*M0**2/f0**2*H0 )  
 Lx = 4*2/kmax  *pi
 Ly = 4*2/kmax  *pi * 2
 Lz = H0
 jet_scale = 0.1*Ly/2
 
 !enable_nonlinear = .false.
 
 
 dx=Lx/nx;dy=Ly/ny;dz=Lz/nz

 eps_AB = 0.01
 dt    = 200./fac

 snapint = 86400.*5
 tsmonint = snapint
 
 runlen =  1e12
 
 enable_vertical_boundaries = .true.
 enable_diag_snap = .true.
 enable_AB_3_order = .true.

 !enable_diag_balance        = .true.

 enable_diag_opt_balance = .true.
 opt_balance_period         = 10./f0
 opt_balance_max_Itts       = 3
 opt_balance_tol            = 1d-12
 
 enable_diag_opt_time_ave = .true.
 opt_balance_average = 2*86400! 9./f0
 opt_balance_average_times = 5
 
end subroutine set_parameter 




subroutine set_initial_conditions
 use main_module  
 use config_module
 use module_diag_opt_balance
 use module_diag_balance
 implicit none
 integer :: i,j,k
 real(real_type) :: xt(nx),yt(ny),zt(nz)
 
 do i=1,nx
   xt(i)=(i-0.5)*dx
 enddo 
 do i=1,ny
   yt(i)=(i-0.5)*dy
 enddo 
 do i=1,nz
   zt(i)=(i-0.5)*dz
 enddo 
 
 do k=ks_pe,ke_pe
   do j=js_pe,je_pe
    do i=is_pe,ie_pe  
     u(i,j,k) = U0_loc*( exp(-(yt(j)-1*Ly/4)**2/jet_scale**2) &
                    -exp(-(yt(j)-3*Ly/4)**2/jet_scale**2 ) &
              + 0.05*sin(xt(i)/Lx*4*pi)*sin(yt(j)/Ly*2*pi)  ) *cos(pi*zt(k)/H0)      
   enddo
  enddo
 enddo  
 call border_exchg_3D(u)

 do k=ks_pe,ke_pe
   b(:,:,k) = -dy*(u(:,:,k+1)-u(:,:,k-1))/(2*dz)*f0
 enddo
 if (my_blk_k==1) b(:,:,1) = b(:,:,2) 
 if (my_blk_k==n_pes_k) b(:,:,nz) = b(:,:,nz-1)
 
 call cumsum_in_y(b)
 call border_exchg_3D(b)

 if (enable_diag_balance .and. .true.) then
  call diag_balance_zero_order
  call diag_balance_first_order
  call diag_balance_second_order 
  u(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) = u0(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) + & 
                                           u1(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) + &
                                           u2(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)
  v(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) = v0(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) + & 
                                           v1(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) + & 
                                           v2(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)
  w(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) = w0(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) + & 
                                           w1(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) + & 
                                           w2(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)
  b(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) = b0(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) + & 
                                           b1(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) + &
                                           b2(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)

  call border_exchg_3D(u)
  call border_exchg_3D(v)
  call border_exchg_3D(w)
  call border_exchg_3D(b)
 endif

 if (enable_diag_opt_balance.and. .true. ) then
 
  call diag_opt_balance()
  !call write_diag_opt_balance
  u = u_bal;  v = v_bal
  w = w_bal;  b = b_bal    
  
 endif
  
end subroutine set_initial_conditions




subroutine set_forcing 
end subroutine set_forcing
