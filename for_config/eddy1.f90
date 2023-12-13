
! idealized eddy

module config_module
 implicit none
 integer :: fac = 1
 real*8 :: Rad,cspeed,mode
end module config_module

subroutine set_parameter
 use main_module 
 use config_module
 implicit none
 
 nx=120*fac; ny=120*fac; nz = 10*fac
 Ro = 1.; dsqr = 1. 
 f0 = 1e-4
 N0 = f0*30
      
 Lz = 2000.;dz=Lz/nz
 
 mode = 1.
 cspeed = N0*Lz/(mode*pi) 
 Rad = 2*cspeed/f0
 Lx = 20*cspeed/f0
 Ly = 20*cspeed/f0
 
 dx=Lx/nx;dy=Ly/ny

 dt    = 400./fac
 runlen =  1e12
 enable_AB_3_order = .true.
 enable_vertical_boundaries = .true.
 
 enable_diag_snap = .true.
 snapint = 86400./2.
 tsmonint = snapint/10.
 
 Ahbi =  dx**4/(30*86400.) ! A = dx**4/T   
 Khbi = Ahbi
  
 enable_diag_balance        = .true.
end subroutine set_parameter 




subroutine set_initial_conditions
 use main_module  
 use config_module
 use module_diag_balance
 implicit none
 integer :: i,j,k
 real(real_type) :: fxa,phiz,x,y
 
 if (.not. enable_diag_balance ) then
  print*,' need balance diagnostic here '
  call halt_stop(' in eddy1.f90')
 endif

 ! u and v from Pablo's idealized eddy 
 do k=ks_pe,ke_pe
  phiz = cos(mode*(k-0.5)*dz*pi/Lz) *sqrt(2*N0) /(Lz*N0) *sqrt(Lz/2)
  do j=js_pe,je_pe 
   do i=is_pe,ie_pe
     x = (i+0.5)*dx + dx - 0.5*Lx ; y = (j+0.5)*dy - Ly/2.
     fxa = sqrt(x**2+y**2) 
     if (  fxa > Rad )  fxa = Rad*exp( - (fxa-Rad)/(Rad/2))
     u(i,j,k) = -fxa/Rad *sin( atan2(y,x) ) *phiz
   enddo
  enddo
  do j=js_pe,je_pe 
   do i=is_pe,ie_pe
     x = (i+0.5)*dx - 0.5*Lx ; y = (j+0.5)*dy + dy - Ly/2.
     fxa = sqrt(x**2+y**2) 
     if (  fxa > Rad )  fxa = Rad*exp( - (fxa-Rad)/(Rad/2))
     v(i,j,k) =  fxa/Rad *cos( atan2(y,x) ) *phiz
   enddo
  enddo 
 enddo 
 call border_exchg_3D(u)
 call border_exchg_3D(v)


 ! b from approximate geostr. balance : p_y = - fu,  b_y = -f u_z
 ! b(j+1) = b(j-1) - f dy u_z 
 do k=ks_pe,ke_pe
   b(:,:,k) = -dy*(u(:,:,k+1)-u(:,:,k-1))/(2*dz)*f0
 enddo
 if (my_blk_k==1) b(:,:,1) = b(:,:,2) 
 if (my_blk_k==n_pes_k) b(:,:,nz) = b(:,:,nz-1)
 call cumsum_in_y(b)
 call border_exchg_3D(b)


 ! scale amplitudes 
 call diag_balance_zero_order
 fxa = maxval( abs( u0(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) ))
 call global_max(fxa)
 fxa = 1./fxa 
 u(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) = fxa*u0(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)
 v(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) = fxa*v0(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)
 w(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) = fxa*w0(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)
 b(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) = fxa*b0(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)
 call border_exchg_3D(u)
 call border_exchg_3D(v)
 call border_exchg_3D(w)
 call border_exchg_3D(b)
  
 ! balance the eddy to 3. order
 call diag_balance_zero_order
 call diag_balance_first_order
 call diag_balance_second_order 
 call diag_balance_third_order 
 u(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) = u0(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) + & 
                                          u1(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) + &
                                          u2(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) + &
                                          u3(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) 
 v(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) = v0(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) + & 
                                          v1(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) + & 
                                          v2(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) + &
                                          v3(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)
 w(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) = w0(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) + & 
                                          w1(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) + & 
                                          w2(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) + &
                                          w3(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)
 b(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) = b0(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) + & 
                                          b1(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) + &
                                          b2(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) + &
                                          b3(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)
 call border_exchg_3D(u)
 call border_exchg_3D(v)
 call border_exchg_3D(w)
 call border_exchg_3D(b)
 
end subroutine set_initial_conditions




subroutine set_forcing 
end subroutine set_forcing
