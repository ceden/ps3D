
! tidal beam hits eddy


module config_module
 implicit none
 real*8, allocatable :: u_r(:,:,:),v_r(:,:,:),w_r(:,:,:),b_r(:,:,:)
 real*8, allocatable :: u_i(:,:,:),v_i(:,:,:),w_i(:,:,:),b_i(:,:,:)
 real*8 :: omega_forc
 integer :: fac = 1
 real*8 :: wavelength = 20e3
end module config_module



subroutine set_parameter
 use main_module 
 use config_module
 implicit none
 
 nx=126*fac; ny=126*fac; nz = 10*fac
 Ro = 1.; dsqr = 1. 
 f0 = 1e-4
 N0 = f0*30  
 
 Lx = 500e3
 Ly = 500e3
 Lz = 2000.  

 dx=Lx/nx;dy=Ly/ny;dz=Lz/nz

 dt    = 400./fac
 runlen =  15*86400.
 enable_AB_3_order = .true.
 enable_vertical_boundaries = .true.
 
 enable_diag_snap = .true.
 snapint = 86400.
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
 real(real_type) :: om2,fxa,phiz,x,y,Rad,fxb
 complex(real_type) :: g0,qx,qy,qz,qb,px,py,pz,pb
 complex(real_type),parameter  :: im = (0d0,1d0)
 
 
 if (.not. enable_diag_balance ) then
  print*,' need balance diagnostic here '
  call halt_stop(' in tide1.f90')
 endif
 
 allocate( u_r(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); u_r=0
 allocate( v_r(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); v_r=0
 allocate( w_r(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); w_r=0
 allocate( b_r(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); b_r=0
 allocate( u_i(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); u_i=0
 allocate( v_i(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); v_i=0
 allocate( w_i(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); w_i=0
 allocate( b_i(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); b_i=0

if (.true.) then
 
 Rad = 2*N0*Lz/(pi) /f0
 ! u and v from Pablo's idealized eddy 
 do k=ks_pe,ke_pe
  phiz = cos((k-0.5)*dz*pi/Lz) *sqrt(2*N0) /(Lz*N0) 
  do j=js_pe,je_pe 
   do i=is_pe,ie_pe
     x = (i+0.5)*dx + dx - 0.5*Lx ; y = (j+0.5)*dy - Ly/2.
     fxa = sqrt(x**2+y**2) 
     if (  fxa > Rad )  fxa = Rad*exp( - (fxa-Rad)/(Rad/2))
     u(i,j,k) = -fxa/Rad *sin( atan2(y,x) ) *phiz*sqrt(Lz/2)
   enddo
  enddo
  do j=js_pe,je_pe 
   do i=is_pe,ie_pe
     x = (i+0.5)*dx - 0.5*Lx ; y = (j+0.5)*dy + dy - Ly/2.
     fxa = sqrt(x**2+y**2) 
     if (  fxa > Rad )  fxa = Rad*exp( - (fxa-Rad)/(Rad/2))
     v(i,j,k) =  fxa/Rad *cos( atan2(y,x) ) *phiz*sqrt(Lz/2)
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
 fxa = 0.2/fxa 
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

endif


 ! wave mode amplitude: construct a single wave as forcing term 
 j=1 
 i = minloc(  (kx - 2*pi/wavelength )**2,1 )
 k = minloc(  (kz - pi/Lz )**2,1 ) ! first vertical mode
 
 call omega2_discrete(kx(i),ky(j),kz(k),dx,dy,dz,f0,N0,dsqr,om2) 
 if (my_pe==0) then
  print*,' beam has kx = ',kx(i),' 1/m','  -> 2pi/kx=',2*pi/kx(i)/1e3,' km'
  print*,'          ky = ',ky(j),' 1/m','  -> 2pi/ky=',2*pi/ky(j)/1e3,' km'
  print*,'          kz = ',kz(k),' 1/m','  -> 2pi/kz=',2*pi/kz(k)/1e3,' km'
  print*,'     omega/f = ',sqrt(om2)/f0
  !print*,'  omega_M2/f = ',2*pi/( 12.*60*60 +  25.2 *60 ) /f0
 endif
 
 US=0;VS=0;WS=0;BS=0
 omega_forc = -sqrt(om2) ! with minus wave propagates in neg. x direction
 call vec_q_discrete(kx(i),ky(j),kz(k),dx,dy,dz,1d0,f0,N0,dsqr,qx,qy,qz,qb)   
 if (fstart(1)<=i.and.fend(1)>=i .and. fstart(2)<=j.and.fend(2) >=j .and.fstart(3)<=k.and.fend(3)>=k)  then
  US(i,j,k) = qx
  VS(i,j,k) = qy
  WS(i,j,k) = qz
  BS(i,j,k) = qb 
 endif 
 call backward_transform(US,VS,WS,BS,u_r,v_r,w_r,b_r)
 if (fstart(1)<=i.and.fend(1)>=i .and. fstart(2)<=j.and.fend(2) >=j .and.fstart(3)<=k.and.fend(3)>=k)  then
  US(i,j,k) = -im*qx
  VS(i,j,k) = -im*qy
  WS(i,j,k) = -im*qz
  BS(i,j,k) = -im*qb  
 endif 
 call backward_transform(US,VS,WS,BS,u_i,v_i,w_i,b_i)
  
 ! mask forcing terms in physical space and scale wavemaker forcing
 ! u_t = ... + cos(omega*t)*u_r , u(t) = sin(omega*t) u_r /omega
 fxb = maxval( abs(u_r(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)))
 fxb = max( fxb, maxval( abs(u_i(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))))
 call global_max(fxb)
 fxb = 0.5e-6*(abs(omega_forc)/f0) /fxb 
 
 do j=js_pe,je_pe
  do i=is_pe,ie_pe     
    fxa =   fxb*exp( - ( (i-0.5)*dx-0.9*Lx)**2/(2*wavelength)**2 ) &
               *exp( - ( (j-0.5)*dy-0.5*Ly)**2/(2*wavelength)**2 ) 
    u_r(i,j,:) = u_r(i,j,:)*fxa
    v_r(i,j,:) = v_r(i,j,:)*fxa
    w_r(i,j,:) = w_r(i,j,:)*fxa
    b_r(i,j,:) = b_r(i,j,:)*fxa
    u_i(i,j,:) = u_i(i,j,:)*fxa
    v_i(i,j,:) = v_i(i,j,:)*fxa
    w_i(i,j,:) = w_i(i,j,:)*fxa
    b_i(i,j,:) = b_i(i,j,:)*fxa
  enddo
 enddo
 
 ! project masked forcing on wave mode again
 call forward_transform(u_r,v_r,w_r,b_r,US,VS,WS,BS)
 do k=fstart(3),fend(3)
  do j=fstart(2),fend(2)
   do i=fstart(1),fend(1)  
      call vec_q_discrete(kx(i),ky(j),kz(k),dx,dy,dz,1d0,f0,N0,dsqr,qx,qy,qz,qb)
      call vec_p_discrete(kx(i),ky(j),kz(k),dx,dy,dz,1d0,f0,N0,dsqr,px,py,pz,pb)
      g0 = px*US(i,j,k) + py*VS(i,j,k) + pz*WS(i,j,k) + pb*BS(i,j,k)
      US(i,j,k) = g0*qx
      VS(i,j,k) = g0*qy
      WS(i,j,k) = g0*qz
      BS(i,j,k) = g0*qb      
    enddo
   enddo
 enddo 
 call backward_transform(US,VS,WS,BS,u_r,v_r,w_r,b_r) 

 call forward_transform(u_i,v_i,w_i,b_i,US,VS,WS,BS)
 do k=fstart(3),fend(3)
  do j=fstart(2),fend(2)
   do i=fstart(1),fend(1)  
      call vec_q_discrete(kx(i),ky(j),kz(k),dx,dy,dz,1d0,f0,N0,dsqr,qx,qy,qz,qb)
      call vec_p_discrete(kx(i),ky(j),kz(k),dx,dy,dz,1d0,f0,N0,dsqr,px,py,pz,pb)
      g0 = px*US(i,j,k) + py*VS(i,j,k) + pz*WS(i,j,k) + pb*BS(i,j,k)
      US(i,j,k) = g0*qx
      VS(i,j,k) = g0*qy
      WS(i,j,k) = g0*qz
      BS(i,j,k) = g0*qb      
    enddo
   enddo
 enddo 
 call backward_transform(US,VS,WS,BS,u_i,v_i,w_i,b_i) 

 call border_exchg_3D(u_r)
 call border_exchg_3D(v_r)
 call border_exchg_3D(w_r)
 call border_exchg_3D(b_r)
 call border_exchg_3D(u_i)
 call border_exchg_3D(v_i)
 call border_exchg_3D(w_i)
 call border_exchg_3D(b_i) 
end subroutine set_initial_conditions




subroutine set_forcing 
 use main_module  
 use config_module
 implicit none
 real(real_type) :: t
 ! add wavemaker forcing:   u = u + dt*Re( exp(-i omega t) (u_r + i u_i) ) 
 !                            = u + dt*Re( (cos(-i omega t) + i sin(..) ) (u_r + i u_i) )
 !                            = u + dt*( cos(omega t) u_r + sin(omega t ) u_i  ) 
 t = dt*itt
 u =  u + dt*(cos(omega_forc*t)*u_r + sin(omega_forc*t)*u_i )
 v =  v + dt*(cos(omega_forc*t)*v_r + sin(omega_forc*t)*v_i )
 w =  w + dt*(cos(omega_forc*t)*w_r + sin(omega_forc*t)*w_i )
 b =  b + dt*(cos(omega_forc*t)*b_r + sin(omega_forc*t)*b_i )
end subroutine set_forcing
