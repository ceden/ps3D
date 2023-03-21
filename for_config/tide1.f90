
! tidal beam hits eddy


module config_module
 implicit none
 real*8, allocatable :: u_r(:,:,:),v_r(:,:,:),w_r(:,:,:),b_r(:,:,:)
 real*8, allocatable :: u_i(:,:,:),v_i(:,:,:),w_i(:,:,:),b_i(:,:,:)
 real*8, allocatable :: u_b(:,:,:),v_b(:,:,:),w_b(:,:,:),b_b(:,:,:)
 real*8 :: omega_forc, rim=20e3
end module config_module



subroutine set_parameter
 use main_module 
! use config_module
 implicit none
 integer :: fac = 2
 nx=120*fac; ny=120*fac; nz = 10
 
 Ro = 1.; dsqr = 1.
 f0 = 1e-4
 N0 = f0*20
   
 Lx = 400e3
 Ly = 400e3
 Lz = 1000. 
 
 dx=Lx/nx;dy=Ly/ny;dz=Lz/nz


 dt    = 800./fac
 runlen =  1e12
 enable_AB_3_order = .true.
 enable_vertical_boundaries = .true.
 
 enable_diag_snap = .true.
 snapint = 86400.
 tsmonint = snapint
 
 enable_diag_balance        = .true.
end subroutine set_parameter 




subroutine set_initial_conditions
 use main_module  
 use config_module
 use module_diag_balance
 implicit none
 integer :: i,j,k
 real*8 :: kx_d,ky_d,om2,fxa,phiz,x,y
 complex(p3dfft_type) :: g0,qx,qy,qz,qb,px,py,pz,pb
 complex(p3dfft_type),parameter  :: im = (0d0,1d0)
 
 allocate( u_r(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); u_r=0
 allocate( v_r(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); v_r=0
 allocate( w_r(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); w_r=0
 allocate( b_r(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); b_r=0
 allocate( u_i(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); u_i=0
 allocate( v_i(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); v_i=0
 allocate( w_i(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); w_i=0
 allocate( b_i(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); b_i=0

 ! u and v from Pablo's idealized eddy 
 do k=ks_pe,ke_pe
  phiz = cos((k-0.5)*dz*pi/Lz) *sqrt(2*N0) /(Lz*N0) 
  do j=js_pe,je_pe 
   do i=is_pe,ie_pe
     x = (i+0.5)*dx + dx - 0.5*Lx ; y = (j+0.5)*dy - Ly/2.
     fxa = sqrt(x**2+y**2) 
     if (  fxa > rim )  fxa = rim*exp( - (fxa-rim)/(rim/2))
     u(i,j,k) = -0.5*fxa/rim *sin( atan2(y,x) ) *phiz*sqrt(Lz/2)
   enddo
  enddo
  do j=js_pe,je_pe 
   do i=is_pe,ie_pe
     x = (i+0.5)*dx - 0.5*Lx ; y = (j+0.5)*dy + dy - Ly/2.
     fxa = sqrt(x**2+y**2) 
     if (  fxa > rim )  fxa = rim*exp( - (fxa-rim)/(rim/2))
     v(i,j,k) =  0.5*fxa/rim *cos( atan2(y,x) ) *phiz*sqrt(Lz/2)
   enddo
  enddo 
 enddo 
 call border_exchg_3D(u)
 call border_exchg_3D(v)

 ! b from approximate geostr. balance
 do k=ks_pe,ke_pe
   b(:,:,k) = -dy*(u(:,:,k+1)-u(:,:,k-1))/(2*dz)*f0
 enddo
 if (my_blk_k==1) b(:,:,1) = b(:,:,2) 
 if (my_blk_k==n_pes_k) b(:,:,nz) = b(:,:,nz-1)
 call cumsum_in_y(b)
 call border_exchg_3D(b)

 ! balance the eddy
 if (enable_diag_balance ) then
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
 endif
 call border_exchg_3D(u)
 call border_exchg_3D(v)
 call border_exchg_3D(w)
 call border_exchg_3D(b)
 

 ! wave mode amplitude
 g0 = 1.25e-08 * Lx*Ly*Lz
 
 ! construct a single wave as forcing term
 i = 1+28;j = 1;k = 2
 
 kx_d = sin(kx(i)*dx)/dx 
 ky_d = sin(ky(j)*dy)/dy 
 call omega2_discrete(kx(i),ky(j),kz(k),dx,dy,dz,f0,N0,dsqr,om2)
 omega_forc = -sqrt(om2) ! with minus wave propagates in neg. x direction
 if (my_pe==0) print*,' omega/f0 = ',omega_forc/f0
 call vec_q_discrete(kx(i),ky(j),kz(k),dx,dy,dz,1d0,f0,N0,dsqr,qx,qy,qz,qb)   
 if (fstart(1)<=i.and.fend(1)>=i .and. fstart(2)<=j.and.fend(2) >=j .and.fstart(3)<=k.and.fend(3)>=k)  then
  US(i,j,k) = g0*qx
  VS(i,j,k) = g0*qy
  WS(i,j,k) = g0*qz
  BS(i,j,k) = g0*qb 
 endif 
 call backward_transform(US,VS,WS,BS,u_r,v_r,w_r,b_r)
 if (fstart(1)<=i.and.fend(1)>=i .and. fstart(2)<=j.and.fend(2) >=j .and.fstart(3)<=k.and.fend(3)>=k)  then
  US(i,j,k) = -im*g0*qx
  VS(i,j,k) = -im*g0*qy
  WS(i,j,k) = -im*g0*qz
  BS(i,j,k) = -im*g0*qb  
 endif 
 call backward_transform(US,VS,WS,BS,u_i,v_i,w_i,b_i)
 
 
 ! mask forcing terms in physical space
 do j=js_pe,je_pe
  do i=is_pe,ie_pe    
    fxa =      exp( - ( (i-0.5)*dx-0.8*Lx)**2/(20e3)**2 ) &
              *exp( - ( (j-0.5)*dy-0.5*Ly)**2/(20e3)**2 ) 
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
 real*8 :: t
 ! add wavemaker forcing:   u = u + dt*Re( exp(-i omega t) (u_r + i u_i) ) 
 !                            = u + dt*Re( (cos(-i omega t) + i sin(..) ) (u_r + i u_i) )
 !                            = u + dt*( cos(omega t) u_r + sin(omega t ) u_i  ) 
 t = dt*itt
 u =  u + dt*(cos(omega_forc*t)*u_r + sin(omega_forc*t)*u_i )
 v =  v + dt*(cos(omega_forc*t)*v_r + sin(omega_forc*t)*v_i )
 w =  w + dt*(cos(omega_forc*t)*w_r + sin(omega_forc*t)*w_i )
 b =  b + dt*(cos(omega_forc*t)*b_r + sin(omega_forc*t)*b_i )
end subroutine set_forcing
