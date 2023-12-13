
! internal wave beam in background current hits critical layer
! writes file forcing.cdf with background state

module config_module
 use main_module  
 implicit none
 real(real_type), allocatable :: u_r(:,:,:),v_r(:,:,:),w_r(:,:,:),b_r(:,:,:)
 real(real_type), allocatable :: u_i(:,:,:),v_i(:,:,:),w_i(:,:,:),b_i(:,:,:)
 real(real_type), allocatable :: u_b(:,:,:),v_b(:,:,:),w_b(:,:,:),b_b(:,:,:)
 real(real_type), allocatable :: u_back(:,:,:),b_back(:,:,:),p_back(:,:,:)
 real(real_type) :: omega_i,u_back0,omega_e
 real(real_type) :: wavelength_x = 0.005, wavelength_z = 0.02
 real(real_type) :: wave_amplitude = 0.01, jet_amplitude = 0.6, jet_width = 0.3
end module config_module



subroutine set_parameter 
 use config_module
 implicit none
 real(real_type) :: fac=2.
 
 nx= int(64*fac); ny=int(64*fac); nz = int(128*fac)
 dt=10e-3/fac 
 Lx=0.05; Ly=Lx; Lz = 2
 f0=1.; N0=1.; Ro=0.05; dsqr = (0.02)**2

 dx = Lx/nx; dy = Ly/ny; dz = Lz/nz 
 
 tsmonint = 0.1
 snapint = 0.1; 
 runlen =  1e12
 enable_diag_snap = .true.
 enable_AB_3_order = .true.
 
end subroutine set_parameter



subroutine set_initial_conditions
 use config_module
 implicit none
 integer :: i,j,k
 real(real_type) :: kx_d,ky_d,om2,fxa,z,y,fxb
 complex(real_type),allocatable :: US(:,:,:),VS(:,:,:),WS(:,:,:),BS(:,:,:)
 complex(real_type) :: qx,qy,qz,qb,px,py,pz,pb,g0
 complex(real_type),parameter  :: im = (0d0,1d0)
 
 allocate(u_back(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) )
 allocate(b_back(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) )
 allocate(p_back(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) )
 
 allocate( US(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); US=0
 allocate( VS(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); VS=0
 allocate( WS(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); WS=0
 allocate( BS(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); BS=0
 allocate( u_r(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); u_r=0
 allocate( v_r(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); v_r=0
 allocate( w_r(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); w_r=0
 allocate( b_r(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); b_r=0
 allocate( u_i(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); u_i=0
 allocate( v_i(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); v_i=0
 allocate( w_i(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); w_i=0
 allocate( b_i(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); b_i=0

 j=1 
 i = minloc(  (kx - 2*pi/wavelength_x )**2,1 )
 k = minloc(  (kz - 2*pi/wavelength_z )**2,1 )  

 call omega2_discrete(kx(i),ky(j),kz(k),dx,dy,dz,f0,N0,dsqr,om2)
 omega_i = sqrt(om2)

 kx_d = sin(kx(i)*dx)/dx 
 ky_d = sin(ky(j)*dy)/dy 
 u_back0 = omega_i/sqrt(kx_d**2+ky_d**2) /Ro  
 omega_e = omega_i - Ro*u_back0*sqrt(kx_d**2+ky_d**2) 

 jet_amplitude = 1.1*(omega_i-f0)/kx_d /Ro

 if (my_pe==0) then
  print*,' beam has kx = ',kx(i),' 1/m','  -> 2pi/kx=',2*pi/kx(i)/1e3,' km'
  print*,'          ky = ',ky(j),' 1/m','  -> 2pi/ky=',2*pi/ky(j)/1e3,' km'
  print*,'          kz = ',kz(k),' 1/m','  -> 2pi/kz=',2*pi/kz(k)/1e3,' km'
  print*,'     omega/f = ',omega_i/f0
  print*,'  beam angle = ', acos(sqrt(kx(i)**2+ky(j)**2) &
              /sqrt(kx(i)**2+ky(j)**2+kz(k)**2))/pi*180,' deg' 
  print*,'     u_back0 = ',u_back0
  print*,' critical layer for Delta U = ',(omega_i-f0)/kx_d /Ro
  print*,'    jet amplt.  = ',jet_amplitude
  print*,'    wave amplt. = ',wave_amplitude
 endif
 
 ! construct a single wave as forcing term with unit amplitude
 US=0.;VS=0.;WS=0.;BS=0.
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
 ! u_t = ... + cos(omega*t)*u_r ->  u_r \sim f u 
 
 fxb = maxval( abs(u_r(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)))
 fxb = max( fxb, maxval( abs(u_i(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))))
 call global_max(fxb)
 fxb =  wave_amplitude*f0/fxb 
 
 do k=ks_pe,ke_pe
  do j=js_pe,je_pe
   do i=is_pe,ie_pe                     
    fxa =  fxb*exp( - ( (i-0.5)*dx-0.5*Lx)**2/(Lx/10.)**2 ) &
              *exp( - ( (j-0.5)*dy-0.5*Ly)**2/(Lx/10.)**2 ) &
              *exp( - ( (k-0.5)*dz-0.25)**2/(Lz/20.)**2 )   
    u_r(i,j,k) = u_r(i,j,k)*fxa
    v_r(i,j,k) = v_r(i,j,k)*fxa
    w_r(i,j,k) = w_r(i,j,k)*fxa
    b_r(i,j,k) = b_r(i,j,k)*fxa
    u_i(i,j,k) = u_i(i,j,k)*fxa
    v_i(i,j,k) = v_i(i,j,k)*fxa
    w_i(i,j,k) = w_i(i,j,k)*fxa
    b_i(i,j,k) = b_i(i,j,k)*fxa
   enddo
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
 
 ! construct background state U(y,z), B(y,z), balanced to zero order
 
 do k=ks_pe,ke_pe
  do j=js_pe,je_pe
   z = (k-0.5)*dz
   y = (j-0.5)*dy
   p_back(:,j,k) = exp(-(z-0.75*Lz)**2/jet_width**2)*sin(2*pi*y/Ly)
  enddo 
 enddo
 call border_exchg_3D(p_back)

 do j=js_pe,je_pe 
  u_back(:,j,:) = (P_back(:,j+1,:)-P_back(:,j,:))/(f0*dy)
 enddo
 call border_exchg_3D(u_back)
 
 ! project on geostr. mode again, to balance exactly to zero order
 v_r=0d0; w_r=0d0; b_back=0.
 call forward_transform(u_back,v_r,w_r,b_back,US,VS,WS,BS)
 do k=fstart(3),fend(3)
  do j=fstart(2),fend(2)
   do i=fstart(1),fend(1)  
      call vec_q_discrete(kx(i),ky(j),kz(k),dx,dy,dz,0d0,f0,N0,dsqr,qx,qy,qz,qb)
      call vec_p_discrete(kx(i),ky(j),kz(k),dx,dy,dz,0d0,f0,N0,dsqr,px,py,pz,pb)
      g0 = px*US(i,j,k) + py*VS(i,j,k) + pz*WS(i,j,k) + pb*BS(i,j,k)
      US(i,j,k) = g0*qx
      VS(i,j,k) = g0*qy
      WS(i,j,k) = g0*qz
      BS(i,j,k) = g0*qb      
    enddo
   enddo
 enddo 
 call backward_transform(US,VS,WS,BS,u_back,v_r,w_r,b_back) 
 
 fxa = maxval( abs( u_back(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) ))
 call global_max(fxa)
 fxa = jet_amplitude/fxa 
 u_back(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) = fxa*u_back(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)
 b_back(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) = fxa*b_back(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)
 call border_exchg_3D(u_back)
 call border_exchg_3D(b_back) 

 u_back = u_back + u_back0
 deallocate(US,VS,WS,BS)
 
 
 call write_forcing()
 
end subroutine set_initial_conditions



subroutine set_forcing
 use config_module
 implicit none
 integer :: i,j,k
 real(real_type) :: t

 ! add mean flow (U(x,y,z),0,0)  and B=B(x,y,z) 
 ! Du/Dt - f v + p_x =  - \vn \v U u - \vn \v u  U  = -(U u)_x - (u U)_x - (v U)_y - (w U)_z
 ! Dv/Dt + f u + p_y =  - \vn \v U v                = -(U v)_x 
 ! Dw/Dt - b   + p_z =  - \vn \v U w                = -(U w)_x 
 ! Db/Dt + w N^2     =  - \vn \v U  b - \vn \v u  B = -(U b)_x - (u B)_x - (v B)_y - (w B)_z
 ! .. tbc

 do i=is_pe-1,ie_pe 
  flux_east(i,:,:) = 2* 0.25*(u_back(i,:,:)+u_back(i+1,:,:))*(u(i+1,:,:)+u(i,:,:))                    
 enddo 
 do j=js_pe-1,je_pe
  do i=is_pe,ie_pe
    flux_north(i,j,:) = 0.25*(u_back(i,j,:)+u_back(i,j+1,:))*(v(i+1,j,:)+v(i,j,:))
  enddo
 enddo
 do k=ks_pe-1,ke_pe
  do i=is_pe,ie_pe
    flux_top(i,:,k) = 0.25*(u_back(i,:,k+1)+u_back(i,:,k))*(w(i,:,k)+w(i+1,:,k))
  enddo
 enddo
 if (enable_vertical_boundaries) then
   if (my_blk_k == 1      ) flux_top(:,:,1-onx:0)   = 0d0
   if (my_blk_k == n_pes_k) flux_top(:,:,nz:nz+onx) = 0d0
 endif
 do k=ks_pe,ke_pe
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
    du(i,j,k,tau) = du(i,j,k,tau) - Ro*( ( flux_east(i,j,k) - flux_east(i-1,j,k) )/dx &
                                        +(flux_north(i,j,k) - flux_north(i,j-1,k) )/dy &
                                        +(  flux_top(i,j,k) - flux_top  (i,j,k-1) )/dz )
   enddo
  enddo
 enddo

 do j=js_pe,je_pe
  do i=is_pe-1,ie_pe
   flux_east(i,j,:) = 0.25*(v(i,j,:)+v(i+1,j,:))*(u_back(i,j+1,:)+u_back(i,j,:))
  enddo
 enddo
 do i=is_pe,ie_pe
  dv(i,:,:,tau) = dv(i,:,:,tau) - Ro* ( flux_east(i,:,:) - flux_east(i-1,:,:) )/dx 
 enddo

 do k=ks_pe,ke_pe
  do i=is_pe-1,ie_pe
   flux_east(i,:,k) = 0.25*(w(i,:,k)+w(i+1,:,k))*(u_back(i,:,k)+u_back(i,:,k+1))
  enddo
 enddo
 do i=is_pe,ie_pe
   dw(i,:,:,tau) = dw(i,:,:,tau) -dsqr*Ro*( flux_east(i,:,:) - flux_east(i-1,:,:) )/dx
 enddo

 do i=is_pe-1,ie_pe
  flux_east(i,:,:)= 0.5*(b(i,:,:) + b(i+1,:,:) )*u_back(i,:,:)  &
                   +0.5*(b_back(i,:,:) + b_back(i+1,:,:) )*u(i,:,:)
 enddo
 do j=js_pe-1,je_pe
  flux_north(:,j,:)=0.5*(b_back(:,j,:) + b_back(:,j+1,:) )*v(:,j,:)
 enddo
 do k=ks_pe-1,ke_pe
  flux_top(:,:,k)=0.5*(b_back(:,:,k) + b_back(:,:,k+1) )*w(:,:,k)
 enddo
 if (enable_vertical_boundaries) then
   if (my_blk_k == 1      ) flux_top(:,:,1-onx:0) = 0d0
   if (my_blk_k == n_pes_k) flux_top(:,:,nz:nz+onx) = 0d0
 endif  
 do k=ks_pe,ke_pe
   do j=js_pe,je_pe
    do i=is_pe,ie_pe
      db(i,j,k,tau)= db(i,j,k,tau)-Ro*( (flux_east(i,j,k)  -  flux_east(i-1,j,k))/dx &
                                       +(flux_north(i,j,k) - flux_north(i,j-1,k))/dy &
                                       +(flux_top(i,j,k)   -   flux_top(i,j,k-1))/dz )
    enddo
   enddo
 enddo


 ! add wavemaker forcing
 t = dt*itt
 !du(:,:,:,tau) =  du(:,:,:,tau) + (cos(-omega_e*t)*u_r(:,:,:) + sin(-omega_e*t)*u_i(:,:,:) )
 !dv(:,:,:,tau) =  dv(:,:,:,tau) + (cos(-omega_e*t)*v_r(:,:,:) + sin(-omega_e*t)*v_i(:,:,:) )
 !dw(:,:,:,tau) =  dw(:,:,:,tau) + (cos(-omega_e*t)*w_r(:,:,:) + sin(-omega_e*t)*w_i(:,:,:) )
 !db(:,:,:,tau) =  db(:,:,:,tau) + (cos(-omega_e*t)*b_r(:,:,:) + sin(-omega_e*t)*b_i(:,:,:) )
 u(:,:,:) =  u(:,:,:) + dt*(cos(-omega_e*t)*u_r(:,:,:) + sin(-omega_e*t)*u_i(:,:,:) )
 v(:,:,:) =  v(:,:,:) + dt*(cos(-omega_e*t)*v_r(:,:,:) + sin(-omega_e*t)*v_i(:,:,:) )
 w(:,:,:) =  w(:,:,:) + dt*(cos(-omega_e*t)*w_r(:,:,:) + sin(-omega_e*t)*w_i(:,:,:) )
 b(:,:,:) =  b(:,:,:) + dt*(cos(-omega_e*t)*b_r(:,:,:) + sin(-omega_e*t)*b_i(:,:,:) )

 
 call border_exchg_3D(u)
 call border_exchg_3D(v)
 call border_exchg_3D(w)
 call border_exchg_3D(b)
 
end subroutine set_forcing





subroutine write_forcing
 use config_module
 implicit none
 include "netcdf.inc"
 include "mpif.h"
 integer :: ncid,iret,xdim,ydim,zdim,id,n
 real(real_type), allocatable :: a(:,:,:)
 integer :: tag=1,ist(3),isz(3),ien(3)
 integer, dimension(MPI_STATUS_SIZE) :: Status
 
 
 if (my_pe==0) then
 iret = nf_create ('forcing.cdf', NF_CLOBBER, ncid)
 if (iret.ne.0) print*,nf_strerror(iret)
 if (iret.ne.0) stop
 iret = nf_put_att_int(ncid,NF_GLOBAL,'nx',nf_int,1,nx)
 iret = nf_put_att_int(ncid,NF_GLOBAL,'ny',nf_int,1,ny)
 iret = nf_put_att_int(ncid,NF_GLOBAL,'nz',nf_int,1,nz)
 
 iret = nf_put_att_double(ncid,NF_GLOBAL,'Lx',nf_double,1,Lx)
 iret = nf_put_att_double(ncid,NF_GLOBAL,'Ly',nf_double,1,Ly)
 iret = nf_put_att_double(ncid,NF_GLOBAL,'Lz',nf_double,1,Lz)
 iret = nf_put_att_double(ncid,NF_GLOBAL,'dx',nf_double,1,dx)
 iret = nf_put_att_double(ncid,NF_GLOBAL,'dy',nf_double,1,dy)
 iret = nf_put_att_double(ncid,NF_GLOBAL,'dz',nf_double,1,dz)
 iret = nf_put_att_double(ncid,NF_GLOBAL,'f0',nf_double,1,f0)
 iret = nf_put_att_double(ncid,NF_GLOBAL,'N0',nf_double,1,N0)
 iret = nf_put_att_double(ncid,NF_GLOBAL,'dt',nf_double,1,dt)

 iret = nf_put_att_double(ncid,NF_GLOBAL,'Ahbi',nf_double,1,0d0)
 iret = nf_put_att_double(ncid,NF_GLOBAL,'Avbi',nf_double,1,0d0)
 iret = nf_put_att_double(ncid,NF_GLOBAL,'Khbi',nf_double,1,0d0)
 iret = nf_put_att_double(ncid,NF_GLOBAL,'Kvbi',nf_double,1,0d0)
 iret = nf_put_att_double(ncid,NF_GLOBAL,'Ah',nf_double,1,Ah)
 iret = nf_put_att_double(ncid,NF_GLOBAL,'Av',nf_double,1,Av)
 iret = nf_put_att_double(ncid,NF_GLOBAL,'Kh',nf_double,1,Kh)
 iret = nf_put_att_double(ncid,NF_GLOBAL,'Kv',nf_double,1,Kv)
 iret = nf_put_att_double(ncid,NF_GLOBAL,'snapint',nf_double,1,snapint)
 iret = nf_put_att_double(ncid,NF_GLOBAL,'runlen',nf_double,1,runlen)
 iret = nf_put_att_double(ncid,NF_GLOBAL,'epsilon',nf_double,1,1d-8)
 iret = nf_put_att_double(ncid,NF_GLOBAL,'epsilon_2D',nf_double,1,1d-8)
 iret = nf_put_att_double(ncid,NF_GLOBAL,'AB_eps',nf_double,1,eps_ab)
 
 xdim  = ncddef(ncid, 'x', nx, iret)
 ydim  = ncddef(ncid, 'y', ny, iret)
 zdim  = ncddef(ncid, 'z', nz, iret)
 id  = ncvdef (ncid,'u',NF_DOUBLE,3,(/xdim, ydim,zdim/),iret)
 id  = ncvdef (ncid,'b',NF_DOUBLE,3,(/xdim, ydim,zdim/),iret)
 iret= nf_enddef(ncid) 
 
 iret=nf_inq_varid(ncid,'u',id)
 iret= nf_put_vara_double(ncid,id,(/istart(1),istart(2),istart(3)/), &
                            (/isize(1),isize(2),isize(3)/),u_back(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))
 iret=nf_inq_varid(ncid,'b',id)
 iret= nf_put_vara_double(ncid,id,(/istart(1),istart(2),istart(3)/), &
                            (/isize(1),isize(2),isize(3)/),b_back(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))

 endif
  
 do n=1,n_pes-1
  call mpi_barrier(MPI_COMM_WORLD, iret)
  if (my_pe==n) then
        call mpi_send(istart,3,mpi_integer,0,tag,mpi_comm_world,iret)
        call mpi_send(iend  ,3,mpi_integer,0,tag,mpi_comm_world,iret)
        call mpi_send(isize ,3,mpi_integer,0,tag,mpi_comm_world,iret)
        call mpi_send(u_back(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe),isize(1)*isize(2)*isize(3),mpi_real8,0,tag,mpi_comm_world,iret)
        call mpi_send(b_back(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe),isize(1)*isize(2)*isize(3),mpi_real8,0,tag,mpi_comm_world,iret)
  else if (my_pe==0) then
        call mpi_recv(ist,3,mpi_integer,n,tag,mpi_comm_world,status,iret)
        call mpi_recv(ien,3,mpi_integer,n,tag,mpi_comm_world,status,iret)
        call mpi_recv(isz,3,mpi_integer,n,tag,mpi_comm_world,status,iret)
        allocate(a(ist(1):ien(1),ist(2):ien(2),ist(3):ien(3)) )
        
        call mpi_recv(a,isz(1)*isz(2)*isz(3),mpi_real8,n,tag,mpi_comm_world,status,iret)
        iret=nf_inq_varid(ncid,'u',id)
        iret= nf_put_vara_double(ncid,id,(/ist(1),ist(2),ist(3)/),(/isz(1),isz(2),isz(3)/),a)
        
        call mpi_recv(a,isz(1)*isz(2)*isz(3),mpi_real8,n,tag,mpi_comm_world,status,iret)
        iret=nf_inq_varid(ncid,'b',id)
        iret= nf_put_vara_double(ncid,id,(/ist(1),ist(2),ist(3)/),(/isz(1),isz(2),isz(3)/),a)
           
        deallocate(a)
  endif
 enddo  

 if (my_pe == 0) iret= nf_close(ncid) 
 
end subroutine write_forcing





