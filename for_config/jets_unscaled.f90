! two jets in unscaled equations

subroutine set_parameter
 use main_module   
 implicit none

 enable_vertical_boundaries = .true.

 !Lx = 4.0; Ly=Lx; Lz = 1.0
 !Ro = 0.1; dsqr = (0.2)**2
 !f0 = 1.; N0 = 1.; dt = 0.1*dsqr ! dt = 0.002 

 Ro = 1.; dsqr = 1. 
 Lx = 500e3; Ly=Lx; Lz = 1000.0
 f0 = 1e-4; N0 = 50*f0; dt = 0.5/N0! 4./N0
  
 nx=100;ny=100;nz=20 
 if (enable_vertical_boundaries) then
   nz=nz/2
   Lz=Lz/2.
 endif
 
 dx=Lx/nx;dy=Ly/ny;dz=Lz/nz 
 snapint =  500*dt
 tsmonint = snapint
 
 runlen =  1e12 
 enable_diag_snap = .true.
 enable_AB_3_order = .true.

end subroutine set_parameter 




subroutine set_initial_conditions
 use main_module  
 implicit none
 integer :: i,j,k
 real(real_type) :: x(nx),y(ny),z(nz)
 
 do i=1,nx
   x(i)=(i-1)*dx
 enddo 
 do i=1,ny
   y(i)=(i-1)*dy
 enddo 
 do i=1,nz
   z(i)=(i-1)*dz
 enddo 
 
 if (enable_vertical_boundaries) then
 
   do k=ks_pe,ke_pe
    do j=js_pe,je_pe
     do i=is_pe,ie_pe  
      u(i,j,k) = .1*( ( exp( -(y(j)-3*Ly/4.)**2/(Lx*0.02)**2 ) &
                       -exp( -(y(j)-1*Ly/4.)**2/(Lx*0.02)**2 )  )*cos(z(k)/Lz*pi) &
           + 0.02*sin(x(i)/Lx*10*pi)*sin(y(j)/Ly*2*pi)*cos(z(k)/Lz*pi) )     
     enddo
    enddo
   enddo 
   
 else
  do k=ks_pe,ke_pe
   do j=js_pe,je_pe
    do i=is_pe,ie_pe  
     u(i,j,k) = .1*( ( exp( -(y(j)-3*Ly/4.)**2/(Lx*0.02)**2 ) &
                      -exp( -(y(j)-1*Ly/4.)**2/(Lx*0.02)**2 )  )*cos(z(k)/Lz*2*pi) &
          + 0.02*sin(x(i)/Lx*10*pi)*sin(y(j)/Ly*2*pi)*cos(z(k)/Lz*2*pi) )     
    enddo
   enddo
  enddo 
 endif
 
 call border_exchg_3D(u)
 call border_exchg_3D(v)

  ! 0 = -p_z + b, 0 = - p_y - fu,  b_y = -fu_z,    
  if (enable_vertical_boundaries) then
    do k=ks_pe,ke_pe
     if (k==1) then
      b(:,:,k) = -dy*(u(:,:,k+1)-u(:,:,k))/(dz)*f0
     elseif (k==nz) then
      b(:,:,k) = -dy*(u(:,:,k)-u(:,:,k-1))/(dz)*f0
     else
      b(:,:,k) = -dy*(u(:,:,k+1)-u(:,:,k-1))/(2*dz)*f0 
     endif  
    enddo 
   else      
    do k=ks_pe,ke_pe
      b(:,:,k) = -dy*(u(:,:,k+1)-u(:,:,k-1))/(2*dz)*f0 
    enddo
  endif 
  call cumsum_in_y(b)

 call border_exchg_3D(b)
end subroutine set_initial_conditions




subroutine set_forcing 
end subroutine set_forcing


