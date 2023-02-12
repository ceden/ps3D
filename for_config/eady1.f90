
subroutine set_parameter
 use main_module   
 implicit none

 nx=100;ny=100;nz=20
 Lx = 4.0; Ly=Lx; Lz = 1.
 dx=Lx/nx;dy=Ly/ny;dz=Lz/nz
 
 Ro = 0.1
 f0 = 1.
 N0 = 1.
 dt = 0.002
 eps_ab = 0.01
 dsqr = (0.2)**2
 tsmonint = 10*dt
 snapint = 50*dt
 runlen =  1e12
 
 enable_vertical_boundaries = .true.
 enable_diag_snap = .true.
 enable_AB_3_order = .true.

end subroutine set_parameter 




subroutine set_initial_conditions
 use main_module  

 implicit none
 integer :: i,j,k
 real*8 :: x(nx),y(ny),z(nz)
 
 do i=1,nx
   x(i)=(i-1)*dx
 enddo 
 do i=1,ny
   y(i)=(i-1)*dy
 enddo 
 do i=1,nz
   z(i)=(i-1)*dz
 enddo 
 
 do k=ks_pe,ke_pe
   do j=js_pe,je_pe
    do i=is_pe,ie_pe  
     u(i,j,k) = 2*( ( exp( -(y(j)-3*Ly/4.)**2/(Lx*0.04)**2 ) &
                     -exp( -(y(j)-1*Ly/4.)**2/(Lx*0.04)**2 )  )*cos(z(k)/Lz*pi) &
          + 0.1*sin(x(i)/Lx*4*pi)*sin(y(j)/Ly*2*pi)*cos(z(k)/Lz*pi) )
   enddo
  enddo
 enddo 
 call border_exchg_3D(u)

 do k=ks_pe,ke_pe
   b(:,:,k) = -dy*(u(:,:,k+1)-u(:,:,k-1))/(2*dz)*f0
 enddo
 if (my_blk_k==1) b(:,:,1) = 0.
 if (my_blk_k==n_pes_k) b(:,:,nz) = 0.
 
 call cumsum_in_y(b)
 call border_exchg_3D(b)
end subroutine set_initial_conditions




subroutine set_forcing 
end subroutine set_forcing

