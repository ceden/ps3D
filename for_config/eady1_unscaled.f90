
module config_module
 ! use this module only locally in this file
 implicit none
 integer,parameter :: fac=1 
 real*8,parameter :: Ri  = 2000    ! Richardson number 
 real*8,parameter :: delta = 0.02 ! aspect ratio
 real*8,parameter :: H0   = 1000.0       ! total depth
 real*8 :: U0,M0,kmax,jet_scale
end module config_module



subroutine set_parameter
 use main_module 
 use config_module
 implicit none
 
 nx=30*fac; ny=30*2*fac; nz = 10*fac
 
 Ro = 1.; dsqr = 1.
 f0 = 1e-4
 N0 = f0/delta
 U0 = sqrt(1./Ri)*N0*H0   
 M0 = sqrt(f0*U0/H0)  
 kmax   = 1./ ( sqrt(2./5.)*sqrt(1+Ri)*M0**2/f0**2*H0 )  
 Lx = 4*2/kmax  *pi
 Ly = 4*2/kmax  *pi * 2
 Lz = H0
 jet_scale = 0.1*Ly/2
 
 
 dx=Lx/nx;dy=Ly/ny;dz=Lz/nz

 eps_AB = 0.01
 dt    = 800./fac

 tsmonint = 10*dt
 snapint = 50*dt
 runlen =  1e12
 
 enable_vertical_boundaries = .true.
 enable_diag_snap = .true.
 enable_AB_3_order = .true.

end subroutine set_parameter 




subroutine set_initial_conditions
 use main_module  
 use config_module
 implicit none
 integer :: i,j,k
 real*8 :: xt(nx),yt(ny),zt(nz)
 
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
     u(i,j,k) = U0*( exp(-(yt(j)-1*Ly/4)**2/jet_scale**2) &
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

end subroutine set_initial_conditions




subroutine set_forcing 
end subroutine set_forcing
