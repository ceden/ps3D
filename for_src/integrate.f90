
subroutine integrate
 use main_module
 use timing_module   
 implicit none
 integer :: i,j,k
 
 call tic('integrate')
 !---------------------------------------------------------------------------------
 ! buoyancy advection
 !--------------------------------------------------------------------------------- 
 if (enable_nonlinear) then  
  call buoyancy_advection_2nd(db(:,:,:,tau))
 else  
  db(:,:,:,tau)=0d0
 endif 
 
 !---------------------------------------------------------------------------------
 ! advection of background stratification
 !---------------------------------------------------------------------------------
 do k=ks_pe,ke_pe
  db(:,:,k,tau) = db(:,:,k,tau) - 0.5*(w(:,:,k)+w(:,:,k-1))*N0**2
 enddo 
 
 !---------------------------------------------------------------------------------
 ! momentum advection
 !---------------------------------------------------------------------------------
 if (enable_nonlinear) then
   call momentum_advection_2nd(du(:,:,:,tau),dv(:,:,:,tau),dw(:,:,:,tau))
 else  
   du(:,:,:,tau)=0d0;dv(:,:,:,tau)=0d0;dw(:,:,:,tau)=0d0
 endif
 
 !---------------------------------------------------------------------------------
 !  add time tendency due to Coriolis force 
 !---------------------------------------------------------------------------------
 do j=js_pe,je_pe
  do i=is_pe,ie_pe
   du(i,j,:,tau) = du(i,j,:,tau) + (v(i,j,:)+v(i,j-1,:) +v(i+1,j,:)+v(i+1,j-1,:))*0.25*f0
   dv(i,j,:,tau) = dv(i,j,:,tau) - (u(i-1,j,:)+u(i,j,:) +u(i-1,j+1,:)+u(i,j+1,:))*0.25*f0
  enddo
 enddo

 !---------------------------------------------------------------------------------
 ! add time tendency due to buoyancy force
 !---------------------------------------------------------------------------------
 do k=ks_pe,ke_pe
   dw(:,:,k,tau) = dw(:,:,k,tau)  + 0.5*(b(:,:,k)+b(:,:,k+1))
 enddo

 if (enable_vertical_boundaries) then
   if (my_blk_k == 1      ) dw(:,:,1-onx:0,tau) = 0d0
   if (my_blk_k == n_pes_k) dw(:,:,nz:nz+onx,tau) = 0d0
 endif

 !---------------------------------------------------------------------------------
 ! intermediate result due to friction
 !---------------------------------------------------------------------------------
 if (Ah>0. .or. Av>0.) call harmonic(u,Ah,Av)
 if (Ah>0. .or. Av>0.) call harmonic(v,Ah,Av)
 if (Ah>0. .or. Av>0.) then 
   call harmonic(w,Ah,Av) 
   if (enable_vertical_boundaries .and. my_blk_k == n_pes_k)  w(:,:,nz:nz+onx) = 0d0
 endif
 
 if (Kh>0. .or. Kv>0.) call harmonic(b,Kh,Kv)
 
 if (Ahbi>0.) call biharmonic_h(u,Ahbi)
 if (Ahbi>0.) call biharmonic_h(v,Ahbi)
 if (Ahbi>0.) call biharmonic_h(w,Ahbi)
 if (Avbi>0.) call biharmonic_v(u,Avbi)
 if (Avbi>0.) call biharmonic_v(v,Avbi)
 if (Avbi>0.) then 
   call biharmonic_v(w,Avbi) 
   if (enable_vertical_boundaries .and. my_blk_k == n_pes_k)  w(:,:,nz:nz+onx) = 0d0
 endif  
 if (Khbi>0.) call biharmonic_h(b,Khbi)
 if (Kvbi>0.) call biharmonic_v(b,Kvbi)
 
 
 ! allow for user defined forcing
 call set_forcing()
 !---------------------------------------------------------------------------------
 ! intermediate result with Adams Bashforth time stepping
 !---------------------------------------------------------------------------------
 
 if (itt==0) then
  u = u +      dt*du(:,:,:,tau) 
  v = v +      dt*dv(:,:,:,tau) 
  w = w + dt/dsqr*dw(:,:,:,tau) 
  b = b +      dt*db(:,:,:,tau) 
 elseif (itt==1 .or. .not. enable_AB_3_order) then     
  u = u +      dt*( (1.5+eps_ab)*du(:,:,:,tau) - ( 0.5+eps_ab)*du(:,:,:,taum1)) 
  v = v +      dt*( (1.5+eps_ab)*dv(:,:,:,tau) - ( 0.5+eps_ab)*dv(:,:,:,taum1)) 
  w = w + dt/dsqr*( (1.5+eps_ab)*dw(:,:,:,tau) - ( 0.5+eps_ab)*dw(:,:,:,taum1)) 
  b = b +      dt*( (1.5+eps_ab)*db(:,:,:,tau) - ( 0.5+eps_ab)*db(:,:,:,taum1)) 
 else
  u = u +      dt*( AB3_a*du(:,:,:,tau) + AB3_b*du(:,:,:,taum1) + AB3_c*du(:,:,:,taum2)) 
  v = v +      dt*( AB3_a*dv(:,:,:,tau) + AB3_b*dv(:,:,:,taum1) + AB3_c*dv(:,:,:,taum2)) 
  w = w + dt/dsqr*( AB3_a*dw(:,:,:,tau) + AB3_b*dw(:,:,:,taum1) + AB3_c*dw(:,:,:,taum2)) 
  b = b +      dt*( AB3_a*db(:,:,:,tau) + AB3_b*db(:,:,:,taum1) + AB3_c*db(:,:,:,taum2)) 
 endif 
 
 call toc('integrate')

 call tic('boundary')
 call border_exchg_3D(u) 
 call border_exchg_3D(v) 
 call border_exchg_3D(w) 
 call border_exchg_3D(b) 
 call toc('boundary')

 call pressure_solver
 
end subroutine integrate

