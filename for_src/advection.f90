
subroutine buoyancy_advection_2nd(db_loc)
!=======================================================================
! Advection of buoyancy with second order 
!=======================================================================
  use main_module   
  implicit none
  integer :: i,j,k
  real(real_type), dimension(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx, &
                             ks_pe-onx:ke_pe+onx) :: db_loc
 
   do k=ks_pe,ke_pe
   do j=js_pe,je_pe
    do i=is_pe-1,ie_pe
      flux_east(i,j,k)=0.5*(b(i,j,k) + b(i+1,j,k) )*u(i,j,k)
    enddo
   enddo
  enddo
  do k=ks_pe,ke_pe
   do j=js_pe-1,je_pe
    do i=is_pe,ie_pe
      flux_north(i,j,k)=0.5*( b(i,j,k) + b(i,j+1,k) )*v(i,j,k)
    enddo
   enddo
  enddo
  do k=ks_pe-1,ke_pe
   do j=js_pe,je_pe
    do i=is_pe,ie_pe
      flux_top(i,j,k)=0.5*( b(i,j,k) + b(i,j,k+1) )*w(i,j,k)
    enddo
   enddo
  enddo

  if (enable_vertical_boundaries) then
   if (my_blk_k == 1      ) flux_top(:,:,1-onx:0) = 0d0
   if (my_blk_k == n_pes_k) flux_top(:,:,nz:nz+onx) = 0d0
  endif
  
  do k=ks_pe,ke_pe
   do j=js_pe,je_pe
    do i=is_pe,ie_pe
      db_loc(i,j,k)= -Ro*(  (flux_east(i,j,k)  -  flux_east(i-1,j,k))/dx &
                           +(flux_north(i,j,k) - flux_north(i,j-1,k))/dy &
                           +(flux_top(i,j,k)   - flux_top(i,j,k-1))/dz )
    enddo
   enddo
  enddo
end subroutine buoyancy_advection_2nd 
 



subroutine momentum_advection_2nd(du_loc,dv_loc,dw_loc)
!=======================================================================
! Advection of momentum with second order which is energy conserving
!=======================================================================
  use main_module   
  implicit none
  integer :: i,j,k
  real(real_type), dimension(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx, &
                             ks_pe-onx:ke_pe+onx) :: du_loc, dv_loc, dw_loc
 !---------------------------------------------------------------------------------
 ! for zonal momentum
 !---------------------------------------------------------------------------------
 do k=ks_pe,ke_pe
  do j=js_pe,je_pe
   do i=is_pe-1,ie_pe 
    flux_east(i,j,k) = 0.25*(u(i,j,k)+u(i+1,j,k))*(u(i+1,j,k)+u(i,j,k))
   enddo
  enddo
 enddo
 do k=ks_pe,ke_pe
  do j=js_pe-1,je_pe
   do i=is_pe,ie_pe
     flux_north(i,j,k) = 0.25*(u(i,j,k)+u(i,j+1,k))*(v(i+1,j,k)+v(i,j,k))
   enddo
  enddo
 enddo
 do k=ks_pe-1,ke_pe
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
      flux_top(i,j,k) = 0.25*(u(i,j,k+1)+u(i,j,k))*(w(i,j,k)+w(i+1,j,k))
   enddo
  enddo
 enddo

 if (enable_vertical_boundaries) then
   if (my_blk_k == 1      ) flux_top(:,:,1-onx:0)   = 0d0
   if (my_blk_k == n_pes_k) flux_top(:,:,nz:nz+onx) = 0d0
 endif

 do k=ks_pe,ke_pe
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
    du_loc(i,j,k) = -Ro*( ( flux_east(i,j,k) - flux_east(i-1,j,k) )/dx &
                          +(flux_north(i,j,k) - flux_north(i,j-1,k) )/dy &
                          +(  flux_top(i,j,k) - flux_top  (i,j,k-1) )/dz )
   enddo
  enddo
 enddo
 !---------------------------------------------------------------------------------
 ! for meridional momentum
 !---------------------------------------------------------------------------------
 do k=ks_pe,ke_pe
  do j=js_pe,je_pe
   do i=is_pe-1,ie_pe
    flux_east(i,j,k) = 0.25*(v(i,j,k)+v(i+1,j,k))*(u(i,j+1,k)+u(i,j,k))
   enddo
  enddo
 enddo
 do k=ks_pe,ke_pe
  do j=js_pe-1,je_pe
   do i=is_pe,ie_pe
     flux_north(i,j,k) = 0.25*(v(i,j,k)+v(i,j+1,k))*(v(i,j+1,k)+v(i,j,k))
   enddo
  enddo
 enddo
 do k=ks_pe-1,ke_pe
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
      flux_top(i,j,k) = 0.25*(v(i,j,k+1)+v(i,j,k))*(w(i,j,k)+w(i,j+1,k))
   enddo
  enddo
 enddo

 if (enable_vertical_boundaries) then
   if (my_blk_k == 1      ) flux_top(:,:,1-onx:0)   = 0d0
   if (my_blk_k == n_pes_k) flux_top(:,:,nz:nz+onx) = 0d0
 endif


 do k=ks_pe,ke_pe
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
    dv_loc(i,j,k) = -Ro*( ( flux_east(i,j,k) - flux_east(i-1,j,k) )/dx  &
                          +(flux_north(i,j,k) - flux_north(i,j-1,k) )/dy  &
                          +(flux_top(i,j,k)   - flux_top(i,j,k-1) )/dz )
   enddo
  enddo
 enddo

 !---------------------------------------------------------------------------------
 ! for vertical momentum
 !---------------------------------------------------------------------------------
 do k=ks_pe,ke_pe
  do j=js_pe,je_pe
   do i=is_pe-1,ie_pe
    flux_east(i,j,k) = 0.25*(w(i,j,k)+w(i+1,j,k))*(u(i,j,k)+u(i,j,k+1))
   enddo
  enddo
 enddo
 do k=ks_pe,ke_pe
  do j=js_pe-1,je_pe
   do i=is_pe,ie_pe
     flux_north(i,j,k) = 0.25*(w(i,j,k)+w(i,j+1,k))*(v(i,j,k)+v(i,j,k+1))
   enddo
  enddo
 enddo
 do k=ks_pe-1,ke_pe
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
      flux_top(i,j,k) = 0.25*(w(i,j,k+1)+w(i,j,k))*(w(i,j,k)+w(i,j,k+1))
   enddo
  enddo
 enddo
 
 if (enable_vertical_boundaries) then
   if (my_blk_k == 1      ) flux_top(:,:,1-onx:0)     = 0d0
   if (my_blk_k == n_pes_k) flux_top(:,:,nz-1:nz+onx) = 0d0
 endif


 do k=ks_pe,ke_pe
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
     dw_loc(i,j,k) = -dsqr*Ro*( ( flux_east(i,j,k) - flux_east(i-1,j,k) )/dx &
                               + (flux_north(i,j,k) - flux_north(i,j-1,k) )/dy &
                               + (flux_top(i,j,k)   - flux_top(i,j,k-1) )/dz )
   enddo
  enddo
 enddo


end subroutine momentum_advection_2nd







