


subroutine biharmonic(p1,Hmix,Vmix)
!---------------------------------------------------------------------------------
! biharmonic mixing/friction
!---------------------------------------------------------------------------------
 use main_module   
 implicit none
 integer :: i,j,k
 real(real_type) :: p1(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx),Hmix,Vmix
 real(real_type) :: del2(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx)

 do k=ks_pe-1,ke_pe+1
  do j=js_pe-1,je_pe+1
   do i=is_pe-2,ie_pe+1
    flux_east(i,j,k) = -Hmix*(p1(i+1,j,k)-p1(i,j,k))/dx
   enddo
  enddo
 enddo

 do k=ks_pe-1,ke_pe+1
  do j=js_pe-2,je_pe+1
   do i=is_pe-1,ie_pe+1
    flux_north(i,j,k) = -Hmix*(p1(i,j+1,k)-p1(i,j,k))/dy
   enddo
  enddo
 enddo

 do k=ks_pe-2,ke_pe+1
  do j=js_pe-1,je_pe+1
   do i=is_pe-1,ie_pe+1
    flux_top(i,j,k) = -Vmix*(p1(i,j,k+1)-p1(i,j,k))/dz
   enddo
  enddo
 enddo
 if (enable_vertical_boundaries) then
   if (my_blk_k == 1      ) flux_top(:,:,1-onx:0)   = 0d0
   if (my_blk_k == n_pes_k) flux_top(:,:,nz:nz+onx) = 0d0
 endif

 do k=ks_pe-1,ke_pe+1
  do j=js_pe-1,je_pe+1
   do i=is_pe-1,ie_pe+1
    del2(i,j,k)= (flux_east(i,j,k) - flux_east(i-1,j,k))/dx &
                +(flux_north(i,j,k)- flux_north(i,j-1,k))/dy &
                +(flux_top(i,j,k)- flux_top(i,j,k-1))/dz 
   enddo
  enddo
 enddo

 do k=ks_pe,ke_pe
  do j=js_pe,je_pe
   do i=is_pe-1,ie_pe
    flux_east(i,j,k) = (del2(i+1,j,k)-del2(i,j,k))/dx
   enddo
  enddo
 enddo

 do k=ks_pe,ke_pe
  do j=js_pe-1,je_pe
   do i=is_pe,ie_pe
    flux_north(i,j,k) = (del2(i,j+1,k)-del2(i,j,k))/dy
   enddo
  enddo
 enddo

 do k=ks_pe-1,ke_pe
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
    flux_top(i,j,k) = (del2(i,j,k+1)-del2(i,j,k))/dz
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
      p1(i,j,k)=p1(i,j,k)+dt*( &
             (flux_east(i,j,k)  - flux_east (i-1,j,k))/dx &
            +(flux_north(i,j,k) - flux_north(i,j-1,k))/dy &
            +(flux_top(i,j,k)   - flux_top  (i,j,k-1))/dz   )
   enddo
  enddo
 enddo
end subroutine biharmonic





subroutine biharmonic_h(p1,Hmix)
!---------------------------------------------------------------------------------
! horizontal biharmonic mixing/friction
!---------------------------------------------------------------------------------
 use main_module   
 implicit none
 integer :: i,j,k
 real(real_type) :: p1(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx),Hmix
 real(real_type) :: del2(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx)

 do k=ks_pe-1,ke_pe+1
  do j=js_pe-1,je_pe+1
   do i=is_pe-2,ie_pe+1
    flux_east(i,j,k) = -Hmix*(p1(i+1,j,k)-p1(i,j,k))/dx
   enddo
  enddo
 enddo

 do k=ks_pe-1,ke_pe+1
  do j=js_pe-2,je_pe+1
   do i=is_pe-1,ie_pe+1
    flux_north(i,j,k) = -Hmix*(p1(i,j+1,k)-p1(i,j,k))/dy
   enddo
  enddo
 enddo


 do k=ks_pe-1,ke_pe+1
  do j=js_pe-1,je_pe+1
   do i=is_pe-1,ie_pe+1
    del2(i,j,k)= (flux_east(i,j,k) - flux_east(i-1,j,k))/dx &
                +(flux_north(i,j,k)- flux_north(i,j-1,k))/dy 
   enddo
  enddo
 enddo

 do k=ks_pe,ke_pe
  do j=js_pe,je_pe
   do i=is_pe-1,ie_pe
    flux_east(i,j,k) = (del2(i+1,j,k)-del2(i,j,k))/dx
   enddo
  enddo
 enddo

 do k=ks_pe,ke_pe
  do j=js_pe-1,je_pe
   do i=is_pe,ie_pe
    flux_north(i,j,k) = (del2(i,j+1,k)-del2(i,j,k))/dy
   enddo
  enddo
 enddo


 do k=ks_pe,ke_pe
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
      p1(i,j,k)=p1(i,j,k)+dt*( &
             (flux_east(i,j,k)  - flux_east (i-1,j,k))/dx &
            +(flux_north(i,j,k) - flux_north(i,j-1,k))/dy  )
   enddo
  enddo
 enddo
end subroutine biharmonic_h





subroutine biharmonic_v(p1,Vmix)
!---------------------------------------------------------------------------------
! vertical biharmonic mixing/friction
!---------------------------------------------------------------------------------
 use main_module   
 implicit none
 integer :: i,j,k
 real(real_type) :: p1(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx),Vmix
 real(real_type) :: del2(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx)


 do k=ks_pe-2,ke_pe+1
  do j=js_pe-1,je_pe+1
   do i=is_pe-1,ie_pe+1
    flux_top(i,j,k) = -Vmix*(p1(i,j,k+1)-p1(i,j,k))/dz
   enddo
  enddo
 enddo
 if (enable_vertical_boundaries) then
   if (my_blk_k == 1      ) flux_top(:,:,1-onx:0)   = 0d0
   if (my_blk_k == n_pes_k) flux_top(:,:,nz:nz+onx) = 0d0
 endif

 do k=ks_pe-1,ke_pe+1
  do j=js_pe-1,je_pe+1
   do i=is_pe-1,ie_pe+1
    del2(i,j,k)= (flux_top(i,j,k)- flux_top(i,j,k-1))/dz 
   enddo
  enddo
 enddo

 do k=ks_pe-1,ke_pe
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
    flux_top(i,j,k) = (del2(i,j,k+1)-del2(i,j,k))/dz
   enddo
  enddo
 enddo
 
 if (enable_vertical_boundaries) then
   if (my_blk_k == 1      ) flux_top(:,:,1-onx:0)   = 0d0
   if (my_blk_k == n_pes_k) flux_top(:,:,nz:nz+onx) = 0d0
   
   !flux_top = 0.
 endif

 do k=ks_pe,ke_pe
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
      p1(i,j,k)=p1(i,j,k)+dt*( (flux_top(i,j,k)   - flux_top  (i,j,k-1))/dz   )
   enddo
  enddo
 enddo
end subroutine biharmonic_v







subroutine harmonic(p1,Hmix,Vmix)
!---------------------------------------------------------------------------------
! harmonic mixing/friction, time level taup1 of p1 is updated 
!---------------------------------------------------------------------------------
 use main_module   
 implicit none
 integer :: i,j,k
 real(real_type) :: p1(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx),Hmix,Vmix


 do k=ks_pe-1,ke_pe+1
  do j=js_pe-1,je_pe+1
   do i=is_pe-2,ie_pe+1
    flux_east(i,j,k) = Hmix*(p1(i+1,j,k)-p1(i,j,k))/dx
   enddo
  enddo
 enddo

 do k=ks_pe-1,ke_pe+1
  do j=js_pe-2,je_pe+1
   do i=is_pe-1,ie_pe+1
    flux_north(i,j,k) = Hmix*(p1(i,j+1,k)-p1(i,j,k))/dy
   enddo
  enddo
 enddo

 do k=ks_pe-2,ke_pe+1
  do j=js_pe-1,je_pe+1
   do i=is_pe-1,ie_pe+1
    flux_top(i,j,k) = Vmix*(p1(i,j,k+1)-p1(i,j,k))/dz
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
      p1(i,j,k)=p1(i,j,k)+dt*( &
             (flux_east(i,j,k)  - flux_east (i-1,j,k))/dx &
            +(flux_north(i,j,k) - flux_north(i,j-1,k))/dy &
            +(flux_top(i,j,k)   - flux_top  (i,j,k-1))/dz   )
   enddo
  enddo
 enddo
end subroutine harmonic


