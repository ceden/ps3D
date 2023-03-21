

module module_diag_opt_balance
  use main_module
  implicit none
  real*8 ,allocatable :: u_bal(:,:,:),v_bal(:,:,:),w_bal(:,:,:),b_bal(:,:,:)
  complex(p3dfft_type),allocatable :: US(:,:,:),VS(:,:,:),WS(:,:,:),BS(:,:,:),g0_base(:,:,:)
  integer :: n_end,n_end_ave
  real*8 :: Ro_loc, dt_loc
  real*8 :: norm_diff
  real*8 ,allocatable :: u_base(:,:,:),v_base(:,:,:),w_base(:,:,:),b_base(:,:,:)
  real*8 ,allocatable :: u_save(:,:,:),v_save(:,:,:),w_save(:,:,:),b_save(:,:,:)
  real*8 ,allocatable :: u_ave(:,:,:), v_ave(:,:,:), w_ave(:,:,:), b_ave(:,:,:)
 
end module module_diag_opt_balance



subroutine allocate_module_diag_opt_balance
  use main_module
  use module_diag_opt_balance

  allocate( u_bal(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); u_bal=0
  allocate( v_bal(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); v_bal=0
  allocate( w_bal(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); w_bal=0
  allocate( b_bal(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); b_bal=0
  
  allocate( US(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); US=0
  allocate( VS(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); VS=0
  allocate( WS(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); WS=0
  allocate( BS(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); BS=0

  allocate( u_base(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); u_base=0
  allocate( v_base(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); v_base=0
  allocate( w_base(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); w_base=0
  allocate( b_base(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); b_base=0 
 
  if (enable_diag_opt_time_ave) then
   allocate( u_save(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); u_save=0
   allocate( v_save(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); v_save=0
   allocate( w_save(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); w_save=0
   allocate( b_save(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); b_save=0  
  
   allocate( u_ave(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); u_ave=0
   allocate( v_ave(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); v_ave=0
   allocate( w_ave(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); w_ave=0
   allocate( b_ave(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); b_ave=0  
  else
   allocate( g0_base(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); g0_base=0
  endif
end subroutine allocate_module_diag_opt_balance



subroutine diag_opt_balance()
  use main_module
  use module_diag_opt_balance
  implicit none
 
  
  real*8,dimension(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx,3):: du_loc,dv_loc,dw_loc,db_loc
  real*8,dimension(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx):: u_loc, v_loc, w_loc, b_loc, p_loc
  integer :: tau_loc,taum1_loc,taum2_loc,itt_loc
  integer :: n,m
  logical,save  ::  first = .true. 

  if (my_pe==0) print*,' '  
  if (my_pe==0) print*,'Entering optimal balance procedure '  
    
  if (first.and.my_pe==0) then
   print*,'opt_balance_max_Itts = ',opt_balance_max_Itts 
   print*,'opt_balance_period   = ',opt_balance_period/86400., ' d'
   print*,'opt_balance_tol      = ',opt_balance_tol
   print*,'enable_diag_opt_time_ave = ',enable_diag_opt_time_ave
   if (enable_diag_opt_time_ave) &
      print*,'opt_balance_average  = ',opt_balance_average/86400.,' d'
      print*,'opt_balance_average_time  = ',opt_balance_average_times
  endif
    
  ! Save model state
  u_loc = u; v_loc = v; w_loc = w; b_loc = b; p_loc = p
  du_loc = du; dv_loc = dv; dw_loc = dw; db_loc = db
  Ro_loc = Ro; dt_loc = dt
  tau_loc = tau; taum1_loc = taum1; taum2_loc = taum2; itt_loc = itt
  u_bal = u; v_bal = v; w_bal = w; b_bal = b

  ! Calculate number of iteration steps for turning off/on the non linear terms
  n_end = int(opt_balance_period / dt) + 1 
  n_end_ave = int(opt_balance_average/dt) + 1
  norm_diff = 1d0

  ! Save base point
  if (my_pe==0)  print*,'calculating base point'      
  if (enable_diag_opt_time_ave) then
     if (my_pe==0)  print '(a,i5,a,i5,a)',' integrating linear model for ',&
                       n_end_ave*opt_balance_average_times,' timesteps in ', &
                       opt_balance_average_times,' chunks'
   
    do m=1,opt_balance_average_times
     call opt_time_ave(u_base,v_base,w_base,b_base) 
     u = u_base; v = v_base; w = w_base; b = b_base
    enddo
    
    u = u_bal; v = v_bal; w = w_bal; b = b_bal
  else
    call geostrophic_base_point
  endif
  

  ! Start optimal balance iteration
  outer: do n=1,opt_balance_max_Itts
    if (my_pe==0)  print '(a,i5)',' start optimal balance iteration:',n

    ! integrate backward to the linear end
    if (my_pe==0)  print '(a,i5,a,e8.2,a,f5.2,a)',' integrating backward for ', &
       n_end,' timesteps or ', n_end*dt,' s = ',n_end*dt/86400.,' d'
    call opt_backward_integration

    ! project on geostrophic mode at linear end
    if (enable_diag_opt_time_ave) then
     if (my_pe==0)  print '(a,i5,a,i5,a)',' integrating linear model for ',&
                       n_end_ave*opt_balance_average_times,' timesteps in ', &
                       opt_balance_average_times,' chunks'
     
     do m=1,opt_balance_average_times
      call opt_time_ave(u_ave,v_ave,w_ave,b_ave)
      u = u_ave; v = v_ave; w = w_ave; b = b_ave
     enddo
    else
     call normal_mode_projection
    endif
    
    ! integrate forward to the non linear end
    if (my_pe==0)  print '(a,i5,a)',' integrating forward for ',n_end,' timesteps '
    call opt_forward_integration

    ! Apply boundary condition at non linear end
    if (enable_diag_opt_time_ave) then
     u_save = u; v_save = v; w_save = w; b_save = b 
     if (my_pe==0)  print '(a,i5,a,i5,a)',' integrating linear model for ',&
                       n_end_ave*opt_balance_average_times,' timesteps in ', &
                       opt_balance_average_times,' chunks'
     
     do m=1,opt_balance_average_times
      call opt_time_ave(u_ave,v_ave,w_ave,b_ave)
      u = u_ave; v = v_ave; w = w_ave; b = b_ave
     enddo
     
     u = u_save - u_ave + u_base; 
     v = v_save - v_ave + v_base; 
     w = w_save - w_ave + w_base; 
     b = b_save - b_ave + b_base; 
    else
     call exchange_geostrophic_base_point
    endif
    
    call opt_balance_norm()

    u_bal = u; v_bal = v; w_bal = w; b_bal = b

    ! Check tolerance criterion
    if (n > 1) then
      if (my_pe==0)  print '(a,i5,a,ES15.3E3,a,ES15.3E3)', & 
              ' norm of difference to n= ',n - 1,' is ', norm_diff,' / ',opt_balance_tol
      if (norm_diff < opt_balance_tol) exit outer
    end if
  end do outer



  if (my_pe==0) print*,'Finishing optimal balance procedure '
  if (my_pe==0) print*,' ' 
    
  ! Recover model state
  u = u_loc;v = v_loc ; w = w_loc; b = b_loc ; p = p_loc
  du = du_loc;dv = dv_loc ; dw = dw_loc; db = db_loc
  Ro = Ro_loc; dt = dt_loc
  tau = tau_loc; taum1 = taum1_loc; taum2 = taum2_loc; itt = itt_loc
  
  first  = .false.
end subroutine diag_opt_balance



subroutine opt_backward_integration
  use main_module
  use module_diag_opt_balance
  implicit none
  real*8 :: rho, ramp, Ahbi_loc, Avbi_loc, Khbi_loc, Kvbi_loc
  integer :: n

  ! Disable friction for backward integration
  Ahbi_loc = Ahbi; Avbi_loc = Avbi; Khbi_loc = Khbi; Kvbi_loc = Kvbi
  Ahbi = 0; Avbi = 0; Khbi = 0; Kvbi = 0

  ! Use the fact, that in this model, nonlinear terms are scaled with Rossby number
  itt = 0 ; taum2 = 1 ; taum1 = 2 ; tau = 3
  dt = - abs(dt)
  do n=0,n_end-1
    rho = 1. - ramp(n * dt_loc / opt_balance_period)
    Ro = rho * Ro_loc
    call integrate
    tau   = mod(tau,3)+1
    taum1 = mod(taum1,3)+1
    taum2 = mod(taum2,3)+1
  enddo

  ! Recover friction coefficients
  Ahbi = Ahbi_loc; Avbi = Avbi_loc; Khbi = Khbi_loc; Kvbi = Kvbi_loc

end subroutine opt_backward_integration



subroutine opt_forward_integration
  use main_module
  use module_diag_opt_balance
  implicit none
  real*8 :: rho, ramp
  integer :: n

  itt = 0 ; taum2 = 1; taum1 = 2 ; tau = 3
  dt = + abs(dt)
  do n=0,n_end-1
    rho = ramp(n * dt_loc / opt_balance_period)
    Ro = rho * Ro_loc
    call integrate
    tau   = mod(tau,3)+1
    taum1 = mod(taum1,3)+1
    taum2 = mod(taum2,3)+1
  enddo

end subroutine opt_forward_integration


subroutine opt_time_ave(u_loc,v_loc,w_loc,b_loc)
  use main_module
  use module_diag_opt_balance
  implicit none
  real*8, dimension(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) :: u_loc,v_loc,w_loc,b_loc
  integer :: n
  
  itt = 0 ; taum2 = 1; taum1 = 2 ; tau = 3
  dt = + abs(dt)
  Ro = 0. 
  u_loc = 0; v_loc = 0; w_loc = 0; b_loc = 0
  do n=0,n_end_ave-1     
    call integrate
    tau   = mod(tau,3)+1
    taum1 = mod(taum1,3)+1
    taum2 = mod(taum2,3)+1
    u_loc = u_loc + u
    v_loc = v_loc + v
    w_loc = w_loc + w
    b_loc = b_loc + b
  enddo
  u_loc = u_loc/n_end_ave; v_loc = v_loc/n_end_ave 
  w_loc = w_loc/n_end_ave; b_loc = b_loc/n_end_ave

end subroutine opt_time_ave




subroutine geostrophic_base_point
  use main_module
  use module_diag_opt_balance
  implicit none
  integer :: i,j,k
  complex*16 :: qx,qy,qz,qb,px,py,pz,pb

  call forward_transform(u,v,w,b,US,VS,WS,BS)    ! dfft für aktuellen chunk
  do k=fstart(3),fend(3)         ! Loop through spectral array
    do j=fstart(2),fend(2)
     do i=fstart(1),fend(1)
       call vec_q(kx(i),ky(j),kz(k),0d0,qx,qy,qz,qb)
       call vec_p(kx(i),ky(j),kz(k),0d0,px,py,pz,pb) 
       g0_base(i,j,k) = px*US(i,j,k) + py*VS(i,j,k) + pz*WS(i,j,k) + pb*BS(i,j,k)
       US(i,j,k) = g0_base(i,j,k)*qx
       VS(i,j,k) = g0_base(i,j,k)*qy
       WS(i,j,k) = g0_base(i,j,k)*qz
       BS(i,j,k) = g0_base(i,j,k)*qb
     enddo
    enddo
  enddo
  call backward_transform(US,VS,WS,BS,u_base,v_base,w_base,b_base)
end subroutine geostrophic_base_point


subroutine normal_mode_projection
  use main_module
  use module_diag_opt_balance
  implicit none
  integer :: i,j,k
  complex*16 :: qx,qy,qz,qb,px,py,pz,pb,g0

  call forward_transform(u,v,w,b,US,VS,WS,BS)    ! dfft für aktuellen chunk
  do k=fstart(3),fend(3)         ! Loop through spectral array
    do j=fstart(2),fend(2)
     do i=fstart(1),fend(1)
       call vec_q(kx(i),ky(j),kz(k),0d0,qx,qy,qz,qb)
       call vec_p(kx(i),ky(j),kz(k),0d0,px,py,pz,pb)
       g0 = px*US(i,j,k) + py*VS(i,j,k) + pz*WS(i,j,k) + pb*BS(i,j,k)
       US(i,j,k) = g0*qx
       VS(i,j,k) = g0*qy
       WS(i,j,k) = g0*qz
       BS(i,j,k) = g0*qb
     enddo
    enddo
  enddo
  call backward_transform(US,VS,WS,BS,u,v,w,b)
end subroutine normal_mode_projection


subroutine exchange_geostrophic_base_point
  use main_module
  use module_diag_opt_balance
  implicit none
  integer :: i,j,k
  complex*16 :: qx,qy,qz,qb,px,py,pz,pb,g0
  
  call forward_transform(u,v,w,b,US,VS,WS,BS)    ! dfft für aktuellen chunk
  do k=fstart(3),fend(3)         ! Loop through spectral array
    do j=fstart(2),fend(2)
     do i=fstart(1),fend(1)
       call vec_q(kx(i),ky(j),kz(k),0d0,qx,qy,qz,qb)
       call vec_p(kx(i),ky(j),kz(k),0d0,px,py,pz,pb) 
       g0 = px*US(i,j,k) + py*VS(i,j,k) + pz*WS(i,j,k) + pb*BS(i,j,k)
       US(i,j,k) = US(i,j,k)  + (g0_base(i,j,k)-g0)*qx
       VS(i,j,k) = VS(i,j,k)  + (g0_base(i,j,k)-g0)*qy 
       WS(i,j,k) = WS(i,j,k)  + (g0_base(i,j,k)-g0)*qz 
       BS(i,j,k) = BS(i,j,k)  + (g0_base(i,j,k)-g0)*qb 
     enddo
    enddo
  enddo
  call backward_transform(US,VS,WS,BS,u,v,w,b)
end subroutine exchange_geostrophic_base_point



subroutine opt_balance_norm()
  use main_module
  use module_diag_opt_balance
  implicit none
  real*8 :: e1, e2, d
  integer :: i, j, k

  e1=0d0; e2=0d0; d=0d0
  do k=ks_pe,ke_pe
    do j=js_pe,je_pe
      do i=is_pe,ie_pe
        e1 = e1 + u(i,j,k)**2 + v(i,j,k)**2 + dsqr*w(i,j,k)**2 + b(i,j,k)**2/N0**2
        e2 = e2 + u_bal(i,j,k)**2 + v_bal(i,j,k)**2 &
              & + dsqr*w_bal(i,j,k)**2 + b_bal(i,j,k)**2/N0**2
        d = d + (u(i,j,k) - u_bal(i,j,k))**2 + (v(i,j,k) - v_bal(i,j,k))**2 &
            & + dsqr*(w(i,j,k) - w_bal(i,j,k))**2 + (b(i,j,k) - b_bal(i,j,k))**2/N0**2
      end do
    end do
  end do
  call global_sum(e1)
  call global_sum(e2)
  call global_sum(d) 
  norm_diff = 2d0 * sqrt(d) / (sqrt(e1) + sqrt(e2)) 
end subroutine opt_balance_norm




function ramp(theta)
  implicit none 
  real*8, intent(in) :: theta
  real*8 :: ramp
  !real*8 :: pi = 3.14159265358979323846264338327950588
  real*8 :: t1,t2
  !ramp = (1d0-cos(pi*theta))*0.5d0
  t1 = 1./max(1d-32,theta )
  t2 = 1./max(1d-32,1d0-theta )
  ramp = exp(-t1)/(exp(-t1)+exp(-t2) )  
end function ramp





