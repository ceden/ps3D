

module module_diag_balance
 use main_module
 implicit none
 complex(p3dfft_type),allocatable :: US(:,:,:),VS(:,:,:),WS(:,:,:),BS(:,:,:)
 complex*16, allocatable :: f1p(:,:,:),f1m(:,:,:)
 complex*16, allocatable :: I_00(:,:,:),I_0p(:,:,:),I_0m(:,:,:)
 complex*16, allocatable :: I_f1p(:,:,:),I_f1m(:,:,:)
 complex*16, allocatable :: I_0f1p(:,:,:),I_0f1m(:,:,:)
 complex*16, allocatable :: I_pT0p(:,:,:),I_pT0m(:,:,:)
 real*8 ,allocatable :: u0(:,:,:),v0(:,:,:),w0(:,:,:),b0(:,:,:)
 real*8 ,allocatable :: pTu0(:,:,:),pTv0(:,:,:),pTw0(:,:,:),pTb0(:,:,:)
 real*8 ,allocatable :: u1(:,:,:),v1(:,:,:),w1(:,:,:),b1(:,:,:)
 real*8 ,allocatable :: u2(:,:,:),v2(:,:,:),w2(:,:,:),b2(:,:,:)
 real*8 ,allocatable :: u3(:,:,:),v3(:,:,:),w3(:,:,:),b3(:,:,:)
 
 complex*16, allocatable :: f2p(:,:,:),f2m(:,:,:)
 complex*16, allocatable :: I_f2p(:,:,:),I_f2m(:,:,:)
 complex*16, allocatable :: I_0f2p(:,:,:),I_0f2m(:,:,:)
 complex*16, allocatable :: pTf2p(:,:,:),pTf2m(:,:,:)

 character*80 :: balance_name
 
end module module_diag_balance 



subroutine allocate_module_diag_balance
 use main_module
 use module_diag_balance
 allocate( US(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); US=0
 allocate( VS(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); VS=0
 allocate( WS(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); WS=0
 allocate( BS(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); BS=0
 
 allocate( f1p(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); f1p=0
 allocate( f1m(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); f1m=0
 
 allocate( I_00(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); I_00=0
 allocate( I_0m(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); I_0m=0
 allocate( I_0p(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); I_0p=0
 allocate( I_f1p(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); I_f1p=0
 allocate( I_f1m(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); I_f1m=0
 
 allocate( I_0f1p(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); I_0f1p=0
 allocate( I_0f1m(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); I_0f1m=0
 
 allocate( I_pT0p(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); I_pT0p=0
 allocate( I_pT0m(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); I_pT0m=0
 
 allocate( U0(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); U0=0
 allocate( V0(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); V0=0
 allocate( W0(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); W0=0
 allocate( B0(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); B0=0
 
 allocate( U1(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); U1=0
 allocate( V1(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); V1=0
 allocate( W1(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); W1=0
 allocate( B1(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); B1=0
 
 allocate( U2(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); U2=0
 allocate( V2(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); V2=0
 allocate( W2(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); W2=0
 allocate( B2(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); B2=0
 
 allocate( U3(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); U3=0
 allocate( V3(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); V3=0
 allocate( W3(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); W3=0
 allocate( B3(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); B3=0
 
 allocate( pTU0(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); pTU0=0
 allocate( pTV0(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); pTV0=0
 allocate( pTW0(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); pTW0=0
 allocate( pTB0(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx) ); pTB0=0

 allocate( f2p(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); f2p=0
 allocate( f2m(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); f2m=0
 allocate( pTf2p(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); pTf2p=0
 allocate( pTf2m(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); pTf2m=0
 
 allocate( I_f2p(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); I_f2p=0
 allocate( I_f2m(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); I_f2m=0
 allocate( I_0f2p(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); I_0f2p=0
 allocate( I_0f2m(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) ); I_0f2m=0
 
end subroutine allocate_module_diag_balance


subroutine diag_balance
 use main_module
 call diag_balance_zero_order
 call diag_balance_first_order
 call diag_balance_second_order
 call diag_balance_third_order
 
 if (enable_diag_balance_chunks) then
   call write_diag_balance_chunks
 else
   call write_diag_balance
 endif 
end subroutine diag_balance


subroutine diag_balance_zero_order
 use main_module
 use module_diag_balance
 implicit none
 integer :: i,j,k
 complex*16 :: qx,qy,qz,qb,px,py,pz,pb,g0
 complex(p3dfft_type),parameter  :: im = (0d0,1d0)
 
 call forward_transform(u,v,w,b,US,VS,WS,BS)
 do k=fstart(3),fend(3)
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
 call backward_transform(US,VS,WS,BS,u0,v0,w0,b0)
end subroutine diag_balance_zero_order


subroutine diag_balance_first_order
 use main_module
 use module_diag_balance
 implicit none
 integer :: i,j,k
 complex*16 :: px0,py0,pz0,pb0
 complex*16 :: qxm,qym,qzm,qbm,pxm,pym,pzm,pbm
 complex*16 :: qxp,qyp,qzp,qbp,pxp,pyp,pzp,pbp
 real*8 :: om2,om
 complex(p3dfft_type),parameter  :: im = (0d0,1d0)
 
 call model_nonlin(u0,v0,w0,b0,US,VS,WS,BS)
 
 do k=fstart(3),fend(3)
  do j=fstart(2),fend(2)
   do i=fstart(1),fend(1)            
      call vec_p(kx(i),ky(j),kz(k), 0d0,px0,py0,pz0,pb0)
      call vec_p(kx(i),ky(j),kz(k),+1d0,pxp,pyp,pzp,pbp)
      call vec_p(kx(i),ky(j),kz(k),-1d0,pxm,pym,pzm,pbm)
      call vec_q(kx(i),ky(j),kz(k),+1d0,qxp,qyp,qzp,qbp)
      call vec_q(kx(i),ky(j),kz(k),-1d0,qxm,qym,qzm,qbm)
      call omega2(kx(i),ky(j),kz(k),om2)
      om = sqrt(om2)  
      I_00(i,j,k) = im*(px0*US(i,j,k) + py0*VS(i,j,k) + pz0*WS(i,j,k) + pb0*BS(i,j,k))
      I_0p(i,j,k) = im*(pxp*US(i,j,k) + pyp*VS(i,j,k) + pzp*WS(i,j,k) + pbp*BS(i,j,k) )
      I_0m(i,j,k) = im*(pxm*US(i,j,k) + pym*VS(i,j,k) + pzm*WS(i,j,k) + pbm*BS(i,j,k))        
      if (om>0.) then
       f1p(i,j,k) =  I_0p(i,j,k)/om
       f1m(i,j,k) = -I_0m(i,j,k)/om
      else
       f1p(i,j,k) = 0.
       f1m(i,j,k) = 0.
      endif          
      US(i,j,k) = Ro*( f1p(i,j,k)*qxp + f1m(i,j,k)*qxm)
      VS(i,j,k) = Ro*( f1p(i,j,k)*qyp + f1m(i,j,k)*qym)
      WS(i,j,k) = Ro*( f1p(i,j,k)*qzp + f1m(i,j,k)*qzm)
      BS(i,j,k) = Ro*( f1p(i,j,k)*qbp + f1m(i,j,k)*qbm)       
   enddo
  enddo 
 enddo 
 call backward_transform(US,VS,WS,BS,u1,v1,w1,b1)
end subroutine diag_balance_first_order



subroutine diag_balance_second_order
 use main_module
 use module_diag_balance
 implicit none
 integer :: i,j,k
 complex*16 :: qx0,qy0,qz0,qb0
 complex*16 :: qxm,qym,qzm,qbm,pxm,pym,pzm,pbm,pTf1p,pTf1m
 complex*16 :: qxp,qyp,qzp,qbp,pxp,pyp,pzp,pbp,pTg0,I_0pT0p,I_0pT0m
 real*8 :: om2,om
 complex(p3dfft_type),parameter  :: im = (0d0,1d0)

 call model_nonlin(u1/Ro,v1/Ro,w1/Ro,b1/Ro,US,VS,WS,BS)
 
 do k=fstart(3),fend(3)
  do j=fstart(2),fend(2)
   do i=fstart(1),fend(1)            
      call vec_p(kx(i),ky(j),kz(k),+1d0,pxp,pyp,pzp,pbp)
      call vec_p(kx(i),ky(j),kz(k),-1d0,pxm,pym,pzm,pbm)     
      I_f1p(i,j,k) = im*(pxp*US(i,j,k) + pyp*VS(i,j,k) + pzp*WS(i,j,k) + pbp*BS(i,j,k) )
      I_f1m(i,j,k) = im*(pxm*US(i,j,k) + pym*VS(i,j,k) + pzm*WS(i,j,k) + pbm*BS(i,j,k))      
   enddo
  enddo 
 enddo 
 
 call model_nonlin(u0+u1/Ro,v0+v1/Ro,w0+w1/Ro,b0+b1/Ro,US,VS,WS,BS)
 
 do k=fstart(3),fend(3)
  do j=fstart(2),fend(2)
   do i=fstart(1),fend(1)  
      call vec_p(kx(i),ky(j),kz(k),+1d0,pxp,pyp,pzp,pbp)
      call vec_p(kx(i),ky(j),kz(k),-1d0,pxm,pym,pzm,pbm)
      call vec_q(kx(i),ky(j),kz(k), 0d0,qx0,qy0,qz0,qb0)
      I_0f1p(i,j,k) = im*(pxp*US(i,j,k) + pyp*VS(i,j,k) + pzp*WS(i,j,k) + pbp*BS(i,j,k) )
      I_0f1m(i,j,k) = im*(pxm*US(i,j,k) + pym*VS(i,j,k) + pzm*WS(i,j,k) + pbm*BS(i,j,k))     
      pTg0 = -im*I_00(i,j,k)
      US(i,j,k) = pTg0*qx0
      VS(i,j,k) = pTg0*qy0
      WS(i,j,k) = pTg0*qz0
      BS(i,j,k) = pTg0*qb0       
   enddo
  enddo 
 enddo 
 
 call backward_transform(US,VS,WS,BS,pTu0,pTv0,pTw0,pTb0)   

 call model_nonlin(pTu0,pTv0,pTw0,pTb0,US,VS,WS,BS)

 do k=fstart(3),fend(3)
  do j=fstart(2),fend(2)
   do i=fstart(1),fend(1)  
      call vec_p(kx(i),ky(j),kz(k),+1d0,pxp,pyp,pzp,pbp)
      call vec_p(kx(i),ky(j),kz(k),-1d0,pxm,pym,pzm,pbm)
      I_pT0p(i,j,k) = im*(pxp*US(i,j,k) + pyp*VS(i,j,k) + pzp*WS(i,j,k) + pbp*BS(i,j,k) )
      I_pT0m(i,j,k) = im*(pxm*US(i,j,k) + pym*VS(i,j,k) + pzm*WS(i,j,k) + pbm*BS(i,j,k))
   enddo
  enddo 
 enddo 
 
 call model_nonlin(u0+pTu0,v0+pTv0,w0+pTw0,b0+pTb0,US,VS,WS,BS)

 do k=fstart(3),fend(3)
  do j=fstart(2),fend(2)
   do i=fstart(1),fend(1)  
      call vec_p(kx(i),ky(j),kz(k),+1d0,pxp,pyp,pzp,pbp)
      call vec_p(kx(i),ky(j),kz(k),-1d0,pxm,pym,pzm,pbm)  
      call vec_q(kx(i),ky(j),kz(k),+1d0,qxp,qyp,qzp,qbp)
      call vec_q(kx(i),ky(j),kz(k),-1d0,qxm,qym,qzm,qbm)
      call omega2(kx(i),ky(j),kz(k),om2)
      om = sqrt(om2)
      
      I_0pT0p = im*(pxp*US(i,j,k) + pyp*VS(i,j,k) + pzp*WS(i,j,k) + pbp*BS(i,j,k) )
      I_0pT0m = im*(pxm*US(i,j,k) + pym*VS(i,j,k) + pzm*WS(i,j,k) + pbm*BS(i,j,k))      
      if (om>0) then
       pTf1p = (I_0pT0p-I_0p(i,j,k)-I_pT0p(i,j,k))/om
       pTf1m =-(I_0pT0m-I_0m(i,j,k)-I_pT0m(i,j,k))/om
       f2p(i,j,k) =  ( I_0f1p(i,j,k) - I_0p(i,j,k) -I_f1p(i,j,k) - im*pTf1p )/om
       f2m(i,j,k) = -( I_0f1m(i,j,k) - I_0m(i,j,k) -I_f1m(i,j,k) - im*pTf1m )/om 
      else
       f2p(i,j,k)=0.
       f2m(i,j,k)=0.
      endif 
           
      US(i,j,k) = Ro**2*( f2p(i,j,k)*qxp + f2m(i,j,k)*qxm)
      VS(i,j,k) = Ro**2*( f2p(i,j,k)*qyp + f2m(i,j,k)*qym)
      WS(i,j,k) = Ro**2*( f2p(i,j,k)*qzp + f2m(i,j,k)*qzm)
      BS(i,j,k) = Ro**2*( f2p(i,j,k)*qbp + f2m(i,j,k)*qbm)     
   enddo
  enddo 
 enddo 

 call backward_transform(US,VS,WS,BS,u2,v2,w2,b2)

end subroutine diag_balance_second_order 


subroutine diag_balance_third_order
 use main_module
 use module_diag_balance
 implicit none
 integer :: i,j,k
 complex*16 :: qx0,qy0,qz0,qb0,px0,py0,pz0,pb0,f3p,f3m,g0
 complex*16 :: qxm,qym,qzm,qbm,pxm,pym,pzm,pbm
 complex*16 :: qxp,qyp,qzp,qbp,pxp,pyp,pzp,pbp
 real*8 :: om2,om
 complex(p3dfft_type),parameter  :: im = (0d0,1d0)
    
  call model_nonlin(u2/Ro**2,v2/Ro**2,w2/Ro**2,b2/Ro**2,US,VS,WS,BS)
  do k=fstart(3),fend(3)
   do j=fstart(2),fend(2)
    do i=fstart(1),fend(1)            
      call vec_p(kx(i),ky(j),kz(k),+1d0,pxp,pyp,pzp,pbp)
      call vec_p(kx(i),ky(j),kz(k),-1d0,pxm,pym,pzm,pbm)     
      I_f2p(i,j,k) = im*(pxp*US(i,j,k) + pyp*VS(i,j,k) + pzp*WS(i,j,k) + pbp*BS(i,j,k) )
      I_f2m(i,j,k) = im*(pxm*US(i,j,k) + pym*VS(i,j,k) + pzm*WS(i,j,k) + pbm*BS(i,j,k))      
    enddo
   enddo 
  enddo  
    
  call model_nonlin(u0+u2/Ro**2,v0+v2/Ro**2,w0+w2/Ro**2,b0+b2/Ro**2,US,VS,WS,BS)
  do k=fstart(3),fend(3)
   do j=fstart(2),fend(2)
    do i=fstart(1),fend(1)            
      call vec_p(kx(i),ky(j),kz(k),+1d0,pxp,pyp,pzp,pbp)
      call vec_p(kx(i),ky(j),kz(k),-1d0,pxm,pym,pzm,pbm)     
      I_0f2p(i,j,k) = im*(pxp*US(i,j,k) + pyp*VS(i,j,k) + pzp*WS(i,j,k) + pbp*BS(i,j,k) )
      I_0f2m(i,j,k) = im*(pxm*US(i,j,k) + pym*VS(i,j,k) + pzm*WS(i,j,k) + pbm*BS(i,j,k))      
    enddo
   enddo 
  enddo 

  call model_run(u0+u1+u2,v0+v1+v2,w0+w1+w2,b0+b1+b2,u3,v3,w3,b3,10)
    
  call forward_transform(u3,v3,w3,b3,US,VS,WS,BS)
  do k=fstart(3),fend(3)
   do j=fstart(2),fend(2)
    do i=fstart(1),fend(1)  
      call vec_q(kx(i),ky(j),kz(k),0d0,qx0,qy0,qz0,qb0)
      call vec_p(kx(i),ky(j),kz(k),0d0,px0,py0,pz0,pb0)
      g0 = px0*US(i,j,k) + py0*VS(i,j,k) + pz0*WS(i,j,k) + pb0*BS(i,j,k)
      US(i,j,k) = g0*qx0
      VS(i,j,k) = g0*qy0
      WS(i,j,k) = g0*qz0
      BS(i,j,k) = g0*qb0      
    enddo
   enddo
 enddo 
 call backward_transform(US,VS,WS,BS,u0,v0,w0,b0)
 
 pTf2p = f2p; pTf2m = f2m
 call diag_balance_first_order
 call diag_balance_second_order
  
 pTf2p = (f2p-pTf2p)/(10.*dt)/Ro
 pTf2m = (f2m-pTf2m)/(10.*dt)/Ro
 
 call diag_balance_zero_order
 call diag_balance_first_order
 call diag_balance_second_order  
 
 do k=fstart(3),fend(3)
  do j=fstart(2),fend(2)
   do i=fstart(1),fend(1)  
      call vec_q(kx(i),ky(j),kz(k),+1d0,qxp,qyp,qzp,qbp)
      call vec_q(kx(i),ky(j),kz(k),-1d0,qxm,qym,qzm,qbm)
      call omega2(kx(i),ky(j),kz(k),om2)
      om = sqrt(om2)
      if (om>0) then
       f3p =  (I_0f2p(i,j,k) + I_f1p(i,j,k) - I_f2p(i,j,k) - I_0p(i,j,k) - im*pTf2p(i,j,k) )/om
       f3m = -(I_0f2m(i,j,k) + I_f1m(i,j,k) - I_f2m(i,j,k) - I_0m(i,j,k) - im*pTf2m(i,j,k) )/om 
      else
       f3p=0.
       f3m=0.
      endif       
      US(i,j,k) = Ro**3*( f3p*qxp + f3m*qxm)
      VS(i,j,k) = Ro**3*( f3p*qyp + f3m*qym)
      WS(i,j,k) = Ro**3*( f3p*qzp + f3m*qzm)
      BS(i,j,k) = Ro**3*( f3p*qbp + f3m*qbm)     
   enddo
  enddo 
 enddo 

 call backward_transform(US,VS,WS,BS,u3,v3,w3,b3)
end subroutine diag_balance_third_order



subroutine vec_q(kx_loc,ky_loc,kz_loc,s,qx,qy,qz,qb)
 use main_module
 implicit none
 real*8, intent(in) :: kx_loc,ky_loc,kz_loc,s
 complex*16, intent(out) :: qx,qy,qz,qb
 if (enable_vertical_boundaries) then
       call vec_q_vbc(kx_loc,ky_loc,kz_loc,dx,dy,dz,s,f0,N0,dsqr,qx,qy,qz,qb)
 else
       call vec_q_discrete(kx_loc,ky_loc,kz_loc,dx,dy,dz,s,f0,N0,dsqr,qx,qy,qz,qb)
 endif 
end subroutine vec_q


subroutine vec_p(kx_loc,ky_loc,kz_loc,s,qx,qy,qz,qb)
 use main_module
 implicit none
 real*8, intent(in) :: kx_loc,ky_loc,kz_loc,s
 complex*16, intent(out) :: qx,qy,qz,qb
 if (enable_vertical_boundaries) then
       call vec_p_vbc(kx_loc,ky_loc,kz_loc,dx,dy,dz,s,f0,N0,dsqr,qx,qy,qz,qb)
 else
       call vec_p_discrete(kx_loc,ky_loc,kz_loc,dx,dy,dz,s,f0,N0,dsqr,qx,qy,qz,qb)
 endif 
end subroutine vec_p


subroutine omega2(kx_loc,ky_loc,kz_loc,om2)
 use main_module
 implicit none
 real*8, intent(in) :: kx_loc,ky_loc,kz_loc
 real*8, intent(out) :: om2

 if (enable_vertical_boundaries) then
       call omega2_vbc(kx_loc,ky_loc,kz_loc,dx,dy,dz,f0,N0,dsqr,om2)
 else
       call omega2_discrete(kx_loc,ky_loc,kz_loc,dx,dy,dz,f0,N0,dsqr,om2)
 endif 
end subroutine omega2




subroutine coeff_p(kx_loc,ky_loc,kz_loc,s,p_loc) 
 use main_module
 implicit none
 real*8, intent(in) :: kx_loc,ky_loc,kz_loc,s
 real*8, intent(out) :: p_loc
  if (enable_vertical_boundaries) then
      call coeff_p_vbc(kx_loc,ky_loc,kz_loc,dx,dy,dz,s,f0,N0,dsqr,p_loc)
  else
      call coeff_p_discrete(kx_loc,ky_loc,kz_loc,dx,dy,dz,s,f0,N0,dsqr,p_loc)      
  endif 
end subroutine coeff_p


subroutine forward_transform(u_in,v_in,w_in,b_in,u_out,v_out,w_out,b_out)
 use main_module
 implicit none
 real*8,dimension(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx)::u_in,v_in,w_in,b_in
 complex*16,dimension(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3))::u_out,v_out,w_out,b_out  

 if (enable_vertical_boundaries) then
  call my_fft_vert_mode (u_in(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe),u_out,phi)
  call my_fft_vert_mode (v_in(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe),v_out,phi)
  call my_fft_vert_mode (w_in(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe),w_out,psi)
  call my_fft_vert_mode (b_in(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe),b_out,psib)
 else
  call my_fft_3D(u_in(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe),u_out)
  call my_fft_3D(v_in(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe),v_out)
  call my_fft_3D(w_in(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe),w_out)
  call my_fft_3D(b_in(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe),b_out)
 endif
end subroutine forward_transform


subroutine backward_transform(u_in,v_in,w_in,b_in,u_out,v_out,w_out,b_out)
 use main_module
 implicit none
 complex*16,dimension(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3))::u_in,v_in,w_in,b_in
 real*8,dimension(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx)::u_out,v_out,w_out,b_out  

 if (enable_vertical_boundaries) then
  call my_ifft_vert_mode(u_in,u_out(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe),phi)
  call my_ifft_vert_mode(v_in,v_out(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe),phi)
  call my_ifft_vert_mode(w_in,w_out(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe),psi) 
  call my_ifft_vert_mode(b_in,b_out(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe),psib)
 else
  call my_ifft_3D(u_in,u_out(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))
  call my_ifft_3D(v_in,v_out(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)) 
  call my_ifft_3D(w_in,w_out(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)) 
  call my_ifft_3D(b_in,b_out(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))  
 endif
 call border_exchg_3D(u_out) 
 call border_exchg_3D(v_out) 
 call border_exchg_3D(w_out) 
 call border_exchg_3D(b_out) 
end subroutine backward_transform


subroutine model_run(u_in,v_in,w_in,b_in,u_out,v_out,w_out,b_out,max_steps)
 use main_module
 implicit none
 real*8,dimension(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx):: u_in,v_in
 real*8,dimension(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx):: w_in,b_in
 real*8,dimension(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx):: u_out,v_out
 real*8,dimension(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx):: w_out,b_out
 integer :: n,max_steps,tau_loc,taum1_loc,taum2_loc,itt_loc
 real*8,dimension(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx):: u_loc,v_loc
 real*8,dimension(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx):: w_loc,b_loc,p_loc
 real*8,dimension(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx,2):: du_loc,dv_loc
 real*8,dimension(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx,2):: dw_loc,db_loc
 
 ! save model state
 u_loc = u; v_loc = v ; w_loc = w; b_loc = b; p_loc = p
 du_loc = du; dv_loc = dv ; dw_loc = dw; db_loc = db
 tau_loc = tau; taum1_loc = taum1; taum2_loc=taum2; itt_loc = itt
  
 u = u_in; v = v_in; w = w_in; b = b_in
 itt=0; taum2 = 1; taum1 = 2 ; tau = 3
 do n=1,max_steps
   call integrate     
   tau   = mod(tau,3)+1
   taum1 = mod(taum1,3)+1 
   taum2 = mod(taum2,3)+1 
 enddo
 u_out = u; v_out = v; w_out = w; b_out = b
 
 ! restore model state
 u = u_loc;v = v_loc ; w = w_loc; b = b_loc ; p = p_loc
 du = du_loc;dv = dv_loc ; dw = dw_loc; db = db_loc
 tau = tau_loc; taum1 = taum1_loc; taum2 = taum2_loc; itt = itt_loc
end subroutine model_run


subroutine model_nonlin(u_in,v_in,w_in,b_in,u_out,v_out,w_out,b_out)
 use main_module
 implicit none
 real*8,dimension(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx):: u_in,v_in,w_in,b_in
 complex*16,dimension(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3))::u_out,v_out,w_out,b_out
 integer :: i,j,k
 real*8,dimension(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx):: u_loc,v_loc
 real*8,dimension(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx):: w_loc,b_loc
 real*8,dimension(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx):: du_loc,dv_loc
 real*8,dimension(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,ks_pe-onx:ke_pe+onx):: dw_loc,db_loc
 real*8 :: fxa
 
 u_loc = u; v_loc = v ; w_loc = w; b_loc = b
 u = u_in; v = v_in; w = w_in; b = b_in
 call buoyancy_advection_2nd(db_loc)
 call momentum_advection_2nd(du_loc,dv_loc,dw_loc)
 db_loc = db_loc/Ro; du_loc = du_loc/Ro; dv_loc = dv_loc/Ro; dw_loc = dw_loc/Ro/dsqr
 u = u_loc;v = v_loc ; w = w_loc; b = b_loc 
 call forward_transform(du_loc,dv_loc,dw_loc,db_loc,u_out,v_out,w_out,b_out)
 
 if (enable_diag_balance_filter) then
  fxa = abs(kz(nz/2+1+diag_balance_filter_width)) !was 5
  do k=fstart(3),fend(3)
   do j=fstart(2),fend(2)
    do i=fstart(1),fend(1)  
     if (kz(k)**2>= fxa**2 ) then
       u_out(i,j,k) = 0.
       v_out(i,j,k) = 0.
       w_out(i,j,k) = 0.
       b_out(i,j,k) = 0.
     endif
    enddo
   enddo
  enddo
 endif 
end subroutine model_nonlin

