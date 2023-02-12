
!---------------------------------------------------------------------------------
! eigenvectors and values for vertical modes
!---------------------------------------------------------------------------------

subroutine omega2_vbc(kx,ky,kz,dx,dy,dz,f0,N0,dsqr,om2)
 implicit none
 real*8, intent(in) :: kx,ky,kz,dx,dy,dz,f0,N0,dsqr
 real*8, intent(out) :: om2
 real*8 :: hat_ksqr,hat_epm,ksqr,k2,m2,N2_hat,f2_hat,arg

 !arg =  mode_nr*pi/(2.*N)  ,  kz =  mode_nr*pi/(Nz*dz) , arg = kz/2*dz
 arg = kz/2*dz
 k2 = hat_ksqr(kx,dx) + hat_ksqr(ky,dy)
 m2 =  4*sin(arg)**2 / dz**2 
 N2_hat =  N0**2 * cos(arg)**2
 f2_hat =  f0**2 * hat_epm(kx,dx)*hat_epm(ky,dy)
 Ksqr = dsqr*k2 + m2
 if (Ksqr>0d0) then
   om2 = (N2_hat*k2+f2_hat*m2)/Ksqr
 else
  om2=0
 endif
end subroutine omega2_vbc


subroutine coeff_p_vbc(kx,ky,kz,dx,dy,dz,s,f0,N0,dsqr,p)
 implicit none
 real*8, intent(in) :: kx,ky,kz,dx,dy,dz,s,f0,N0,dsqr
 real*8, intent(out) :: p
 real*8 :: om2,hat_ksqr,hat_epm,d,k2,f2_hat,N2_hat,m2,arg
 
 arg = kz/2*dz
 call omega2_vbc(kx,ky,kz,dx,dy,dz,f0,N0,dsqr,om2)
 k2 = hat_ksqr(kx,dx) + hat_ksqr(ky,dy)
 m2 =  4*sin(arg)**2 / dz**2 
 f2_hat = f0**2*hat_epm(kx,dx)*hat_epm(ky,dy)
 N2_hat =  N0**2 * cos(arg)**2
 d = k2*N2_hat + f2_hat*m2
 if (s==0d0 .and. d>0d0) then
    p = f2_hat*N2_hat/d
 else if (d>0d0) then
    P =  (N2_hat-dsqr*om2)*(om2-f2_hat)/(2*d) 
 else
    p = 0d0
 endif  
end subroutine coeff_p_vbc


subroutine vec_q_vbc(kx,ky,kz,dx,dy,dz,s,f0,N0,dsqr,qx,qy,qz,qb)
 implicit none
 real*8, intent(in) :: kx,ky,kz,dx,dy,dz,s,f0,N0,dsqr
 complex*16, intent(out) :: qx,qy,qz,qb
 complex*16,parameter :: i = cmplx(0.,1.,kind=8)
 real*8 :: om2,om,f2_hat,hat_epm,N2_hat,fxa,lam,m2,arg
 complex*16 :: hat_ep,hat_em,hat_kp
 
 arg = kz/2*dz
 call omega2_vbc(kx,ky,kz,dx,dy,dz,f0,N0,dsqr,om2)
 om = sqrt(om2)
 f2_hat =  f0**2 * hat_epm(kx,dx)*hat_epm(ky,dy)
 N2_hat =  N0**2 * cos(arg)**2
 lam = 2*sin(arg) / dz
 m2 = lam**2
 
 if ( kx**2+ky**2>0d0) then
  fxa = (f2_hat-s**2*om2)
  qx = (-i*f0*hat_ep(kx,dx)*hat_em(ky,dy)*hat_kp(ky,dy) + s*om*hat_kp(kx,dx))/fxa
  qy = (+i*f0*hat_em(kx,dx)*hat_ep(ky,dy)*hat_kp(kx,dx) + s*om*hat_kp(ky,dy))/fxa
  if (kz**2>0d0) then
   fxa = N2_hat-dsqr*s**2*om2 
   qz = i*lam*s*om/fxa 
   qb = -lam*cos(arg)*N0**2/fxa 
  else
   qz=0d0;qb=0d0
  endif 
 else
  qx = -i*s
  qy = s**2
  qz =  0d0
  qb =  0d0
 endif 
 if (s==0d0 .and. (f2_hat==0d0 .or. N2_hat==0d0) ) then
  qx=0d0;qy=0d0;qz=0d0;qb=0d0
 endif
 if (s/=0d0 .and. om2 == 0d0) then
  qx=0d0;qy=0d0;qz=0d0;qb=0d0
 endif 
end subroutine vec_q_vbc


subroutine vec_p_vbc(kx,ky,kz,dx,dy,dz,s,f0,N0,dsqr,px,py,pz,pb)
 implicit none
 real*8, intent(in) :: kx,ky,kz,dx,dy,dz,s,f0,N0,dsqr
 complex*16, intent(out) :: px,py,pz,pb
 real*8 :: p

 call coeff_p_vbc(kx,ky,kz,dx,dy,dz,s,f0,N0,dsqr,p)
 call vec_q_vbc  (kx,ky,kz,dx,dy,dz,s,f0,N0,dsqr,px,py,pz,pb)
 px = p*conjg(px)
 py = p*conjg(py)
 pz = p*conjg(pz)*dsqr
 pb = p*conjg(pb)/N0**2
end subroutine vec_p_vbc

!---------------------------------------------------------------------------------
! eigenvectors and values for freely propagating waves
!---------------------------------------------------------------------------------

subroutine omega2_discrete(kx,ky,kz,dx,dy,dz,f0,N0,dsqr,om2)
 implicit none
 real*8, intent(in) :: kx,ky,kz,dx,dy,dz,f0,N0,dsqr
 real*8, intent(out) :: om2
 real*8 :: hat_ksqr,hat_epm,ksqr,k2
 
 k2 = hat_ksqr(kx,dx) + hat_ksqr(ky,dy)
 Ksqr = dsqr*k2 + hat_ksqr(kz,dz)
 if (Ksqr>0d0) then
   om2 = (hat_epm(kz,dz)*N0**2*k2+f0**2*hat_ksqr(kz,dz)*hat_epm(kx,dx)*hat_epm(ky,dy))/ksqr
 else
  om2=0
 endif
end subroutine omega2_discrete

subroutine coeff_p_discrete(kx,ky,kz,dx,dy,dz,s,f0,N0,dsqr,p)
 implicit none
 real*8, intent(in) :: kx,ky,kz,dx,dy,dz,s,f0,N0,dsqr
 real*8, intent(out) :: p
 real*8 :: om2,hat_ksqr,hat_epm,d,k2,f2_hat,N2_hat
 
 call omega2_discrete(kx,ky,kz,dx,dy,dz,f0,N0,dsqr,om2)
 k2 = hat_ksqr(kx,dx) + hat_ksqr(ky,dy)
 f2_hat = f0**2*hat_epm(kx,dx)*hat_epm(ky,dy)
 N2_hat =  N0**2 * hat_epm(kz,dz) 
 d = k2*N2_hat + f2_hat*hat_ksqr(kz,dz)
 if (s==0d0 .and. d>0d0) then
    p = f2_hat*N2_hat/d
 else if (d>0d0) then
    P =  (N2_hat-dsqr*om2)*(om2-f2_hat)/(2*d) 
 else
    p = 0d0
 endif 
end subroutine coeff_p_discrete


subroutine vec_q_discrete(kx,ky,kz,dx,dy,dz,s,f0,N0,dsqr,qx,qy,qz,qb)
 implicit none
 real*8, intent(in) :: kx,ky,kz,dx,dy,dz,s,f0,N0,dsqr
 complex*16, intent(out) :: qx,qy,qz,qb
 complex*16,parameter :: i = cmplx(0.,1.,kind=8)
 real*8 :: om2,om,f2_hat,hat_epm,N2_hat,fxa
 complex*16 :: hat_ep,hat_em,hat_kp

 call omega2_discrete(kx,ky,kz,dx,dy,dz,f0,N0,dsqr,om2)
 om = sqrt(om2)
 f2_hat =  f0**2 * hat_epm(kx,dx)*hat_epm(ky,dy)
 N2_hat =  N0**2 * hat_epm(kz,dz)

 if ( kx**2+ky**2>0d0) then
  fxa = (f2_hat-s**2*om2)
  qx = (-i*f0*hat_ep(kx,dx)*hat_em(ky,dy)*hat_kp(ky,dy) + s*om*hat_kp(kx,dx))/fxa
  qy = (+i*f0*hat_em(kx,dx)*hat_ep(ky,dy)*hat_kp(kx,dx) + s*om*hat_kp(ky,dy))/fxa
  if (kz**2>0d0) then
   fxa = (N2_hat-dsqr*s**2*om2 ) 
   qz =  hat_kp(kz,dz)*s*om/fxa 
   qb =  i*hat_em(kz,dz)*hat_kp(kz,dz)*N0**2/fxa
  else
   qz=0d0;qb=0d0
  endif 
 else
  qx = -i*s
  qy = s**2
  qz =  0d0
  qb =  0d0
 endif
 
 if (s==0d0 .and. (f2_hat==0d0 .or. N2_hat==0d0) ) then
  qx=0d0;qy=0d0;qz=0d0;qb=0d0
 endif
 if (s/=0d0 .and. om2 == 0d0) then
  qx=0d0;qy=0d0;qz=0d0;qb=0d0
 endif 
 
 ! Silvano's hack
 if ( (kx**2+ky**2>0d0) .and. (kz**2==0d0) .and. (s/=0d0) ) then
   qx = 0d0
   qy = 0d0
   qz = i*s
   qb = s**2*N0*sqrt(dsqr)
 end if 
end subroutine vec_q_discrete


subroutine vec_p_discrete(kx,ky,kz,dx,dy,dz,s,f0,N0,dsqr,px,py,pz,pb)
 implicit none
 real*8, intent(in) :: kx,ky,kz,dx,dy,dz,s,f0,N0,dsqr
 complex*16, intent(out) :: px,py,pz,pb
 real*8 :: p
 complex*16,parameter :: i = cmplx(0.,1.,kind=8)
 
 call coeff_p_discrete(kx,ky,kz,dx,dy,dz,s,f0,N0,dsqr,p)
 call vec_q_discrete  (kx,ky,kz,dx,dy,dz,s,f0,N0,dsqr,px,py,pz,pb)
 px = p*conjg(px)
 py = p*conjg(py)
 pz = p*conjg(pz)*dsqr
 pb = p*conjg(pb)/N0**2
 
 ! Silvano's hack 
 if ( (kx**2 + ky**2 == 0d0) .and. (kz**2 > 0d0) ) then
   px=i*s/2d0; py=s**2/2d0
 end if
 if ( (kx**2+ky**2>0d0) .and. (kz**2==0d0) .and. (s/=0d0) ) then
   px = 0d0
   py = 0d0
   pz = -i*s/2d0
   pb = s**2/(2d0*N0*sqrt(dsqr))
 end if 
end subroutine vec_p_discrete

!---------------------------------------------------------------------------------
! helper functions
!---------------------------------------------------------------------------------

function hat_ksqr(k,dx)
 implicit none
 real*8, intent(in) :: k,dx
 real*8 :: hat_ksqr
 hat_ksqr = 2d0*(1d0-cos(k*dx))/dx**2
end function hat_ksqr

function hat_epm(k,dx)
 implicit none
 real*8, intent(in) :: k,dx
 real*8 :: hat_epm
 hat_epm =  0.5d0*(1+cos(k*dx))
end function hat_epm

function hat_kp(k,dx)
 implicit none
 real*8, intent(in) :: k,dx
 complex*16,parameter :: i = cmplx(0.,1.,kind=8)
 complex*16 :: hat_kp
 hat_kp =  -i*(exp(i*k*dx) -1.)/dx
end function hat_kp

function hat_km(k,dx)
 implicit none
 real*8, intent(in) :: k,dx
 complex*16,parameter :: i = cmplx(0.,1.,kind=8)
 complex*16 :: hat_km
 hat_km =  -i*(1-exp(-i*k*dx) )/dx
end function hat_km

function hat_ep(k,dx)
 implicit none
 real*8, intent(in) :: k,dx
 complex*16,parameter :: i = cmplx(0.,1.,kind=8)
 complex*16 :: hat_ep
 hat_ep =  (exp(i*k*dx) +1.)*0.5d0
end function hat_ep

function hat_em(k,dx)
 implicit none
 real*8, intent(in) :: k,dx
 complex*16,parameter :: i = cmplx(0.,1.,kind=8)
 complex*16 :: hat_em
 hat_em =  (exp(-i*k*dx) +1.)*0.5d0
end function hat_em



