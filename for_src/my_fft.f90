

subroutine my_ifft_2D(QS,QP)
 use main_module
 implicit none
 real(p3dfft_type)    :: QP(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3))
 complex(p3dfft_type) :: QS(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3))
 call p3dfft_btran_c2r(QS,QP,'nff')
 QP=QP /(Lx*Ly) 
end subroutine my_ifft_2D


subroutine my_fft_2D(QP,QS)
 use main_module
 implicit none
 real(p3dfft_type)    :: QP(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3))
 complex(p3dfft_type) :: QS(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3))
 call p3dfft_ftran_r2c(QP,QS,'ffn')
 QS = QS*dx*dy
end subroutine my_fft_2D





subroutine my_ifft_3D(QS,QP)
 use main_module
 implicit none
 real(p3dfft_type)    :: QP(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3))
 complex(p3dfft_type) :: QS(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3))
 call p3dfft_btran_c2r(QS,QP,'fft')
 QP=QP /(Lx*Ly*Lz) 
end subroutine my_ifft_3D




subroutine my_fft_3D(QP,QS)
 use main_module
 implicit none
 real(p3dfft_type)    :: QP(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3))
 complex(p3dfft_type) :: QS(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3))
 call p3dfft_ftran_r2c(QP,QS,'fft')
 QS = QS*dx*dy*dz 
end subroutine my_fft_3D





subroutine my_ifft_vert_mode(QS,QP,mode)
 use main_module
 implicit none
 real(p3dfft_type)    :: QP(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3))
 complex(p3dfft_type) :: QS(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3))
 real*8 :: mode(nz,nz)
 complex(p3dfft_type) :: Qloc(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) 
 integer :: k1,k2
   
 Qloc = 0
 do k1=1,nz ! mode
    do k2=1,nz ! z 
       Qloc(:,:,k2) = Qloc(:,:,k2) + QS(:,:,k1)*mode(k2,k1)
    enddo     
 enddo    
 call p3dfft_btran_c2r(Qloc,QP,'nff')
 QP=QP /(Lx*Ly) 
end subroutine my_ifft_vert_mode




subroutine my_fft_vert_mode(QP,QS,mode)
 use main_module
 implicit none
 real(p3dfft_type)    :: QP(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3))
 complex(p3dfft_type) :: QS(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3))
 real*8 :: mode(nz,nz)
 complex(p3dfft_type) :: Qloc(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) 
 integer :: k1,k2
 
 call p3dfft_ftran_r2c(QP,Qloc,'ffn')
 Qloc = Qloc*dx*dy
 QS = 0
 do k1=1,nz ! mode
    do k2=1,nz ! z 
       QS(:,:,k1) = QS(:,:,k1) + Qloc(:,:,k2)*mode(k2,k1)*dz
    enddo     
 enddo  
end subroutine my_fft_vert_mode


function hat_ksqr(k,dx)
 implicit none
 real*8, intent(in) :: k,dx
 real*8 :: hat_ksqr
 hat_ksqr = 2d0*(1d0-cos(k*dx))/dx**2
end function hat_ksqr


