

subroutine init_diag_opt_balance
  use main_module
  use module_diag_opt_balance
  implicit none
  include "netcdf.inc"
  integer :: ncid,iret,i,xdim,ydim,zdim,timeid,timedim,id
  real*8 :: x(nx),y(ny),z(nz)

  call allocate_module_diag_opt_balance

  if (my_pe==0) then

    print*,'creating file diag_opt_balance.cdf'
    iret = NF_CREATE ('diag_opt_balance.cdf',or(NF_CLASSIC_MODEL,nf_64bit_offset ),ncid)
    if (iret.ne.0) print*,nf_strerror(iret)

    iret=nf_set_fill(ncid, NF_NOFILL, iret)
    xdim  = ncddef(ncid, 'x', nx, iret)
    ydim  = ncddef(ncid, 'y', ny, iret)
    zdim  = ncddef(ncid, 'z', nz, iret)

    Timedim = ncddef(ncid,'Time', nf_unlimited, iret)
    timeid  = ncvdef (ncid,'Time', NF_DOUBLE,1,(/timedim/),iret)
    iret = nf_put_att_text(ncid,timeid,'long_name',len_trim('time'),'time')
    iret = nf_put_att_text(ncid,timeid,'unit',len_trim('seconds'),'seconds')

    id  = ncvdef (ncid,'x',NF_DOUBLE,1,(/xdim/),iret)
    iret = nf_put_att_text(ncid,id,'long_name',len_trim('x'),'x')
    iret = nf_put_att_text(ncid,id,'unit',len_trim('m'),'m')

    id  = ncvdef (ncid,'y',NF_DOUBLE,1,(/ydim/),iret)
    iret = nf_put_att_text(ncid,id,'long_name',len_trim('y'),'y')
    iret = nf_put_att_text(ncid,id,'unit',len_trim('m'),'m')

    id  = ncvdef (ncid,'z',NF_DOUBLE,1,(/zdim/),iret)
    iret = nf_put_att_text(ncid,id,'long_name',len_trim('z'),'z')
    iret = nf_put_att_text(ncid,id,'unit',len_trim('m'),'m')

    id  = ncvdef (ncid,'u_bal',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
    iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
    iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')

    id  = ncvdef (ncid,'v_bal',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
    iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
    iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')

    id  = ncvdef (ncid,'w_bal',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
    iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
    iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')

    id  = ncvdef (ncid,'b_bal',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
    iret = nf_put_att_text(ncid,id,'long_name',len_trim('buoyancy'),'buoyancy')
    iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s^2'),'m/s^2')

    iret= nf_enddef(ncid)

    do i=1,nx
     x(i)=(i-1)*dx
    enddo
    iret=nf_inq_varid(ncid,'x',id)
    iret= nf_put_vara_double(ncid,id,(/1/),(/nx/),x)

    do i=1,ny
     y(i)=(i-1)*dy
    enddo
    iret=nf_inq_varid(ncid,'y',id)
    iret= nf_put_vara_double(ncid,id,(/1/),(/ny/),y)

    do i=1,nz
     z(i)=(i-1)*dz
    enddo
    iret=nf_inq_varid(ncid,'z',id)
    iret= nf_put_vara_double(ncid,id,(/1/),(/nz/),z)

    iret= nf_close(ncid)

   endif

end subroutine init_diag_opt_balance



subroutine write_diag_opt_balance
  use main_module
  use module_diag_opt_balance
  implicit none
  include "netcdf.inc"
  include "mpif.h"
  integer :: ncid,iret,n
  integer :: tdimid,ilen,timeid,id
  integer :: tag=1,ist(3),isz(3),ien(3)
  integer, dimension(MPI_STATUS_SIZE) :: Status
  real*8, allocatable :: a(:,:,:)

  if (my_pe==0) then
    print*,'writing to file diag_opt_balance.cdf at t=',time
    iret=nf_open('diag_opt_balance.cdf',NF_WRITE,ncid)
    iret=nf_set_fill(ncid, NF_NOFILL, iret)
    iret=nf_inq_dimid(ncid,'Time',tdimid)
    iret=nf_inq_dimlen(ncid, tdimid,ilen)
    iret=nf_inq_varid(ncid,'Time',timeid)
    ilen=ilen+1
    iret= nf_put_vara_double(ncid,timeid,(/ilen/),(/1/),(/time/))

    iret=nf_inq_varid(ncid,'u_bal',id)
    iret= nf_put_vara_double(ncid,id,(/istart(1),istart(2),istart(3),ilen/), &
                             (/isize(1),isize(2),isize(3),1/),&
                             u_bal(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))
    iret=nf_inq_varid(ncid,'v_bal',id)
    iret= nf_put_vara_double(ncid,id,(/istart(1),istart(2),istart(3),ilen/), &
                             (/isize(1),isize(2),isize(3),1/),&
                             v_bal(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))
    iret=nf_inq_varid(ncid,'w_bal',id)
    iret= nf_put_vara_double(ncid,id,(/istart(1),istart(2),istart(3),ilen/), &
                             (/isize(1),isize(2),isize(3),1/),&
                             w_bal(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))
    iret=nf_inq_varid(ncid,'b_bal',id)
    iret= nf_put_vara_double(ncid,id,(/istart(1),istart(2),istart(3),ilen/), &
                             (/isize(1),isize(2),isize(3),1/),&
                             b_bal(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))

  end if

  do n=1,n_pes-1
   call mpi_barrier(MPI_COMM_WORLD, iret)
   if (my_pe==n) then
         call mpi_send(istart,3,mpi_integer,0,tag,mpi_comm_world,iret)
         call mpi_send(iend  ,3,mpi_integer,0,tag,mpi_comm_world,iret)
         call mpi_send(isize ,3,mpi_integer,0,tag,mpi_comm_world,iret)

         call mpi_send(u_bal(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)&
                      ,isize(1)*isize(2)*isize(3),mpi_real8,0,tag,mpi_comm_world,iret)
         call mpi_send(v_bal(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) &
                       ,isize(1)*isize(2)*isize(3),mpi_real8,0,tag,mpi_comm_world,iret)
         call mpi_send(w_bal(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) &
                      ,isize(1)*isize(2)*isize(3),mpi_real8,0,tag,mpi_comm_world,iret)
         call mpi_send(b_bal(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) &
                      ,isize(1)*isize(2)*isize(3),mpi_real8,0,tag,mpi_comm_world,iret)

    else if (my_pe==0) then
          call mpi_recv(ist,3,mpi_integer,n,tag,mpi_comm_world,status,iret)
          call mpi_recv(ien,3,mpi_integer,n,tag,mpi_comm_world,status,iret)
          call mpi_recv(isz,3,mpi_integer,n,tag,mpi_comm_world,status,iret)
          allocate(a(ist(1):ien(1),ist(2):ien(2),ist(3):ien(3)) )

          call mpi_recv(a,isz(1)*isz(2)*isz(3),mpi_real8,n,tag,mpi_comm_world,status,iret)
          iret=nf_inq_varid(ncid,'u_bal',id)
          iret= nf_put_vara_double(ncid,id,(/ist(1),ist(2),ist(3),ilen/),(/isz(1),isz(2),isz(3),1/),a)

          call mpi_recv(a,isz(1)*isz(2)*isz(3),mpi_real8,n,tag,mpi_comm_world,status,iret)
          iret=nf_inq_varid(ncid,'v_bal',id)
          iret= nf_put_vara_double(ncid,id,(/ist(1),ist(2),ist(3),ilen/),(/isz(1),isz(2),isz(3),1/),a)

          call mpi_recv(a,isz(1)*isz(2)*isz(3),mpi_real8,n,tag,mpi_comm_world,status,iret)
          iret=nf_inq_varid(ncid,'w_bal',id)
          iret= nf_put_vara_double(ncid,id,(/ist(1),ist(2),ist(3),ilen/),(/isz(1),isz(2),isz(3),1/),a)

          call mpi_recv(a,isz(1)*isz(2)*isz(3),mpi_real8,n,tag,mpi_comm_world,status,iret)
          iret=nf_inq_varid(ncid,'b_bal',id)
          iret= nf_put_vara_double(ncid,id,(/ist(1),ist(2),ist(3),ilen/),(/isz(1),isz(2),isz(3),1/),a)

          deallocate(a)
    end if
  end do

  if (my_pe == 0) iret= nf_close(ncid)
end subroutine write_diag_opt_balance
