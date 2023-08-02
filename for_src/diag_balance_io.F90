#ifdef with_netcdf_parallel 

subroutine init_diag_balance
 use main_module
 use module_diag_balance
 implicit none
 include "netcdf.inc"
 include "mpif.h"
 integer :: ncid,iret,xdim,ydim,zdim,tdim,id,i,j,k,nc_mode
 real*8 :: x(nx),y(ny),z(nz)
 
 call allocate_module_diag_balance
 
 if (enable_diag_balance_chunks) then
  if (my_pe==0) print*,' output in chunks not available for parallel netcdf'
  call halt_stop(' in diag_balance_io.F90')
 endif
  
 if (my_pe==0)  print*,'preparing file diag_balance.cdf'

 nc_mode = IOR(nf_clobber,nf_mpiio)
 nc_mode = IOR(nc_mode,nf_netcdf4 )
 nc_mode = IOR(nc_mode,nf_classic_model )
 
 iret = NF_CREATE_par('diag_balance.cdf',nc_mode,MPI_COMM_WORLD,MPI_INFO_NULL, ncid) 
 if (my_pe==0.and.iret/=0) print*,'at opening:',nf_strerror(iret)
 iret=nf_set_fill(ncid, NF_NOFILL, iret)

 iret = nf_def_dim(ncid,'x', nx,xdim)
 iret = nf_def_dim(ncid,'y', ny,ydim)
 iret = nf_def_dim(ncid,'z', nz,zdim)  
 iret = nf_def_dim(ncid,'Time', nf_unlimited,tdim)
 
 iret = nf_def_var(ncid, 'x',NF_DOUBLE,1,xdim,id) 
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('x'),'x')
 iret = nf_put_att_text(ncid,id,'units',len_trim('m'),'m')
 
 iret = nf_def_var(ncid, 'y',NF_DOUBLE,1,ydim,id) 
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('y'),'y')
 iret = nf_put_att_text(ncid,id,'units',len_trim('m'),'m') 
 
 iret = nf_def_var(ncid, 'z',NF_DOUBLE,1,zdim,id) 
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('z'),'z')
 iret = nf_put_att_text(ncid,id,'units',len_trim('m'),'m') 

 iret = nf_def_var(ncid, 'Time',NF_DOUBLE,1,tdim,id)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('time'),'time')
 iret = nf_put_att_text(ncid,id,'unit',len_trim('seconds'),'seconds')

 id  = ncvdef (ncid,'u0',NF_DOUBLE,4,(/xdim, ydim,zdim,tdim/),iret)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
 iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
 id  = ncvdef (ncid,'u1',NF_DOUBLE,4,(/xdim, ydim,zdim,tdim/),iret)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
 iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
 id  = ncvdef (ncid,'u2',NF_DOUBLE,4,(/xdim, ydim,zdim,tdim/),iret)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
 iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
 id  = ncvdef (ncid,'u3',NF_DOUBLE,4,(/xdim, ydim,zdim,tdim/),iret)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
 iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')

 id  = ncvdef (ncid,'v0',NF_DOUBLE,4,(/xdim, ydim,zdim,tdim/),iret)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
 iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
 id  = ncvdef (ncid,'v1',NF_DOUBLE,4,(/xdim, ydim,zdim,tdim/),iret)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
 iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
 id  = ncvdef (ncid,'v2',NF_DOUBLE,4,(/xdim, ydim,zdim,tdim/),iret)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
 iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
 id  = ncvdef (ncid,'v3',NF_DOUBLE,4,(/xdim, ydim,zdim,tdim/),iret)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
 iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
  
 id  = ncvdef (ncid,'w0',NF_DOUBLE,4,(/xdim, ydim,zdim,tdim/),iret)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
 iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
 id  = ncvdef (ncid,'w1',NF_DOUBLE,4,(/xdim, ydim,zdim,tdim/),iret)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
 iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
 id  = ncvdef (ncid,'w2',NF_DOUBLE,4,(/xdim, ydim,zdim,tdim/),iret)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
 iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
 id  = ncvdef (ncid,'w3',NF_DOUBLE,4,(/xdim, ydim,zdim,tdim/),iret)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
 iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
  
 id  = ncvdef (ncid,'b0',NF_DOUBLE,4,(/xdim, ydim,zdim,tdim/),iret)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('buoyancy'),'buoyancy')
 iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s^2'),'m/s^2')
 id  = ncvdef (ncid,'b1',NF_DOUBLE,4,(/xdim, ydim,zdim,tdim/),iret)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('buoyancy'),'buoyancy')
 iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s^2'),'m/s^2')
 id  = ncvdef (ncid,'b2',NF_DOUBLE,4,(/xdim, ydim,zdim,tdim/),iret)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('buoyancy'),'buoyancy')
 iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s^2'),'m/s^2')
 id  = ncvdef (ncid,'b3',NF_DOUBLE,4,(/xdim, ydim,zdim,tdim/),iret)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('buoyancy'),'buoyancy')
 iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s^2'),'m/s^2')
   
 iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'dx',NF_DOUBLE,1,dx)
 iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'dy',NF_DOUBLE,1,dy)
 iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'dz',NF_DOUBLE,1,dz)
 iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'N0',NF_DOUBLE,1,N0)
 iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'f0',NF_DOUBLE,1,f0)
 iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'Ro',NF_DOUBLE,1,Ro)
 iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'dsqr',NF_DOUBLE,1,dsqr)

 iret= nf_enddef(ncid)
 
 do i=1,nx
   x(i)=(i-1)*dx
 enddo     
 iret=nf_inq_varid(ncid,'x',id)
 iret= nf_put_vara_double(ncid,id,is_pe,ie_pe-is_pe+1,x(is_pe:ie_pe)) 

 do j=1,ny
   y(j)=(j-1)*dy
 enddo 
 iret=nf_inq_varid(ncid,'y',id)
 iret= nf_put_vara_double(ncid,id,js_pe,je_pe-js_pe+1,y(js_pe:je_pe)) 

 do k=1,nz
   z(k)=(k-1)*dz
 enddo   
 iret=nf_inq_varid(ncid,'z',id)
 iret= nf_put_vara_double(ncid,id,ks_pe,ke_pe-ks_pe+1,z(ks_pe:ke_pe)) 

 iret= nf_close(ncid) 
 if (my_pe==0) print*,' done'   
   
end subroutine init_diag_balance


subroutine write_diag_balance 
 use main_module
 use module_diag_balance
 implicit none
 include "netcdf.inc"
 include "mpif.h"
 integer :: ncid,iret,id,nc_mode,start(4),count(4),ilen,k
 real*8 :: aloc(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)
 
 if (my_pe==0)  print*,'writing to file diag_balance.cdf'

 nc_mode = IOR(nf_write,nf_mpiio)
 nc_mode = IOR(nc_mode,nf_netcdf4 )
 nc_mode = IOR(nc_mode,nf_classic_model )
 
 iret = NF_open_par('diag_balance.cdf',nc_mode,MPI_COMM_WORLD,MPI_INFO_NULL, ncid) 
 if (my_pe==0.and.iret/=0) print*,'at opening:',nf_strerror(iret)
 iret=nf_set_fill(ncid, NF_NOFILL, iret)
 
 iret=nf_inq_dimid(ncid,'Time',id)
 iret=nf_inq_dimlen(ncid,id,ilen)
 ilen=ilen+1
 
 iret=nf_inq_varid(ncid,'Time',id) 
 iret = nf_var_par_access(ncid, id, nf_collective)
 iret= nf_put_vara_double(ncid,id,ilen,1,time) 
 
 start = (/is_pe,        js_pe,        ks_pe        ,ilen/)
 count = (/ie_pe-is_pe+1,je_pe-js_pe+1,ke_pe-ks_pe+1,1/)
 
 iret=nf_inq_varid(ncid,'u0',id)
 iret = nf_var_par_access(ncid, id, nf_collective)
 aloc = u0(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)
 iret= nf_put_vara_double(ncid,id,start,count,aloc)
 iret=nf_inq_varid(ncid,'u1',id)
 iret = nf_var_par_access(ncid, id, nf_collective)
 aloc = u1(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)
 iret= nf_put_vara_double(ncid,id,start,count,aloc)
 iret=nf_inq_varid(ncid,'u2',id)
 iret = nf_var_par_access(ncid, id, nf_collective)
 aloc = u2(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)
 iret= nf_put_vara_double(ncid,id,start,count,aloc)
 iret=nf_inq_varid(ncid,'u3',id)
 iret = nf_var_par_access(ncid, id, nf_collective)
 aloc = u3(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)
 iret= nf_put_vara_double(ncid,id,start,count,aloc)

 iret=nf_inq_varid(ncid,'v0',id)
 iret = nf_var_par_access(ncid, id, nf_collective)
 aloc = v0(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)
 iret= nf_put_vara_double(ncid,id,start,count,aloc)
 iret=nf_inq_varid(ncid,'v1',id)
 iret = nf_var_par_access(ncid, id, nf_collective)
 aloc = v1(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)
 iret= nf_put_vara_double(ncid,id,start,count,aloc)
 iret=nf_inq_varid(ncid,'v2',id)
 iret = nf_var_par_access(ncid, id, nf_collective)
 aloc = v2(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)
 iret= nf_put_vara_double(ncid,id,start,count,aloc)
 iret=nf_inq_varid(ncid,'v3',id)
 iret = nf_var_par_access(ncid, id, nf_collective)
 aloc = v3(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)
 iret= nf_put_vara_double(ncid,id,start,count,aloc)

 iret=nf_inq_varid(ncid,'w0',id)
 iret = nf_var_par_access(ncid, id, nf_collective)
 aloc = w0(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)
 iret= nf_put_vara_double(ncid,id,start,count,aloc)
 iret=nf_inq_varid(ncid,'w1',id)
 iret = nf_var_par_access(ncid, id, nf_collective)
 aloc = w1(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)
 iret= nf_put_vara_double(ncid,id,start,count,aloc)
 iret=nf_inq_varid(ncid,'w2',id)
 iret = nf_var_par_access(ncid, id, nf_collective)
 aloc = w2(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)
 iret= nf_put_vara_double(ncid,id,start,count,aloc)
 iret=nf_inq_varid(ncid,'w3',id)
 iret = nf_var_par_access(ncid, id, nf_collective)
 aloc = w3(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)
 iret= nf_put_vara_double(ncid,id,start,count,aloc)

 iret=nf_inq_varid(ncid,'b0',id)
 iret = nf_var_par_access(ncid, id, nf_collective)
 aloc = b0(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)
 iret= nf_put_vara_double(ncid,id,start,count,aloc)
 iret=nf_inq_varid(ncid,'b1',id)
 iret = nf_var_par_access(ncid, id, nf_collective)
 aloc = b1(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)
 iret= nf_put_vara_double(ncid,id,start,count,aloc)
 iret=nf_inq_varid(ncid,'b2',id)
 iret = nf_var_par_access(ncid, id, nf_collective)
 aloc = b2(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)
 iret= nf_put_vara_double(ncid,id,start,count,aloc)
 iret=nf_inq_varid(ncid,'b3',id)
 iret = nf_var_par_access(ncid, id, nf_collective)
 aloc = b3(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)
 iret= nf_put_vara_double(ncid,id,start,count,aloc)
     
 iret= nf_close(ncid) 
  
end subroutine write_diag_balance 


subroutine write_diag_balance_chunks
 use main_module
 implicit none
 if (my_pe==0) print*,' output in chunks not available for parallel netcdf'
 call halt_stop(' in diag_balance_io.F90')
end subroutine write_diag_balance_chunks

#else

subroutine init_diag_balance
 use main_module
 use module_diag_balance
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret,i,xdim,ydim,zdim,timeid,timedim,id
 real*8 :: x(nx),y(ny),z(nz)
 
 call allocate_module_diag_balance
 
 if (enable_diag_balance_chunks) then

  write(balance_name,'(a,i5,a,i5,a,i5,a,i5,a,i5,a,i5,a)')  &
  'diag_balance_i=',is_pe,':',ie_pe,'_j=',js_pe,':',je_pe,'_k=',ks_pe,':',ke_pe,'.cdf'
  call replace_space_zero(balance_name)
 
  if (my_pe==0) print*,'creating file ',balance_name(1:len_trim(balance_name))
  
  iret = NF_CREATE (balance_name,or(NF_CLASSIC_MODEL,nf_64bit_offset ),ncid) 
  if (iret.ne.0) print*,nf_strerror(iret)
  
  iret=nf_set_fill(ncid, NF_NOFILL, iret)
  xdim  = ncddef(ncid, 'x', ie_pe-is_pe+1, iret)
  ydim  = ncddef(ncid, 'y', je_pe-js_pe+1, iret)
  zdim  = ncddef(ncid, 'z', ke_pe-ks_pe+1, iret)
 
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
  
  id  = ncvdef (ncid,'u0',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
  
  id  = ncvdef (ncid,'v0',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
  
  id  = ncvdef (ncid,'w0',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
  
  id  = ncvdef (ncid,'b0',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('buoyancy'),'buoyancy')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s^2'),'m/s^2')
  
  id  = ncvdef (ncid,'u1',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
  
  id  = ncvdef (ncid,'v1',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
  
  id  = ncvdef (ncid,'w1',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
  
  id  = ncvdef (ncid,'b1',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('buoyancy'),'buoyancy')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s^2'),'m/s^2')

  id  = ncvdef (ncid,'u2',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
  
  id  = ncvdef (ncid,'v2',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
  
  id  = ncvdef (ncid,'w2',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
  
  id  = ncvdef (ncid,'b2',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('buoyancy'),'buoyancy')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s^2'),'m/s^2')

  id  = ncvdef (ncid,'u3',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
  
  id  = ncvdef (ncid,'v3',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
  
  id  = ncvdef (ncid,'w3',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
  
  id  = ncvdef (ncid,'b3',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('buoyancy'),'buoyancy')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s^2'),'m/s^2')
 
  iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'dx',NF_DOUBLE,1,dx)
  iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'dy',NF_DOUBLE,1,dy)
  iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'dz',NF_DOUBLE,1,dz)
  iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'N0',NF_DOUBLE,1,N0)
  iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'f0',NF_DOUBLE,1,f0)
  iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'Ro',NF_DOUBLE,1,Ro)
  iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'dsqr',NF_DOUBLE,1,dsqr)

  iret= nf_enddef(ncid) 
  
  do i=is_pe,ie_pe
   x(i)=(i-1)*dx
  enddo 
  iret=nf_inq_varid(ncid,'x',id)
  iret= nf_put_vara_double(ncid,id,(/1/),(/ie_pe-is_pe+1/),x(is_pe:ie_pe))  
  
  do i=js_pe,je_pe
   y(i)=(i-1)*dy
  enddo 
  iret=nf_inq_varid(ncid,'y',id)
  iret= nf_put_vara_double(ncid,id,(/1/),(/je_pe-js_pe+1/),y(js_pe:je_pe))  
  
  do i=ks_pe,ke_pe
   z(i)=(i-1)*dz
  enddo 
  iret=nf_inq_varid(ncid,'z',id)
  iret= nf_put_vara_double(ncid,id,(/1/),(/ke_pe-ks_pe+1/),z(ks_pe:ke_pe))   
  iret= nf_close(ncid)   

 else
 
 if (my_pe==0) then
 
  if (my_pe==0) print*,'creating file diag_balance.cdf'

  iret = NF_CREATE ('diag_balance.cdf',or(NF_CLASSIC_MODEL,nf_64bit_offset ),ncid) 
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
  
  id  = ncvdef (ncid,'u0',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
  
  id  = ncvdef (ncid,'v0',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
  
  id  = ncvdef (ncid,'w0',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
  
  id  = ncvdef (ncid,'b0',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('buoyancy'),'buoyancy')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s^2'),'m/s^2')

  id  = ncvdef (ncid,'u1',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
  
  id  = ncvdef (ncid,'v1',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
  
  id  = ncvdef (ncid,'w1',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
  
  id  = ncvdef (ncid,'b1',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('buoyancy'),'buoyancy')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s^2'),'m/s^2')

  id  = ncvdef (ncid,'u2',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
  
  id  = ncvdef (ncid,'v2',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
  
  id  = ncvdef (ncid,'w2',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
  
  id  = ncvdef (ncid,'b2',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('buoyancy'),'buoyancy')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s^2'),'m/s^2')

  id  = ncvdef (ncid,'u3',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
  
  id  = ncvdef (ncid,'v3',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
  
  id  = ncvdef (ncid,'w3',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
  
  id  = ncvdef (ncid,'b3',NF_DOUBLE,4,(/xdim, ydim,zdim,Timedim/),iret)
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
 
 endif
end subroutine init_diag_balance



subroutine write_diag_balance 
 use main_module
 use module_diag_balance
 implicit none
 include "netcdf.inc"
 include "mpif.h"
 integer :: ncid,iret,n
 integer :: tdimid,ilen,timeid,id
 integer :: tag=1,ist(3),isz(3),ien(3)
 integer, dimension(MPI_STATUS_SIZE) :: Status
 real*8, allocatable :: a(:,:,:)

 if (my_pe==0) then
   print*,'writing to file diag_balance.cdf at t=',time
   iret=nf_open('diag_balance.cdf',NF_WRITE,ncid)
   iret=nf_set_fill(ncid, NF_NOFILL, iret)
   iret=nf_inq_dimid(ncid,'Time',tdimid)
   iret=nf_inq_dimlen(ncid, tdimid,ilen)
   iret=nf_inq_varid(ncid,'Time',timeid)
   ilen=ilen+1
   iret= nf_put_vara_double(ncid,timeid,(/ilen/),(/1/),(/time/))
   
   iret=nf_inq_varid(ncid,'u0',id)
   iret= nf_put_vara_double(ncid,id,(/istart(1),istart(2),istart(3),ilen/), &
                            (/isize(1),isize(2),isize(3),1/),&
                            u0(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))
   iret=nf_inq_varid(ncid,'v0',id)
   iret= nf_put_vara_double(ncid,id,(/istart(1),istart(2),istart(3),ilen/), &
                            (/isize(1),isize(2),isize(3),1/),&
                            v0(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))
   iret=nf_inq_varid(ncid,'w0',id)
   iret= nf_put_vara_double(ncid,id,(/istart(1),istart(2),istart(3),ilen/), &
                            (/isize(1),isize(2),isize(3),1/),&
                            w0(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))
   iret=nf_inq_varid(ncid,'b0',id)
   iret= nf_put_vara_double(ncid,id,(/istart(1),istart(2),istart(3),ilen/), &
                            (/isize(1),isize(2),isize(3),1/),&
                            b0(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))   

   iret=nf_inq_varid(ncid,'u1',id)
   iret= nf_put_vara_double(ncid,id,(/istart(1),istart(2),istart(3),ilen/), &
                            (/isize(1),isize(2),isize(3),1/),&
                            u1(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))
   iret=nf_inq_varid(ncid,'v1',id)
   iret= nf_put_vara_double(ncid,id,(/istart(1),istart(2),istart(3),ilen/), &
                            (/isize(1),isize(2),isize(3),1/),&
                            v1(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))
   iret=nf_inq_varid(ncid,'w1',id)
   iret= nf_put_vara_double(ncid,id,(/istart(1),istart(2),istart(3),ilen/), &
                            (/isize(1),isize(2),isize(3),1/),&
                            w1(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))
   iret=nf_inq_varid(ncid,'b1',id)
   iret= nf_put_vara_double(ncid,id,(/istart(1),istart(2),istart(3),ilen/), &
                            (/isize(1),isize(2),isize(3),1/),&
                            b1(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))   

   iret=nf_inq_varid(ncid,'u2',id)
   iret= nf_put_vara_double(ncid,id,(/istart(1),istart(2),istart(3),ilen/), &
                            (/isize(1),isize(2),isize(3),1/),&
                            u2(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))
   iret=nf_inq_varid(ncid,'v2',id)
   iret= nf_put_vara_double(ncid,id,(/istart(1),istart(2),istart(3),ilen/), &
                            (/isize(1),isize(2),isize(3),1/),&
                            v2(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))
   iret=nf_inq_varid(ncid,'w2',id)
   iret= nf_put_vara_double(ncid,id,(/istart(1),istart(2),istart(3),ilen/), &
                            (/isize(1),isize(2),isize(3),1/),&
                            w2(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))
   iret=nf_inq_varid(ncid,'b2',id)
   iret= nf_put_vara_double(ncid,id,(/istart(1),istart(2),istart(3),ilen/), &
                            (/isize(1),isize(2),isize(3),1/),&
                            b2(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))   

   iret=nf_inq_varid(ncid,'u3',id)
   iret= nf_put_vara_double(ncid,id,(/istart(1),istart(2),istart(3),ilen/), &
                            (/isize(1),isize(2),isize(3),1/),&
                            u3(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))
   iret=nf_inq_varid(ncid,'v3',id)
   iret= nf_put_vara_double(ncid,id,(/istart(1),istart(2),istart(3),ilen/), &
                            (/isize(1),isize(2),isize(3),1/),&
                            v3(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))
   iret=nf_inq_varid(ncid,'w3',id)
   iret= nf_put_vara_double(ncid,id,(/istart(1),istart(2),istart(3),ilen/), &
                            (/isize(1),isize(2),isize(3),1/),&
                            w3(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))
   iret=nf_inq_varid(ncid,'b3',id)
   iret= nf_put_vara_double(ncid,id,(/istart(1),istart(2),istart(3),ilen/), &
                            (/isize(1),isize(2),isize(3),1/),&
                            b3(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))   
 

 endif

 do n=1,n_pes-1
  call mpi_barrier(MPI_COMM_WORLD, iret)
  if (my_pe==n) then
        call mpi_send(istart,3,mpi_integer,0,tag,mpi_comm_world,iret)
        call mpi_send(iend  ,3,mpi_integer,0,tag,mpi_comm_world,iret)
        call mpi_send(isize ,3,mpi_integer,0,tag,mpi_comm_world,iret)
        
        call mpi_send(u0(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)&
                     ,isize(1)*isize(2)*isize(3),mpi_real8,0,tag,mpi_comm_world,iret)
        call mpi_send(v0(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) &
                      ,isize(1)*isize(2)*isize(3),mpi_real8,0,tag,mpi_comm_world,iret)
        call mpi_send(w0(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) &
                     ,isize(1)*isize(2)*isize(3),mpi_real8,0,tag,mpi_comm_world,iret)
        call mpi_send(b0(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) &
                     ,isize(1)*isize(2)*isize(3),mpi_real8,0,tag,mpi_comm_world,iret)

        call mpi_send(u1(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)&
                     ,isize(1)*isize(2)*isize(3),mpi_real8,0,tag,mpi_comm_world,iret)
        call mpi_send(v1(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) &
                      ,isize(1)*isize(2)*isize(3),mpi_real8,0,tag,mpi_comm_world,iret)
        call mpi_send(w1(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) &
                     ,isize(1)*isize(2)*isize(3),mpi_real8,0,tag,mpi_comm_world,iret)
        call mpi_send(b1(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) &
                     ,isize(1)*isize(2)*isize(3),mpi_real8,0,tag,mpi_comm_world,iret)
       
        call mpi_send(u2(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)&
                     ,isize(1)*isize(2)*isize(3),mpi_real8,0,tag,mpi_comm_world,iret)
        call mpi_send(v2(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) &
                      ,isize(1)*isize(2)*isize(3),mpi_real8,0,tag,mpi_comm_world,iret)
        call mpi_send(w2(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) &
                     ,isize(1)*isize(2)*isize(3),mpi_real8,0,tag,mpi_comm_world,iret)
        call mpi_send(b2(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) &
                     ,isize(1)*isize(2)*isize(3),mpi_real8,0,tag,mpi_comm_world,iret)

        call mpi_send(u3(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe)&
                     ,isize(1)*isize(2)*isize(3),mpi_real8,0,tag,mpi_comm_world,iret)
        call mpi_send(v3(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) &
                      ,isize(1)*isize(2)*isize(3),mpi_real8,0,tag,mpi_comm_world,iret)
        call mpi_send(w3(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) &
                     ,isize(1)*isize(2)*isize(3),mpi_real8,0,tag,mpi_comm_world,iret)
        call mpi_send(b3(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe) &
                     ,isize(1)*isize(2)*isize(3),mpi_real8,0,tag,mpi_comm_world,iret)
  


  else if (my_pe==0) then
        call mpi_recv(ist,3,mpi_integer,n,tag,mpi_comm_world,status,iret)
        call mpi_recv(ien,3,mpi_integer,n,tag,mpi_comm_world,status,iret)
        call mpi_recv(isz,3,mpi_integer,n,tag,mpi_comm_world,status,iret)
        allocate(a(ist(1):ien(1),ist(2):ien(2),ist(3):ien(3)) )
        
        call mpi_recv(a,isz(1)*isz(2)*isz(3),mpi_real8,n,tag,mpi_comm_world,status,iret)
        iret=nf_inq_varid(ncid,'u0',id)
        iret= nf_put_vara_double(ncid,id,(/ist(1),ist(2),ist(3),ilen/),(/isz(1),isz(2),isz(3),1/),a)
        
        call mpi_recv(a,isz(1)*isz(2)*isz(3),mpi_real8,n,tag,mpi_comm_world,status,iret)
        iret=nf_inq_varid(ncid,'v0',id)
        iret= nf_put_vara_double(ncid,id,(/ist(1),ist(2),ist(3),ilen/),(/isz(1),isz(2),isz(3),1/),a)
        
        call mpi_recv(a,isz(1)*isz(2)*isz(3),mpi_real8,n,tag,mpi_comm_world,status,iret)
        iret=nf_inq_varid(ncid,'w0',id)
        iret= nf_put_vara_double(ncid,id,(/ist(1),ist(2),ist(3),ilen/),(/isz(1),isz(2),isz(3),1/),a)
        
        call mpi_recv(a,isz(1)*isz(2)*isz(3),mpi_real8,n,tag,mpi_comm_world,status,iret)
        iret=nf_inq_varid(ncid,'b0',id)
        iret= nf_put_vara_double(ncid,id,(/ist(1),ist(2),ist(3),ilen/),(/isz(1),isz(2),isz(3),1/),a)
 
        call mpi_recv(a,isz(1)*isz(2)*isz(3),mpi_real8,n,tag,mpi_comm_world,status,iret)
        iret=nf_inq_varid(ncid,'u1',id)
        iret= nf_put_vara_double(ncid,id,(/ist(1),ist(2),ist(3),ilen/),(/isz(1),isz(2),isz(3),1/),a)
        
        call mpi_recv(a,isz(1)*isz(2)*isz(3),mpi_real8,n,tag,mpi_comm_world,status,iret)
        iret=nf_inq_varid(ncid,'v1',id)
        iret= nf_put_vara_double(ncid,id,(/ist(1),ist(2),ist(3),ilen/),(/isz(1),isz(2),isz(3),1/),a)
        
        call mpi_recv(a,isz(1)*isz(2)*isz(3),mpi_real8,n,tag,mpi_comm_world,status,iret)
        iret=nf_inq_varid(ncid,'w1',id)
        iret= nf_put_vara_double(ncid,id,(/ist(1),ist(2),ist(3),ilen/),(/isz(1),isz(2),isz(3),1/),a)
        
        call mpi_recv(a,isz(1)*isz(2)*isz(3),mpi_real8,n,tag,mpi_comm_world,status,iret)
        iret=nf_inq_varid(ncid,'b1',id)
        iret= nf_put_vara_double(ncid,id,(/ist(1),ist(2),ist(3),ilen/),(/isz(1),isz(2),isz(3),1/),a)
  
  
        call mpi_recv(a,isz(1)*isz(2)*isz(3),mpi_real8,n,tag,mpi_comm_world,status,iret)
        iret=nf_inq_varid(ncid,'u2',id)
        iret= nf_put_vara_double(ncid,id,(/ist(1),ist(2),ist(3),ilen/),(/isz(1),isz(2),isz(3),1/),a)
        
        call mpi_recv(a,isz(1)*isz(2)*isz(3),mpi_real8,n,tag,mpi_comm_world,status,iret)
        iret=nf_inq_varid(ncid,'v2',id)
        iret= nf_put_vara_double(ncid,id,(/ist(1),ist(2),ist(3),ilen/),(/isz(1),isz(2),isz(3),1/),a)
        
        call mpi_recv(a,isz(1)*isz(2)*isz(3),mpi_real8,n,tag,mpi_comm_world,status,iret)
        iret=nf_inq_varid(ncid,'w2',id)
        iret= nf_put_vara_double(ncid,id,(/ist(1),ist(2),ist(3),ilen/),(/isz(1),isz(2),isz(3),1/),a)
        
        call mpi_recv(a,isz(1)*isz(2)*isz(3),mpi_real8,n,tag,mpi_comm_world,status,iret)
        iret=nf_inq_varid(ncid,'b2',id)
        iret= nf_put_vara_double(ncid,id,(/ist(1),ist(2),ist(3),ilen/),(/isz(1),isz(2),isz(3),1/),a)
 
  
        call mpi_recv(a,isz(1)*isz(2)*isz(3),mpi_real8,n,tag,mpi_comm_world,status,iret)
        iret=nf_inq_varid(ncid,'u3',id)
        iret= nf_put_vara_double(ncid,id,(/ist(1),ist(2),ist(3),ilen/),(/isz(1),isz(2),isz(3),1/),a)
        
        call mpi_recv(a,isz(1)*isz(2)*isz(3),mpi_real8,n,tag,mpi_comm_world,status,iret)
        iret=nf_inq_varid(ncid,'v3',id)
        iret= nf_put_vara_double(ncid,id,(/ist(1),ist(2),ist(3),ilen/),(/isz(1),isz(2),isz(3),1/),a)
        
        call mpi_recv(a,isz(1)*isz(2)*isz(3),mpi_real8,n,tag,mpi_comm_world,status,iret)
        iret=nf_inq_varid(ncid,'w3',id)
        iret= nf_put_vara_double(ncid,id,(/ist(1),ist(2),ist(3),ilen/),(/isz(1),isz(2),isz(3),1/),a)
        
        call mpi_recv(a,isz(1)*isz(2)*isz(3),mpi_real8,n,tag,mpi_comm_world,status,iret)
        iret=nf_inq_varid(ncid,'b3',id)
        iret= nf_put_vara_double(ncid,id,(/ist(1),ist(2),ist(3),ilen/),(/isz(1),isz(2),isz(3),1/),a)
  
  
        deallocate(a)
  endif
 enddo  


 if (my_pe == 0) iret= nf_close(ncid) 
end subroutine write_diag_balance 





subroutine write_diag_balance_chunks
 use main_module
 use module_diag_balance
 implicit none
 include "netcdf.inc" 
 integer :: ncid,iret,tdimid,ilen,timeid,id,start(4),count(4)

 if (my_pe==0) print*,'writing to file ',balance_name(1:len_trim(balance_name)),' at t=',time
 iret=nf_open(balance_name,NF_WRITE,ncid)
 iret=nf_set_fill(ncid, NF_NOFILL, iret)
 iret=nf_inq_dimid(ncid,'Time',tdimid)
 iret=nf_inq_dimlen(ncid, tdimid,ilen)
 iret=nf_inq_varid(ncid,'Time',timeid)
 ilen=ilen+1
 iret= nf_put_vara_double(ncid,timeid,(/ilen/),(/1/),(/time/))
   
 start = (/1,1,1,ilen/)
 count = (/ie_pe-is_pe+1,je_pe-js_pe+1,ke_pe-ks_pe+1,1/)
 
 iret=nf_inq_varid(ncid,'u0',id)
 iret= nf_put_vara_double(ncid,id,start,count,u0(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))
 iret=nf_inq_varid(ncid,'v0',id)
 iret= nf_put_vara_double(ncid,id,start,count,v0(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))
 iret=nf_inq_varid(ncid,'w0',id)
 iret= nf_put_vara_double(ncid,id,start,count,w0(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))
 iret=nf_inq_varid(ncid,'b0',id)
 iret= nf_put_vara_double(ncid,id,start,count,b0(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))   
 
 iret=nf_inq_varid(ncid,'u1',id)
 iret= nf_put_vara_double(ncid,id,start,count,u1(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))
 iret=nf_inq_varid(ncid,'v1',id)
 iret= nf_put_vara_double(ncid,id,start,count,v1(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))
 iret=nf_inq_varid(ncid,'w1',id)
 iret= nf_put_vara_double(ncid,id,start,count,w1(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))
 iret=nf_inq_varid(ncid,'b1',id)
 iret= nf_put_vara_double(ncid,id,start,count,b1(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))   

 iret=nf_inq_varid(ncid,'u2',id)
 iret= nf_put_vara_double(ncid,id,start,count,u2(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))
 iret=nf_inq_varid(ncid,'v2',id)
 iret= nf_put_vara_double(ncid,id,start,count,v2(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))
 iret=nf_inq_varid(ncid,'w2',id)
 iret= nf_put_vara_double(ncid,id,start,count,w2(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))
 iret=nf_inq_varid(ncid,'b2',id)
 iret= nf_put_vara_double(ncid,id,start,count,b2(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))   

 iret=nf_inq_varid(ncid,'u3',id)
 iret= nf_put_vara_double(ncid,id,start,count,u3(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))
 iret=nf_inq_varid(ncid,'v3',id)
 iret= nf_put_vara_double(ncid,id,start,count,v3(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))
 iret=nf_inq_varid(ncid,'w3',id)
 iret= nf_put_vara_double(ncid,id,start,count,w3(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))
 iret=nf_inq_varid(ncid,'b3',id)
 iret= nf_put_vara_double(ncid,id,start,count,b3(is_pe:ie_pe,js_pe:je_pe,ks_pe:ke_pe))   

 iret= nf_close(ncid) 
end subroutine write_diag_balance_chunks

#endif
