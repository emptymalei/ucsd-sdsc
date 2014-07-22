!-*-f90-*-
subroutine restart_output_h5

#if 0
  use GR1D_module
  use hdf5
#if HAVE_LEAK_ROS
  use leakage_rosswog, only : leak_tau, have_old_tau
#endif
  implicit none

  integer i,j
  double precision x1_16,D_16,S_16,T_16,press_16,X_16,Y_16,temp_16
  integer sw1,sw2,sw3,sw4,sw7,sw8,sw9,sw10,b1
  character(len=1024) filename
  character(len=1024) filenamelist
  character(len=256) basename
  character(len=1024) formatname
  real*8 test_ye

  integer error,rank,cerror
  integer(HID_T) file_id,dset_id,dspace_id,aspace_id,attr_id
  integer(HID_T) atype_id
  integer(SIZE_T) attrlen
  
  integer(HSIZE_T) dims1(1), dims2(2), dims3(3)


  !This routine outputs all the necessary variables needed to restart from any given point in time in HDF5 format
  sw1 = 0
  sw2 = 0
  sw3 = 0
  b1 = 0
  sw4 = 0
  sw7 = 0
  sw8 = 0
  sw9 = 0
  sw10 = 0
  cerror = 0

  if (switch1) sw1 = 1
  if (switch2) sw2 = 1
  if (switch3) sw3 = 1
  if (bounce) b1 = 1
  if (do_nupress) sw4 = 1
  if (GR) sw7 = 1
  if (fake_neutrinos) sw8 = 1
#if HAVE_LEAK_ROS
  if (have_old_tau) sw9 = 1
#endif
  if (do_rotation) sw10 = 1

  basename = "restart"
  call generate_filename(basename,outdir,&
       time,nt,"h5",filename)

  !append filename to list for post-processing
  filenamelist = trim(adjustl(outdir))//"/"//trim(adjustl("list_of_restartfiles.txt"))
  open(unit=473,file=trim(adjustl(filenamelist)),status="unknown", &
       form='formatted',position="append")
  write(473,*) trim(adjustl(filename))
  close(473)

  !open HDF5 file
  call h5open_f(error)
  cerror = cerror + error

  call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,error)
  cerror = cerror + error

  !write scalars
  rank = 1
  dims1(1) = 1
  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "time", H5T_NATIVE_DOUBLE, dspace_id,&
       & dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, time, dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)  
  cerror = cerror + error

  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "time_c", H5T_NATIVE_DOUBLE, dspace_id,&
       & dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, time_c, dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)  
  cerror = cerror + error

  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "dt", H5T_NATIVE_DOUBLE, dspace_id,&
       & dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dt, dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)
  cerror = cerror + error  

  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "nt", H5T_NATIVE_INTEGER, dspace_id,&
       & dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, nt, dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)  
  cerror = cerror + error

  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "tdump", H5T_NATIVE_DOUBLE, dspace_id,&
       & dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tdump, dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)  
  cerror = cerror + error

  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "tdump_restart", H5T_NATIVE_DOUBLE,&
       & dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tdump_restart, dims1,&
       & error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)  
  cerror = cerror + error

  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "tdump_scalar", H5T_NATIVE_DOUBLE,&
       & dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tdump_scalar, dims1,&
       & error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)  
  cerror = cerror + error

  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "t_bounce", H5T_NATIVE_DOUBLE, dspace_id,&
       & dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, t_bounce, dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)  
  cerror = cerror + error

  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "shock_radius", H5T_NATIVE_DOUBLE, dspace_id,&
       & dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, shock_radius, dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)  
  cerror = cerror + error

  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "ishock", H5T_NATIVE_INTEGER, dspace_id,&
       & dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, ishock(1), dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)  
  cerror = cerror + error

  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "sw1", H5T_NATIVE_INTEGER, dspace_id,&
       & dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, sw1, dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)  
  cerror = cerror + error

  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "sw2", H5T_NATIVE_INTEGER, dspace_id,&
       & dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, sw2, dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)  
  cerror = cerror + error

  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "sw3", H5T_NATIVE_INTEGER, dspace_id,&
       & dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, sw3, dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)  
  cerror = cerror + error

  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "sw4", H5T_NATIVE_INTEGER, dspace_id,&
       & dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, sw4, dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)  
  cerror = cerror + error

  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "sw7", H5T_NATIVE_INTEGER, dspace_id,&
       & dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, sw7, dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)  
  cerror = cerror + error

  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "sw8", H5T_NATIVE_INTEGER, dspace_id,&
       & dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, sw8, dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)  
  cerror = cerror + error

#if HAVE_LEAK_ROS
  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "sw9", H5T_NATIVE_INTEGER, dspace_id,&
       & dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, sw9, dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)  
  cerror = cerror + error
#endif

  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "sw10", H5T_NATIVE_INTEGER, dspace_id,&
       & dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, sw10, dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)  
  cerror = cerror + error

  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "eoskey", H5T_NATIVE_INTEGER, dspace_id,&
       & dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, eoskey, dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)  
  cerror = cerror + error

  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "n1", H5T_NATIVE_INTEGER, dspace_id,&
       & dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, n1, dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)  
  cerror = cerror + error

  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "b1", H5T_NATIVE_INTEGER, dspace_id,&
       & dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, b1, dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)  
  cerror = cerror + error
  
!1D arrays
  rank = 1
  dims1(1) = n1

  if (GR) then
     ! only GR
     call h5screate_simple_f(rank, dims1, dspace_id, error)
     call h5dcreate_f(file_id, "X", H5T_NATIVE_DOUBLE, dspace_id,&
          & dset_id, error)
     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, X, dims1, error)
     call h5dclose_f(dset_id, error)
     call h5sclose_f(dspace_id, error)
     cerror = cerror + error
     
     call h5screate_simple_f(rank, dims1, dspace_id, error)
     call h5dcreate_f(file_id, "Xp", H5T_NATIVE_DOUBLE, dspace_id,&
          & dset_id, error)
     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, Xp, dims1, error)
     call h5dclose_f(dset_id, error)
     call h5sclose_f(dspace_id, error)
     cerror = cerror + error
     
     call h5screate_simple_f(rank, dims1, dspace_id, error)
     call h5dcreate_f(file_id, "Xm", H5T_NATIVE_DOUBLE, dspace_id,&
          & dset_id, error)
     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, Xm, dims1, error)
     call h5dclose_f(dset_id, error)
     call h5sclose_f(dspace_id, error)
     cerror = cerror + error
     
     call h5screate_simple_f(rank, dims1, dspace_id, error)
     call h5dcreate_f(file_id, "alp", H5T_NATIVE_DOUBLE, dspace_id,&
          & dset_id, error)
     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, alp, dims1, error)
     call h5dclose_f(dset_id, error)
     call h5sclose_f(dspace_id, error)
     cerror = cerror + error

     call h5screate_simple_f(rank, dims1, dspace_id, error)
     call h5dcreate_f(file_id, "alpp", H5T_NATIVE_DOUBLE, dspace_id,&
          & dset_id, error)
     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, alpp, dims1, error)
     call h5dclose_f(dset_id, error)
     call h5sclose_f(dspace_id, error)
     cerror = cerror + error

     call h5screate_simple_f(rank, dims1, dspace_id, error)
     call h5dcreate_f(file_id, "alpm", H5T_NATIVE_DOUBLE, dspace_id,&
          & dset_id, error)
     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, alpm, dims1, error)
     call h5dclose_f(dset_id, error)
     call h5sclose_f(dspace_id, error)
     cerror = cerror + error

     call h5screate_simple_f(rank, dims1, dspace_id, error)
     call h5dcreate_f(file_id, "mgrav", H5T_NATIVE_DOUBLE, dspace_id,&
          & dset_id, error)
     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mgrav, dims1, error)
     call h5dclose_f(dset_id, error)
     call h5sclose_f(dspace_id, error)
     cerror = cerror + error

     call h5screate_simple_f(rank, dims1, dspace_id, error)
     call h5dcreate_f(file_id, "mgravi", H5T_NATIVE_DOUBLE, dspace_id&
          &, dset_id, error)
     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mgravi, dims1, error)
     call h5dclose_f(dset_id, error)
     call h5sclose_f(dspace_id, error)
     cerror = cerror + error

     call h5screate_simple_f(rank, dims1, dspace_id, error)
     call h5dcreate_f(file_id, "W", H5T_NATIVE_DOUBLE, dspace_id,&
          & dset_id, error)
     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, W, dims1, error)
     call h5dclose_f(dset_id, error)
     call h5sclose_f(dspace_id, error)
     cerror = cerror + error

     call h5screate_simple_f(rank, dims1, dspace_id, error)
     call h5dcreate_f(file_id, "v", H5T_NATIVE_DOUBLE, dspace_id,&
          & dset_id, error)
     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, v, dims1, error)
     call h5dclose_f(dset_id, error)
     call h5sclose_f(dspace_id, error)
     cerror = cerror + error

     if (do_rotation) then
        call h5screate_simple_f(rank, dims1, dspace_id, error)
        call h5dcreate_f(file_id, "vphi", H5T_NATIVE_DOUBLE, dspace_id,&
             & dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vphi, dims1, error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)
        cerror = cerror + error
     endif

  else
     if (do_rotation) then
        call h5screate_simple_f(rank, dims1, dspace_id, error)
        call h5dcreate_f(file_id, "vphi1", H5T_NATIVE_DOUBLE, dspace_id,&
             & dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vphi1, dims1, error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)
        cerror = cerror + error
     endif

  endif

  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "v1", H5T_NATIVE_DOUBLE, dspace_id,&
       & dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, v1, dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)
  cerror = cerror + error
  
  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "mass", H5T_NATIVE_DOUBLE, dspace_id,&
       & dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mass, dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)
  cerror = cerror + error

  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "temperature", H5T_NATIVE_DOUBLE,&
       & dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, temp, dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)
  cerror = cerror + error

  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "rho", H5T_NATIVE_DOUBLE, dspace_id,&
       & dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, rho, dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)
  cerror = cerror + error

  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "ye", H5T_NATIVE_DOUBLE, dspace_id,&
       & dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, ye, dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)
  cerror = cerror + error

  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "press", H5T_NATIVE_DOUBLE, dspace_id,&
       & dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, press, dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)
  cerror = cerror + error

  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "eps", H5T_NATIVE_DOUBLE, dspace_id,&
       & dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, eps, dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)
  cerror = cerror + error

#if HAVE_LEAK_ROS
  call h5screate_simple_f(rank, dims1, dspace_id, error)
  call h5dcreate_f(file_id, "leak_tau", H5T_NATIVE_DOUBLE, dspace_id,&
       & dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, leak_tau, dims1, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)
  cerror = cerror + error
#endif
     
  if (do_nupress) then
     call h5screate_simple_f(rank, dims1, dspace_id, error)
     call h5dcreate_f(file_id, "press_nu", H5T_NATIVE_DOUBLE,&
          & dspace_id, dset_id, error)
     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, press_nu, dims1,&
          & error)
     call h5dclose_f(dset_id, error)
     call h5sclose_f(dspace_id, error)
     cerror = cerror + error
     
     call h5screate_simple_f(rank, dims1, dspace_id, error)
     call h5dcreate_f(file_id, "dnupdr", H5T_NATIVE_DOUBLE, dspace_id&
          &, dset_id, error)
     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dnupdr, dims1, error)
     call h5dclose_f(dset_id, error)
     call h5sclose_f(dspace_id, error)
     cerror = cerror + error
  endif

  if (cerror.ne.0) then
     write(*,*) "We have errors on writing HDF5 restart file", cerror
     stop
  endif

  call h5fclose_f(file_id,error)
  call h5close_f(error)


contains
     subroutine generate_filename(varname,outdir,time,nt,suffix,fname)

      implicit none


      real*8 time
      integer nt
      character(*) varname
      character(len=128) outdir
      character*(*) suffix
      character*(*) fname
      character*(400) aa
      character(len=100) outtime
      character(len=20) outnt
      integer i,ii

      aa=" "
      fname=" "
      write(outnt,"(i10.10)") nt
      fname = trim(adjustl(outdir))//"/"//trim(adjustl(varname))&
           &//"_nt_"//outnt
      write(outtime,"(f9.7)") time
      fname = trim(adjustl(fname))//"_time_"//trim(adjustl(outtime))&
           &//"."//trim(adjustl(suffix))

    end subroutine generate_filename

#endif
end subroutine restart_output_h5


subroutine restart_init_h5
#if 0
  use GR1D_module
  use hdf5
#if HAVE_LEAK_ROS
  use leakage_rosswog, only : leak_tau,have_old_tau
#endif
  implicit none

  real*8 temp_time
  integer sw1,sw2,sw3,sw4,sw7,sw8,sw9,sw10,b1,ibuffer,lbuffer
  real*8 x1_16,D_16,S_16,T_16,press_16,X_16,Y_16,temp_16

  integer :: tint, dtint, tdumpint, tmsint
  real*8 :: ttemp
  
  integer error,rank,cerror
  integer(HID_T) file_id,dset_id,dspace_id,aspace_id,attr_id
  integer(HID_T) atype_id,i,j
  integer(SIZE_T) attrlen
  
  integer(HSIZE_T) dims1(1), dims2(2), dims3(3)

  cerror = 0

  sw1 = 0
  sw2 = 0
  sw3 = 0
  sw4 = 0
  sw7 = 0
  sw8 = 0
  sw9 = 0
  sw10 = 0

  !This routine reloads all the variables
  call h5open_f(error)
  if (error.ne.0) then
     stop "Error reading in restart file"
  endif

  call h5fopen_f(trim(adjustl(restart_file_name)),H5F_ACC_RDONLY_F&
       &,file_id,error)
  if (error.ne.0) then
     stop "Error reading in restart file"
  endif

  !read scalars
  rank = 1
  dims1(1) = 1
  call h5dopen_f(file_id,"time",dset_id,error)
  call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,time,dims1,error)
  call h5dclose_f(dset_id,error)
  cerror = cerror + error

  call h5dopen_f(file_id,"time_c",dset_id,error)
  call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,time_c,dims1,error)
  call h5dclose_f(dset_id,error)
  cerror = cerror + error

  call h5dopen_f(file_id,"dt",dset_id,error)
  call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,dt,dims1,error)
  call h5dclose_f(dset_id,error)
  cerror = cerror + error

  call h5dopen_f(file_id,"nt",dset_id,error)
  call h5dread_f(dset_id,H5T_NATIVE_INTEGER,nt,dims1,error)
  call h5dclose_f(dset_id,error)
  cerror = cerror + error

  call h5dopen_f(file_id,"sw1",dset_id,error)
  call h5dread_f(dset_id,H5T_NATIVE_INTEGER,sw1,dims1,error)
  call h5dclose_f(dset_id,error)
  cerror = cerror + error

  call h5dopen_f(file_id,"sw2",dset_id,error)
  call h5dread_f(dset_id,H5T_NATIVE_INTEGER,sw2,dims1,error)
  call h5dclose_f(dset_id,error)
  cerror = cerror + error

  call h5dopen_f(file_id,"sw3",dset_id,error)
  call h5dread_f(dset_id,H5T_NATIVE_INTEGER,sw3,dims1,error)
  call h5dclose_f(dset_id,error)
  cerror = cerror + error
  
  call h5dopen_f(file_id,"b1",dset_id,error)
  call h5dread_f(dset_id,H5T_NATIVE_INTEGER,b1,dims1,error)
  call h5dclose_f(dset_id,error)
  cerror = cerror + error

  switch1 = .false. 
  switch2 = .false.
  switch3 = .false.
  bounce = .false.

  if (sw1.eq.1) switch1 = .true.
  if (sw2.eq.1) switch2 = .true.
  if (sw3.eq.1) switch3 = .true.
  if (b1.eq.1) bounce = .true.  
  
  call h5dopen_f(file_id,"sw4",dset_id,error)
  call h5dread_f(dset_id,H5T_NATIVE_INTEGER,sw4,dims1,error)
  call h5dclose_f(dset_id,error)
  cerror = cerror + error
  
  call h5dopen_f(file_id,"sw7",dset_id,error)
  call h5dread_f(dset_id,H5T_NATIVE_INTEGER,sw7,dims1,error)
  call h5dclose_f(dset_id,error)
  cerror = cerror + error
  
  call h5dopen_f(file_id,"sw8",dset_id,error)
  call h5dread_f(dset_id,H5T_NATIVE_INTEGER,sw8,dims1,error)
  call h5dclose_f(dset_id,error)
  cerror = cerror + error

#if HAVE_LEAK_ROS
  call h5dopen_f(file_id,"sw9",dset_id,error)
  call h5dread_f(dset_id,H5T_NATIVE_INTEGER,sw9,dims1,error)
  call h5dclose_f(dset_id,error)
  cerror = cerror + error

  have_old_tau = .false.
  if (sw9.eq.1) have_old_tau = .true.
#endif

  call h5dopen_f(file_id,"sw10",dset_id,error)
  call h5dread_f(dset_id,H5T_NATIVE_INTEGER,sw10,dims1,error)
  call h5dclose_f(dset_id,error)
  cerror = cerror + error
    
  lbuffer = 0
  if (do_nupress) lbuffer = 1
  if (sw4.ne.lbuffer) then
     write(*,*) "In Restart: old do_nupress =",sw4," current do_nupress =",lbuffer
     stop
  endif
  
  lbuffer = 0
  if (GR) lbuffer = 1
  if (sw7.ne.lbuffer) then
     write(*,*) "In Restart: old GR =",sw7," current GR =",lbuffer
     stop
  endif

  lbuffer = 0
  if (fake_neutrinos) lbuffer = 1
  if (sw8.ne.lbuffer) then
     write(*,*) "In Restart: old fake_neutrinos =",sw8," current fake_neutrinos =",lbuffer
     stop
  endif

  lbuffer = 0
  if (do_rotation) lbuffer = 1
  if (sw10.ne.lbuffer) then
     write(*,*) "In Restart: old do_rotation =",sw10," current do_rotation =",lbuffer
     stop
  endif

  call h5dopen_f(file_id,"eoskey",dset_id,error)
  call h5dread_f(dset_id,H5T_NATIVE_INTEGER,ibuffer,dims1,error)
  call h5dclose_f(dset_id,error)
  cerror = cerror + error

  if (ibuffer.ne.eoskey) then
     write(*,*) "In Restart: old eoskey =",ibuffer," current eoskey =",eoskey
     stop
  endif

  call h5dopen_f(file_id,"n1",dset_id,error)
  call h5dread_f(dset_id,H5T_NATIVE_INTEGER,ibuffer,dims1,error)
  call h5dclose_f(dset_id,error)
  cerror = cerror + error

  if (ibuffer.ne.n1) then
     write(*,*) "In Restart: old n1 =",ibuffer," current n1 =",n1
     stop
  endif

  call h5dopen_f(file_id,"tdump",dset_id,error)
  call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,tdump,dims1,error)
  call h5dclose_f(dset_id,error)
  cerror = cerror + error

  call h5dopen_f(file_id,"tdump_restart",dset_id,error)
  call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,tdump_restart,dims1,error)
  call h5dclose_f(dset_id,error)
  cerror = cerror + error

  call h5dopen_f(file_id,"tdump_scalar",dset_id,error)
  call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,tdump_scalar,dims1,error)
  call h5dclose_f(dset_id,error)
  cerror = cerror + error

  call h5dopen_f(file_id,"t_bounce",dset_id,error)
  call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,t_bounce,dims1,error)
  call h5dclose_f(dset_id,error)
  cerror = cerror + error

  call h5dopen_f(file_id,"shock_radius",dset_id,error)
  call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,shock_radius,dims1,error)
  call h5dclose_f(dset_id,error)
  cerror = cerror + error

  call h5dopen_f(file_id,"ishock",dset_id,error)
  call h5dread_f(dset_id,H5T_NATIVE_INTEGER,ishock(1),dims1,error)
  call h5dclose_f(dset_id,error)
  cerror = cerror + error

  !1D arrays
  rank=1
  dims1(1) = n1

  if (GR) then
     !GR
     call h5dopen_f(file_id, "X", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, X, dims1, error)
     call h5dclose_f(dset_id,error)  
     cerror = cerror + error

     call h5dopen_f(file_id, "Xp", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, Xp, dims1, error)
     call h5dclose_f(dset_id,error)  
     cerror = cerror + error

     call h5dopen_f(file_id, "Xm", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, Xm, dims1, error)
     call h5dclose_f(dset_id,error)  
     cerror = cerror + error

     call h5dopen_f(file_id, "alp", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alp, dims1, error)
     call h5dclose_f(dset_id,error)  
     cerror = cerror + error

     call h5dopen_f(file_id, "alpp", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alpp, dims1, error)
     call h5dclose_f(dset_id,error)  
     cerror = cerror + error

     call h5dopen_f(file_id, "alpm", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alpm, dims1, error)
     call h5dclose_f(dset_id,error)  
     cerror = cerror + error

     call h5dopen_f(file_id, "mgrav", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, mgrav, dims1, error)
     call h5dclose_f(dset_id,error) 
     cerror = cerror + error

     call h5dopen_f(file_id, "mgravi", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, mgravi, dims1, error)
     call h5dclose_f(dset_id,error) 
     cerror = cerror + error

     call h5dopen_f(file_id, "W", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, W, dims1, error)
     call h5dclose_f(dset_id,error) 
     cerror = cerror + error

     call h5dopen_f(file_id, "v", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, v, dims1, error)
     call h5dclose_f(dset_id,error) 
     cerror = cerror + error

     if (do_rotation) then
        call h5dopen_f(file_id, "vphi", dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, vphi, dims1, error)
        call h5dclose_f(dset_id,error) 
        cerror = cerror + error
     endif

  else
     if (do_rotation) then
        call h5dopen_f(file_id, "vphi1", dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, vphi1, dims1, error)
        call h5dclose_f(dset_id,error) 
        cerror = cerror + error
     endif

  endif

  call h5dopen_f(file_id, "v1", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, v1, dims1, error)
  call h5dclose_f(dset_id,error) 
  cerror = cerror + error
  
  call h5dopen_f(file_id, "mass", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, mass, dims1, error)
  call h5dclose_f(dset_id,error) 
  cerror = cerror + error
  
  call h5dopen_f(file_id, "temperature", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, temp, dims1, error)
  call h5dclose_f(dset_id,error) 
  cerror = cerror + error

  call h5dopen_f(file_id, "rho", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, rho, dims1, error)
  call h5dclose_f(dset_id,error) 
  cerror = cerror + error

  call h5dopen_f(file_id, "ye", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, ye, dims1, error)
  call h5dclose_f(dset_id,error) 
  cerror = cerror + error

  call h5dopen_f(file_id, "press", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, press, dims1, error)
  call h5dclose_f(dset_id,error) 
  cerror = cerror + error

  call h5dopen_f(file_id, "eps", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, eps, dims1, error)
  call h5dclose_f(dset_id,error) 
  cerror = cerror + error

#if HAVE_LEAK_ROS
  call h5dopen_f(file_id, "leak_tau", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, leak_tau, dims1, error)
  call h5dclose_f(dset_id,error) 
  cerror = cerror + error
#endif

  if (do_nupress) then
     call h5dopen_f(file_id, "press_nu", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, press_nu, dims1, error)
     call h5dclose_f(dset_id,error) 
     cerror = cerror + error
     
     call h5dopen_f(file_id, "dnupdr", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dnupdr, dims1, error)
     call h5dclose_f(dset_id,error) 
     cerror = cerror + error

  endif

  if (cerror.ne.0) then
     write(*,*) "Error(s) reading in HDF5 restart file", cerror
     stop
  endif

  call h5fclose_f(file_id,error)
  call h5close_f(error)
#endif
end subroutine restart_init_h5
