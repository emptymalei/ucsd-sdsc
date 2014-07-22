! -*-f90-*-
subroutine start

  use GR1D_module
  use ppm
#if HAVE_LEAK_ROS
  use leakage_rosswog,only: initialize_leakage_rosswog
#endif
  implicit none

  character(len=128) cpstring
  character(len=128) rmstring

  write(*,*) "Parsing input file..."
  !this time to get details of arrays sizes etc
  call input_parser

  if(eoskey.eq.3) then
#if HAVE_NUC_EOS
     ! Ott EOS routines: read table
     call readtable(eos_table_name)
#else
     stop "start.F90: NUC_EOS code not present but want eoskey=3"
#endif
  endif

  !total zones
  n1 = radial_zones+ghosts1*2

  if(do_rotation) then
     n_cons = 5
  else
     n_cons = 4
  endif

  !allocate & initialize variables
  call allocate_vars
  call initialize_vars
  !this time to set all variables requested values
  call input_parser
  write(*,*) "Your job name: ", trim(adjustl(jobname))
  write(*,*) "Using ", radial_zones, " radial zones and ", ghosts1, " ghosts zones"
  write(*,*) "Using ",trim(adjustl(reconstruction_method))," reconstruction"
  write(*,*) "Using ",trim(adjustl(flux_type)), " for hydro"

  ! wipe output directory
  rmstring="rm -rf "//trim(adjustl(outdir))//"/*"
  call system(rmstring)
  ! copy parameter file
  cpstring="cp parameters "//trim(adjustl(outdir))
  call system(cpstring)

  !setting up initial data
  call problem

  !Collapse specific setups
  if(initial_data.eq."Collapse") then
     write(*,*) "Finding envelope binding energy"
     call map_envelope_binding_energy(profile_name)
     write(*,*) "Finding accretion radii"
     call findaccretionradii
  endif

  !setup PPM coefficients
  if(reconstruction_method.eq.'ppm') then
     write(*,*) "Setting up PPM coefficients"
     call ppm_coefficients
  endif

  !setup dynamic output control
  if(dynamic_output_control) then
     call output_control
  endif

  !if leaking, initialize variables etc.
#if HAVE_LEAK_ROS
  call initialize_leakage_rosswog
#endif

end subroutine start
