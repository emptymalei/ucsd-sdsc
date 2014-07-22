! -*-f90-*-
subroutine driver

  use GR1D_module
  use timers
  use atmos
  implicit none

  integer i, gi

  real*8 discrim

  logical :: OutputFlag = .false.
  logical :: CheckFlag = .false.
  logical :: OutputFlagScalar  = .false.
  logical :: OutputFlagRestart = .false.

  if(dynamic_output_control) then
     call output_control
  endif

  call start_timers

  tdump_scalar = 0.0d0 !output scalars right away
  tdump = 0.0d0 !output singles right away
  tdump_restart = 0.0d0 !make a restart file at the beginning
  time = 0.0d0
  time_c = 0.0d0
  nt = 0

  OutputFlag = .false.
  OutputFlagScalar = .false.

  !if doing restart here is where we call to read in file.
#if 0
  if (do_restart) then
     write(*,*) "Doing restart: ", trim(adjustl(restart_file_name))
     call restart_init_h5
  endif
#endif

  write(*,*) "Begin time integration loop:"
  IntegrationLoop: do 

     call SetTimeStep

     if(dynamic_output_control) then
        call output_control
     endif
     if(mod(nt,ntinfo).eq. 0) then
        call outinfo
     endif 

     if (nt.ge.ntmax) then
        write(*,*) "Done! :-) ntmax reached"
        call output_all(1)
        call output_all(2)
        call output_timers
!        call restart_output_h5
        exit
     endif

     if (time.eq.tend) then
        write(*,*) "Done! :-) tend reached"
        call output_all(1)
        call output_all(2)
        call output_timers
!        call restart_output_h5
        exit
     endif
!!   Output/Checking

     if ( mod(nt,ntout) .eq. 0 .and. ntout.ne.-1) OutputFlag = .true.

     if ( mod(nt,ntout_restart) .eq. 0 .and. &
          ntout_restart.ne.-1) OutputFlagRestart = .true.

     if ( mod(nt,ntout_scalar).eq.0 .and. ntout_scalar.ne.-1) &
          OutputFlagScalar = .true.

     if ( time.ge.tdump) then
        tdump=tdump+dtout
        OutputFlag = .true.
     endif

     if ( time.ge.tdump_restart) then
        OutputFlagRestart = .true.
     endif

     if ( time.ge.tdump_scalar) then
        tdump_scalar=tdump_scalar+dtout_scalar
        OutputFlagScalar = .true.
     endif

     if(nt.eq.0) then
        OutputFlag = .true.
        OutputFlagScalar = .true.
     endif

     if (OutputFlag) then
        call output_all(1)
        call output_timers
        OutputFlag = .false.
     endif

     if (OutputFlagRestart) then
        call restart_output_h5
        tdump_restart=tdump_restart+dtout_restart
        OutputFlagRestart = .false.
     endif

     if (OutputFlagScalar) then
        call output_all(2)
        OutputFlagScalar = .false.
     endif

     if ((time+dt/time_gf).gt.tend) dt = (tend-time)*time_gf

!!   Integrate
     call Step(dt)
     
     if (initial_data.eq."Collapse") then
        call get_shock_radius
        call mass_analysis    
        call get_binding_energy
        if (do_rotation) then
           call rotation_analysis
        endif
     endif

     nt=nt+1
     time = time+dt/time_gf
     time_c = time_c + dt/time_gf*alp(ghosts1+1)

  enddo IntegrationLoop

contains

  subroutine SetTimeStep

    use GR1D_module
    implicit none


    real*8 sound,dtnew
    integer keytemp,eosflag,keyerr
    integer i
    logical nan,inf

    if (GR) then
       call findnaninf(v,n1,nan,inf)
    else
       call findnaninf(v1,n1,nan,inf)
    endif

    if(nan) stop "NaNs found :-/"

    ! sets time step
    ! Note: so far only cfl-stability conditions is implemented
    ! Other conditions may be added

    ! get the speed of sound from the eos
    if (eoskey.eq.3) then
       do i=1,n1
          keytemp = 0 ! always should come in with eps to maintain accuracy
          eosflag = 6 ! we want cs2 to be reset
          keyerr = 0
          call eos(i,rho(i),temp(i),ye(i),eps(i),cs2(i), &
               keytemp,keyerr,eosflag,eoskey,eos_rf_prec)
          if(keyerr.ne.0) then
             stop "problem in eos: cs2 at SetTimeStep"
          endif
       enddo
    else
       do i=1,n1
          keytemp = 0 ! not coming in with temperature (that would reset the energy), needs to be zero for the hybrid/poly/ideal EOS
          eosflag = 6 ! we want cs2 to be reset
          keyerr = 0
          call eos(i,rho(i),temp(i),ye(i),eps(i),cs2(i), &
               keytemp,keyerr,eosflag,eoskey,eos_rf_prec)
          if(keyerr.ne.0) then
             stop "problem in eos: cs2 at SetTimeStep"
          endif
       enddo
    endif
    
    dtp = dt
    
    dtnew = 1.0d0*time_gf
    do i=ghosts1+1, n1-ghosts1
       sound = sqrt(cs2(i))
       if (GR) then 
          if(.not.do_rotation) then
             dtnew = min(dtnew,  (x1(i+1)-x1(i)) /  & 
                  max(abs(v(i) + sound),abs(v(i) - sound)))
          else 
             dtnew = min(dtnew,  (x1(i+1)-x1(i)) /  & 
                  max(max(abs(v(i) + sound),abs(v(i) - sound)),&
                  max(abs(vphi(i)+sound),abs(vphi(i)-sound))))
          endif
       else
          if(.not.do_rotation) then
             dtnew = min(dtnew,  (x1(i+1)-x1(i)) /  & 
                  max(abs(v1(i) + sound),abs(v1(i) - sound)))
          else
             dtnew = min(dtnew,  (x1(i+1)-x1(i)) /  & 
                  max(max(abs(v1(i) + sound),abs(v1(i) - sound)),&
                  max(abs(vphi1(i)+sound),abs(vphi1(i)-sound))))
          endif

       endif
    enddo
       
    dt = min(cffac*dtnew,1.05d0*dtp)
   
  end subroutine SetTimeStep

end subroutine driver


