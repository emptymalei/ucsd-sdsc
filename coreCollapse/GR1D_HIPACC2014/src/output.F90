!-*-f90-*-
subroutine output_all(modeflag)

  use GR1D_module
  implicit none

  character*1024 filename
  character*256 basename
  integer modeflag

  integer i,tempghosts
  real*8, allocatable :: cs(:)
  real*8 alp_ana(n1),rho_ana(n1),vel_ana(n1),X_ana(n1)
  real*8 maxr
  
  integer, parameter :: nscalars0 = 64
  real*8 scalars(nscalars0)
  integer nscalars

  integer ishock_radius(1)

  if(modeflag.eq.0) then
     
  else if(modeflag.eq.1) then
     
     filename = trim(adjustl(outdir))//"/v1.xg"
     if (.not.small_output) call output_single(v1*clite,filename)
     
     if(do_rotation) then
        filename = trim(adjustl(outdir))//"/omega.xg"
        call output_single(omega*time_gf,filename)
        filename = trim(adjustl(outdir))//"/ToverW.xg"
        call output_single(ToverW,filename)        
        
        if(GR) then
           filename = trim(adjustl(outdir))//"/vphi.xg"
           call output_single(vphi,filename)
        else
           filename = trim(adjustl(outdir))//"/vphi1.xg"
           call output_single(vphi1,filename)
        endif
     endif
     
     filename = trim(adjustl(outdir))//"/rho.xg"
     call output_single(rho/rho_gf,filename)
     
     if (do_nupress) then
        filename = trim(adjustl(outdir))//"/nuchem.xg"
        if (.not.small_output) call output_single(nuchem,filename)
        filename = trim(adjustl(outdir))//"/press_nu.xg"
        if (.not.small_output) call output_single(press_nu/press_gf,filename)
        filename = trim(adjustl(outdir))//"/dnupdr.xg"
        if (.not.small_output) call output_single(dnupdr/press_gf*length_gf,filename)
     endif

     filename = trim(adjustl(outdir))//"/ye.xg"
     call output_single(ye,filename)
     
     filename = trim(adjustl(outdir))//"/press.xg"
     call output_single(press/press_gf,filename)
     
     filename = trim(adjustl(outdir))//"/eps.xg"
     if (.not.small_output) call output_single(eps/eps_gf,filename)
     
     if (GR) then
        filename = trim(adjustl(outdir))//"/mass_grav.xg"
        call output_single(mgrav/mass_gf,filename)
        filename = trim(adjustl(outdir))//"/mass_bary.xg"
        call output_single(mass/mass_gf,filename)
     else
        filename = trim(adjustl(outdir))//"/mass_bary.xg"
        call output_single(mass/mass_gf,filename)
     endif
     
     if (eoskey.eq.3) then
        filename = trim(adjustl(outdir))//"/xn.xg"
        if (.not.small_output) call output_single(massfrac_n,filename)
        
        filename = trim(adjustl(outdir))//"/xp.xg"
        if (.not.small_output) call output_single(massfrac_p,filename)
        
        filename = trim(adjustl(outdir))//"/xa.xg"
        if (.not.small_output) call output_single(massfrac_a,filename)
        
        filename = trim(adjustl(outdir))//"/xh.xg"
        if (.not.small_output) call output_single(massfrac_h,filename)
        
        filename = trim(adjustl(outdir))//"/xabar.xg"
        if (.not.small_output) call output_single(massfrac_abar,filename)
        
        filename = trim(adjustl(outdir))//"/xzbar.xg"
        if (.not.small_output) call output_single(massfrac_zbar,filename)
     endif

     if (eoskey.eq.1) then
        filename = trim(adjustl(outdir))//"/pressth.xg"
        if (.not.small_output) call output_single(pressth/press_gf,filename)
     endif

     filename = trim(adjustl(outdir))//"/eps_kin.xg"
     if (.not.small_output) call output_single(eps_kin/eps_gf,filename)

     filename = trim(adjustl(outdir))//"/cs.xg"
     allocate(cs(n1))
     cs(:) = sqrt(cs2(:))*clite
     if (.not.small_output) call output_single(cs,filename)
     deallocate(cs)
     
     if(GR) then
        if(initial_data.eq."OSC") then
           call analytic_OSC_alpha(time*time_gf,10.d0,1.0d0, &
                alp_ana,rho_ana,vel_ana,X_ana,maxr)
           filename = trim(adjustl(outdir))//"/alpha_analytic.xg"
           call output_singlemod(alp_ana,filename,maxr)
           filename = trim(adjustl(outdir))//"/rho_analytic.xg"
           call output_single(rho_ana/rho_gf,filename)
           filename = trim(adjustl(outdir))//"/vel_analytic.xg"
           call output_singlemod(vel_ana*clite,filename,maxr)
           filename = trim(adjustl(outdir))//"/alphamod.xg"
           call output_singlemod(alp,filename,maxr)
           filename = trim(adjustl(outdir))//"/X_analytic.xg"
           call output_single(X_ana,filename)
        endif
        filename = trim(adjustl(outdir))//"/alpha.xg"
        if (.not.small_output) call output_single(alp,filename)
        filename = trim(adjustl(outdir))//"/X.xg"
        if (.not.small_output) call output_single(X,filename)
        filename = trim(adjustl(outdir))//"/W.xg"
        if (.not.small_output) call output_single(W,filename)
        filename = trim(adjustl(outdir))//"/v.xg"
        call output_single(v*clite,filename)
     endif
     
     if(eoskey.eq.3) then
        filename = trim(adjustl(outdir))//"/entropy.xg"
        if (.not.small_output) call output_single(entropy,filename)
        filename = trim(adjustl(outdir))//"/temperature.xg"
        call output_single(temp,filename)
     endif
     
  else if(modeflag.eq.2) then
     
     if (initial_data.eq.'Collapse') then
        !Shock radius
        if (bounce) then
           filename = trim(adjustl(outdir))//"/shock_radius_t.dat"
           call output_scalar(shock_radius/length_gf,filename)
           filename = trim(adjustl(outdir))//"/binding_energy_total.dat"
           call output_scalar(binding_energy_total/energy_gf,filename)
        endif
        
        ! Mass values
        filename = trim(adjustl(outdir))//"/mgrav_Xmax.dat"
        call output_scalar(mgravX,filename)        
        
        filename = trim(adjustl(outdir))//"/mbary_shock.dat"
        call output_scalar(mass(ishock(1)),filename)       
        
        filename = trim(adjustl(outdir))//"/mgrav_shock.dat"
        call output_scalar(mgrav(ishock(1)),filename)                
        
        filename = trim(adjustl(outdir))//"/mgrav_rho1e12.dat"
        call output_scalar(mgrav12,filename)      
        
        filename = trim(adjustl(outdir))//"/mbary_Xmax.dat"
        call output_scalar(mbaryX,filename)      
        
        filename = trim(adjustl(outdir))//"/mbary_rho1e12.dat"
        call output_scalar(mbary12,filename)     
        
        filename = trim(adjustl(outdir))//"/r_Xmax.dat"
        call output_scalar(rXmax/length_gf,filename)      
        
        filename = trim(adjustl(outdir))//"/r_rho1e12.dat"
        call output_scalar(r12max/length_gf,filename)    
        
        filename = trim(adjustl(outdir))//"/M_innercore.dat"
        call output_scalar(mass_inner_core,filename)
        
        !rotation scalars
        if (do_rotation) then
           filename = trim(adjustl(outdir))//"/total_angular_momentum.dat"
           call output_scalar(angular_momentum/mass_gf*time_gf/length_gf**2,filename)
           filename = trim(adjustl(outdir))//"/ToverW_edge.dat"
           call output_scalar(ToverW(n1-ghosts1-1),filename)
        endif
     endif
        
     ! central values
     filename = trim(adjustl(outdir))//"/rho_c_t.dat"
     call output_central(rho/rho_gf,filename)
        
     filename = trim(adjustl(outdir))//"/ye_c_t.dat"
     call output_central(ye,filename)
        
     filename = trim(adjustl(outdir))//"/csound_c_t.dat"
     call output_central(sqrt(cs2),filename)
     
     if(eoskey.eq.3) then
        filename = trim(adjustl(outdir))//"/entropy_c_t.dat"
        call output_central(entropy,filename)
        filename = trim(adjustl(outdir))//"/temperature_c_t.dat"
        call output_central(temp,filename)
     endif
     
     filename = trim(adjustl(outdir))//"/totalmass.dat"
     call output_scalar(totalmass/mass_gf,filename)
     
     if(GR) then
        if(initial_data.eq."OSC") then
           call analytic_OSC_alpha(time*time_gf,10.d0,1.0d0, &
                alp_ana,rho_ana,vel_ana,X_ana,maxr)
           filename = trim(adjustl(outdir))//"/alpha_analytic_c_t.dat"
           call output_central(alp_ana,filename)
           filename = trim(adjustl(outdir))//"/rho_analytic_c_t.dat"
           call output_central(rho_ana/rho_gf,filename)
           filename = trim(adjustl(outdir))//"/vel_analytic_c_t.dat"
           call output_central(vel_ana*clite,filename)
        endif
        filename = trim(adjustl(outdir))//"/alpha_c_t.dat"
        call output_central(alp,filename)
        filename = trim(adjustl(outdir))//"/time_c.dat"
        call output_scalar(time_c,filename)
     endif
     
     if(initial_data.eq."Sedov") then
        ishock_radius = maxloc(abs(v1))
        filename = trim(adjustl(outdir))//"/Sedov_radius.dat"
        call output_scalar(x1(ishock_radius(1)),filename)
        filename = trim(adjustl(outdir))//"/Sedov_velocity.dat"
        call output_scalar(v1(ishock_radius(1)),filename)
        filename = trim(adjustl(outdir))//"/Sedov_press.dat"
        call output_scalar(press(ishock_radius(1)),filename)
        filename = trim(adjustl(outdir))//"/Sedov_density.dat"
        call output_scalar(rho(ishock_radius(1)),filename)
        
     endif
     
     if (initial_data.eq.'Collapse') then
        filename = trim(adjustl(outdir))//"/accretion_rates.dat"
        call output_accretion(accretion_rates,filename)
        filename = trim(adjustl(outdir))//"/accreted_mass.dat"
        call output_accretion(accreted_mass,filename)
     endif
     
  endif
  
end subroutine output_all

! *******************************************************************
subroutine output_accretion(var,filename)
      
  use GR1D_module, only: time
  
  implicit none
  real*8 var(*)
  character(*) filename
  
  open(unit=666,file=trim(adjustl(filename)),status="unknown",&
       form='formatted',position="append")
  
  write(666,"(1P20E18.9)") time, var(1), var(2), var(3), var(4), &
       var(5), var(6), var(7), var(8), var(9), var(10), var(11)
  
  close(666)
  
end subroutine output_accretion
! ******************************************************************
subroutine output_single(var,filename)
  
  use GR1D_module, only: vs_mass,x1,mass,time,n1,length_gf
  
  implicit none
  real*8 var(n1)
  character(len=100) filename
  integer nt
  integer i
  
  open(unit=666,file=trim(adjustl(filename)),status="unknown",&
       form='formatted',position="append")
  write(666,*) '"Time = ',time
  
  if(vs_mass) then
     do i=1,n1
        write(666,"(1P20E18.9)") mass(i),var(i)
     enddo
  else
     do i=1,n1
        write(666,"(1P20E18.9)") x1(i)/length_gf,var(i)
     enddo
  endif
  write(666,*) " "
  write(666,*) " "
  close(666)
  
end subroutine output_single

! ******************************************************************

subroutine output_central(var,filename)
  
  use GR1D_module, only: x1,time,n1,ghosts1
  
  implicit none
  real*8 var(n1)
  character(len=100) filename
  integer nt
  integer i
  
  open(unit=666,file=filename,status="unknown",form='formatted',position="append")
  
  write(666,"(1P20E18.9)") time,var(ghosts1+1)
  
  close(666)
  
end subroutine output_central

! *******************************************************************
subroutine output_scalar(var,filename)
  
  use GR1D_module, only: time
  implicit none
  real*8 var
  character(len=100) filename
  integer nt
  integer i
  
  open(unit=666,file=filename,status="unknown",form='formatted',position="append")
  
  write(666,"(1P20E18.9)") time,var
  
  close(666)

end subroutine output_scalar
! *******************************************************************
subroutine output_many_scalars(var,n0,n,filename)
  
  use GR1D_module, only: time
  implicit none
  integer n,n0
  real*8 var(n0)
  character(len=100) filename
  integer nt
  integer i
  
  open(unit=666,file=filename,status="unknown",form='formatted',position="append")
  
  write(666,"(1P64E18.9)") time,var(1:n)
  
  close(666)
  
end subroutine output_many_scalars
! *******************************************************************
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
  
  fname = trim(adjustl(outdir))//"/"//trim(adjustl(varname))//"_nt_"//outnt
  write(outtime,"(f11.7)") time
  
  fname = trim(adjustl(fname))//"_time_"//trim(adjustl(outtime))//".dat"
  
end subroutine generate_filename
! ******************************************************************
subroutine output_singlemod(var,filename,maxr)
  
  use GR1D_module, only: x1,time,n1,length_gf,ghosts1
  
  implicit none
  real*8 var(*)
  character(*) filename
  integer nt
  real*8 maxr
  integer i
  
  open(unit=666,file=trim(adjustl(filename)),status="unknown",&
       form='formatted',position="append")
  write(666,*) '"Time = ',time
  
  do i=ghosts1,n1
     if (x1(i).lt.maxr) then
        write(666,"(1P20E18.9)") x1(i),var(i)
     endif
  enddo
  write(666,*) " "
  write(666,*) " "
  close(666)
  
end subroutine output_singlemod

! ******************************************************************
