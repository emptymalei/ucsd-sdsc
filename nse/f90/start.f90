
       subroutine start


!-----linkages.

!     called by - [subroutine] driver
!     calls     - [subroutine] rate1, bessel, rate0, calcst


!-----remarks.

!     sets initial conditions.


!------modules.

       use bbnvar
       use bessel


       implicit none


!------local variables.

       integer, parameter :: phspbins=100   !number of intervals for the phase-space integral
       integer i,n,m        !indicies
       real(dl) hubcst0       !defined by hubcst=(8/3*pi*g*rho+lambda/3)^0.5
       real(dl) ggs           !initial statistical weight in entropy
       real(dl) z             !defined by z = m(electron)*c**2/k*t9.
       real(dl) parfact	     !parity factor
       real(dl), dimension(phspbins) :: fermifactphsp
       real(dl) ediv
       real(dl) enu
       real(dl) ee
       real(dl) phspint
       real(dl) coulcorr
       character char1*200,char2*12,char3*12


!------procedure.

!10----initialize flags and counters.

       is = 1 !first iteration coming up.
       ip = inc !set to maximum allowed # of iteration.
       it = 0 !no accumulation yet.
       mbad = 0 !no computational errors.

!20----settings.

!------computational settings.

       t9 = t9i                    !initial temperature.
       t9mev = 8.6173324e-02_dl*t9i
       tcm = t9i                    !initial co-moving temperature.
       tcmev = 8.6173324e-02_dl*t9i
       ti = 1._dl/(const1*t9)**2       !initial time (ref 1).
       t = ti
       dt = dt1                    !initial time step.

!------model settings.

       hubcst0 = 608.5_dl     !initial expansion rate (s^-1) at t9~348?
       hubcst = hubcst0
       ggs = 11._dl/2._dl        !initial entropic statistical weight

       !Convert rho_{rad} into proper units at start of computation
       rhors = 232011._dl*rhors*tcmev**4*(11._dl/4._dl)**(4._dl/3._dl)

!------find constant in front of n<->p rates

       ediv = (deltamnp - xmelec)/float(phspbins)
       fermifactphsp(:) = 0.98
       phspint = 0.0
       do i=1,phspbins
         enu = (2.0*float(i)-1.0)*ediv/2.0
         ee = deltamnp - enu
         coulcorr = fermifactphsp(i)
         phspint = phspint + coulcorr*enu**2*ee*sqrt(ee**2 - xmelec**2)*ediv
       end do
!       cnorm=1.0/tau/phspint*1.106/1.130        !has units of 1/s/mev^5
       cnorm=1.0/tau/phspint       !has units of 1/s/mev^5

!30----compute initial abundances for neutron and proton

       y0(1) = 1.0/(exp(15.011/t9+xi(1))+1.) !initial n abundance (ref 3).
       y0(2) = 1.0/(exp(-15.011/t9-xi(1))+1.) !initial p abundance (ref 3).
       y(1) = y0(1)
       y(2) = y0(2)
       
!40----find ratio of baryon density to temperature cubed

       z = 5.930/t9            !inverse of temperature.
       call bessel_eval(z)
       hv = 3.3683e+04*eta1*2.75!(ref 4 but with final eta).
       seninit = 121310.0*ggs/2.0/hv
       sen = seninit
       dsendt = 0.0
       phiei = hv*(1.784e-5*y(2))/ &
              (.5*z**3*(blz(1)-2.*blz(2)+3.*blz(3)-4.*blz(4)+5.*blz(5)))
       phie = phiei !chemical potential of electron (ref 5).
                 
       rhob0 = hv*t9**3 !baryon density.
 
!50--------set abundances for rest of nuclides-------------------

       y0(3) = y(1)*y(2)*rhob0*exp(25.82/t9)/(.471e+10*t9**1.5)!(ref 7).
       y(3) = y0(3)
       do i = 4,isize
         y0(i)  = ytmin     !set rest to minimum abundance.
         y(i) = y0(i)      !init abundances at beginning of iteration.
       end do

       call ratedecay          !compute weak decay rates.
       
       return


!----------references-----------------------------------------------
!     1) wagoner, r.v., fowler, w.a., and hoyle, f. 1967, ap. j. 148,
!        page 44, equation a15.
!     2) coulomb correction obtained by dividing by correction factor fp(t9)
!        fp(t9) = 1 - 0.5(pi/(137<v>/c))
!        wagoner, r.v. 1973, ap. j. 179, page 358.
!     3) for the nondegenerate case:
!        wagoner, r.v., fowler, w.a., and hoyle, f. 1967, ap. j. 148,
!        page 4, equation 3.
!        for the case with neutrino degeneracy:
!        beaudet,g. and goret,p., 1976, astron. & astrophys., 49,
!        page 417, equation 9.
!     4) wagoner, r.v. 1969, ap j. suppl. no. 162, 18, page 250, equation 4.
!        3.3683e+4 = mu(ng/t9**3) with mu the atomic mass, ng the
!        photon density.  2.75 is for the 11/4 factor difference
!        between the initial and final values of eta.
!     5) kawano, l., 1992, fermilab preprint fermilab-pub-92/04-a,
!        kellogg radiation lab preprint oap-714.
!        equation d.2.
!     6) wagoner, r.v., fowler, w.a., and hoyle, f. 1967, ap. j. 148,
!        page 43, equation a4.
!        7.366 is used instead of 14.73 as the latter is the sum total
!        for 2 neutrino species.
!     7) initial deuterium abundance from nuclear statistical equilibrium
!        wagoner, r.v., fowler, w.a., and hoyle, f. 1967, ap. j. 148,
!        page 19, equation 17.


       end subroutine start
