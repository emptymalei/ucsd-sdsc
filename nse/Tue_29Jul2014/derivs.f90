
       subroutine derivs(loop)

!------linkages.

!      called by - [subroutine] driver
!      calls     - [subroutine] therm, ratenp, ratenuc, sol


!-----remarks.

!     computes derivatives of
!       - temperature
!       - hv
!       - chemical potential
!       - abundances


!------modules.

       use bbnvar
       use nucsolver


       implicit none


       save


!------throughput variables.

       integer, intent(in) :: loop


!------local variables.

       integer i        !index
       real(dl) sumy         !sum of abundances.
       real(dl) sumzy        !sum of charge*abundances.
       real(dl) sumdy        !sum of abundance flows.
       real(dl) summdy       !sum of (mass excess)*(abundance flows).
       real(dl) sumzdy       !sum of (charge)*(abundance flows).
       real(dl) dphdt9       !d(phi e)/d(t9).
       real(dl) dphdln       !d(phi e)/d(h).
       real(dl) dphdzy       !d(phi e)/d(sumzy).
       real(dl) bar          !baryon density and pressure terms.


!------procedure.

!10----compute derivatives for abundances-------------

       rnb = hv*t9*t9*t9/rhob0 !baryon mass density (ratio to init value).

!------various thermodynamic quantities.
       call therm(loop)
       hubcst = sqrt(8._dl/3._dl*pi*g*thm(10) + cosmo/3._dl) !expansion rate.
       rhob = thm(9) !baryon mass density.

!------compute reaction rate coefficients.
       call ratenp !n <-> p rates.
       call ratenuc !forward rate for reactions with a < 10.

!------solve coupled differential equations.
       call sol(loop)
       if (mbad.gt.0) return !abort in case matrix not invertible.

!20----compute derivatives for temperature, hv, and chemical potential

!------accumulate to get sum.
       sumy = 0._dl
       sumzy = 0._dl
       sumdy = 0._dl
       summdy = 0._dl
       sumzdy = 0._dl
       do i = 1,isize
         sumy = sumy + y(i)           !sum of abundance.
         sumzy = sumzy + zm(i)*y(i)     !sum of charge*abundance.
         sumdy = sumdy + dydt(i)        !sum of abundance flow.
         summdy = summdy + dm(i)*dydt(i) !sum of (mass excess)*(abundance flow).
         sumzdy = sumzdy + zm(i)*dydt(i) !sum of (charge)*(abundance flow).
       end do

!------changes in temperature, hv, and chemical potential.
       dphdt9 = thm(12)*(-1.070e-4*hv*sumzy/t9 - thm(11))
       dphdln = -thm(12)*3.568e-5*hv*sumzy
       dphdzy = thm(12)*3.568e-5*hv
       bar = 9.25e-5*t9*sumy + 1.388e-4*t9*sumdy/(3.*hubcst) + &
             summdy/(3.*hubcst)
       dt9 = -(3.*hubcst)*(thm(1) + thm(3) + thm(4) + thm(7) + &
             thm(9)*bar + thm(6)*(dphdln + dphdzy*sumzdy/(3.*hubcst)))/ &
             (thm(2) + thm(5) + thm(6)*dphdt9 + thm(9)*1.388e-4*sumy)   !(ref 1).
       dlt9dt = dt9/t9
       dhv = -hv*((3.*hubcst) + 3.*dlt9dt)                    !(ref 2).
       dphie = dphdt9*dt9 + dphdln*(3.*hubcst) + dphdzy*sumzdy  !(ref 3).

       return

!----------references-----------------------------------------------
!     1)  kawano, l., 1992, fermilab preprint fermilab-pub-92/04-a,
!         kellogg radiation lab preprint oap-714,
!         equation d.35.
!     2)  kawano, l., 1992, preprint, equation d.19.
!     3)  kawano, l., 1992, preprint, equation d.20.


       end subroutine derivs
