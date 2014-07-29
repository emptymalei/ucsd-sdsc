
       subroutine accum


!------linkages.

!      called by - [subroutine] driver
!      calls     - [subroutine] nse


!------modules.

       use bbnvar
       use nse
       use ye


       implicit none


!-----remarks.

!     output accumulator.


!----------local variables.

       integer i,j        !indicies

       real(dl)	alln !quantities to determine n/p ratio
       real(dl)	allp !quantities to determine n/p ratio
       real(dl)	allnp !quantities to determine n/p ratio
       character char1*200,char2*12


!------procedure.

       it = it + 1                  !set up accumulation counter.


!40----set up output variables.

!------divide number fraction by that of proton.
       do i = 1,isize
         xout(it,i) = y(i)/y(2)
       end do
       xout(it,2) = y(2)*am(2)      !exception for proton.
       xout(it,6) = y(6)*am(6)      !exception for helium.

!------calculate n/p ratio for all n and p
       alln = 0._dl
       allp = 0._dl
       do i=1,isize
         alln = alln + (am(i) - zm(i))*y(i)
         allp = allp + zm(i)*y(i)
       end do
       allnp = alln/allp

!..........relabel temperature, time, thermodynamic variables, etc.
       t9out(it) = t9            !temperature.
       tcmevout(it) = tcmev            !co-moving temperature.
       tout(it) = t             !time.
       thmout(it,1) = thm(1)       !rho photon.
       thmout(it,2) = thm(4)        !rho electron.
       thmout(it,3) = thm(8)        !rho neutrino.
       thmout(it,4) = thm(9)        !rho baryon.
       thmout(it,5) = phie          !chemical potential.
       thmout(it,6) = thm(10)       !rho total.
       thmout(it,7) = allnp !n/p ratio
       thmout(it,8)   = hv/(3.3683d+4)!baryon to photon ratio.
       dtout(it)    = dt            !time step.
       senout(it)   = sen           !entropy ber baryon
       hubout(it)   = hubcst        !expansion rate.

!------NSE abundances.

       if (equilflag) then
         call equil(t9mev,y(1),y(2),rhob0*rnb)
       end if
       if (yeflag) then
         call ye_det(t9mev,alln,allp)
       end if

       return


       end subroutine accum
