
       subroutine driver

!------linkages.

!      called by - [subroutine] run
!      calls     - [subroutine] start, derivs, accum


!------remarks.

!      runge-kutta computational routine


!------remarks.

       use bbnvar


       implicit none


       save


!------local variables.

       integer loop	   !rk loop number
       integer i,m,n        !indicies
       real(dl) dtlim,dtmin,dtl        !time step from limitation on abund changes.
       real(dl) t90          !temperature (in units of 10**9 k) at beginning of iteration.
       real(dl) hv0          !defined by hv = m(atomic)n(baryon)/t9**3 at beginning of iteration.
       real(dl) phie0        !chemical potential for electron at beginning of iteration.
       real(dl) dt90                  !change in temperature at beginning of iteration.
       real(dl) dhv0                  !change in hv at beginning of iteration.
       real(dl) dphie0                !change in chemical potential at beginning of iteration.
       real(dl), dimension(nnuc) :: dydt0 !change in rel number abundances at beginning of iteration.


!10----input initialization information, relabel-------------

       idc = 1
       call start           !input initialization information.

       dt = 1.e-10_dl

!20----loop one-------------------------------------------------

       do                     !begin runge-kutta looping.

       loop = 1                     !loop indicator.

!------compute derivatives of variables to be evolved.
       call derivs(loop)

!------accumulate.
       if ((t9.le.t9f).or.(dt.lt.abs(dtlow/dlt9dt))) then
         call accum
         exit
       end if

       if (ip.eq.inc) then
         call accum
         if (it.eq.itmax) exit
       end if                 
     
!     the above refer to low temp, small dt, enough iterations.
     
       dtlim = 0.1_dl*dt

!------adjust time step.
       if (is.gt.3) then           
        !adjust time step after 3 iterations.
         dtmin = abs(1._dl/dlt9dt)*ct  
         !trial value for minimum time step (ref 1).
         do i = 1,isize             
         !go through all abundance changes.
           if ((dydt(i).ne.0.).and.(y(i).gt.ytmin)) then
             dtl = abs(y(i)/dydt(i))*cy &
                   *(1.+(dlog10(y(i))/dlog10(ytmin))**2)  !(ref 2).
             if (dtl.lt.dtmin) dtmin = dtl        
              !find smallest time step.
           end if
         end do
         if (dtmin.gt.1.5*dt) then
           dtmin = 1.5*dt
          !limit change in time step.
         end if
         if ((dtmin.lt.dt).and.(idc.gt.3)) then
           dtmin = dt
         end if
         dt = dtmin                 !set new time step.
       end if

       t = t + dt                   !increment time.

!------store and increment values (ref 3).
       t90 = t9
       dt90 = dt9
       t9 = t90 + dt90*dt
       hv0 = hv
       dhv0 = dhv
       hv = hv0 + dhv0*dt
       phie0 = phie
       dphie0 = dphie
       phie = phie0 + dphie0*dt
       do i = 1,nnuc
         y0(i)    = y(i)
         dydt0(i) = dydt(i)
         y(i)     = y0(i) + dydt0(i)*dt
         if ((y(i).lt.ytmin)) y(i) = ytmin  
       end do
       sen0 = sen
       dsendt0 = dsendt
       sen = sen0 + dsendt0*dt

!30----loop two----------------------------------------------

       loop = 2                     !step up loop counter.

!------compute derivatives of variables to be evolved.
       call derivs(loop)

!------increment values.
       t9 = t90 + 0.5*(dt9 + dt90)*dt
       hv = hv0 + 0.5*(dhv + dhv0)*dt
       phie = phie0 + 0.5*(dphie + dphie0)*dt
       do i = 1,nnuc
         y(i) = y0(i) + 0.5*(dydt(i) + dydt0(i))*dt
         if ((y(i).lt.ytmin)) y(i) = ytmin  
         !set at minimum value.
       end do
       sen = sen0 + 0.5*(dsendt + dsendt0)*dt

!------reset counters.
       if (ip.eq.inc) then          !reset iteration counters.
         ip = 0
       end if
       ip = ip + 1
       is = is + 1
       idc = idc + 1

       end do

       
       return


!------references
!     1)  constraint on dt from the requirement that
!         (d(t9)/dt)*(dt/t9) < ct
!         wagoner, r.v. 1969, ap j. suppl. no. 162, 18, page 293,
!         equation c6.
!     2)  constraint on dt from
!         dtl < y/(dy/dt)*cy*(1+(log(y)/log(ytmin))**2)
!         wagoner, r.v. 1969, page 293, equation c7 but with log 
!         term squared.
!     3)  wagoner, r.v. 1969, page 292, equations c1, c2.


       end subroutine driver
