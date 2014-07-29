
       module nucsolver


!------modules.

       use bbnvar


!------components of matrix equation.

       double precision, dimension(nnuc,nnuc) :: a !relates y(t-dt) to y(t).
       real(dl), dimension(nnuc) :: b              !contains y0 in inverse order.
       real(dl), dimension(nnuc) :: yx            !yy in reverse order.


!------number of nuclides in reaction types 1-13.

       real(dl), dimension(13) :: si = & !# of nuclide i
                                  (/1.,1.,1.,1.,1.,2.,3.,2.,1.,1.,2.,1.,1./)

       real(dl), dimension(13) :: sj = & !# of nuclide j
                                  (/0.,1.,1.,0.,1.,0.,0.,1.,1.,1.,0.,0.,0./)

       real(dl), dimension(13) :: sk = & !# of nuclide k
                                  (/0.,0.,1.,0.,0.,1.,0.,0.,1.,0.,2.,1.,1./)

       real(dl), dimension(13) :: sl = & !# of nuclide l
                                  (/1.,1.,1.,2.,2.,1.,1.,1.,2.,3.,1.,1.,2./)


       contains


       subroutine sol(loop)

!------linkages.

!      called by - [subroutine] derivs
!      calls     - [subroutine] eqslin


!------modules.

       use bbnvar


       implicit none


!------remarks.

!      computes reverse strong and electromagnetic reaction rates.
!      fills and solves matrix equation for dydt(i).


!------input variables.

       integer loop	   !rk loop number


!------local variables.

       integer ierror               !element which does not converge.
       integer i,j,l,k,n         !indicies
       integer ind                  !equate to iform.
       integer isize1               !isize1 = isize - 1
       integer i1,j1                !i1 = isize1 - i, j1 = isize1 - j
       real(dl) ri,rj,rk,rl          !equate to si,sj,sk,sl.
       real(dl) ci,cj,ck,cl          !coefficients of rate equation.
       real(dl), dimension(nnuc) :: yy            !abundances at end of iteration.
       real(dl) bdln                 !(10**(-5))*volume expansion rate.         


!------procedure.

!10----temperature factors and initial values-------------------

!------matrix size.
       isize1 = isize + 1

!------initialize a-matrix.
       do i = 1,isize
         do j = 1,isize
           a(j,i) = 0.d0            !set a-matrix to zero.
         end do
       end do

!20----compute factors for the a-matrix-------------------------

       do n = 1,jsize

!------equate variables to arrays.
         ind = iform(n)             !type of reaction.
         i = ii(n)                  !id # of incoming nuclide i.
         j = jj(n)                  !id # of incoming nuclide j.
         k = kk(n)                  !id # of outgoing nuclide k.
         l = ll(n)                  !id # of outgoing nuclide l.

         if ((ind.ne.0).and.(i.le.isize).and.(l.le.isize)) then  
           ri = si(ind)             !# of incoming nuclide i.
           rj = sj(ind)             !# of incoming nuclide j.
           rk = sk(ind)             !# of outgoing nuclide k.
           rl = sl(ind)             !# of outgoing nuclide l.

!------compute different reaction rates.

           if (ind.eq.1) then
             ci = f(n)              !(ref 1).
             cj = 0._dl
             ck = 0._dl
             cl = r(n)
           else if (ind.eq.2) then !1-1-0-1 configuration.
             r(n) = rev(n)*1.e+10_dl*t9**1.5*exp(-q9(n)/t9)*f(n)  !(ref 2).         
             f(n) = rhob*f(n)
             ci = y(j)*f(n)/2._dl
             cj = y(i)*f(n)/2._dl
             ck = 0._dl
             cl = r(n)
           else if (ind.eq.3) then !1-1-1-1 configuration.
             f(n) = rhob*f(n)
             r(n) = rev(n)*exp(-q9(n)/t9)*f(n)  !(ref 3).
             ci = y(j)*f(n)/2._dl
             cj = y(i)*f(n)/2._dl
             ck = y(l)*r(n)/2._dl
             cl = y(k)*r(n)/2._dl
           else if (ind.eq.4) then !1-0-0-2 configuration.
             ci = f(n)
             cj = 0._dl
             ck = 0._dl
             cl = y(l)*r(n)/2._dl
           else if (ind.eq.5) then !1-1-0-2 configuration.
             f(n) = rhob*f(n)
             r(n) = rev(n)*exp(-q9(n)/t9)*f(n)  !(ref 3).
             ci = y(j)*f(n)/2._dl
             cj = y(i)*f(n)/2._dl
             ck = 0._dl
             cl = y(l)*r(n)/2._dl
           else if (ind.eq.6) then !2-0-1-1 configuration.
             f(n) = rhob*f(n)
             r(n) = rev(n)*exp(-q9(n)/t9)*f(n)  !(ref 3).
             ci = y(i)*f(n)/2._dl
             cj = 0._dl
             ck = y(l)*r(n)/2._dl
             cl = y(k)*r(n)/2._dl
           else if (ind.eq.7) then !3-0-0-1 configuration.
             r(n) = rev(n)*1.e+20_dl*t9**3.0*exp(-q9(n)/t9)*f(n)  
             !(ref 4).
             f(n) = rhob*rhob*f(n)
             ci = y(i)*y(i)*f(n)/6._dl
             cj = 0._dl
             ck = 0._dl
             cl = r(n)
           else if (ind.eq.8) then !2-1-0-1 configuration.
             r(n) = rev(n)*1.e+20_dl*t9**3.0*exp(-q9(n)/t9)*f(n)  
             !(ref 4).
             f(n) = rhob*rhob*f(n)
             ci = y(j)*y(i)*f(n)/3._dl
             cj = y(i)*y(i)*f(n)/6._dl
             ck = 0._dl
             cl = r(n)
           else if (ind.eq.9) then !1-1-1-2 configuration.
             f(n) = rhob*f(n)
             r(n) = rev(n)*1.e-10_dl*t9**(-1.5)*rhob*exp(-q9(n)/t9)*f(n)  
             !(ref 5).
             ci = y(j)*f(n)/2._dl
             cj = y(i)*f(n)/2._dl
             ck = y(l)*y(l)*r(n)/6._dl
             cl = y(k)*y(l)*r(n)/3._dl
           else if (ind.eq.10) then !1-1-0-3 configuration.
             f(n) = rhob*f(n)
             r(n) = rev(n)*1.e-10_dl*t9**(-1.5)*rhob*exp(-q9(n)/t9)*f(n)  
             !(ref 5).
             ci = y(j)*f(n)/2._dl
             cj = y(i)*f(n)/2._dl
             ck = 0._dl
             cl = y(l)*y(l)*r(n)/6._dl
           else if (ind.eq.11) then !2-0-2-1 configuration.
             f(n) = rhob*f(n)
             r(n) = rev(n)*1.e-10_dl*t9**(-1.5)*rhob*exp(-q9(n)/t9)*f(n) 
             !(ref 5).
             ci = y(i)*f(n)/2._dl
             cj = 0._dl
             ck = y(l)*y(k)*r(n)/3._dl
             cl = y(k)*y(k)*r(n)/6._dl
           else if (ind.eq.12) then !1-0-1-1 configuration.
             ci = f(n)
             cj = 0._dl
             ck = 0._dl
             cl = 0._dl
           else if (ind.eq.13) then !1-0-1-2 configuration.
             ci = f(n)
             cj = 0._dl
             ck = 0._dl
             cl = 0._dl
           end if


!30----construct the a-matrix-----------------------------------

           i = isize1 - i           !invert i index.
           j = isize1 - j           !invert j index.
           k = isize1 - k           !invert k index.
           l = isize1 - l           !invert l index.

!------fill i nuclide column.
           if (j.le.isize) a(j,i) = a(j,i) +  rj*ci
           if (k.le.isize) a(k,i) = a(k,i) -  rk*ci
           a(i,i) = a(i,i) +  ri*ci
           a(l,i) = a(l,i) -  rl*ci

!------fill j nuclide column.
           if (j.le.isize) then
             a(j,j) = a(j,j) +  rj*cj
             if (k.le.isize) a(k,j) = a(k,j) -  rk*cj
             a(i,j) = a(i,j) +  ri*cj
             a(l,j) = a(l,j) -  rl*cj
           end if

!------fill k nuclide column.
           if (k.le.isize) then
             if (j.le.isize) a(j,k) = a(j,k) -  rj*ck
             a(k,k) = a(k,k) +  rk*ck
             a(i,k) = a(i,k) -  ri*ck
             a(l,k) = a(l,k) +  rl*ck
           end if

!------fill l nuclide column.
           if (j.le.isize) a(j,l) = a(j,l) -  rj*cl
           if (k.le.isize) a(k,l) = a(k,l) +  rk*cl
           a(i,l) = a(i,l) -  ri*cl
           a(l,l) = a(l,l) +  rl*cl

         end if !((ind.ne.0).and.(i.le.isize).and.(l.le.isize))

       end do !n = 1,jsize

!40----put a-matrix and b-vector in final form of matrix equation--

       bdln   = 1.0e-05_dl*(3._dl*hubcst)   !(10**(-5))*(expansion rate).

       do i = 1,isize
         i1 = isize1 - i            !invert the rows.
         do j = 1,isize
           j1 = isize1 - j          !invert the columns.
           if (dabs(a(j,i)).lt.bdln*y0(j1)/y0(i1)) then
             a(j,i) = 0.d0          !set 0 if tiny.
           else
             a(j,i) = a(j,i)*dt     !bring dt over to other side.
           end if
         end do
         a(i,i) = 1.d0 + a(i,i)     !add identity matrix to a-matrix.
         b(i1)  = y0(i)             !initial abundances.
       end do

!50----solve equations to get derivative--------------------------

!------set monitor flag and solve by gaussian elimination.
       if (loop.eq.1) then
         call eqslin(ip,ierror)
       else
         call eqslin(0,ierror)
       end if

!------obtain derivative.
       do i = 1,isize
         yy(i)   = yx(isize1 - i)     !abundance at t+dt.
         dydt(i) = (yy(i) - y0(i))/dt         !take derivative.
       end do

!60----possible error messages and exit------------------------

       if (mbad.ne.0) then          !problem in gaussian elimination.
         if (mbad.eq.-1) write (iw,6000) ierror !error message.
         if (mbad.ge. 1) write (iw,6002) mbad   !error message.
6000     format (' ','** y(', i2, ') fails to converge **')
6002     format (' ','** ', i2, ' th diagonal term equals zero **')
       end if

       return

!------references----------------------------------------------
!     1) the coefficients are given in general as:
!             ci = ri*(y(j)**rj)*(y(i)**(ri-1)*f(n)/
!                  ((ri+rj)*fac(ri)*fac(rj))
!             cj = rj*(y(i)**ri)*(y(j)**(rj-1)*f(n)/
!                  ((ri+rj)*fac(ri)*fac(rj))
!             ck = rk*(y(l)**rl)*(y(k)**(rk-1)*f(n)/
!                  ((rk+rl)*fac(rk)*fac(rl))
!             cl = rl*(y(k)**rk)*(y(l)**(rl-1)*f(n)/
!                  ((rk+rl)*fac(rk)*fac(rl))
!        in which fac(x) is the factorial of x.
!     2) form of reverse rate given in
!        wagoner, r.v.,1969, ap. j. suppl. no. 162, 18, 247,
!          tables 1b, 4b, 7b.
!     3) form of reverse rate given in
!        wagoner, r.v.,1969, ap. j. suppl. no. 162, 18, 247,
!          tables 2b, 3b, 5b, 6b, 8b, 9b, 10b.
!     4) form of reverse rate given in
!        wagoner, r.v.,1969, ap. j. suppl. no. 162, 18, 247,
!          table 11b.
!     5) form of reverse rate given in
!        wagoner, r.v.,1969, ap. j. suppl. no. 162, 18, 247,
!          tables 12b, 13b.


       end subroutine sol


!-----------------------------------------------------------------------


       subroutine eqslin(icnvm,ierror)

!------linkages.

!      called by - [subroutine] sol
!      calls     - none


!------remarks.

!      solves for new abundances using gaussian elimination with back substitution, no pivoting.


!------remarks.

       use bbnvar


       implicit none


!------throughput variables.

       integer, intent(in) :: icnvm                !convergence monitor.
       integer, intent(out) :: ierror               !element which does not converge.


!------local variables.

       integer i,j,k                !indicies
       integer nord                 !order of correction.
       real(dl) rerror               !find error in right-hand vector
       double precision, dimension(nnuc,nnuc) :: a0 !coefficient array w/o manipulation.       
       double precision, dimension(nnuc) :: xx     !right-hand vector.
       double precision cx     !scaling factor in triangularization.
       double precision sum    !sum for backsubstitution.
       real(dl) xdy              !relative error.


!------procedure.

!10----initialize vector--------------------------------------

!------set counters to zero.
       nord = 0 !no corrections yet.
       mbad = 0 !no errors yet.

!------set right-hand and solution vectors to initial values.
       do i = 1,isize
         xx(i) = b(i) !right-hand vector.
         yx(i) = 0._dl !solution vector.
       end do

!------save matrix.
       if (icnvm.eq.inc) then       !monitor convergence.
         do i = 1,isize
           do j = 1,isize
             a0(j,i) = a(j,i)       !initial value of coefficient array.        
           end do
         end do
       end if

!20--------triangularize matrix and save operator-----------------------

!------check to see that there are no zeroes at pivot points.

       do i = 1,isize-1
 
        if (a(i,i).eq.0.d0) then !don't want to divide by zero.
           mbad = i !position of zero coefficient.
           return !terminate matrix evaluation.
         end if

!------triangularize matrix.

         do j = i+1,isize

           if (a(j,i).ne.0.d0) then 
           !progress diagonally down the column.
             cx = a(j,i)/a(i,i)     !scaling factor down the column.
             do k = i+1,isize       !progress diagonally along row.
               a(j,k) = a(j,k) - cx*a(i,k)  
               !subtract scaled coeff along row.
             end do
             a(j,i) = cx            !scaled coefficient.
!-------operate on right-hand vector.
             xx(j) = xx(j) - cx*xx(i)  !subtract off scaled coefficient.
           end if

         end do

       end do

!30--------do back substitution-----------------------------------

300    continue
       xx(isize) = xx(isize)/a(isize,isize)   
       !solution for ultimate position.
       yx(isize) = yx(isize) + xx(isize)
       do i = isize-1,1,-1          !from i = penultimate to i = 1.
         sum = 0.d0
         do j = i+1,isize
           sum = sum + a(i,j)*xx(j)  !sum up all previous terms.
         end do
         xx(i) = (xx(i) - sum)/a(i,i)
         yx(i) = yx(i) + xx(i)         !add difference to initial value.
       end do

!40----tests and exits---------------------------------------

       if (icnvm.eq.inc) then

         do i = 1,isize

           if (yx(i).ne.0._dl) then
             xdy = dabs(xx(i)/yx(i))  !relative error.

             if (xdy.gt.epstol) then

               if (nord.lt.mord) then !continue to higher orders.
                 nord = nord + 1

!------find error in right-hand vector.

                 do j = 1,isize
                   rerror = 0.d0         !initialize r.
                   do k = 1,isize
                     rerror = rerror + a0(j,k)*yx(k) 
                     !left side with approximate solution.
                   end do
                   xx(j) = b(j) - rerror  
                   !subtract difference from right side.
                 end do

!------operate on right-hand vector.

                 do j = 1,isize-1
                   do k = j+1,isize
                     xx(k) = xx(k) - a(k,j)*xx(j)  
                     !subtract off scaled coefficient.
                   end do
                  
                 end do

                 go to 300          !go for another iteration.

               else
!------not enough convergence.
                 mbad = -1          !signal error problem.
                 ierror = i         !ith nuclide for which x/y checked.
                 return

               end if !(nord.lt.mord)

             end if !(xdy.gt.epstol)

           end if !(y(i).ne.0)

         end do !i = 1,isize

       end if !(icnvm.eq.inc)

       return               !no more iterations & relative error small.


       end subroutine eqslin


       end module nucsolver
