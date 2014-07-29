
       module bessel


!-----linkages.

!     called by - [subroutine] driver
!     calls     - [subroutine] rate1, bessel, rate0, calcst


!-----remarks.

!     sets initial conditions.


!------modules.

       use bbnvar


       implicit none


!------evaluation of functions bl,bm,bn.

       real(dl) bk0,bk1,bk2,bk3,bk4 !values k0,k1,k2,k3,k4 as function of z.
       real(dl), dimension(5) :: blz !array containing values from function bl.
       real(dl), dimension(5) :: bmz !array containing values from function bm.
       real(dl), dimension(5) :: bnz !array containing values from function bn.


!------expansion coefficients.

       real(dl), dimension(7) ::  ci0 = (/ 1., & !expansion coefficients for i0 (z.le.2).
                  3.5156229,      3.0899424,      1.2067492, &
                  0.2659732,      0.0360768,      0.0045813/)

       real(dl), dimension(7) ::  ci1 = (/ 0.5, & !expansion coefficients for i1 (z.le.2).
                  0.87890594,     0.51498869,     0.15084934, &
                  0.02658733,     0.00301532,     0.00032411/)

       real(dl), dimension(7) ::  ck0 = (/-0.57721566, & !expansion coefficients for k0 (z.le.2).
                  0.42278420,     0.23069756,     0.03488590, &
                  0.00262698,     0.00010750,     0.00000740/)

       real(dl), dimension(7) ::  ck1 = (/ 1., & !expansion coefficients for k1 (z.le.2).
                  0.15443144,    -0.67278579,    -0.18156897, &
                 -0.01919402,    -0.00110404,    -0.00004686/)

       real(dl), dimension(7) ::  c0  = (/ 1.25331414, & !expansion coefficients for k0 (z.gt.2).
                 -0.07832358,     0.02189568,    -0.01062446, &
                  0.00587872,    -0.00251540,     0.00053208/)

       real(dl), dimension(7) ::  c1  = (/ 1.25331414, & !expansion coefficients for k1 (z.gt.2).
                  0.23498619,    -0.03655620,     0.01504268, &
                 -0.00780353,     0.00325614,    -0.00068245/)       

       contains

       subroutine bessel_eval(z)


!------linkages.

!      called by - [subroutine] start, therm
!      calls     - [subroutine] knux

!------modules.

!       use bbnvar


       implicit none


!------remarks.

!      evaluates functions bl(z), bm(z), and bn(z) using solutions to modified bessel functions.


!-----throughput variables.

      real(dl), intent(in) :: z             !defined by z = m(electron)*c**2/k*t9.


!------local variables.

       integer i         !indicies



!------procedure.

!20----calculate for 1 thru 5 z.

       do i=1,5
         call knux(float(i)*z)       !get k0(r),k1(r),k2(r),k3(r),k4(r),k(5).
         blz(i) = bk2/float(i)/z     !put value from function bl into array.
         bmz(i) = 0.25*(3.0*bk3 + bk1)/float(i)/z     !put value from function bm into array.
         bnz(i) = 0.5*(bk4 + bk2)/float(i)/z     !put value from function bn into array.
       end do

       return


       end subroutine bessel_eval


!-----------------------------------------------------------------------


       subroutine knux(z)

!------linkages.
!      called by - [subroutine] bessel
!      calls     - [function] exp

!------modules.

!       use bbnvar


       implicit none


!------remarks.

!     a subroutine for modified bessel functions of the second kind k-nu(z).


!------throughput variables.
       real(dl) z             !defined by z = m(electron)*c**2/k*t9.


!------local variables.
       integer i         !indicies
       real(dl) yz              !expansion variable = z/2.
       real(dl) tz              !expansion variable = z/3.75.
       real(dl) coeff          !logrithmic or exponential coefficient.
       real(dl) bi0 !values i0(z).
       real(dl) bi1 !values i1(z).


!------procedure.

!10----compute k0 and k1----------------------------------------

       if (z.le.2.) then            !(ref. 1).
!------compute factors.
         tz = (z/3.75)
         yz = (z/2.0)
         coeff = dlog(yz)

!------values for i0(z) and i1(z).
         bi0 = ci0(1)
         bi1 = ci1(1)
         bk0 = ck0(1)
         bk1 = ck1(1)
         do i = 2,7
           bi0 = bi0 + ci0(i)*tz**(2*(i-1))
           bi1 = bi1 + ci1(i)*tz**(2*(i-1))
           bk0 = bk0 + ck0(i)*yz**(2*(i-1))
           bk1 = bk1 + ck1(i)*yz**(2*(i-1))
         end do

!------values for k0(z) and k1(z).
         bk0 = -coeff*bi0 + bk0
         bk1 = coeff*bi1*z + bk1/z

       else !(z.le.2.)               !(ref. 2).

!------compute factors.
         yz = (2.0/z)
         coeff = (exp(-z)/sqrt(z))

!------values for k0(z) and k1(z).
         bk0 = c0(1)
         bk1 = c1(1)
         do i = 2,7
           bk0 = bk0 + c0(i)*yz**(i-1)
           bk1 = bk1 + c1(i)*yz**(i-1)
         end do
         bk0 = coeff*bk0
         bk1 = coeff*bk1

       end if !(z.le.2.)

!20----find k2, k3, and k4 by iteration (ref. 3)----------------

       bk2 = 2.0*(bk1/z) + bk0       !k2(z).
       bk3 = 4.0*(bk2/z) + bk1       !k3(z).
       bk4 = 6.0*(bk3/z) + bk2       !k4(z).

       return

!----------references-----------------------------------------------
!     handbook of mathematical functions (abramowitz and stegun),
!       dover publications, inc., new york
!       1) polynomial approximations for z.le.2
!         page 378, equations 9.8.1 and 9.8.3.
!         page 379, equations 9.8.5 and 9.8.7.
!       2) polynomial approximations for z > 2
!         page 379, equations 9.8.6 and 9.8.8.
!       3) recursion relation from 1st line of 9.6.26, page 376.


       end subroutine knux


       end module bessel
