
       subroutine therm(loop)

!------linkages.

!      called by - [subroutine] derivs
!      calls     - [subroutine] bessel, nudens, tevolst

!------remarks.

!      computes various temperature dependent thermodynamic quantities.


!------modules.

       use bbnvar
       use bessel


       implicit none


!------throughput variables.

       integer, intent(in) :: loop !rk loop number


!------local variables.

       integer i,m,n         !indicies
       real(dl) z             !defined by z = m(electron)*c**2/k*t9.
       real(dl) cosh1,cosh2,cosh3,cosh4,cosh5
       real(dl) sinh1,sinh2,sinh3,sinh4,sinh5
       real(dl) rhor


!------procedure division.

!10--------compute factors---------------------------------------

       z = 5.930/t9                 !z = m(electron)c**2/k(t9).
       t9mev = 8.6173324e-02*t9
       tcm = rnb**(1._dl/3._dl)*t9i   !co-moving temperature.
       tcmev = 8.6173324e-02*tcm
       
!..........trignometric function values.
       if (phie.le.17._dl) then        !no chance of overflow.
         cosh1 = cosh(phie)
         cosh2 = cosh(2._dl*phie)
         cosh3 = cosh(3._dl*phie)
         cosh4 = cosh(4._dl*phie)
         cosh5 = cosh(5._dl*phie)
         sinh1 = sinh(phie)
         sinh2 = sinh(2._dl*phie)
         sinh3 = sinh(3._dl*phie)
         sinh4 = sinh(4._dl*phie)
         sinh5 = sinh(5._dl*phie)
       else
         cosh1 = 0._dl
         cosh2 = 0._dl
         cosh3 = 0._dl
         cosh4 = 0._dl
         cosh5 = 0._dl
         sinh1 = 0._dl
         sinh2 = 0._dl
         sinh3 = 0._dl
         sinh4 = 0._dl
         sinh5 = 0._dl
       end if

       call bessel_eval(z)

!20--------compute thermodynamic variables--------------------

       thm(1)  = 8.418*t9*t9*t9*t9                        !(ref 1).
       thm(2)  = 4.0*thm(1)/t9                             !(ref 2).
       thm(3)  = thm(1)/3.0                              !(ref 3).
       thm(4)  = 3206.0*(bmz(1)*cosh1 - bmz(2)*cosh2 + bmz(3)*cosh3 &
                  - bmz(4)*cosh4 + bmz(5)*cosh5)                !(ref 4).
       thm(5)  = 3206.0*(z/t9)*(bnz(1)*cosh1 - 2.0*bnz(2)*cosh2 &
            + 3.0*bnz(3)*cosh3 - 4.0*bnz(4)*cosh4 + 5.0*bnz(5)*cosh5) !(ref 5).
       thm(6)  = 3206.0*(bmz(1)*sinh1 - 2.0*bmz(2)*sinh2 + 3.0*bmz(3)*sinh3 &
            - 4.0*bmz(4)*sinh4 + 5.0*bmz(5)*sinh5)                !(ref 6).
       thm(7)  = 3206.0*(blz(1)*cosh1/z - blz(2)*cosh2/(2.0*z) &
                 + blz(3)*cosh3/(3.0*z) - blz(4)*cosh4/(4.0*z) &
                 + blz(5)*cosh5/(5.0*z))                      !(ref 7).
       rhonu = 7._dl/8._dl*pi**2/30._dl*2._dl*tcmev**4*3._dl !units of mev^4
       thm(8) = 232011.0*rhonu        !units of g/cm^3
       thm(9)  = rhob0*rnb
       rhor = rhors*rnb**(4._dl/3._dl)
       !(ref 9).
       thm(10) = thm(1) + thm(4) + thm(8) + thm(9) &
                 + rhor
       !(ref 10).
       thm(11) = -(z**3/t9)*(sinh1*(3.0*blz(1) - z*bmz(1)) - sinh2*(3.0*blz(2) - &
                 2.0*z*bmz(2)) + sinh3*(3.0*blz(3) - 3.0*z*bmz(3)) - sinh4* &
                 (3.0*blz(4) - 4.0*z*bmz(4)) + sinh5*(3.0*blz(5) - 5.0*z*bmz(5)))   !(ref 11).
       thm(12) = z**3*(cosh1*blz(1) - 2.0*cosh2*blz(2) + &
                 3.0*cosh3*blz(3) - 4.0*cosh4*blz(4) + 5.0*cosh5*blz(5))   !(ref 12). 
       if (thm(12).ne.0._dl) thm(12) = 1._dl/thm(12)
       thm(13) = 1.000 + 0.565/z - 6.382/z**2 + 11.108/z**3 &!this is now an artifact         
                  + 36.492/z**4 + 27.512/z**5                  !(ref 13).
       thm(14) = (5.252/z - 16.229/z**2 + 18.059/z**3 + 34.181/z**4 &!this is now an artifact 
                   + 27.617/z**5)*exp(-q*z)

       return

!----------references and notes----------------------------------------
!     1)  thm(1)  = rho photon
!         (wagoner, r.v., fowler, w.a., and hoyle, f. 1967, ap. j. 148,
!          page 43, equation a2.)
!     2)  thm(2)  = d(rho photon)/d(t9)
!     3)  thm(3)  = (p photon)/c**2
!         (wagoner, r.v., fowler, w.a., and hoyle, f. 1967,
!          page 43, equation a3.)
!     4)  thm(4)  = rho electron+positron
!         (fowler, w.a. and hoyle, f., 1964, ap. j. suppl. no. 91, 9,
!          page 281, equation b44.)
!     5)  thm(5)  = d(rho electron+positron)/d(t9)
!     6)  thm(6)  = d(rho electron+positron)/d(phi e)
!     7)  thm(7)  = (p electron+positron)/c**2
!         (fowler, w.a. and hoyle, f., 1964, ap. j. suppl. no. 91, 9,
!          page 279, equation b27.)
!     8)  thm(8)  = rho neutrino
!                 = # neutrino species x rho electron neutrino (nondegenerate)
!                 = rho nu(e) + rho nu(m) + rho nu(t)          (degenerate)
!     9)  thm(9)  = rho baryon
!     10) thm(10) = rho total
!                 = rho photon + rho electron+positron + rho neutrino
!                              + rho baryon
!     11) thm(11) = d     /pi**2(hbar*c)**3(ne- - ne+)*z**3\
!                   d(t9) \  2  (mc**2)**3                 /
!     12) thm(12) = d        /pi**2(hbar*c)**3(ne- - ne+)*z**3\
!                   d(phi e) \  2  (mc**2)**3                 /
!     13) thm(13) = rate for n->p
!     14) thm(14) = rate for p->n


       end subroutine therm
