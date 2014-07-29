
       subroutine ratedecay


!------linkages.
!      called by - [subroutine] start
!      calls     - none


!------remarks.
!      generates weak decay rates.


!------modules.

       use bbnvar


       implicit none


!------procedure.

!10----set decay rate coefficients

!------H3 -> e- + v + He3  Tilly-Weller-Hasan 1987
       f(2) = 1.79e-09_dl

!------Li8 -> e- + v + 2He4  Ajzenberg-Selove 1988
       f(3) = 8.27e-01_dl

!------B12 -> e- + v + C12  Ajzenberg-Selove 1990
       f(4) = 3.43e+01_dl

!------C14 -> e- + v + N14  Aajzenberg-Selove 1986
       f(5) = 3.834e-12_dl

!------B8 -> e+ + v + 2He4  Ajzenberg-Selove 1988
       f(6) = 9.00e-01_dl

!------C11 -> e+ + v + B11  Ajzenberg-Selove 1990
       f(7) = 5.668e-04_dl

!------N12 -> e+ + v + C12  Ajzenberg-Selove 1990
       f(8) = 6.301e+01_dl

!------N13 -> e+ + v + C13  Ajzenberg-Selove 1986
       f(9) = 1.159e-03_dl

!------O14 -> e+ + v + N14  Ajzenberg-Selove 1986
       f(10) = 9.8171e-03_dl

!------O15 -> e+ + v + N15  Ajzenberg-Selove 1986
       f(11) = 5.6704e-03_dl

       return      


!------references 
!      Ajzenberg-Selove, F., 1990, nucl. phys. a506, 1.
!      Ajzenberg-Selove, F., 1988, nucl. phys. a490, 1.
!      Ajzenberg-Selove, F., 1986, nucl. phys. a449, 1.
!      Tilley, D.R., Weller, H.R., and Hasan, H.H., 1987, nucl. phys. a474, 1.


       end subroutine ratedecay
