
       subroutine ratenuc

!------linkages.

!      called by - [subroutine] derivs         

                                      
!------modules.

       use bbnvar


       implicit none


!------remarks.

!      generates rate coefficients for reactions involving nuclides up to a = 9.


!------local variables

       real(dl) t913        !t9**(1/3)
       real(dl) t923        !t9**(2/3)
       real(dl) t943        !t9**(4/3)
       real(dl) t953        !t9**(5/3)
       real(dl) t912        !t9**(1/2)
       real(dl) t932        !t9**(3/2)
       real(dl) t9m1        !t9**(-1)
       real(dl) t9m23       !t9**(-2/3)
       real(dl) t9m32       !t9**(-3/2)
       real(dl) t9a         !for reaction 17.
       real(dl) t9a32       !t9a**(3/2)
       real(dl) t9b         !for reaction 18.
       real(dl) t9b32       !t9b**(3/2)
       real(dl) t9c         !for reaction 22
       real(dl) t9c13       !t9c**(1/3)
       real(dl) t9c56       !t9c**(5/6)
       real(dl) t9d         !for reaction 24.
       real(dl) t9d13       !t9d**(1/3)
       real(dl) t9d56       !t9d**(5/6)
       real(dl) t9e         !for reaction 26.
       real(dl) t9e13       !t9e**(1/3)
       real(dl) t9e56       !t9e**(5/6)
       real(dl) t9g         !for reaction 27.
       real(dl) t9g13       !t9g**(1/3)
       real(dl) t9g56       !t9g**(5/6)


!------procedure.

!10----temperature factors

       t913  = t9**(.33333333)      !t9**(1/3)
       t923  = t913*t913            !t9**(2/3)
       t943  = t923*t923            !t9**(4/3)
       t953  = t9*t923              !t9**(5/3)
       t912  = sqrt(t9)             !t9**(1/2)
       t932  = t9*t912              !t9**(3/2)
       t9m1  = 1/t9                 !t9**(-1)
       t9m23 = 1.0/t923             !t9**(-2/3)
       t9m32 = 1.0/t932             !t9**(-3/2)
       t9a   = t9/(1.0+13.076*t9)   !for reaction 17.
       t9a32 = t9a**(1.5)           !t9a**(3/2)
       t9b   = t9/(1.+49.18*t9)     !for reaction 18.
       t9b32 = t9b**(1.5)           !t9b**(3/2)
       if (t9.gt.10.) then          !for reaction 22.
         t9c = 1.
       else
         t9c = t9/(1.-9.69e-2*t9+2.84e-2*t953/(1.-9.69e-2*t9)**(2./3.))
       end if
       t9c13 = t9c**(.3333333)      !t9c**(1/3)
       t9c56 = t9c**(.8333333)      !t9c**(5/6)
       t9d   = t9/(1.+0.759*t9)     !for reaction 24.
       t9d13 = t9d**(.3333333)      !t9d**(1/3)
       t9d56 = t9d**(.8333333)      !t9d**(5/6)
       t9e   = t9/(1.+0.1378*t9)    !for reaction 26.
       t9e13 = t9e**(.3333333)      !t9e**(1/3)
       t9e56 = t9e**(.8333333)      !t9e**(5/6)
       t9g   = t9/(1.+0.1071*t9)    !for reaction 27.
       t9g13 = t9g**(.3333333)      !t9g**(1/3)
       t9g56 = t9g**(.8333333)      !t9g**(5/6)


!20----neutron, photon reactions

!------H(n,g)H2  Smith-Kawano-Malaney 1992
       f(12) = 4.742e+4*(1.-.8504*t912+.4895*t9-.09623*t932 + &
               8.471e-3*t9*t9-2.80e-4*t9*t932)

!------H2(n,g)H3  Wagoner 1969
       f(13) = 6.62e+1*(1.+18.9*t9)

!------He3(n,g)He4  Wagoner 1969
       f(14) = 6.62e+0*(1.+905.*t9)

!------Li6(n,g)Li7  Malaney-Fowler 1989
       f(15) = 5.10e+3


!30----neutron, proton reactions

!------He3(n,p)H3  Smith-Kawano-Malaney 1992
       f(16) = 7.21e+8*(1.-.508*t912+.228*t9)

!------Be7(n,p)Li7  Smith-Kawano-Malaney 1992
       f(17) = 2.675e+9*(1.-.560*t912+.179*t9-.0283*t932 + &
               2.214e-3*t9*t9-6.851e-5*t9*t932) + &
               9.391e+8*t9a32*t9m32 + &
               4.467e+7*t9m32*exp(-0.07486/t9)


!40----neutron, alpha reactions

!------Li6(n,a)H3  Caughlan-Fowler 1988
       f(18) = 2.54e+9*t9m32*exp(-2.39/t9) + 1.68e+8*(1.-.261*t9b32/t932)

!------Be7(n,a)He4  Wagoner 1969
       f(19) = 2.05e+4*(1.+3760.*t9)


!50----proton, photon reactions

!------H2(p,g)He3  Smith-Kawano-Malaney 1992
       f(20) = 2.65e+3*t9m23*exp(-3.720/t913)* &
               (1.+.112*t913+1.99*t923+1.56*t9+.162*t943+.324*t953)

!------H3(p,g)He4  Caughlan-Fowler 1988
       f(21) = 2.20e+4*t9m23*exp(-3.869/t913)* &
               (1.+.108*t913+1.68*t923+1.26*t9+.551*t943+1.06*t953)

!------Li6(p,g)Be7  Caughlan-Fowler 1988
       f(22) = 6.69e+5*t9c56*t9m32*exp(-8.413/t9c13)


!60----proton, alpha reactions

!------Li6(p,a)He3  Caughlan-Fowler 1988
       f(23) = 3.73e+10*t9m23*exp(-8.413/t913-(t9/5.50)**2)* &
               (1.+.050*t913-.061*t923-.021*t9+.006*t943+.005*t953) + &
               1.33e+10*t9m32*exp(-17.763/t9) + &
               1.29e+09*t9m1*exp(-21.820/t9)

!------Li7(p,a)He4  Smith-Kawano-Malaney 1992
       f(24) = 1.096e+9*t9m23*exp(-8.472/t913) - &
               4.830e+8*t9d56*t9m32*exp(-8.472/t9d13) + &
               1.06e+10*t9m32*exp(-30.442/t9) + &
               1.56e+5*t9m23*exp((-8.472/t913)-(t9/1.696)**2)* &
               (1.+.049*t913-2.498*t923+.860*t9+3.518*t943+3.08*t953) + &
               1.55e+6*t9m32*exp(-4.478/t9)


!70----alpha, photon reactions

!------H2(a,g)Li6  Caughlan-Fowler 1988
       f(25) = 3.01e+01*t9m23*exp(-7.423/t913)* &
               (1.+.056*t913-4.85*t923+8.85*t9-.585*t943-.584*t953) + &
               8.55e+1*t9m32*exp(-8.228/t9)

!------H3(a,g)Li7  Smith-Kawano-Malaney 1992
       f(26) = 3.032e+5*t9m23*exp(-8.090/t913)* &
               (1.+.0516*t913+.0229*t923+8.28e-3*t9 - &
               3.28e-4*t943-3.01e-4*t953) + &
               5.109e+5*t9e56*t9m32*exp(-8.068/t9e13)

!------He3(a,g)Be7  Smith-Kawano-Malaney 1992
       f(27) = 4.817e+6*t9m23*exp(-14.964/t913)* &
               (1.+.0325*t913-1.04e-3*t923-2.37e-4*t9 - &
               8.11e-5*t943-4.69e-5*t953) + &
               5.938e+6*t9g56*t9m32*exp(-12.859/t9g13)


!80----deuterium, neutron and deuterium, proton reactions

!------H2(d,n)He3  Smith-Kawano-Malaney 1992
       f(28) = 3.95e+8*t9m23*exp(-4.259/t913)* &
               (1.+.098*t913+.765*t923+.525*t9+9.61e-3*t943+.0167*t953)

!------H2(d,p)H3  Smith-Kawano-Malaney 1992
       f(29) = 4.17e+8*t9m23*exp(-4.258/t913)* &
               (1.+.098*t913+.518*t923+.355*t9-.010*t943-.018*t953)

!------H3(d,n)He4  Smith-Kawano-Malaney 1992
       f(30) = 1.063e+11*t9m23*exp(-4.559/t913-(t9/.0754)**2)* &
               (1.+.092*t913-.375*t923-.242*t9+33.82*t943+55.42*t953) + &
               8.047e+8*t9m23*exp(-0.4857/t9)

!------He3(d,p)He4  Smith-Kawano-Malaney 1992
       f(31) = 5.021e+10*t9m23*exp(-7.144/t913-(t9/.270)**2)* &
               (1.+.058*t913+.603*t923+.245*t9+6.97*t943+7.19*t953) + &
               5.212e+8/t912*exp(-1.762/t9)


!90-----three particle reactions

!------He3(he3,2p)He4  Caughlan-Fowler 1988
       f(32) = 6.04e+10*t9m23*exp(-12.276/t913)* &
               (1.+.034*t913-.522*t923-.124*t9+.353*t943+.213*t953)         

!------Li7(d,na)He4  Caughlan-Fowler 1988
       f(33) = 2.92e+11*t9m23*exp(-10.259/t913)

!------Be7(d,pa)He4  Caughlan-Fowler 1988
       f(34) = 1.07e+12*t9m23*exp(-12.428/t913)

       return


!------references 

! smith, m., kawano, l.h., and malaney, r.a., 1992, submitted to ap. j.        
! malaney, r.a., and fowler, w.a., 1989, astrophys. j., 345, l5.
! caughlan, g.r., and fowler, w.a., 1988, atomic data and nuclear data
! tables, 40, 283.
! wagoner, r.v.,1969, ap. j. suppl. no. 162, 18, 247.


       end subroutine ratenuc
