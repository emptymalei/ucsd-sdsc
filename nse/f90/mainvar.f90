
       module mainvar


!------linkages.

!      used by - [program] main


!-----remarks.

!     Contains a list of variables and related functions for main


       implicit none


!------computational numbers.

       integer, parameter :: dl = kind(1.d0)


!------numerical constants.

       real(dl), parameter :: pi = 3.1415926536_dl
       real(dl), parameter :: boltzmann = 8.6173324e+04_dl !eV/(10^9 K)


!------flags for output files.

       logical bbnflag !flag for bbn.dat
       logical equilflag !flag for equil.dat
       logical yeflag !flag for ye.dat
       integer runind


       end module mainvar
