
       program main


!------linkages.

!      called by - none
!      calls     - [subroutine] bbn


!-----remarks.

!     Program to call both bbn.


!------modules.

       use mainvar


       implicit none


       save


!------local variables.

       integer barind, counter, i

       integer barprec
       real(dl) barlow
       real(dl) barhigh
       logical barflag

       real(dl) obh2 !ob*h**2
       real(dl) :: obh2p = 0.022068_dl !ob*h**2 (from Planck)

       real(dl) ypri !4He Mass fraction
       real(dl) yd !D relative abundance
       real(dl) yhe3 !3He relative abundance
       real(dl) yli7 !7Li relative abundance
       real(dl), dimension(:,:), allocatable :: runarray

       logical filefound

       real(dl) time1, time2


!------procedure.


!------set values from main_params.ini.

       open (unit=10, file='main_params.ini', status='unknown')

       read (10,*) barprec
       read (10,*) barlow
       read (10,*) barhigh
       read (10,*) bbnflag
       read (10,*) equilflag
       read (10,*) yeflag
       read (10,*) barflag

       close (unit=10)

       allocate(runarray(5,barprec))


!------open output files.

       if (bbnflag) then

       inquire(file='bbn.dat', exist=filefound)
       if (filefound) then
         open(unit=21, file='bbn.dat', form='unformatted')
         close(21, status='delete')
       end if
       
       open (unit=21, file='bbn.dat', status='new')  !output file.

       end if !bbnflag


       if (equilflag) then

       inquire(file='equil.dat', exist=filefound)
       if (filefound) then
         open(unit=22, file='equil.dat', form='unformatted')
         close(22, status='delete')
       end if
       
       open (unit=22, file='equil.dat', status='new')  !output file.

       end if !equilflag


       if (yeflag) then

       inquire(file='ye.dat', exist=filefound)
       if (filefound) then
         open(unit=23, file='ye.dat', form='unformatted')
         close(23, status='delete')
       end if
       
       open (unit=23, file='ye.dat', status='new')  !output file.

       end if !yeflag


       if (barflag) then

       inquire(file='bar.dat', exist=filefound)
       if (filefound) then
         open(unit=24, file='bar.dat', form='unformatted')
         close(24, status='delete')
       end if
       
       open (unit=24, file='bar.dat', status='new')  !output file.


       end if !barflag


!------get time at start of computation.

       call cpu_time(time1) 

!------run through cases with different baryon numbers.

       runarray = 0._dl

       runind = 1
       do barind=1,barprec

         !Determine \Omega_b*h^2:
         if (barprec.eq.1) then
           obh2 = barlow
         else
           obh2 = real(barind-1, kind=dl)*(barhigh - barlow)/real(barprec - 1, kind=dl) &
                  + barlow
         end if

         if (equilflag)  write (22,10) runind, obh2
10       format (/,'run: ',i6,'; obh2 = ',1pe9.3,/,4x,'Temp',10x,'D/H',9x, &
                 'T/H',8x,'3He/H',7x,'4He/H',/,60('-'))
         if (yeflag)  write (23,12) runind, obh2
12       format (/,'run: ',i6,'; obh2 = ',1pe9.3,/,5x,'Temp',8x,'ye_c',8x, &
                 'ye_a',/,36('-'))

         call bbn(obh2,ypri,yd,yhe3,yli7)

         runarray(1,runind) = obh2
         runarray(2,runind) = ypri
         runarray(3,runind) = yd
         runarray(4,runind) = yhe3
         runarray(5,runind) = yli7

         runind = runind + 1
       end do !barind


       write (24,14) barprec
14     format ('it = ',i6,//,5x,'obh2',9x,'Yp',10x,'D/H',8x,'3He/H' &
               ,7x,'7Li/H',/,60('-'))
       if (barflag) then
         do barind=1,barprec
           write (24,15) (runarray(i,barind),i=1,5)
         end do
       end if
15     format (1p10e12.3)

!------close data files.

       if (bbnflag) close(unit=21, status = 'keep')
       if (equilflag) close(unit=22, status = 'keep')
       if (yeflag) close(unit=23, status = 'keep')
       if (barflag) close(unit=24, status = 'keep')

!------get time at end of computation and print duration.

       call cpu_time(time2) 

       write (*,50) time2 - time1, barprec

50     format('Duration (in s): ',1pe12.3,/,'Number of runs:',i6)

       stop


       end program main
