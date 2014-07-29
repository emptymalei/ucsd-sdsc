
       subroutine bbn(obh2,ypri,yd,yhe3,yli7)
 
!------linkages.

!      called by - [subroutine] main
!      calls     - [subroutine] driver


!-----remarks.

!     Computes light-element abundances using a baryon number and extra radiation energy density.


!------modules.

       use bbnvar


       implicit none


       save


!------throughput variables.

       real(dl), intent(in) :: obh2
       real(dl), intent(out) :: ypri !4He Mass Fraction
       real(dl), intent(out) :: yd !D relative abundance
       real(dl), intent(out) :: yhe3 !3He relative abundance
       real(dl), intent(out) :: yli7 !7Li relative abundance


!------local variables.

       integer i,j        !index
       integer inum                 !selection number.


!------procedure.


!20--------input initialization information and pause-----------------

       do i  = 1,nrec
!..........read in reaction parameters.
         iform(i) = int(reacpr(i,2))!reaction type.
         ii(i)    = int(reacpr(i,3))!incoming nuclide type.
         jj(i)    = int(reacpr(i,4))!incoming nuclide type.
         kk(i)    = int(reacpr(i,5))!outgoing nuclide type.
         ll(i)    = int(reacpr(i,6))!outgoing nuclide type.
         rev(i)   = reacpr(i,7)     !reverse reaction coefficient.
         q9(i)    = reacpr(i,8)     !energy released.
!..........initialize reaction rates.
         f(i) = 0._dl                 !forward rate coeff.
         r(i) = 0._dl                 !reverse rate coeff.
!..........set run options to default.
       end do
       irun       = 1               !do full run.
       isize      = nnuc
       jsize      = nrec
!..........set output option to default.
       nout    = 0                  !no output requests.

!------set values from bbn_params.ini.

       open (unit=20, file='bbn_params.ini', status='unknown')

       read (20,*) cy
       read (20,*) ct
       read (20,*) t9i
       read (20,*) t9f
       read (20,*) ytmin
       read (20,*) inc
       read (20,*) tau
       read (20,*) cosmo
       read (20,*) xi(1)
       read (20,*) xi(2)
       read (20,*) xi(3)
       read (20,*) dt1
       read (20,*) dtlow
       read (20,*) epstol

       close (unit=20)

!------set values from throughput variables.

       eta1 = 3._dl/8._dl/pi*pi**2/1.2020569031_dl/2._dl/(2.725_dl*8.6173324e-11_dl)**3 &
              /931.494061_dl*(1.221e+22_dl)**2/(2997.92458_dl)**2 &
              *(197.3269718_dl*3.24077929e-38_dl)**2*obh2
       rhors = 0._dl

       call driver        !do nucleosynthesis computation.

       
       do i = 1,it            !Temperature in MeV.
         t9out(i) = t9out(i)*.08617
       end do


       if (bbnflag) then

       write (21,2000) runind, obh2
2000   format ('run: ',i6,'; obh2 = ',1pe9.3)

       write (21,2002) cy,ct,t9i,t9f,ytmin
2002   format (' computational parameters:',/, &
               '   cy = ',f5.3,'/  ct = ',f5.3, &
               '/  initial temp = ',1pe8.2, &
               '/  final temp = ',1pe8.2, &
               '/  smallest abundances allowed = ',1pe8.2)

       write (21,2004) g,tau,3.0,cosmo,xi(1),xi(2),xi(3)
2004   format (' model parameters:',/, &
               '   g = ',f5.2,'/  tau = ',f6.2, &
               '/  # nu = ',f5.2,'/  lambda = ',1pe10.3, &
               '/  xi-e = ',e10.3,'/  xi-m = ',e10.3, &
               '/  xi-t = ',e10.3)

       write (21,2005) it
2005   format ('it = ',i6,/)

!..........print headings, abundances for neutron to li8.
       write (21,2006)
2006   format (4x,'temp',8x,'n/H',10x,'p',10x,'D/H',9x,'T/H',8x, &
                  '3He/H',8x,'4He',8x,'6Li/H',7x,'7Li/H',7x, &
                  '7Be/H',/,120('-'))
       do j = 1,it
         write (21,2008) t9out(j),(xout(j,i),i=1,9)
2008     format (1pe10.3,1p10e12.3)
       end do

!..........print thermodynamic quantities.
       write (21,2010)
2010   format (' ',/,4x,'temp',9x,'t',10x,'rhog',8x,'rhoe',7x, &
                        'rhone',8x,'rhob',8x,'phie',9x,'dt',9x, &
                        'sen',10x,'H',9x,'eta',/,132('-'))
       do j = 1,it
         write (21,2012) t9out(j),tout(j),(thmout(j,i),i=1,4), &
                         thmout(j,5),dtout(j),senout(j),hubout(j), &
                         thmout(j,8)
2012     format (1pe10.3,1p10e12.3)
       end do

       write (21,2014)
2014   format (///)

       end if


       ypri = y(6)*4._dl
       yd = y(3)/y(2)
       yhe3 = (y(4) + y(5))/y(2) !use A = 3 isobar
       yli7 = (y(9) + y(8))/y(2) !use A = 7 isobar


       return


       end subroutine bbn
