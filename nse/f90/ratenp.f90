
       subroutine ratenp


!------linkages.

!      called by - [subroutine] start, derivs  
!      calls     - none
                                       

!------modules.

       use bbnvar


       implicit none                                            


!------remarks.

!      generates rate coefficients for weak n->p and p->n reactions.
!      does not use thermal nu distributions, but rather ffact


!------local variables.

       integer, parameter :: nbins = 101 !total number of nu-epsilon bins
       integer i,n         !indicies
       integer	eind1,eind2
       real(dl), dimension(nbins) :: fermifact    !fermi factors for coulomb corrections
       real(dl)	coulcorr               !coulcorr = fermifact(i)
       real(dl), dimension(2,nbins) :: nuoccprob  !Electron neutrino occupation probabilities
       real(dl), dimension(nbins) :: epsmid
       real(dl), dimension(nbins) :: deltaeps
       real(dl), dimension(nbins-1) :: epsparam
       real(dl)	eps      
       real(dl)	snue,snuebar,sep,sem        !occupation fraction for nu's and e^+/-
       real(dl)	deltae,ee          !depsilon and epsilon for integrals
       real(dl)	enu	           !nu epsilon
       real(dl)	etadeg             !e^+/- degeneracy (etadeg = +/-phie)
       real(dl)	parfact             
       real(dl)	intstepfor,intsteprev,lamforint,lamrevint  !integral steps and sums
       real(dl)	lamfor1        !nue on n
       real(dl)	lamfor2        !n decay
       real(dl)	lamfor3        !e^+ on n
       real(dl)	lamrev1        !e^- on p
       real(dl)	lamrev2        !nuebar on p
       real(dl)	lamrev3        !inverse n decay
       real(dl) meps !xmelec/tnmev
       real(dl) mnpeps !deltamnp/tnmev


!------procedure.

!------epsilon values of deltamnp and xmelec.

       mnpeps = deltamnp/tcmev
       meps = xmelec/tcmev

       eps = mnpeps - meps
       if (eps.lt.100._dl) then
         eind1 = int(eps/100._dl*float(nbins-1)) + 1
       else
         eind1 = nbins
       end if

       eps = mnpeps + meps
       if (eps.lt.100._dl) then
         eind2 = int(eps/100._dl*float(nbins-1)) + 1
       else
         eind2 = nbins
       end if
              
!------Fermi correction factors.

       do i=1,nbins
         fermifact(i) = 0.98
       end do
       
!------Fermi correction factors.

       do i=1,nbins
         if (i.lt.nbins) then
           epsmid(i) = (2._dl*float(i) - 1._dl)/2._dl*100._dl/float(nbins-1)
           deltaeps(i) = 100._dl/float(nbins-1)
           epsparam(i) = float(i)
           parfact = 1._dl
           do n=1,2
             nuoccprob(n,i) = 1._dl/(exp(epsmid(i) - parfact*phie) + 1._dl)
             parfact = -1._dl
           end do
         else
           epsmid(i) = 100._dl
           deltaeps(i) = 100._dl
           parfact = 1._dl
           do n=1,2
             nuoccprob(n,i) = 1._dl/(exp(epsmid(i) - parfact*phie) + 1._dl)
             parfact = -1._dl
           end do
         end if
       end do
       

!------for n->p
       lamfor1 = 0.0
       lamfor2 = 0.0
       lamfor3 = 0.0

       lamforint = 0.0        !nue on n
       etadeg = phie
       do i=1,nbins
         enu = epsmid(i)
         ee = enu + mnpeps
         sem = 1.0/(exp(ee*tcmev/t9mev - etadeg) + 1.0)
         deltae = deltaeps(i)           
         snue = nuoccprob(1,i)
         coulcorr = fermifact(i)
         intstepfor = coulcorr*enu**2*ee*squarert(ee**2-meps**2,1.0e-06_dl) &
                      *snue*(1.0-sem)*deltae
         lamforint = lamforint + intstepfor
       end do
       lamfor1 = tcmev**5*cnorm*lamforint
       
       lamforint = 0.0        !n decay
       etadeg = phie
       do i=1,eind1
         if (i.eq.eind1) then
           if (i.ne.1) then
             enu = 0.5*(epsparam(i-1) + mnpeps - meps)
             deltae = mnpeps - meps - epsparam(i-1)
             if (deltae.eq.0.0) cycle
           else
             enu = 0.5*(mnpeps - meps)
             deltae = mnpeps - meps
           end if
         else
           enu = epsmid(i)
           deltae = deltaeps(i)           
         end if
         snuebar = nuoccprob(2,i)
         ee = mnpeps - enu
         sem = 1.0/(exp(ee*tcmev/t9mev - etadeg) + 1.0)
         coulcorr = fermifact(i)
         intstepfor = coulcorr*enu**2*ee*squarert(ee**2-meps**2,1.0e-06_dl) &
                      *(1.0-snuebar)*(1.0-sem)*deltae
         lamforint = lamforint + intstepfor
       end do
       lamfor2 = tcmev**5*cnorm*lamforint

       lamforint = 0.0        !e^+ on n
       etadeg = -phie
       do i=eind2,nbins
         if (i.eq.eind2) then
           if (eind2.ne.nbins) then
             enu = 0.5*(epsparam(i) + mnpeps + meps)
             deltae = epsparam(i) - mnpeps - meps
           else
             enu = mnpeps + meps
             deltae = mnpeps + meps - epsparam(nbins-1)
           end if
           if (deltae.eq.0.0) cycle
         else
           enu = epsmid(i)
           deltae = deltaeps(i)           
         end if
         snuebar = nuoccprob(2,i)
         ee = enu - mnpeps
         sep = 1.0/(exp(ee*tcmev/t9mev - etadeg) + 1.0)
         intstepfor = enu**2*ee*squarert(ee**2-meps**2,1.0e-06_dl) &
                      *(1.0-snuebar)*sep*deltae
         lamforint = lamforint + intstepfor
       end do
       lamfor3 = tcmev**5*cnorm*lamforint

!------for p->n
       lamrev1 = 0.0
       lamrev2 = 0.0
       lamrev3 = 0.0

       lamrevint = 0.0        !e^- on p
       etadeg = phie
       do i=1,nbins
         enu = epsmid(i)
         ee = enu + mnpeps
         sem = 1.0/(exp(ee*tcmev/t9mev - etadeg) + 1.0)
         deltae = deltaeps(i)           
         snue = nuoccprob(1,i)
         coulcorr = fermifact(i)
         intsteprev = coulcorr*enu**2*ee*squarert(ee**2-meps**2,1.0e-06_dl) &
                      *(1.0-snue)*sem*deltae
         lamrevint = lamrevint + intsteprev
       end do
       lamrev1 = tcmev**5*cnorm*lamrevint

       lamrevint = 0.0        !nuebar on p
       etadeg = -phie
       do i=eind2,nbins
         if (i.eq.eind2) then
           if (eind2.ne.nbins) then
             enu = 0.5*(epsparam(i) + mnpeps + meps)
             deltae = epsparam(i) - mnpeps - meps
           else
             enu = mnpeps + meps
             deltae = mnpeps + meps - epsparam(nbins-1)
           end if
           if (deltae.eq.0.0) cycle
         else
           enu = epsmid(i)
           deltae = deltaeps(i)           
         end if
         snuebar = nuoccprob(2,i)
         ee = enu - mnpeps
         sep = 1.0/(exp(ee*tcmev/t9mev - etadeg) + 1.0)
         intsteprev = enu**2*ee*squarert(ee**2-meps**2,1.0e-06_dl) &
                      *snuebar*(1.0-sep)*deltae
         lamrevint = lamrevint + intsteprev
       end do
       lamrev2 = tcmev**5*cnorm*lamrevint

       lamrevint = 0.0        !inverse n decay
       etadeg = phie
       do i=1,eind1
         if (i.eq.eind1) then
           if (i.ne.1) then
             enu = 0.5*(epsparam(i-1) + mnpeps - meps)
             deltae = mnpeps - meps - epsparam(i-1)
             if (deltae.eq.0.0) cycle
           else
             enu = 0.5*(mnpeps - meps)
             deltae = mnpeps - meps
           end if
         else
           enu = epsmid(i)
           deltae = deltaeps(i)           
         end if
         snuebar = nuoccprob(2,i)
         ee = mnpeps - enu
         sem = 1.0/(exp(ee*tcmev/t9mev - etadeg) + 1.0)
         coulcorr = fermifact(i)
         intsteprev = coulcorr*enu**2*ee*squarert(ee**2-meps**2,1.0e-06_dl) &
                      *snuebar*sem*deltae
         lamrevint = lamrevint + intsteprev
       end do
       lamrev3 = tcmev**5*cnorm*lamrevint
       
       f(1) = lamfor1 + lamfor2 + lamfor3        !no unit conversions necessary
       r(1) = lamrev1 + lamrev2 + lamrev3

       return


       end subroutine ratenp
