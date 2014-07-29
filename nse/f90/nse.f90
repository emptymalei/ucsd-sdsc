
       module nse


!------modules.

       use mainvar


       implicit none


!------data.

!    nuclide and corresponding number
!    --------------------------------
!    1) D
!    2) T
!    3) 3He
!    4) 4He

       real(dl), dimension(4) :: z = (/1._dl,1._dl,2._dl,2._dl/) !atomic numbers
       real(dl), dimension(4) :: n = (/1._dl,2._dl,1._dl,2._dl/) !neutron numbers
       real(dl), dimension(4) :: deltaq = (/2.2246_dl,8.4818_dl,7.7180_dl,28.296_dl/) !binding energies (in MeV)
       real(dl), dimension(4) :: m = (/1875.6_dl,2808.9_dl,2809.4_dl,3727.4_dl/) !atomic masses (in MeV)
       real(dl), dimension(4) :: spina =(/3._dl,2._dl,2._dl,1._dl/) !spin partition functions


       contains


       subroutine equil(tpl,yneu,ypro,rb)


!------linkages.

!      called by - [subroutine] accum
!      calls     - none




!------remarks.

!      Calculates NSE abundances.

       real(dl), intent(in) :: tpl !plasma temperature in MeV
       real(dl), intent(in) :: yneu !free neutron abundance
       real(dl), intent(in) :: ypro !free proton abundance
       real(dl), intent(in) :: rb !baryon energy density in cgs

!------local variables.

       integer j
       real(dl) nb !baryon number density in MeV
       real(dl), dimension(4) :: ya !nse abundances
       real(dl), dimension(4) :: equil_y !nse relative abundances

       real(dl) :: amu = 931.494061_dl !atomic mass unit in MeV
       real(dl) :: spinp = 2._dl !proton spin
       real(dl) :: spinn = 2._dl !neutron spin

!------procedure.

       !all units MeV.

       !baryon number density:
       nb = rb/amu/232011.568_dl       

       do j=1,4
         !put nse expression here:
         ya(j) = nb*yneu*ypro*(spinn*spinp/spina(j))*(2*pi*(z(j)+n(j))*m(j)/(tpl))**(-2/3)*exp((m(j)-z(j)*amu-n(j)*amu)/tpl)
       end do

       !Make relative abundances:
       do j=1,4
         equil_y(j) = ya(j)/ypro
         if (equil_y(j).gt.2._dl) then
           equil_y(j) = 2._dl
         end if
       end do

       !Write to file:
       write (22,25) tpl,(equil_y(j),j=1,4)

25     format (1p10e12.3)


       return


       end subroutine equil


       end module nse
