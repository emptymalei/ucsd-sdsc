
       module ye


!------modules.

       use mainvar


       implicit none


!------data.


       contains


       subroutine ye_det(tpl,yntot,yptot)


!------linkages.

!      called by - [subroutine] accum
!      calls     - none


!-----remarks.

!     Calculates the equilibrium and actual electron fractions.


!------throughput variables.

       real(dl), intent(in) :: tpl !plasma temperature in MeV
       real(dl), intent(in) :: yntot !total neutron abundance
       real(dl), intent(in) :: yptot !total proton abundance


!------data.

       real(dl) :: deltamnp = 1.29333217_dl !mass difference in MeV between proton and neutron


!------local variables.

       real(dl) ye_equil !equilibrium ye
       real(dl) ye_actual !actual ye


!------procedure.

       !equilibrium electron fraction:
       ye_equil = 0._dl

       !actual electron fraction:
       ye_actual = 0._dl

       !write to file:
       write (23,35) tpl,ye_equil,ye_actual

35     format (1p10e12.3)


       return


       end subroutine ye_det


       end module ye
