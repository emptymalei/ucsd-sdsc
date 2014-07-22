!-*-f90-*-
program GR1D	
    	  
  use GR1D_module
  implicit none

  !Welcome to GR1D
  write(*,*) "#################################################"
  write(*,*) "#################################################"
  write(*,*) "########### GR1D SPHERICAL HYDRO v1.03 ##########"
  write(*,*) "##################Jan 8th, 2011##################"
  write(*,*) "#################################################"

  ! Call problem setup and allocate/initialize variables 
  call start
  write(*,*) "Done with initial data :-)"

  ! Call driver routine to do evolution
  call driver
    
  write(*,*) "Shutting down!"
  write(*,*) " "

end program GR1D
