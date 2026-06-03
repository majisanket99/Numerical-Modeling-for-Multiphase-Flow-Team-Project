module Module_GlobalVariables

! *****************************************************************************
! 
!     purpose: 		global variables
! 
!     log:           2015 / 03 - s.tschisgale
!                           05 - b.krull
!
! *****************************************************************************

  implicit none
  
  
! *** time integration scheme *****************************
  real              :: t, dt, t_end
  integer           :: n_steps, nt, nt_out

! *** further variables ***********************************
  logical           :: exit_time_loop = .false.
  
! *** parameters ******************************************
  real, parameter   :: pi = 3.14159265359, small = 1e-12
  
 
! *** line for line plot **********************************
  real, allocatable, &
      dimension(:,:) :: plotlineX, plotlineU
  
  integer :: plotlineN
 
contains 
  
  
  subroutine dummy
    
  end subroutine dummy
  
  
end module Module_GlobalVariables
