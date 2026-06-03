program IBM2D

! *****************************************************************************
! 
!     purpose: 		2D Immersed Boundary Method 
!
!                    based on Peskin's matlab solver
!                    see www.math.nyu.edu/faculty/peskin/ib_lecture_notes/
! 
!     log:           2015 / 03 - b. krull / s.tschisgale
!
! *****************************************************************************

  use Module_GlobalVariables
  use Module_ImmersedBoundary
  use Module_Fluid

  implicit none

  write(*,*)
  write(*,*)'-----------------------------------------------------------------'
  write(*,*)'           2D Immersed Boundary Method for education             '
  write(*,*)'-----------------------------------------------------------------'
  write(*,*)'-----------------------------------------------------------------'


  ! initialization --------------------
  write(*,*)
  write(*,*) ' initialization ... '
  write(*,*)
  call initialize
  t = 0.d0


  ! solving... ------------------------
  write(*,*)
  write(*,*) ' solving ... '
  write(*,*)
  
  nt=1
  do while(.not.exit_time_loop)
    
    ! update time --------------
    t = t+dt

    ! time integration ---------
    call solve

    ! write plot files ---------
    if(modulo(nt,nt_out).eq.0) call output

    ! exit time loop ? ---------
    inquire(file='xstop', exist=exit_time_loop)
    if(exit_time_loop) write(*,*) 'xstop found. Stop.'
    if(nt.eq.n_steps)  exit_time_loop=.true.

    nt=nt+1
    
  end do
  
  ! deallocation ----------
  call deallocation
  
  write(*,*)
  write(*,*)'-----------------------------------------------------------------'
  write(*,*)' done ... '
  write(*,*)'-----------------------------------------------------------------'
  write(*,*)


end program IBM2D
