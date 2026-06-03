subroutine deallocation()

! *****************************************************************************
! 
!     purpose: 		deallocate fields
! 
!     log:           2015 / 03 - s.tschisgale
!                           05 - b.krull
!
! *****************************************************************************
    
  use Module_GlobalVariables
  use Module_ImmersedBoundary
  use Module_Fluid
  
  implicit none

  ! deallocation
  deallocate( u, a, IB )
  deallocate( work, wsave )
  deallocate( plotlineX, plotlineU )
  
end subroutine deallocation
