module module_floodfill
  
! *****************************************************************************
! 
!     purpose: 		variables and subroutines for fluid
! 
!     log:           2016 / 06 - s.tschisgale/b.krull 
!
! *****************************************************************************

  use module_GlobalVariables
  use module_Fluid

  implicit none
  save
  
  logical, allocatable :: bool(:,:)     ! logical field
  
  contains

  subroutine setIBborder(X,Nk)
  ! set border
    implicit none
    real    :: X(Nk,2)
    integer :: Nk
    
    integer :: k, indx, indy
    
    do k=1,Nk ! compute border
      indx = int( X(k,1)/dx ) + 1
      indy = int( X(k,2)/dx ) + 1
      bool(indx,indy) = .true.
    enddo
    
  end subroutine setIBborder
  
  recursive subroutine floodFill(i,j)
  ! not safe (stack overfow) but very fast
    implicit none
    integer :: i,j
    
    if(.NOT.bool(i,j)) then ! false / 0.
        bool(i,j) = .true.  ! true / 1.
        call floodFill(i+1,j)
        call floodFill(i-1,j)
        call floodFill(i,j+1)
        call floodFill(i,j-1)
    end if
    
  end subroutine floodFill
  
  subroutine setBoolTrueWithinIB(X,Nk)
  ! fill inner range of immersed object --> z.B. "where (bool) alpha=0."
    implicit none
    real    :: X(Nk,2)
    integer :: Nk
    integer :: indx, indy
    
    call setIBborder(X,Nk)
    indx = int( sum(X(:,1))/Nk/dx ) + 1
    indy = int( sum(X(:,2))/Nk/dx ) + 1

    call floodfill(indx,indy)
    
  end subroutine setBoolTrueWithinIB
  
end module module_floodfill