subroutine solve()

! *****************************************************************************
! 
!     purpose: 	  IBM solver
! 
!     log:          2015 / 03 - s.tschisgale
!                        / 05 - b.krull
!
! *****************************************************************************
  
  use Module_GlobalVariables
  use Module_ImmersedBoundary
  use Module_Fluid
  use Module_floodfill
  
  implicit none

#include "main.def"

  real    :: uu(N,N,2), &  ! preliminary velocity
             ff(N,N,2)     ! force field acting on fluid
             
  integer :: b

ff = 0.

! _force phi on inlet and outlet
phi(1,:) = 0
phi(1:15,ceiling(0.40*N):ceiling(0.60*N)) = 1.

! _transport phi
phi = phi - dt*( uddx(phi,u) + vddy(phi,u) )

! clip to range
phi = min(max(0.,phi),1.)


! calculate immersed boundary forces

  do b=1,Nb ! for each immersed boundary object
    
    ! _get velocity at Lagrangian points from Eulerian velocity field
    IB(b)%U  = interpolation(u, IB(b)%X, IB(b)%Nk)

       ! _prediction step (Runge-Kutta method)
       IB(b)%XX = IB(b)%X + dt/2 * IB(b)%U
       
    ! _calculate Lagrangian force
    IB(b)%F  = forcing(IB(b)%XX, IB(b)%X0, IB(b)%U, IB(b)%UD, IB(b)%Nk)
    
    ! _get Eulerian representation of Lagrangian forces
    ff = ff + spreading(IB(b)%F, IB(b)%XX, IB(b)%Nk)
    
  end do
  
! calculate new velocity field u with Eulerian force field ff
  call fluid(u, uu, ff)

! move immersed boundaries
  do b=1,1
    
    IB(b)%U = interpolation(u, IB(b)%XX, IB(b)%Nk)

    ! transport Lagrangian points with mean velocity
    IB(b)%X(:,1) = IB(b)%X(:,1) + dt*sum( IB(b)%U(:,1) )/IB(b)%Nk
    IB(b)%X(:,2) = IB(b)%X(:,2) + dt*sum( IB(b)%U(:,2) )/IB(b)%Nk
    
    ! keep object in domain
    IB(b)%X(:,1) = IB(b)%X(:,1) - max( 0., maxval(IB(b)%X(:,1)) - L ) - min( 0., minval(IB(b)%X(:,1)) )
    IB(b)%X(:,2) = IB(b)%X(:,2) - max( 0., maxval(IB(b)%X(:,2)) - L ) - min( 0., minval(IB(b)%X(:,2)) )
 
  end do
  
! set phi=0 within the object
 bool = .false.
 call setBoolTrueWithinIB(IB(1)%X,IB(1)%Nk)
 where (bool) phi=0.

end subroutine solve
