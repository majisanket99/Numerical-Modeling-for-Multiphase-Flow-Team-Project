module Module_ImmersedBoundary
  
! *****************************************************************************
! 
!     purpose: 		variables and subroutines for Immersed Boundaries
! 
!     log:           2015 / 03 - s.tschisgale
!                         / 04 - b.krull: forcing() adapted
!
! *****************************************************************************

  use module_GlobalVariables
  use module_Fluid

  implicit none
  save
  
  integer :: speed_mode = 2

    real :: control_gain_slow   = 100.0
    real :: control_gain_medium = 300.0
    real :: control_gain_fast   = 500.0

  integer :: Nb

  type ImmersedBoundary
    
    real :: C(2), D
    ! center position, diameter
    
    real, allocatable, &
      dimension(:,:) :: X, X0, XX, U, UD, F

    integer :: Nk
    
  end type ImmersedBoundary
  
  type(ImmersedBoundary), allocatable :: IB(:)
  

  contains
  
  
  function interpolation(u,X,Nk) result(UX)
  ! interpolation of velocities to forcing points 
  ! with two-point dirac hat function
    real    :: UX(Nk,2)
    real    :: u(N,N,2)
    real    :: X(Nk,2)
    integer :: Nk

    integer :: i, j, k,    &
               id, jd,     &
               ids, ide,   &
               jde, jds,   &
               indx, indy, &
               idso, jdso

    real    :: dh, dirx(2), diry(2), &
               x_indx_c, x_indx_p,   &
               y_indy_c, y_indy_p


    UX = 0.0

    do k=1,Nk
      
      ! _compute support range of dirac function
      indx = int( X(k,1)/dx ) + 1
      indy = int( X(k,2)/dx ) + 1
      
      ids = indx ; ide = indx + 1
      jds = indy ; jde = indy + 1
      
      ! _limit support to local index range
      idso = ids
      jdso = jds
      
      ids  = max(1, ids) ; ide = min(N, ide)
      jds  = max(1, jds) ; jde = min(N, jde)
      
      idso = ids - idso
      jdso = jds - jdso
      
      ! _get grid coordinates (in order to avoid memory corruption)
      x_indx_c = (indx - 1.0) * dx
      x_indx_p = (indx      ) * dx
      y_indy_c = (indy - 1.0) * dx
      y_indy_p = (indy      ) * dx
      
      ! _compute radius and dirac function
      dirx(1) = 1.0 - abs(x_indx_c - X(k,1))/dx
      dirx(2) = 1.0 - abs(x_indx_p - X(k,1))/dx
      diry(1) = 1.0 - abs(y_indy_c - X(k,2))/dx
      diry(2) = 1.0 - abs(y_indy_p - X(k,2))/dx
      
      ! _interpolation => U
      jd=0 
      do j=jds,jde
        
        id=0; jd=jd+1
        
        do i=ids,ide
          
          id=id+1
          
          dh = dirx(id + idso) * diry(jd + jdso)
          
          UX(k,1) = UX(k,1) + u(i,j,1) * dh
          
        end do
      end do
      
      ! _interpolation => V
      jd=0 
      do j=jds,jde
        
        id=0; jd=jd+1
        
        do i=ids,ide
          
          id=id+1
          
          dh = dirx(id + idso) * diry(jd + jdso)
          
          UX(k,2) = UX(k,2) + u(i,j,2) * dh
          
        end do
      end do
      
    end do !Nk

  end function interpolation
  

  function spreading(FX,X,Nk) result(f)
  ! spreading of forces on forcing points to force field
  ! with two-point dirac hat function
    real    :: f(N,N,2)
    real    :: FX(Nk,2), X(Nk,2), DV
    integer :: Nk

    integer :: i, j, k,    &
               id, jd,     &
               ids, ide,   &
               jde, jds,   &
               indx, indy, &
               idso, jdso

    real    :: dh, dirx(2), diry(2), &
               x_indx_c, x_indx_p,   &
               y_indy_c, y_indy_p,   &
               weight(Nk)

    f = 0.0
    
    ! _compute volume ratio Lagrangian / Eulerian volume = DV/(dx*dx)
      weight(1)  = 0.5/dx*( sqrt( (X(Nk  ,1)-X(1 ,1))**2 + (X(Nk  ,2)-X(1 ,2))**2) &
                        +   sqrt( (X(1   ,1)-X(2 ,1))**2 + (X(1   ,2)-X(2 ,2))**2) )
      
      
      ! _A) uniform distribution of Lagrangian points: all weights are the same
      weight     = weight(1)
      
      ! _B) non-uniform distribution of Lagrangian points
      ! ...

      
    do k=1,Nk
      
      ! _compute support range of dirac function
      indx = int( X(k,1)/dx ) + 1
      indy = int( X(k,2)/dx ) + 1
      
      ids = indx ; ide = indx + 1
      jds = indy ; jde = indy + 1
      
      ! _limit support to local index range
      idso = ids
      jdso = jds
      
      ids  = max(1, ids) ; ide = min(N, ide)
      jds  = max(1, jds) ; jde = min(N, jde)
      
      idso = ids - idso
      jdso = jds - jdso
      
      ! _get grid coordinates (in order to avoid memory corruption)
      x_indx_c = (indx - 1.0) * dx
      x_indx_p = (indx      ) * dx
      y_indy_c = (indy - 1.0) * dx
      y_indy_p = (indy      ) * dx
      
      ! _compute radius and dirac function
      dirx(1) = 1.0 - abs(x_indx_c - X(k,1))/dx
      dirx(2) = 1.0 - abs(x_indx_p - X(k,1))/dx
      diry(1) = 1.0 - abs(y_indy_c - X(k,2))/dx
      diry(2) = 1.0 - abs(y_indy_p - X(k,2))/dx
        
      ! _spreading => fx
      jd=0 
      do j=jds,jde
        
        id=0; jd=jd+1
        
        do i=ids,ide
          
          id=id+1
          
          dh = dirx(id + idso) * diry(jd + jdso)
          
          f(i,j,1) = f(i,j,1) + FX(k,1) * dh*weight(k)
          
        end do 
      end do
      
      ! _spreading => fy
      jd=0 
      do j=jds,jde
        
        id=0; jd=jd+1
          
        do i=ids,ide
          
          id=id+1
          
          dh = dirx(id + idso) * diry(jd + jdso)
          
          f(i,j,2) = f(i,j,2) + FX(k,2) * dh*weight(k)
          
        end do
      end do
      
    end do

  end function spreading


  function forcing(X,X0,U,UD,Nk) result(FX)
  ! Lagrangian forces F on forcing points of IB(b)
    real    :: FX(Nk,2)
    real    :: X(Nk,2), X0(Nk,2), U(Nk,2), UD(Nk,2)
    integer :: Nk
    integer :: speed_mode = 2

    real :: control_gain_slow   = 100.0
    real :: control_gain_medium = 300.0
    real :: control_gain_fast   = 500.0
    integer :: k
    logical :: slow_file, medium_file, fast_file
    real    :: C, DIF(Nk,2)
    
    real    :: Xsoll(2)      ! target position

    FX = 0.0
    
!     ! _elastic force
!     !  F = S*( curvature - curvature0 )
!     C = (2*pi/Nk)*(2*pi/Nk)
!     C = 2500.0/C
! 
!     DIF = X-X0 ! displacement
!     
!       FX(1 ,:) = FX(1 ,:) + C*(DIF(2  ,:)+DIF(Nk  ,:)-2*DIF(1 ,:))
!     do k=2,Nk-1
!       FX(k ,:) = FX(k ,:) + C*(DIF(k+1,:)+DIF(k-1 ,:)-2*DIF(k ,:))
!     end do
!       FX(Nk,:) = FX(Nk,:) + C*(DIF(1  ,:)+DIF(Nk-1,:)-2*DIF(Nk,:))
    
    ! artificial force
    Xsoll = getTargetPos()

    ! Check speed control flags
    inquire(file='./controls/slow.flag', exist=slow_file)
    inquire(file='./controls/medium.flag', exist=medium_file)
    inquire(file='./controls/fast.flag', exist=fast_file)

    if (slow_file) then
      speed_mode = 1
      call system("rm -f ./controls/slow.flag")
    end if

    if (medium_file) then
      speed_mode = 2
      call system("rm -f ./controls/medium.flag")
    end if

    if (fast_file) then
      speed_mode = 3
      call system("rm -f ./controls/fast.flag")
    end if

    ! Select control strength
    select case(speed_mode)

    case(1)
      C = control_gain_slow

    case(2)
      C = control_gain_medium

    case(3)
      C = control_gain_fast

    case default
      C = control_gain_medium

    end select
    ! C = 500.0
    FX(:,1) = FX(:,1) + C*sign( 1., Xsoll(1) - sum(X(:,1))/L/Nk )
    FX(:,2) = FX(:,2) + C*sign( 1., Xsoll(2) - sum(X(:,2))/L/Nk )

    ! _direct forcing: impose desired velocity UD
    C = rho/dt
    FX(:,1) = FX(:,1) + C*( UD(:,1)-U(:,1) )
    FX(:,2) = FX(:,2) + C*( UD(:,2)-U(:,2) )

  end function forcing
  
  function getTargetPos() result(X)
  ! get target position from input device
    real    :: X(2)
    real    :: cur(2)
    logical :: lo
    
    ! use xdotool to get mouse coordinates
    call system("eval $(xdotool getmouselocation --shell); echo $X $Y > ./targetPos/targetPos.dat")
        
    inquire(file='./targetPos/targetPos.dat', exist=lo)
    if(lo) then
      open (unit=5, file='./targetPos/targetPos.dat', action='read')
        read(5,*,end=250) cur(1), cur(2)
250   close(unit=5)

      ! relative position, mirror y axis
      ! get display resolution with $ xrandr | grep '*'
      cur(1) =     cur(1)/1920.
      cur(2) = 1 - cur(2)/1080.

      cur(1) = min(max(cur(1), 0.0), 1.0)
      cur(2) = min(max(cur(2), 0.0), 1.0)

      ! transform to physical coordinates
      X(1) =   cur(1)*L
      X(2) =   cur(2)*L   ! y-axis upwards      

    end if !lo  
  
  end function getTargetPos
  
end module module_ImmersedBoundary
