subroutine initialize

! *****************************************************************************
! 
!     purpose: 		initialization
! 
!     log:           2015 / 03 - s.tschisgale
!                           05 - b.krull: direct forcing, wall
!
! *****************************************************************************
  
  use Module_GlobalVariables
  use Module_ImmersedBoundary
  use Module_Fluid
  use Module_floodfill
  
  implicit none
  
#include "main.def"

  real    :: x, s1, s2, ss, dtheta, theta, d

  integer :: i, j, b, k, Nk, ct

  logical :: lo
  
  namelist /general/   n_steps,    &! number of time steps
                        t_end,     &! final time
                        dt,        &! time step
                        nt_out,    &! output-files every nt_bin timesteps
			rho,       &! density
			rho2,      &! density 2
                        nu,        &! kinematic viscosity
                        N           ! number of gridpoints in x & y
                        
  namelist /velocity_bc/  setbc_north, setbc_east, &
                          setbc_south, setbc_west, &
                          u_north,  v_north, &
                          u_south,  v_south, &
                          u_east,   v_east,  &
                          u_west,   v_west
  
  ! read parameters from input.dat --------------------------------------------
  open(unit=1, file='input.dat')
    read(1, general)
    read(1, velocity_bc)
  close(1)
  ! _estimate number of time steps
  if(n_steps.gt.0) &
  n_steps = min(n_steps,ceiling(t_end/dt))

  ! initialize floodfill field ------------------------------------------------
  allocate( bool(N,N) )
  

  ! initialize fluid ----------------------------------------------------------
  ! _derived variables
  L  = 1.0             ! fluid domain == square box with length L=1
  dx = L/(N-1)         ! grid spacing

  ! _initial condition of velocity field
  allocate( u(N,N,2) )

  u = 0.0
!   u(:,:,1) = 0.6
!   do j=0,(N-1)
!     x = j*dx
!     u(j+1,:,2) = sin(2*pi*x/L)
!   end do

  ! _initial condition of scalar field
  allocate( phi(N,N) )
  phi = 0.0
 
!   ! _ellipse
!   d = N/5. ! width
!   do i=1,N
!    do j=1,N
!     if ( ((real(i)-N/2.)/( N/12. ))**2. + ((real(j)-4./5.*N)/( N/12. ))**2. <= 1. ) phi(i,j) = 1.0
! !     if ( ((real(i)-N/2.)/( N/7. ))**2. + ((real(j)-0.65*N)/( N/7. ))**2. <= 1. ) phi(i,j) = 1.0 ! coalescence
! !     if ( ((real(i)-N/2.)/( N/7. ))**2. + ((real(j)-0.35*N)/( N/7. ))**2. <= 1. ) phi(i,j) = 1.0
!    enddo
!   enddo

  ! _square
!   phi(int(2./5.*N):int(3./5.*N),int(2./5.*N):int(3./5.*N)) = 1.
  
!   _bottom
!     phi(:,1:int(N/3.)) = 1.
    
! !   _top
!     phi(:,int(2.*N/3.):N) = 1.

   ! _horizontal block
!   phi(1:N,int(3./5.*N):int(4./5.*N)) = 1.
!   phi(11:N-10,int(3.5/5.*N):int(4.5/5.*N)) = 1.
  
!   
! !   phi(int(1./4.*N),int(3.5/5.*N))        = 0. ! disturb
!   phi(int(0.5*N  ),int(3.5/5.*N))        = 0. ! disturb
! !   phi(int(3./4.*N),int(3.5/5.*N))        = 0. ! disturb

!   ! _sloshing
!   phi(3:int(1./4.*N),3:int(0.9*N)) = 1.

  
  ! _array  a  that is used in the fluid solver
  allocate( a(N,N,2,2) )

  a = 0.0

  do i=0,(N-1)
    do j=0,(N-1)
      a(i+1,j+1,1,1) = cmplx(1.0,0.0)
      a(i+1,j+1,2,2) = cmplx(1.0,0.0)
    end do
  end do

  do i=0,(N-1)
    do j=0,(N-1)
      if(.not.(     (i.eq.0 .or. i.eq.N/2) &
               .and.(j.eq.0 .or. j.eq.N/2))) then
        
        s1 = sin(2*pi*i/N)
        s2 = sin(2*pi*j/N)
        ss = s1*s1 + s2*s2 + 1.E-12
        
        a(i+1,j+1,1,1) = a(i+1,j+1,1,1) - cmplx(s1*s1/ss,0.0)
        a(i+1,j+1,1,2) = a(i+1,j+1,1,2) - cmplx(s1*s2/ss,0.0)
        a(i+1,j+1,2,1) = a(i+1,j+1,2,1) - cmplx(s2*s1/ss,0.0)
        a(i+1,j+1,2,2) = a(i+1,j+1,2,2) - cmplx(s2*s2/ss,0.0)
      end if
    end do
  end do

  do i=0,(N-1)
    do j=0,(N-1)
      
      s1 = sin(2*pi*i/N)
      s2 = sin(2*pi*j/N)
      ss = s1*s1 + s2*s2
      
      a(i+1,j+1,:,:) = a(i+1,j+1,:,:) &
                     /(1+(dt/2)*nu*(4/(dx*dx))*ss)
    end do
  end do

  ! _fast fourier transformation
  lenwrk = 2*N*N
  lensav = 2*N + int(log(real(N,kind=4))/log(2.0E+00)) &
         + 2*N + int(log(real(N,kind=4))/log(2.0E+00)) &
         + 8

  allocate( work (1:lenwrk), &
            wsave(1:lensav)  )

  call cfft2i(N,N,wsave,lensav,ier) ! initializes the transform


  ! initialize plot line ------------------------------------------------------
  plotlineN = N-1                                  ! number of points
  allocate( plotlineX(plotlineN,2), plotlineU(plotlineN,2) )
  do k=1,plotlineN
      plotlineX(k,1) =    real(k-1)/real(plotlineN-1) * L  ! chose x position
      plotlineX(k,2) =    0.5 * L                          ! chose y position
  end do
  
  ! initialize immersed boundaries (IB)----------------------------------------
  ! _specify immersed boundaries
  Nb = 1                  ! number of immersed boundaries
  allocate(IB(Nb))

  ! ----Boundary: sphere----
  
  b = 1

    IB(b)%D  = L/6.            ! circle diameter
    IB(b)%C  = (/0.3*L,0.5*L/) ! circle center position  +dx/3.0

    Nk     = int(2*pi*IB(b)%D/dx)+1
    
    dtheta = 2*pi/Nk

    IB(b)%Nk = Nk
    
    allocate( IB(b)%X (Nk,2), &
              IB(b)%X0(Nk,2), &
              IB(b)%XX(Nk,2), &
              IB(b)%U (Nk,2), &
              IB(b)%UD(Nk,2), &
              IB(b)%F (Nk,2)  )

  ! _positions of forcing points
    do k=1,Nk
    
      theta = (k-1)*dtheta
      IB(b)%X(k,1) = IB(b)%C(1) + (IB(b)%D/2)*cos(theta)
      IB(b)%X(k,2) = IB(b)%C(2) + (IB(b)%D/2)*sin(theta)
      
    end do

    ! _set reference configuration
    IB(b)%X0 = IB(b)%X
  
    IB(b)%UD(:,1) = 0.0
    IB(b)%UD(:,2) = 0.0


    
!       ! ----Boundary: NACA profile ----
!   
!   b = 1
! 
!     IB(b)%D  = L/3.            ! circle diameter
!     IB(b)%C  = (/0.25*L,0.5*L/) ! circle center position  +dx/3.0
! 
!     Nk     = int(2*pi*IB(b)%D/dx)+1
!     
!     Nk     = (Nk+1)/2*2 ! Nk should be even
!     
!     dtheta = 2*pi/Nk
! 
!     IB(b)%Nk = Nk
!     
!     allocate( IB(b)%X (Nk,2), &
!               IB(b)%X0(Nk,2), &
!               IB(b)%XX(Nk,2), &
!               IB(b)%U (Nk,2), &
!               IB(b)%UD(Nk,2), &
!               IB(b)%F (Nk,2)  )
! 
!   ! _positions of forcing points
!     IB(b)%X(1,1)  = 0. !01*IB(b)%D !IB(b)%D/Nk
!     IB(b)%X(Nk,1) = IB(b)%X(1,1)
!     do k=2,Nk/2
!        IB(b)%X(k,1)      = IB(b)%X(k-1,1) + 2.*IB(b)%D/Nk
!        IB(b)%X(Nk-k+1,1) = IB(b)%X(k,1)
!     end do
!     IB(b)%X(:,2)         = 0.2*IB(b)%D * sqrt( IB(b)%X(:,1)/IB(b)%D )*(1.- IB(b)%X(:,1)/IB(b)%D )
!     IB(b)%X(Nk/2+1:Nk,2) = - IB(b)%X(Nk/2+1:Nk,2)
!     
!     IB(b)%X(:,1) = IB(b)%C(1) + IB(b)%X(:,1)
!     IB(b)%X(:,2) = IB(b)%C(2) + IB(b)%X(:,2)
! 
!     ! _set reference configuration
!     IB(b)%X0 = IB(b)%X
!   
!     IB(b)%UD(:,1) = 0.0
!     IB(b)%UD(:,2) = 0.0
    
  ! remove old result files ---------------------------------------------------
  inquire(file='./results/fluid.bin', exist=lo)
  if(lo) then
    open (unit=5, file='./results/fluid.bin', status='old')
    close(unit=5, status='delete')
  end if
  
  inquire(file='./results/ib.dat', exist=lo)
  if(lo) then
    open (unit=5, file='./results/ib.dat', status='old')
    close(unit=5, status='delete')
  end if

  inquire(file='./results/slice.dat', exist=lo)
  if(lo) then
    open (unit=5, file='./results/slice.dat', status='old')
    close(unit=5, status='delete')
  end if
  
  ! remove xstop --------------------------------------------------------------
  inquire(file='./xstop', exist=lo)
  if(lo) then
    open (unit=5, file='./xstop', status='old')
    close(unit=5, status='delete')
  end if

  ! remove cursor file --------------------------------------------------------
    inquire(file='./results/cursor.dat', exist=lo)
    if(lo) then
      open (unit=5, file='./results/cursor.dat', status='old')
      close(unit=5, status='delete')
    end if

  ! write info.dat ------------------------------------------------------------
  open(unit=5, file = './results/info.dat', action = 'write')
    write(5,*) N, 0 !sum(IB(:)%Nk)
  close(5)
  
  
  ! _immersed boundaries
  open(unit=12, file = './results/ib0.dat', access='direct', recl=25, form='formatted')
    ct=0
    do b=1,Nb
      do k=1,IB(b)%Nk
        ct=ct+1
        write(12,'(2E12.4,A)',rec=ct) Ib(b)%X(k,1)/dx, Ib(b)%X(k,2)/dx, char(10)
      end do
    end do
  close(12)
  
end subroutine initialize
