module module_Fluid
  
! *****************************************************************************
! 
!     purpose: 		variables and subroutines for fluid
! 
!     log:           2015 / 03 - s.tschisgale
!                           05 - b.krull: 
!
! *****************************************************************************

  use module_GlobalVariables

  implicit none
  save
  
  ! material parameters
  real              :: rho, rho2, nu
  
  ! velocity boundary conditions
  logical           :: setbc_north, setbc_east, &
                       setbc_south, setbc_west
  real              :: u_north,  v_north, &
                       u_south,  v_south, &
                       u_east,   v_east,  &
                       u_west,   v_west 
  
  ! domain properties
  real              :: L, dx
  integer           :: N

  ! velocity field
  real, allocatable :: u(:,:,:)
  
  ! scalar field
  real, allocatable :: phi(:,:)

  ! fast fourier transformation
  real, allocatable :: work(:)
  real, allocatable :: wsave(:)

  integer*4         :: ier
  integer*4         :: lensav
  integer*4         :: lenwrk

  ! further variables & fields
  complex(kind=4), &
    allocatable     :: a(:,:,:,:)


  contains


  subroutine fluid(u,uu,f)
  ! solves NSE for periodic domain
    real            :: u (N,N,2), &
                       uu(N,N,2), &
                       f (N,N,2)

    real            :: w(N,N,2)
    complex(kind=4) :: c(N,N,2)
    integer :: k
    
    w = u - dt/2. * skew(u) + dt/2. * f

    c = cmplx(w,0.0)

    call cfft2f(N,N,N,c(:,:,1),wsave,lensav,work,lenwrk,ier)
    call cfft2f(N,N,N,c(:,:,2),wsave,lensav,work,lenwrk,ier)

    c(:,:,1) = a(:,:,1,1) * c(:,:,1) + a(:,:,1,2) * c(:,:,2)
    c(:,:,2) = a(:,:,2,1) * c(:,:,1) + a(:,:,2,2) * c(:,:,2)

    call cfft2b(N,N,N,c(:,:,1),wsave,lensav,work,lenwrk,ier)
    call cfft2b(N,N,N,c(:,:,2),wsave,lensav,work,lenwrk,ier)
    uu(:,:,1) = real(c(:,:,1))
    uu(:,:,2) = real(c(:,:,2))

    w = u - dt * skew(uu) + dt * f + dt/2. * nu * laplacian(u)

    c = cmplx(w,0.0)
    call cfft2f(N,N,N,c(:,:,1),wsave,lensav,work,lenwrk,ier)
    call cfft2f(N,N,N,c(:,:,2),wsave,lensav,work,lenwrk,ier)

    c(:,:,1) = a(:,:,1,1) * c(:,:,1) + a(:,:,1,2) * c(:,:,2)
    c(:,:,2) = a(:,:,2,1) * c(:,:,1) + a(:,:,2,2) * c(:,:,2)

    call cfft2b(N,N,N,c(:,:,1),wsave,lensav,work,lenwrk,ier)
    call cfft2b(N,N,N,c(:,:,2),wsave,lensav,work,lenwrk,ier)
    u(:,:,1) = real(c(:,:,1))
    u(:,:,2) = real(c(:,:,2))

    ! _set boundary conditions
    if (setbc_north) then
     u(:,N-1:N,1) = u_north
     u(:,N-1:N,2) = v_north
    endif
    if (setbc_south) then
     u(:,1:2,1)   = u_south
     u(:,1:2,2)   = v_south
    endif
    if (setbc_east) then
     u(N-1:N,:,1) = u_east
     u(N-1:N,:,2) = v_east
    endif
    if (setbc_west) then
      u(1:2,:,1)   = u_west
      u(1:2,:,2)   = v_west
    endif     
    
  end subroutine fluid


  function skew(u) result(w)
  ! skew of velocity field u
    real    :: w(N,N,2)
    real    :: u(N,N,2)

    real    :: C
    integer :: i , j,  &
               ip, im, &
               jp, jm
               
    C = 0.25/dx

    do j=1,N
      
      jp = j+1; if(jp.gt.N) jp = 1 ! for periodicity
      jm = j-1; if(jm.lt.1) jm = N
      
      do i=1,N
        
        ip = i+1; if(ip.gt.N) ip = 1 ! for periodicity
        im = i-1; if(im.lt.1) im = N
        
        w(i,j,:) = C*( (u(ip,j ,1)+u(i,j,1))*u(ip,j ,:) &
                      -(u(im,j ,1)+u(i,j,1))*u(im,j ,:) &
                      +(u(i ,jp,2)+u(i,j,2))*u(i ,jp,:) &
                      -(u(i ,jm,2)+u(i,j,2))*u(i ,jm,:) )
        
      end do
    end do

  end function skew 

  function laplacian(u) result(w)
  ! calculates laplacian of velocity field u
    real    :: w(N,N,2)
    real    :: u(N,N,2)

    real    :: C
    integer :: i , j,  &
               ip, im, &
               jp, jm


    C = 1.0/(dx*dx)

    do j=1,N
      
      jp = j+1; if(jp.gt.N) jp = 1 ! for periodicity
      jm = j-1; if(jm.lt.1) jm = N
      
      do i=1,N
        
        ip = i+1; if(ip.gt.N) ip = 1 ! for periodicity
        im = i-1; if(im.lt.1) im = N
        
        w(i,j,:) = C*(    u(ip,j ,:) &
                      +   u(im,j ,:) &
                      +   u(i ,jp,:) &
                      +   u(i ,jm,:) &
                      - 4*u(i ,j ,:) )
        
      end do
    end do

  end function laplacian
  

! edit krull: transport scalar

  function ddxc(s) result(w) ! central diff.
    real    :: w(N,N)
    real    :: s(N,N)
   
    real    :: C
    
    C = 1.0/dx
    w = 0.0

    w(2:N-1,:) = 0.5*C*( s(3:N,:) - s(1:N-2,:) )
    w(1,:)     = 0.5*C*( s(2  ,:) - s(N    ,:) ) ! periodicity
    w(N,:)     = 0.5*C*( s(1  ,:) - s(N-1  ,:) ) ! periodicity
    
  end function ddxc

  function ddxp(s) result(w) ! forward
    real    :: w(N,N)
    real    :: s(N,N)
   
    real    :: C
    
    C = 1.0/dx
    w = 0.0

    w(1:N-1,:) = C*( s(2:N,:) - s(1:N-1,:) )
    w(N,:)     = C*( s(1  ,:) - s(N    ,:) ) ! periodicity
    
!     ! 2nd order
!     w(1:N-2,:) = 0.5*C*( - s(3:N,:) + 4.* s(2:N-1,:) - 3.*s(1:N-2,:) )
!     w(N-1  ,:) = 0.5*C*( - s(1  ,:) + 4.* s(N    ,:) - 3.*s(N-1  ,:) ) ! periodicity
!     w(N    ,:) = 0.5*C*( - s(2  ,:) + 4.* s(1    ,:) - 3.*s(N    ,:) ) ! periodicity
    
  end function ddxp
  
  function ddxm(s) result(w) ! forward
    real    :: w(N,N)
    real    :: s(N,N)
   
    real    :: C
    
    C = 1.0/dx
    w = 0.0

    w(2:N,:) = C*( s(2:N,:) - s(1:N-1,:) )
    w(1,:)   = C*( s(1  ,:) - s(N    ,:) ) ! periodicity
    
!     ! 2nd order
!     w(3:N,:) = 0.5*C*( 3.*s(3:N,:) - 4.* s(2:N-1,:) + 1.*s(1:N-2,:) )
!     w(2  ,:) = 0.5*C*( 3.*s(2  ,:) - 4.* s(1    ,:) + 1.*s(N    ,:) ) ! periodicity
!     w(1  ,:) = 0.5*C*( 3.*s(1  ,:) - 4.* s(N    ,:) + 1.*s(N-1  ,:) ) ! periodicity
    
  end function ddxm
  
  function uddx(s,u) result(w)
  
    real    :: w(N,N)
    real    :: s(N,N)
    real    :: u(N,N,2)
  
    w = max(u(:,:,1),0.)*ddxm(s) + min(u(:,:,1),0.)*ddxp(s)
  
  end function uddx
  
  function ddxu(s,u) result(w)
  
    real    :: w(N,N)
    real    :: s(N,N)
    real    :: u(N,N,2)
  
    w = max(sign(1.,u(:,:,1)),0.)*ddxm(s) - min(sign(1.,u(:,:,1)),0.)*ddxp(s)
  
  end function ddxu
  
  function ddyc(s) result(w)
    real    :: w(N,N)
    real    :: s(N,N)
   
    real    :: C
    
    C = 1.0/dx
    w = 0.0
    
    w(:,2:N-1) = 0.5*C*( s(:,3:N) - s(:,1:N-2) )
    w(:,1)     = 0.5*C*( s(:,2  ) - s(:,N    ) ) ! periodicity
    w(:,N)     = 0.5*C*( s(:,1  ) - s(:,N-1  ) ) ! periodicity
  
  end function ddyc
  
  function ddyp(s) result(w) ! forward
    real    :: w(N,N)
    real    :: s(N,N)
   
    real    :: C
    
    C = 1.0/dx
    w = 0.0

    w(:,1:N-1) = C*( s(:,2:N) - s(:,1:N-1) )
    w(:,N)     = C*( s(:,1  ) - s(:,N    ) ) ! periodicity
    
!     ! 2nd order
!     w(:,1:N-2) = 0.5*C*( - s(:,3:N) + 4.* s(:,2:N-1) - 3.*s(:,1:N-2) )
!     w(:,N-1  ) = 0.5*C*( - s(:,1  ) + 4.* s(:,N    ) - 3.*s(:,N-1  ) ) ! periodicity
!     w(:,N    ) = 0.5*C*( - s(:,2  ) + 4.* s(:,1    ) - 3.*s(:,N    ) ) ! periodicity
    
  end function ddyp
  
  function ddym(s) result(w) ! forward
    real    :: w(N,N)
    real    :: s(N,N)
   
    real    :: C
    
    C = 1.0/dx
    w = 0.0

    w(:,2:N) = C*( s(:,2:N) - s(:,1:N-1) )
    w(:,1)   = C*( s(:,1  ) - s(:,N    ) ) ! periodicity

!     ! 2nd order
!     w(:,3:N) = 0.5*C*( 3.*s(:,3:N) - 4.* s(:,2:N-1) + 1.*s(:,1:N-2) )
!     w(:,2  ) = 0.5*C*( 3.*s(:,2  ) - 4.* s(:,1    ) + 1.*s(:,N    ) ) ! periodicity
!     w(:,1  ) = 0.5*C*( 3.*s(:,1  ) - 4.* s(:,N    ) + 1.*s(:,N-1  ) ) ! periodicity
    
  end function ddym

  function vddy(s,u) result(w)
  
    real    :: w(N,N)
    real    :: s(N,N)
    real    :: u(N,N,2)
  
    w = max(u(:,:,2),0.)*ddym(s) + min(u(:,:,2),0.)*ddyp(s)
    
  end function vddy
  
  function ddyu(s,u) result(w)
  
    real    :: w(N,N)
    real    :: s(N,N)
    real    :: u(N,N,2)
  
    w = max(sign(1.,u(:,:,2)),0.)*ddym(s) - min(sign(1.,u(:,:,2)),0.)*ddyp(s)
  
  end function ddyu
  
  function smooth(s) result(w) ! test, only inner domain
  
    real :: s(N,N)
    real :: w(N,N)
  
    w(2:N-1,2:N-1) = 0.5*s(2:N-1,2:N-1) + 0.125*( s(2:N-1,1:N-2) + s(2:N-1,3:N) + s(1:N-2,2:N-1) + s(3:N,2:N-1) )
  
  end function smooth
  
  function normal(s,u) result(no)
  
     real :: s(N,N)
     real :: no(N,N,2)
     real :: u(N,N,2)
     
!      no(:,:,1) = ddxu(s,u)
!      no(:,:,2) = ddyu(s,u)
     no(:,:,1) = ddxc(s)
     no(:,:,2) = ddyc(s)
     
     no(:,:,1) = no(:,:,1)/( sqrt( no(:,:,1)**2 + no(:,:,2)**2 ) + small )
     no(:,:,2) = no(:,:,2)/( sqrt( no(:,:,1)**2 + no(:,:,2)**2 ) + small )
     
  end function normal
  
  function div(s,u) result(w)
  
    real :: s(N,N,2)
    real :: u(N,N,2)
    real :: w(N,N)
  
    w = ddxc(s(:,:,1)) + ddyc(s(:,:,2)) ! CDS
!     w = ddxu(s(:,:,1),u) + ddyu(s(:,:,2),u) ! UDS2
  
  end function div
 
! end edit krull
  
end module Module_Fluid
