!***************************************************************************
! parameters
!***************************************************************************
!
!
!***************************************************************************
module parameters
!***************************************************************************

    implicit none
    !simulation parameters:
    integer(kind=4), parameter      :: Nx = 512
    integer(kind=4), parameter      :: Ny = 512
    integer(kind=4), parameter      :: Nz = 512
    integer(kind=4), parameter      :: Nt = 1000000      ! Number of t-steps
    integer(kind=4), parameter      :: iSaveRate1 = 1000 ! save rate for state files
    integer(kind=4), parameter      :: iSaveRate2 = 100  ! save rate for analysis files
    real(kind=4), parameter         :: alpha_x = 1.0d0   ! 2 pi / Lx 
    real(kind=4), parameter         :: alpha_y = 1.0d0   ! 2 pi / Ly
    real(kind=4), parameter         :: alpha_z = 1.0d0   ! 2 pi / Lz
    real(kind=8), parameter         :: dt= 0.00005d0      ! Time step size
    real(kind=8), parameter         :: nu = 3.0d-3       ! Kinematic viscosity 
    real(kind=8), parameter         :: Q = 1.0d0      ! Forcing multiplier
    real(kind=8), parameter         :: tol=0.1d0**10     ! tolerance for t-steps
    real(kind=8), parameter         :: Deltak=1.0d0      ! k-window for shell
                                                         ! averaging
    real(kind=8), parameter         :: Courant= 0.2d0    ! Courant number limit
    logical                         :: initrand = .false.! if true, random 
                                                         ! initial condition    
    real(kind=8), parameter &
    :: pi=3.14159265358979323846264338327950288419716939937510d0
    
!***************************************************************************
 end module parameters
!***************************************************************************
