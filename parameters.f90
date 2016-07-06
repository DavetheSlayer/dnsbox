!***************************************************************************
! parameters
!***************************************************************************
!
! N		n in [1,N] radial r_n
! L		l in [0,L) theta
! M		m in [0,M) phi
! Mp		index m=0,1,2,.. corresponds to m=0,Mp,2Mp,...
!
!***************************************************************************
module parameters
!***************************************************************************

    implicit none
    !simulation parameters:
    integer(kind=4), parameter      :: Nx = 256
    integer(kind=4), parameter      :: Ny = 256
    integer(kind=4), parameter      :: Nz = 256
    integer(kind=4), parameter      :: Nt = 10000        ! Number of t-steps
    integer(kind=4), parameter      :: iSaveRate1 = 1000 ! save rate for state files
    integer(kind=4), parameter      :: iSaveRate2 = 10   ! save rate for analysis files
    real(kind=4), parameter         :: alpha_x = 1.0d0   ! 2 pi / Lx 
    real(kind=4), parameter         :: alpha_y = 1.0d0   ! 2 pi / Ly
    real(kind=4), parameter         :: alpha_z = 1.0d0   ! 2 pi / Lz
    real(kind=8), parameter         :: dt= 0.001d0        ! Time step size
    real(kind=8), parameter         :: nu = 4.491d-3     ! Kinematic viscosity 
    real(kind=8), parameter         :: Q = 0.0667d0      ! Forcing multiplier
    real(kind=8), parameter         :: tol=0.1d0**10     ! tolerance for t-steps

    real(kind=8), parameter &
    :: pi=3.14159265358979323846264338327950288419716939937510d0
    real(kind=8)            :: scalemodes, scalemodessquare, chg, factor

!***************************************************************************
 end module parameters
!***************************************************************************
