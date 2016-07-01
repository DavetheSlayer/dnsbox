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
    integer(kind=4), parameter      :: Lx = 1 
    integer(kind=4), parameter      :: Ly = 1
    integer(kind=4), parameter      :: Lz = 1
    integer(kind=4), parameter      :: Nt = 1000          ! Number of t-steps
    integer(kind=4), parameter      :: iSaveRate1 = 100   ! save rate for state files
    integer(kind=4), parameter      :: iSaveRate2 = 1     ! save rate for analysis files
    real(kind=8), parameter         :: dt= 0.1d-3         ! Time step size
    real(kind=8), parameter         :: nuk = 4.5681d-4    ! Kinematic viscosity 
    real(kind=8), parameter         :: Re= 1.0d0 / nuk    ! 1/nu 
    real(kind=8), parameter         :: Q = 1.667d-1       ! Forcing multiplier
    real(kind=8), parameter         :: tol=0.1d0**10   ! tolerance for t-steps
    real(kind=8), parameter         :: theta=0.0d0     ! dunno

    real(kind=8), parameter &
    :: pi=3.14159265358979323846264338327950288419716939937510d0
    real(kind=8), parameter :: ReInv=1.0d0/real(Re, kind(0d0))
    real(kind=8), parameter :: dtInv=1.0d0/real(dt, kind(0d0))
    real(kind=8)            :: scalemodes, chg, factor


!***************************************************************************
 end module parameters
!***************************************************************************
