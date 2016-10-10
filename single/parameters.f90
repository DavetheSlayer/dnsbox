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
    integer(kind=4), parameter      :: Nx = 32           ! 256
    integer(kind=4), parameter      :: Ny = 32           ! 256
    integer(kind=4), parameter      :: Nz = 32           ! 256
    integer(kind=4), parameter      :: Nh = Nx/2 + 1     ! Real fft, x dimension
    integer(kind=4), parameter      :: Nt = 100000        ! Number of t-steps
    integer(kind=4), parameter      :: iSaveRate1 = 1000 ! save rate for state files
    integer(kind=4), parameter      :: iSaveRate2 = 10  ! save rate for analysis files
    real(kind=4), parameter         :: alpha_x = 1.0d0   ! 2 pi / Lx 
    real(kind=4), parameter         :: alpha_y = 1.0d0   ! 2 pi / Ly
    real(kind=4), parameter         :: alpha_z = 1.0d0   ! 2 pi / Lz
    real(kind=8), parameter         :: nu = 0.055d0       ! Kinematic viscosity 
    real(kind=8), parameter         :: Q = 0.0d0         ! Forcing multiplier ! set to 0 for band-lim
    real(kind=8), parameter         :: kCutOff = 2.5d0   ! 2.5     ! Forcing cut-off frequency
    real(kind=8), parameter         :: Pin = 0.1d0       ! Power input-rate
    
    real(kind=8), parameter         :: Deltak = 1.0d0    ! k-window for shell
                                                         ! averaging
    real(kind=8), parameter         :: kzero = 5.0d0     ! Center frequency for 
                                                         ! initial field generation
    real(kind=8), parameter         :: uzero = 1.0d0     ! rms velocity for 
                                                         ! initial field generation
    real(kind=8), parameter         :: CourantMin= 0.15d0! min Courant number
    real(kind=8), parameter         :: CourantMax= 0.2d0 ! max Courant number
    real(kind=8), parameter         :: tStepMax = 1.0d-2 ! max time step 
    logical, parameter              :: initrand = .true. ! if true, random 
                                                         ! initial condition    
    logical, parameter              :: tStepFix = .false.! fixed time step
    logical, parameter              :: analytic = .false. ! analytic solution 
                                                         ! for checking integrator
    real(kind = 8), parameter       :: c = 0.5d0         ! implicitness param.
    logical, parameter              :: spherical = .false.! spherical trunc. 
                                                         ! cubic trunc. if false
    logical, parameter              :: bandlim = .true.  ! if true, band-limited
                                                         ! forcing applied
    
    real(kind=8), parameter &
    :: pi=3.14159265358979323846264338327950288419716939937510d0
    complex(kind=8), parameter      :: ii = cmplx(0.0d0, 1.0d0)
    
!***************************************************************************
 end module parameters
!***************************************************************************
