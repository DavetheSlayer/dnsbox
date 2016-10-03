module variables
    
    use parameters

    implicit none
    include 'mpif.h'

    !simulation variables:
    real(kind=8)            :: dt= 0.001d0             ! Time step size
    real(kind=8), dimension(:), allocatable     :: x, y, z, time, mychg, allchg
    real(kind=8)            :: scalemodes, chg, factor!, scalemodessquare
    real(kind=8)            :: eps                    ! epsilon to avoid divs
                                                      ! by 0
    
    real(kind=8), dimension(:, :, :), allocatable :: u, v, w, &
                                                     ux, uy, uz, &
                                                     vx, vy, vz, &
                                                     wx, wy, wz, &
                                                     omegax, &
                                                     omegay, &
                                                     omegaz, &
                                                     utemp, vtemp, wtemp,& 
                                                     temp_r 
    real(kind=8), dimension(:), allocatable       :: kx, ky, kz
    complex(kind=8), dimension(:, :, :), allocatable :: uhat, vhat, what,&
                                                        uhattemp, &
                                                        vhattemp, &
                                                        whattemp, &
                                                        uhatold, &
                                                        vhatold, &
                                                        whatold, &
                                                        nonlinuhat, &
                                                        nonlinvhat, &
                                                        nonlinwhat, &
                                                        nonlinuhatold, &
                                                        nonlinvhatold, &
                                                        nonlinwhatold, &
                                                        rhsuhatfix, &
                                                        rhsvhatfix, &
                                                        rhswhatfix, &
                                                        temp_c, &
                                                        intFact, &
                                                        phat
    
    ! extraneous variables, for easy code reading:
    real(kind = 8) :: Lx, Ly, Lz
        
    !counters and logicals:
    integer(kind=4) :: ind, iter, i, j, k, n, t, AllocateStatus
    logical         :: running_exist
      
    !measurement variables
    real(kind = 8) :: Ekin        !total kinetic energy
    real(kind = 8) :: Eband       !total kinetic energy of input band
    real(kind = 8) :: Courant     !Maximum Courant number for time-step control 
                                  !control 
    real(kind = 8) :: kk          !|k|^2
    real(kind = 8) :: phase       !phase for initiation
    real(kind=8), dimension(:), allocatable     :: kSpec      ! |k| array for  
    real(kind=8), dimension(:), allocatable     :: Espec      ! Energy 
                                                              ! spectrum
    real(kind = 8) :: absk        ! |k|
    integer        :: Nspec       ! length of the energy spectrum arrays
    integer        :: nk          ! position of the window corresponding to k
    real(kind = 8) :: Disp        ! Total dissipation
    real(kind = 8) :: divMax      ! Maximum divergence for error control
    real(kind = 8) :: EZero       ! Energy contained in k=0
    
    contains
    
    subroutine var_init()
        
        Lx = 2 * pi / alpha_x
        Ly = 2 * pi / alpha_y
        Lz = 2 * pi / alpha_z
        
        ! Cubic truncation:
        if (spherical) then
            ! Isotropic truncation:
            Nspec = int(((real(Nx, kind=8) / 2.0d0) & 
                        * (2.0d0 * alpha_x / 3.0d0)) / Deltak) - 1
        else
            Nspec = int(sqrt(((real(Nx, kind=8) / 2.0d0) &
                          * (2.0d0 * alpha_x / 3.0d0)) ** 2  + &
                          ((real(Ny, kind=8) / 2.0d0) &
                          * (2.0d0 * alpha_y / 3.0d0)) ** 2  + &
                          ((real(Nz, kind=8) / 2.0d0) &
                          * (2.0d0 * alpha_z / 3.0d0)) ** 2) / Deltak) - 1
        end if
        
         
        print *, 'Nspec = ', Nspec
        
        ! Allocate common arrays:
        allocate (x(1:Nx), y(1:Ny), z(1:Nz), kx(1:Nx), ky(1:Ny), kz(1:Nz), & 
                  time(1:Nt+1), mychg(1:3),allchg(1:3), &
                  kSpec(1:Nspec), myEspec(1:Nspec), Espec(1:Nspec), &
                  stat=AllocateStatus)
                  
        if (AllocateStatus .ne. 0) then
            print *, 'Error ', AllocateStatus, ' allocating common arrays'
            stop
        endif
        
        ! Allocate simulation arrays:
        allocate (u(1:Nx, 1:Ny, 1:Nz), &
                  v(1:Nx, 1:Ny, 1:Nz), &
                  w(1:Nx, 1:Ny, 1:Nz), &
                  ux(1:Nx, 1:Ny, 1:Nz), &
                  uy(1:Nx, 1:Ny, 1:Nz), &
                  uz(1:Nx, 1:Ny, 1:Nz), &
                  vx(1:Nx, 1:Ny, 1:Nz), &
                  vy(1:Nx, 1:Ny, 1:Nz), &
                  vz(1:Nx, 1:Ny, 1:Nz), &
                  wx(1:Nx, 1:Ny, 1:Nz), &
                  wy(1:Nx, 1:Ny, 1:Nz), &
                  wz(1:Nx, 1:Ny, 1:Nz), &
                  omegax(1:Nx, 1:Ny, 1:Nz), &
                  omegay(1:Nx, 1:Ny, 1:Nz), &
                  omegaz(1:Nx, 1:Ny, 1:Nz), &
                  utemp(1:Nx, 1:Ny, 1:Nz), &
                  vtemp(1:Nx, 1:Ny, 1:Nz), &
                  wtemp(1:Nx, 1:Ny, 1:Nz), &
                  temp_r(1:Nx, 1:Ny, 1:Nz), &
                  uhat(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), &
                  vhat(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), &
                  what(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), &
                  uhatold(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), &
                  vhatold(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), &
                  whatold(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), &
                  uhattemp(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), &
                  vhattemp(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), &
                  whattemp(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), &
                  rhsuhatfix(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), &
                  rhsvhatfix(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), &
                  rhswhatfix(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), &
                  nonlinuhat(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), &
                  nonlinvhat(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), &
                  nonlinwhat(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), &
                  nonlinuhatold(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), &
                  nonlinvhatold(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), &
                  nonlinwhatold(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), &
                  temp_c(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), &
                  intFact(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), &
                  phat(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), &
                  stat=AllocateStatus)
                    
        if(AllocateStatus .ne. 0) then
            print *, 'Error ', AllocateStatus, ' allocating simulation arrays'
            stop
        end if
        
        if (proc_id .eq. 0) then
            print *, 'Successfully allocated memory'
        end if
        
        ! Fourier frequencies in x-direction
        do i = 1, Nx/2 + 1
            kx(i) = cmplx(0.0d0, 1.0d0) * real(i - 1, kind(0d0)) * alpha_x
        end do
        kx(1 + Nx/2) = 0.0d0
        
        do i = 1, Nx/2 - 1
            kx(i + 1 + Nx/2) = -kx(1 - i + Nx/2)
        end do

        ind = 1
        do i = -Nx/2, Nx/2 - 1
            x(ind) = real(i, kind(0d0)) * Lx / real(Nx, kind(0d0))
            ind = ind + 1
        end do
        
        ! Fourier frequencies in y-direction
        do i = 1, Ny/2 + 1
            ky(i) = cmplx(0.0d0, 1.0d0) * real(i - 1, kind(0d0)) * alpha_y
        end do
        ky(1 + Ny/2) = 0.0d0
        do i = 1, Ny/2 - 1
            ky(i + 1 + Ny/2) = -ky(1 - i + Ny/2)
        end do
        ind = 1
        do i = -Ny/2, Ny/2 - 1
            y(ind) = real(i, kind(0d0)) * Ly / real(Ny, kind(0d0))
            ind = ind + 1
        end do
            
        ! Fourier frequencies in z-direction
        do i = 1, Nz/2 + 1
            kz(i) = cmplx(0.0d0, 1.0d0) * real(i - 1, kind(0d0)) * alpha_z
        end do
        kz(1 + Nz/2) = 0.0d0
        do i = 1, Nz/2 - 1
            kz(i + 1 + Nz/2) = -kz(1 - i + Nz/2)
        end do
        ind = 1
        do i = -Nz/2, Nz/2 - 1
            z(ind) = real(i, kind(0d0)) * Lz / real(Nz, kind(0d0))
            ind = ind + 1
        end do
        
        do i = 1, Nspec
            kSpec(i) = i * Deltak
        end do
        
        scalemodes = 1.0d0 / real(Nx * Ny * Nz, kind(0d0))
        !scalemodessquare = 1.0d0 / real(Nx * Ny * Nz, kind(0d0))
        
        if (proc_id .eq. 0) then
            print *, 'Setup grid and Fourier frequencies'
        end if    
        
    end subroutine var_init
    
end module variables
