module variables
    
    use parameters
    use p3dfft
    implicit none
    include 'mpif.h'

    !simulation variables:
    real(kind=8)            :: dt= 0.001d0             ! Time step size
    real(kind=8), dimension(:), allocatable     :: x, y, z, time, mychg, allchg
    real(kind=8)            :: scalemodes, chg, factor!, scalemodessquare
    real(kind=8)            :: eps                    ! epsilon to avoid divs
                                                      ! by 0
    
    real(p3dfft_type), dimension(:, :, :), allocatable :: u, v, w, &
                                                          ux, uy, uz, &
                                                          vx, vy, vz, &
                                                          wx, wy, wz, &
                                                          omegax, &
                                                          omegay, &
                                                          omegaz, &
!                                                          uold, vold, wold, &
                                                          utemp, vtemp, wtemp,& 
                                                          temp_r 
    complex(kind=8), dimension(:), allocatable      :: kx, ky, kz, linearterm
    complex(p3dfft_type), dimension(:, :, :), allocatable :: uhat, vhat, what,&
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
  
    !parallelization variables:
    integer ierr, nup, ndim, dims(2), nproc, proc_id ! MPI vars
    integer istart(3), iend(3), isize(3)             ! conf. space indices
    integer fstart(3), fend(3), fsize(3)             ! Fourier space indices
    integer iproc, jproc, nxc, nyc, nzc              ! double-paralellization
    logical iex
	integer memsize(3)
    integer p3_error
    
    !measurement variables
    real(kind = 8) :: Ekin        !total kinetic energy
    real(kind = 8) :: myEkin      !kinetic energy on a cpu
    real(kind = 8) :: Eband       !total kinetic energy of input band
    real(kind = 8) :: myEband     !kinetic energy of input band on a cpu 
    real(kind = 8) :: Courant     !Maximum Courant number for time-step control 
    real(kind = 8) :: myCourant   !Maximum Courant number for time-step 
                                    !control 
    real(kind = 8) :: myuSum      !sum of values on single cpu
    real(kind = 8) :: uSum        !sum all
    real(kind = 8) :: myvSum      !sum of values on single cpu
    real(kind = 8) :: vSum        !sum all
    real(kind = 8) :: mywSum      !sum of values on single cpu
    real(kind = 8) :: wSum        !sum all
    real(kind = 8) :: kk          !|k|^2
    real(kind = 8) :: phase       !phase for initiation
    real(kind=8), dimension(:), allocatable     :: kSpec      ! |k| array for  
    real(kind=8), dimension(:), allocatable     :: myEspec, Espec ! Energy 
                                                                  ! spectrum
    real(kind = 8) :: absk        ! |k|
    integer        :: Nspec       ! length of the energy spectrum arrays
    integer        :: nk          ! position of the window corresponding to k
    real(kind = 8) :: myDisp      ! Dissipation on one cpu
    real(kind = 8) :: Disp        ! Total dissipation
    real(kind = 8) :: mydivMax      ! Maximum divergence for error control
    real(kind = 8) :: divMax      ! Maximum divergence for error control
    real(kind = 8) :: myEZero     ! Energy contained in k=0
    real(kind = 8) :: EZero       ! Energy contained in k=0
    real(kind = 8) :: myEOne     ! Energy contained in k=1
    real(kind = 8) :: EOne       ! Energy contained in k=1
    real(kind = 8) :: myEsqrt2     ! Energy contained in k=sqrt(2)
    real(kind = 8) :: Esqrt2       ! Energy contained in k=sqrt(2)
    
    real(kind = 8) :: myError     ! For tests
    real(kind = 8) :: maxError    ! For tests
    
    contains
    
    subroutine var_init()
        
        Lx = 2 * pi / alpha_x
        Ly = 2 * pi / alpha_y
        Lz = 2 * pi / alpha_z
        
        ! Cubic truncation:
!        Nspec = int(sqrt(((real(Nx, kind=8) / 2.0d0) &
!                         * (2.0d0 * alpha_x / 3.0d0)) ** 2  + &
!                         ((real(Ny, kind=8) / 2.0d0) &
!                         * (2.0d0 * alpha_y / 3.0d0)) ** 2  + &
!                         ((real(Nz, kind=8) / 2.0d0) &
!                         * (2.0d0 * alpha_z / 3.0d0)) ** 2) / Deltak) - 1
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
        
        if(proc_id.eq.0) then
            
            print *, 'Nspec = ', Nspec
            
        end if
        
        !1D decomposition:
        ndim = 1
        dims(1) = 1
        dims(2) = nproc
        
        ! Add 2D decomposition setup here in future, if necessary

        !broadcast ndim
        call mpi_bcast (ndim, 1, MPI_INTEGER, 0, mpi_comm_world, ierr)
            
        iproc = dims(1)
        jproc = dims(2)
        
        if(proc_id .eq. 0) then
            print *, 'Using processor grid ', iproc, ' x ', jproc
        endif
        
        ! Pruned dimensions (see p3dfft user guide)
        nxc = Nx
        nyc = Ny
        nzc = Nz
        
        ! Set up p3dfft:
        call p3dfft_setup (dims, Nx, Ny, Nz, MPI_COMM_WORLD, nxc, nyc, nzc, .true.)
        ! Get dimensions of the current cpu:
        call p3dfft_get_dims (istart, iend, isize, 1)
        call p3dfft_get_dims (fstart, fend, fsize, 2)
        
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
        allocate (u(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  v(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  w(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  ux(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  uy(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  uz(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  vx(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  vy(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  vz(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  wx(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  wy(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  wz(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  omegax(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  omegay(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  omegaz(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
!                  uold(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
!                  vold(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
!                  wold(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  utemp(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  vtemp(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  wtemp(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  temp_r(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
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
