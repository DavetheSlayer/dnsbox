module variables
    
    use parameters
    use p3dfft
    implicit none
    include 'mpif.h'

    !simulation variables:
    real(kind=8), dimension(:), allocatable     :: x, y, z, time, mychg, allchg
    
    real(p3dfft_type), dimension(:, :, :), allocatable :: u, v, w, &
                                                          ux, uy, uz, &
                                                          vx, vy, vz, &
                                                          wx, wy, wz, &
                                                   uold, uxold, uyold, uzold, &
                                                   vold, vxold, vyold, vzold, &
                                                   wold, wxold, wyold, wzold, &
                                                   utemp, vtemp, wtemp, temp_r
    complex(kind=8), dimension(:), allocatable      :: kx, ky, kz
    complex(p3dfft_type), dimension(:, :, :), allocatable :: uhat, vhat, what,&
                                                             rhsuhatfix, &
                                                             rhsvhatfix, &
                                                             rhswhatfix, &
                                                             nonlinuhat, &
                                                             nonlinvhat, &
                                                             nonlinwhat, &
                                                             phat, temp_c
    real(p3dfft_type), dimension(:,:,:), allocatable    :: realtemp

        
    !counters and logicals:
    integer(kind=4) :: count, iol, i, j, k, n, ind, t, AllocateStatus
  
    !parallelization variables:
    integer ierr, nu, ndim, dims(2), nproc, proc_id
    integer istart(3), iend(3), isize(3)
    integer fstart(3), fend(3), fsize(3)
    integer iproc, jproc, nxc, nyc, nzc
    logical iex
	integer memsize(3)
    
    contains
    
    subroutine var_init()
    
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
                  time(1:Nt+1), mychg(1:3),allchg(1:3), stat=AllocateStatus)
                  
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
                  uold(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  vold(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  wold(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  uxold(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  uyold(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  uzold(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  vxold(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  vyold(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  vzold(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  wxold(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  wyold(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  wzold(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  utemp(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  vtemp(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  wtemp(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  temp_r(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  realtemp(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &
                  uhat(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), &
                  vhat(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), &
                  what(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), &
                  rhsuhatfix(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), &
                  rhsvhatfix(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), &
                  rhswhatfix(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), &
                  nonlinuhat(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), &
                  nonlinvhat(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), &
                  nonlinwhat(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), &
                  phat(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), &
                  temp_c(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), &
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
            kx(i) = cmplx(0.0d0, 1.0d0) * real(i - 1, kind(0d0)) / Lx
        end do
        kx(1 + Nx/2) = 0.0d0
        
        do i = 1, Nx/2 - 1
            kx(i + 1 + Nx/2) = -kx(1 - i + Nx/2)
        end do

        ind = 1
        do i = -Nx/2, Nx/2 - 1
            x(ind) = 2.0d0 * pi * real(i, kind(0d0)) * Lx / real(Nx, kind(0d0))
            ind = ind + 1
        end do
        
        ! Fourier frequencies in y-direction
        do i = 1, Ny/2 + 1
            ky(i) = cmplx(0.0d0, 1.0d0) * real(i - 1, kind(0d0)) / Ly
        end do
        ky(1 + Ny/2) = 0.0d0
        do i = 1, Ny/2 - 1
            ky(i + 1 + Ny/2) = -ky(1 - i + Ny/2)
        end do
        ind = 1
        do i = -Ny/2, Ny/2 - 1
            y(ind) = 2.0d0 * pi * real(i, kind(0d0)) * Ly / real(Ny, kind(0d0))
            ind = ind + 1
        end do
            
        ! Fourier frequencies in z-direction
        do i = 1, Nz/2 + 1
            kz(i) = cmplx(0.0d0, 1.0d0) * real(i - 1, kind(0d0)) / Lz
        end do
        kz(1 + Nz/2) = 0.0d0
        do i = 1, Nz/2 - 1
            kz(i + 1 + Nz/2) = -kz(1 - i + Nz/2)
        end do
        ind = 1
        do i = -Nz/2, Nz/2 - 1
            z(ind) = 2.0d0 * pi * real(i, kind(0d0)) * Lz / real(Nz, kind(0d0))
            ind = ind + 1
        end do
        
        scalemodes = 1.0d0 / real(Nx * Ny * Nz, kind(0d0))
        
        if (proc_id .eq. 0) then
            print *, 'Setup grid and Fourier frequencies'
        end if    
        
    end subroutine var_init
    
end module variables
