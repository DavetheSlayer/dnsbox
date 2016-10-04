module variables
    
    use parameters
    
    implicit none

    !simulation variables:
    real(kind=8)            :: dt= 0.001d0             ! Time step size
    real(kind=8), dimension(:), allocatable     :: x, y, z, time, allchg
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
                                                        nonlinuhat, &
                                                        nonlinvhat, &
                                                        nonlinwhat, &
                                                        temp_c, &
                                                        phat
    
    ! extraneous variables, for easy code reading:
    real(kind = 8) :: Lx, Ly, Lz
        
    ! counters and logicals:
    integer(kind=4) :: ind, iter, i, j, k, n, t, AllocateStatus
    logical         :: running_exist
      
    ! measurement variables
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
    real(kind = 8) :: Ezero       ! Energy contained in k=0
    
    ! fftw variables
    integer ( kind = 8 ) plan_forward_u, plan_forward_v, plan_forward_w
    integer ( kind = 8 ) plan_backward_u, plan_backward_v, plan_backward_w
    
    contains
    
    subroutine var_init()
        
        include "fftw3.f"
        
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
                  time(1:Nt+1), allchg(1:3), &
                  kSpec(1:Nspec), Espec(1:Nspec), &
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
                  uhat(1:Nh, 1:Ny, 1:Nz), &
                  vhat(1:Nh, 1:Ny, 1:Nz), &
                  what(1:Nh, 1:Ny, 1:Nz), &
                  uhattemp(1:Nh, 1:Ny, 1:Nz), &
                  vhattemp(1:Nh, 1:Ny, 1:Nz), &
                  whattemp(1:Nh, 1:Ny, 1:Nz), &
                  nonlinuhat(1:Nh, 1:Ny, 1:Nz), &
                  nonlinvhat(1:Nh, 1:Ny, 1:Nz), &
                  nonlinwhat(1:Nh, 1:Ny, 1:Nz), &
                  temp_c(1:Nh, 1:Ny, 1:Nz), &
                  phat(1:Nh, 1:Ny, 1:Nz), &
                  stat=AllocateStatus)
                    
        if(AllocateStatus .ne. 0) then
            print *, 'Error ', AllocateStatus, ' allocating simulation arrays'
            stop
        end if
        

        print *, 'Successfully allocated memory'
        
        ! Fourier frequencies in x-direction
        do i = 1, Nx/2 + 1
            kx(i) = real(i - 1, kind=8) * alpha_x
        end do
        kx(1 + Nx/2) = 0.0d0
        
        do i = 1, Nx/2 - 1
            kx(i + 1 + Nx/2) = -kx(1 - i + Nx/2)
        end do

        ind = 1
        do i = -Nx/2, Nx/2 - 1
            x(ind) = real(i, kind=8) * Lx / real(Nx, kind(0d0))
            ind = ind + 1
        end do
        
        ! Fourier frequencies in y-direction
        do i = 1, Ny/2 + 1
            ky(i) = real(i - 1, kind=8) * alpha_y
        end do
        ky(1 + Ny/2) = 0.0d0
        do i = 1, Ny/2 - 1
            ky(i + 1 + Ny/2) = -ky(1 - i + Ny/2)
        end do
        ind = 1
        do i = -Ny/2, Ny/2 - 1
            y(ind) = real(i, kind=8) * Ly / real(Ny, kind(0d0))
            ind = ind + 1
        end do
            
        ! Fourier frequencies in z-direction
        do i = 1, Nz/2 + 1
            kz(i) = real(i - 1, kind=8) * alpha_z
        end do
        kz(1 + Nz/2) = 0.0d0
        do i = 1, Nz/2 - 1
            kz(i + 1 + Nz/2) = -kz(1 - i + Nz/2)
        end do
        ind = 1
        do i = -Nz/2, Nz/2 - 1
            z(ind) = real(i, kind=8) * Lz / real(Nz, kind(0d0))
            ind = ind + 1
        end do
        
        do i = 1, Nspec
            kSpec(i) = i * Deltak
        end do
        
        scalemodes = 1.0d0 / real(Nx * Ny * Nz, kind=8)
        
        print *, 'Setup grid and Fourier frequencies'
        
        ! Plan ffts:
        
        call dfftw_plan_dft_r2c_3d_(plan_forward_u, nx, ny, nz, & 
                                    utemp, uhattemp, FFTW_ESTIMATE)
        call dfftw_plan_dft_r2c_3d_(plan_forward_v, nx, ny, nz, & 
                                    vtemp, vhattemp, FFTW_ESTIMATE)
        call dfftw_plan_dft_r2c_3d_(plan_forward_w, nx, ny, nz, &
                                    wtemp, whattemp, FFTW_ESTIMATE)

        call dfftw_plan_dft_c2r_3d_(plan_backward_u, nx, ny, nz, &
                                    uhattemp, utemp, FFTW_ESTIMATE)
        call dfftw_plan_dft_c2r_3d_(plan_backward_v, nx, ny, nz, &
                                    vhattemp, vtemp, FFTW_ESTIMATE)
        call dfftw_plan_dft_c2r_3d_(plan_backward_w, nx, ny, nz, &
                                    whattemp, wtemp, FFTW_ESTIMATE)
        
    end subroutine var_init
    
    subroutine var_final()
        
        call dfftw_destroy_plan_(plan_forward_u)
        call dfftw_destroy_plan_(plan_forward_v)
        call dfftw_destroy_plan_(plan_forward_w)
        call dfftw_destroy_plan_(plan_backward_u)
        call dfftw_destroy_plan_(plan_backward_v)
        call dfftw_destroy_plan_(plan_backward_w)

        deallocate(x, y, z, time, allchg, kSpec, Espec, &
                   u, v, w, ux, uy, uz, vx, vy, vz, wx, wy, wz, &
                   omegax, omegay, omegaz, &
                   utemp, vtemp, wtemp,&
                   temp_r, kx, ky, kz, uhat, vhat, what,&
                   uhattemp, vhattemp, whattemp,&
                   nonlinuhat, nonlinvhat, nonlinwhat, temp_c,&
                   phat, stat=AllocateStatus)		
               
        if (AllocateStatus .ne. 0) stop
		print *,'Successfully deallocated memory'
		print *,'Program execution complete'
	
    end subroutine var_final
    
end module variables
