! Program to integrate Navier-Stokes equations in a cubic box

program nsbox
    
    use io
    real(kind=8) :: a, b
        
    ! Initialize:
    call var_init()
    call state_init_rand()
    
    do k=1,Nz; do j=1,Ny; do i=1,Nx
                    
        temp_r(i, j, k) = u(i, j, k) 
        
    end do; end do; end do                    
    
    do n=1,10    
        call state_u2uhat()
        call state_uhat2u()
    end do
    
    do k=1,Nz; do j=1,Ny; do i=1,Nx
                    
        temp_r(i, j, k) = temp_r(i, j, k) - u(i, j, k)
        
    end do; end do; end do                        
    
    chg = maxval(abs(temp_r))
    
    print *, 'maximum error after fwd/bwd fft: ', chg
    
    a = 5.0
    b = 3.0
    
    print *, 'b = ', b
    call copy(a, b)
    
    print *, 'b = ', b
    
    call io_save_state()
    
    call state_u2utemp()
    
    call io_load_state('state0000.nc')
    
    do k=1,Nz; do j=1,Ny; do i=1,Nx
                    
        temp_r(i, j, k) = vtemp(i, j, k) - v(i, j, k)
        
    end do; end do; end do                        
    
    chg = maxval(abs(temp_r))
    
    print *, 'maximum error after save/load: ', chg    
    
    ! call state_project()
    call state_uhat2uhattemp()
    call state_derivatives()
    call state_divergence()
    
    print *, 'maximum divergence: ', divMax
    
    if (kx(1) .eq. 0.0d0) then
        print *, 'kx(1) = 0'
    end if
    
    if (kx(1) .eq. real(0.0d0, kind=8)) then
        print *, 'kx(1) .eq. (0.0d0, kind=8)'
    end if
    
    call var_final()
    
    contains 
    
    subroutine copy(c, d)
        real(kind=8) :: c, d
        d = c
        
    end subroutine copy
    
end program nsbox
