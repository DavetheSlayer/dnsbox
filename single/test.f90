! Program to integrate Navier-Stokes equations in a cubic box

program nsbox
    
    use rhs
    real(kind=8) :: a, b
        
    ! Initialize:
    call var_init()
    call state_init_rand()
    
    do k=1,Nz; do j=1,Ny; do i=1,Nx
                    
        temp_r(i, j, k) = u(i, j, k) 
        
    end do; end do; end do                    
        
    call state_u2uhat()
    call state_uhat2u()
    
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
    
    call var_final()
    
    contains 
    
    subroutine copy(c, d)
        real(kind=8) :: c, d
        d = c
        
    end subroutine copy
    
end program nsbox
