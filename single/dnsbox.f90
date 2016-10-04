! Program to integrate Navier-Stokes equations in a cubic box

program nsbox
    
    use state

    ! Initialize:
    call var_init()
    call state_init_rand()
    
    do k=1,Nz; do j=1,Ny; do i=1,Nx
                    
        temp_r(i, j, k) = u(i, j, k) 
        
    end do; end do; end do                    
        
    call state_u2uhat()
    call state_uhat2u()
    
    
        
    call var_final()
   
end program nsbox
