! Program to integrate Navier-Stokes equations in a cubic box

program nsbox
    
    use rhs

    ! Initialize:
    call var_init()

    if (analytic) then
    
        call state_init_analytic()
        call state_u2uhat()
        call state_dealias()
        
    elseif(initrand) then
        
        call state_init_rand()
    
    else
    
        print *, "Load state has not implemented yet"
        stop
    
    end if
    
    time(1) = 0.0d0
    
    do n = 1, Nt
        
        print *, 't = ', time(n)
        call state_uhat2uhattemp()
        call rhs_nonlinear()
        call state_kinetic()
        print *, 'k = ', Ekin
        
        do i=1,Nh; do j=1,Ny; do k=1,Nz
            
            intfact = exp(-nu * dt * &
                          (kx(i) * kx(i) + ky(j) * ky(j) + kz(k) * kz(k)))
            
            uhattemp(i, j, k) = intfact &
                              * (uhat(i, j, k) + dt * nonlinuhat(i, j, k))
            uhat(i, j, k) = intfact &
                          * (uhat(i, j, k) + dt * nonlinuhat(i, j, k) * 0.5)
            
            vhattemp(i, j, k) = intfact &
                              * (vhat(i, j, k) + dt * nonlinvhat(i, j, k))
            vhat(i, j, k) = intfact &
                          * (vhat(i, j, k) + dt * nonlinvhat(i, j, k) * 0.5)
            
            whattemp(i, j, k) = intfact &
                              * (what(i, j, k) + dt * nonlinwhat(i, j, k))
            what(i, j, k) = intfact &
                          * (what(i, j, k) + dt * nonlinwhat(i, j, k) * 0.5)
            
        end do; end do; end do     
        
        call rhs_nonlinear()
        
        do i=1,Nh; do j=1,Ny; do k=1,Nz
            
            uhat(i, j, k) = uhat(i, j, k) + dt * nonlinuhat(i, j, k) * 0.5
            vhat(i, j, k) = vhat(i, j, k) + dt * nonlinvhat(i, j, k) * 0.5
            what(i, j, k) = what(i, j, k) + dt * nonlinwhat(i, j, k) * 0.5
            
        end do; end do; end do
        
        time(n+1) = n * dt
        
    end do
    
    call state_check_error()
    
    call var_final()
   
end program nsbox
