! Program to integrate Navier-Stokes equations in a cubic box

program nsbox
    
    use rhs
    use io

    ! Initialize:
    call var_init()
    call io_init()

    if (analytic) then
        
        ! Set initial condition from analytic solution
        call state_init_analytic()
        call state_u2uhat()
        call state_dealias()
        
        print *, 'Initial condition set to analytic solution'
        
    elseif(initrand) then
        
        ! Random initial condition
        call state_init_rand()
        
    else
        
        ! Load initial condition:
        call io_load_state('state0000.nc')
        call state_u2uhat()
        call state_dealias()
    
    end if
    
    time(1) = 0.0d0
    n = 0
    call state_uhat2uhattemp()
    
    ! Save initial state:
    call io_save_stats()
    call io_save_spectrum()
    call io_save_state()
    
    open (99, file='RUNNING')
         write(99,*) 'This file indicates a job running in this directory.'
         write(99,*) 'Delete this file to cleanly terminate the process.'
    close(99)
    inquire(file='RUNNING', exist=running_exist)    
    
    do while (running_exist)
    n = n + 1    
            
        ! print *, 't = ', time(n)
        call rhs_nonlinear() ! Compute nonlinear term for initial state
!        call state_kinetic()
!        print *, 'k = ', Ekin 
        
        
        ! Integrating factor and predictor:
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
        
        call rhs_nonlinear() ! Compute nonlinear term for predicted state
        
        ! Corrector
        do i=1,Nh; do j=1,Ny; do k=1,Nz
            
            uhat(i, j, k) = uhat(i, j, k) + dt * nonlinuhat(i, j, k) * 0.5
            vhat(i, j, k) = vhat(i, j, k) + dt * nonlinvhat(i, j, k) * 0.5
            what(i, j, k) = what(i, j, k) + dt * nonlinwhat(i, j, k) * 0.5
            
        end do; end do; end do
        
        call state_project()
                
        time(n+1) = time(n) + dt
        
        if(modulo(n,iSaveRate2)==0) then
            
            call state_uhat2u()             ! Back to the configuration space
            call state_uhat2uhattemp()      ! bfft messes up uhattemp, restore
            call io_save_stats() ! Save statistics
        
        end if
        
        if(modulo(n,iSaveRate1)==0) then
            
            call io_save_spectrum() ! Save spectrum
            call io_save_state()    ! Save state
        
        end if

        call state_uhat2uhattemp()
        
        ! Terminate if maximum time steps are reached or the Courant number 
        ! exceeds the allowed limit:
        if ((tStepFix .eqv. .false.) .and. & 
           ((Courant .gt. CourantMax) .or. &
            (Courant .lt. CourantMin .and. dt .lt. tStepMax))) then
            ! Courant number outside the desired range, change the time step:
            call set_time_step()
            
        else if (tStepFix .and. Courant .gt. CourantMax) then
            ! Courant number exceeded maximum allowed value, terminate the run:
            print *, ' Courant number exceeded the maximum allowed value '
            print *, ' terminating run ... '
            
            inquire(file='RUNNING', exist=running_exist)
            if(running_exist) open(99, file='RUNNING')
            if(running_exist) close(99, status='delete')
            
            ! Save final state and spectrum:
            call io_save_spectrum()
            call io_save_state()
        
        else if (n .eq. Nt) then
            print *, ' Maximum time steps are reached'
            print *, ' terminating run ... '
            
            inquire(file='RUNNING', exist=running_exist)
            if(running_exist) open(99, file='RUNNING')
            if(running_exist) close(99, status='delete')
            
            ! Save final state and spectrum:
            call io_save_spectrum()
            call io_save_state()
        end if
        
        inquire(file='RUNNING', exist=running_exist)        
        
    end do
    
    ! If analytic solution then check the final error:
    if (analytic) call state_check_error()
    
    call var_final()
    
    contains
    
    subroutine set_time_step()
        ! Change the time-step such that the Courant number is set to 
        ! (CourantMin+CourantMax)/2
        ! This subroutine should be called after the stats are computed
    
        dt = ((CourantMax + CourantMin) / (2.0d0 * Courant) ) * dt
    
        if (dt > tStepMax) then
            dt = tStepMax
        end if
    
        call io_Courant()
        print *,'time step set to', dt
        print *, 'new Courant number', Courant
	
        
end subroutine set_time_step
    
end program nsbox
