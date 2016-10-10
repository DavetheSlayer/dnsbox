! Program to integrate Navier-Stokes equations in a cubic box in 
! parallel using p3dfft
! based on: https://github.com/openmichigan/PSNM/blob/master/NavierStokes/
!           Programs/NavierStokes3dFortranMPI/NavierStokes3DfftIMR.f90
!
! part of: https://en.wikibooks.org/wiki/Parallel_Spectral_Numerical_Methods

program nsbox
    
    use rhs
    use state
    use io
    
    ! Initialize:
    call init_mpi()
    call var_init()
    call io_init()
    
    call io_saveInfo()
    
    time(1) = 0.0d0
    
    if(proc_id .eq. 0) then 
        print *, 'Run notes: using rotational form'
    endif
    
    if( initrand ) then 
        
        ! Generate a random initial state
!        call stateInitRand()
        call state_init_rand()
        
    elseif ( analytic ) then 
        
        ! Set the analytic solution as initial condition
        call state_init_analytic()
        
    elseif ( eig ) then 
        
        ! Set the eigenvector as initial condition
        ! call stateInitEig()
        if(proc_id .eq. 0) then 
            print *, 'Eigenvector initiation is not implemented.'
        endif        

    else
        ! Load initial state:
        call io_load_state('state0000.h5')
    end if
    
    ! Transform to Fourier space
    call state_u2uhat()
    
    call state_dealias()
    call state_project()
    
    ! copy uhat -> uhattemp
    call state_uhat2uhattemp()  
    
    call io_save_spectrum()
    call io_save_stats()
    call io_save_state()
    
    if(proc_id .eq. 0) then 
        print *, 'Starting time-stepping'
    endif
    
    
    if (tStepFix) then
        call rhs_int_fact() ! compute integration factor    
    else 
        call set_time_step()
    end if
    
    open (99, file='RUNNING')
         write(99,*) 'This file indicates a job running in this directory.'
         write(99,*) 'Delete this file to cleanly terminate the process.'
    close(99)
    inquire(file='RUNNING', exist=running_exist)
    
    do while (running_exist)
    n = n + 1
        
        ! compute the nonlinear term for uhattemp:
        call rhs_nonlinear()

        ! Predictor
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
            
            ! Predicted next step for corrector calculation:
            uhattemp(i, j, k) = intFact(i, j, k) &
                              * (uhat(i, j, k) + dt * nonlinuhat(i, j, k))
            
            ! Contribution to the final step from the predictor calculation:
            uhat(i, j, k) = intFact(i, j, k) &
                          * (uhat(i, j, k) + dt * nonlinuhat(i, j, k) * 0.5)
            
            ! Predicted next step for corrector calculation:
            vhattemp(i, j, k) = intFact(i, j, k) &
                              * (vhat(i, j, k) + dt * nonlinvhat(i, j, k))
            
            ! Contribution to the final step from the predictor calculation:
            vhat(i, j, k) = intFact(i, j, k) &
                          * (vhat(i, j, k) + dt * nonlinvhat(i, j, k) * 0.5)
                        
            ! Predicted next step for corrector calculation:
            whattemp(i, j, k) = intFact(i, j, k) &
                              * (what(i, j, k) + dt * nonlinwhat(i, j, k))
            
            ! Contribution to the final step from the predictor calculation:
            what(i, j, k) = intFact(i, j, k) &
                          * (what(i, j, k) + dt * nonlinwhat(i, j, k) * 0.5)
                        
        end do; end do; end do    
        
        ! compute the nonlinear term for uhattemp:
        call rhs_nonlinear()
        
        ! Corrector
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
        
            uhat(i, j, k) = uhat(i, j, k) + dt * nonlinuhat(i, j, k) * 0.5
            vhat(i, j, k) = vhat(i, j, k) + dt * nonlinvhat(i, j, k) * 0.5
            what(i, j, k) = what(i, j, k) + dt * nonlinwhat(i, j, k) * 0.5
                        
        end do; end do; end do    
 
        time(n+1) = time(n) + dt
       
        if (proc_id.eq.0) then
            print *, 'time', time(n+1)
        end if
        
        if(modulo(n,iSaveRate1)==0) then
                
            ! Back to the configuration space:
            call state_uhat2u()  ! messes up uhattemp!!!
            
            call io_save_spectrum()
            call io_save_state()
        end if
               
        ! Copy uhat2uhattemp for the next step        
        call state_uhat2uhattemp()
                
        if(modulo(n,iSaveRate2)==0) then
            call io_save_stats()
        end if

                
        ! Terminate if maximum time steps are reached or the Courant number 
        ! exceeds the allowed limit:
        if ((tStepFix .eqv. .false.) .and. & 
           ((Courant .gt. CourantMax) .or. &
            (Courant .lt. CourantMin .and. dt .lt. tStepMax))) then
            ! Courant number outside the desired range, change the time step:
            call set_time_step()
            
        else if (tStepFix .and. Courant .gt. CourantMax) then
            ! Courant number exceeded maximum allowed value, terminate the run:
            if (proc_id.eq.0) then
                print *, ' Courant number exceeded the maximum allowed value '
                print *, ' terminating run ... '
            end if            
            
            inquire(file='RUNNING', exist=running_exist)
            if(running_exist) open(99, file='RUNNING')
            if(running_exist) close(99, status='delete')
            
            ! Save final state and spectrum:
            call io_save_spectrum()
            call io_save_state()
        
        else if (n .eq. Nt) then
            if (proc_id.eq.0) then
                print *, ' Maximum time steps are reached'
                print *, ' terminating run ... '
            end if            
            
            inquire(file='RUNNING', exist=running_exist)
            if(running_exist) open(99, file='RUNNING')
            if(running_exist) close(99, status='delete')
            
            ! Save final state and spectrum:
            call io_save_spectrum()
            call io_save_state()
        end if
        
        inquire(file='RUNNING', exist=running_exist)
            
    end do
    
    if(analytic) then
        call state_check_error()
    end if
    
	deallocate(x, y, z, time, mychg, allchg, kSpec, myEspec, Espec, &
               u, v, w, ux, uy, uz, vx, vy, vz, wx, wy, wz, &
               utemp, vtemp, wtemp,&
               temp_r, kx, ky, kz, uhat, vhat, what,&
               uhattemp, vhattemp, whattemp,&
               uhatold, vhatold, whatold,&
               nonlinuhat, nonlinvhat, nonlinwhat, temp_c,&
               intFact, phat, stat=AllocateStatus)		
	if (AllocateStatus .ne. 0) stop
	if (proc_id.eq.0) then
		print *,'Program execution complete'
	end if
    
    ! Finalize
    call h5close_f(io_error)
    call p3dfft_clean
    call io_finalize()
    call MPI_FINALIZE (ierr)
        
contains

subroutine init_mpi()
    
    !initialize mpi:
    call mpi_init (ierr)
    call mpi_comm_size (mpi_comm_world, nproc, ierr)
    call mpi_comm_rank (mpi_comm_world, proc_id, ierr)

end subroutine init_mpi

subroutine set_time_step()
    ! Change the time-step such that the Courant number is set to 
    ! (CourantMin+CourantMax)/2
    ! This subroutine should be called after the stats are computed
    
    dt = ((CourantMax + CourantMin) / (2.0d0 * Courant) ) * dt
    
    if (dt > tStepMax) then
        dt = tStepMax
    end if
    
    call rhs_int_fact() ! Recompute the integration factor with the new time step
    call io_Courant()
	if (proc_id.eq.0) then
		print *,'time step set to', dt
        print *, 'new Courant number', Courant
	end if
    
        
end subroutine set_time_step

end program nsbox
