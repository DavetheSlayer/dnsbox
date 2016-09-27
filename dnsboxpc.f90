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
    call initMpi()
    call var_init()
    call io_init()
    
    call io_saveInfo()
    
    time(1) = 0.0d0
    
    if( initrand ) then 
        
        ! Generate a random initial state
        call stateInitRand()
        
    elseif ( analytic ) then 
        
        ! Set the analytic solution as initial condition
        call stateInitAnalytic()
        
    else
        ! Load initial state:
        call io_loadState('state0000.h5')
    end if
        
    call p3dfft_ftran_r2c (u, uhat, 'fft')
    call p3dfft_ftran_r2c (v, vhat, 'fft')
    call p3dfft_ftran_r2c (w, what, 'fft')

    call stateDealias()

    ! copy uhat -> uhatold and uhattemp
    do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
        
        uhattemp(i, j, k) = uhat(i, j, k)
        vhattemp(i, j, k) = vhat(i, j, k)
        whattemp(i, j, k) = what(i, j, k)
        
        uhatold(i, j, k) = uhat(i, j, k)
        vhatold(i, j, k) = vhat(i, j, k)
        whatold(i, j, k) = what(i, j, k)
        
    end do; end do; end do
    
    call io_saveSpectrum()
    call io_saveStats()
    call io_saveState()
    
    if(proc_id .eq. 0) then 
        print *, 'Starting time-stepping'
    endif
    
    
    if (tStepFix) then
        call rhsIntFact() ! compute integration factor    
    else 
        call setTimeStep()
    end if
    
    open (99, file='RUNNING')
         write(99,*) 'This file indicates a job running in this directory.'
         write(99,*) 'Delete this file to cleanly terminate the process.'
    close(99)
    inquire(file='RUNNING', exist=running_exist)
    
    do while (running_exist)
    n = n + 1
        
        call rhsNonlinear()
            
            ! Predictor
            do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
                if((kx(k) .eq. cmplx(0.0d0, 0.0d0)) .and. &
                   (ky(j) .eq. cmplx(0.0d0, 0.0d0)) .and. &
                   (kz(i) .eq. cmplx(0.0d0, 0.0d0))) then
                    eps = 0.1d0 ** 13
                else
                    eps = 0.0d0
                end if                
                
                linearterm = (nu * (kx(k) * kx(k) &
                                  + ky(j) * ky(j) &
                                  + kz(i) * kz(i)) + Q)
                ! Predicted next step for corrector calculation:
                uhattemp(i, j, k) = (1.0d0 / dt + (1 - c) * linearterm) &
                                  * uhatold(i, j, k) + nonlinuhat(i, j, k) &
                                  / (1.0d0 / dt - c * (linearterm + eps))
                            
                vhattemp(i, j, k) = (1.0d0 / dt + (1 - c) * linearterm) &
                                  * vhatold(i, j, k) + nonlinvhat(i, j, k) &
                                  / (1.0d0 / dt - c * (linearterm + eps))
                            
                whattemp(i, j, k) = (1.0d0 / dt + (1 - c) * linearterm) &
                                  * whatold(i, j, k) + nonlinwhat(i, j, k) &
                                  / (1.0d0 / dt - c * (linearterm + eps))
                            
            end do; end do; end do    
        
        chg = 1.0d0
        iter = 0.0d0
        do while (chg .gt. tol)
            ! predictor-corrector iterations
            uhat(i, j, k) 
            
        
        end do
        
        call stateProject()
        call stateDealias()
        ! Copy
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
        
            uhattemp(i, j, k) = uhat(i, j, k)
            vhattemp(i, j, k) = vhat(i, j, k)
            whattemp(i, j, k) = what(i, j, k)
            
            if (timestepper .eq. 1) then
            
                uhatold(i, j, k) = uhat(i, j, k) 
                vhatold(i, j, k) = vhat(i, j, k) 
                whatold(i, j, k) = what(i, j, k) 
                
            elseif (timestepper .eq. 2 .and. n .gt. 1) then
                
                nonlinuhatold(i, j, k) = nonlinuhat(i, j, k)
                nonlinvhatold(i, j, k) = nonlinvhat(i, j, k)
                nonlinwhatold(i, j, k) = nonlinwhat(i, j, k)
                
            end if
                
        end do; end do; end do    

        time(n+1) = time(n) + dt
        
        if (proc_id.eq.0) then
            print *, 'time', time(n+1)
        end if
      
        if(modulo(n,iSaveRate1)==0) then
          
            ! Back to the configuration space:
            call p3dfft_btran_c2r (uhat, u, 'tff')
            call p3dfft_btran_c2r (vhat, v, 'tff')
            call p3dfft_btran_c2r (what, w, 'tff')
            
            do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
                
                u(i, j, k) = u(i, j, k) * scalemodes
                v(i, j, k) = v(i, j, k) * scalemodes
                w(i, j, k) = w(i, j, k) * scalemodes
                
            end do; end do; end do
            
            call io_saveSpectrum()
            call io_saveState()
        end if
        
        if(modulo(n,iSaveRate2)==0) then
            call io_saveStats()
        end if
                
        ! Terminate if maximum time steps are reached or the Courant number 
        ! exceeds the allowed limit:
        if ((tStepFix .eqv. .false.) .and. & 
           ((Courant .gt. CourantMax) .or. (Courant .lt. CourantMin))) then
            ! Courant number outside the desired range, change the time step:
            call setTimeStep()
            
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
            call io_saveSpectrum()
            call io_saveState()
        
        else if (n .eq. Nt) then
            if (proc_id.eq.0) then
                print *, ' Maximum time steps are reached'
                print *, ' terminating run ... '
            end if            
            
            inquire(file='RUNNING', exist=running_exist)
            if(running_exist) open(99, file='RUNNING')
            if(running_exist) close(99, status='delete')
            
            ! Save final state and spectrum:
            call io_saveSpectrum()
            call io_saveState()
        end if
        
        inquire(file='RUNNING', exist=running_exist)
            
    end do
    
    if (analytic) then
        call stateCheckError()
    end if 
    
	deallocate(x, y, z, time, mychg, allchg, kSpec, myEspec, Espec, &
               u, v, w, ux, uy, uz, vx, vy, vz, wx, wy, wz, &
               omegax, omegay, omegaz, &
               utemp, vtemp, wtemp,&
               temp_r, kx, ky, kz, uhat, vhat, what,&
               uhattemp, vhattemp, whattemp,&
               uhatold, vhatold, whatold,&
               nonlinuhat, nonlinvhat, nonlinwhat, temp_c,&
               nonlinuhatold, nonlinvhatold, nonlinwhatold, &
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

subroutine initMpi()
    
    !initialize mpi:
    call mpi_init (ierr)
    call mpi_comm_size (mpi_comm_world, nproc, ierr)
    call mpi_comm_rank (mpi_comm_world, proc_id, ierr)

end subroutine initMpi

subroutine setTimeStep()
    ! Change the time-step such that the Courant number is set to 
    ! (CourantMin+CourantMax)/2
    ! This subroutine should be called after the stats are computed
    
    dt = ((CourantMax + CourantMin) / (2.0d0 * Courant) ) * dt
    call rhsIntFact() ! Recompute the integration factor with the new time step
    call io_Courant()
	if (proc_id.eq.0) then
		print *,'time step set to', dt
        print *, 'new Courant number', Courant
	end if
    
        
end subroutine setTimeStep

end program nsbox
