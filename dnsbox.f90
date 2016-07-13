! Program to integrate Navier-Stokes equations in a cubic box in 
! parallel using p3dfft
! based on: https://github.com/openmichigan/PSNM/blob/master/NavierStokes/
!           Programs/NavierStokes3dFortranMPI/NavierStokes3DfftIMR.f90
!
! part of: https://en.wikibooks.org/wiki/Parallel_Spectral_Numerical_Methods

program nsbox
    
    use rhs
    use io
    
    ! Initialize:
    call initMpi()
    call var_init()
    call io_init()
    
    call io_saveInfo()
    
    time(1) = 0.0d0
    
    if( initrand ) then ! Set true to initiate simulation from a random field
        ! Initiate fields in momentum space:
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)    
            ! skip 000 mode:
            if((kx(k) .eq. cmplx(0.0d0, 0.0d0)) .and. &
               (ky(j) .eq. cmplx(0.0d0, 0.0d0)) .and. &
               (kz(i) .eq. cmplx(0.0d0, 0.0d0))) cycle 
            
            kk = real(conjg(kx(k)) * kx(k) &
                    + conjg(ky(j)) * ky(j) &
                    + conjg(kz(i)) * kz(i)) ! |k|^2
            
            call random_number(phase)
            phase = phase * 2 * pi ! a random phase for the Fourier mode
            ! Assign a gaussian distributed amplitude and a random phase to the 
            ! Fourier mode:
            uhat(i, j, k) = (scalemodes**(-1)) * 1.0d-1 &
                          * sqrt(exp(-2.0d0 * kk / 3.0d2)) * exp(cmplx(0.0d0, phase))
            call random_number(phase)
            phase = phase * 2 * pi ! a random phase for the Fourier mode
            ! Assign a gaussian distributed amplitude and a random phase to the 
            ! Fourier mode:
            vhat(i, j, k) = (scalemodes**(-1)) * 1.0d-1 &
                          * sqrt(exp(-2.0d0 * kk / 3.0d2)) * exp(cmplx(0.0d0, phase))
            call random_number(phase)
            phase = phase * 2 * pi ! a random phase for the Fourier mode
            
            ! Assign a gaussian distributed amplitude and a random phase to the 
            ! Fourier mode:
            what(i, j, k) = (scalemodes**(-1)) * 1.0d-1 &
                          * sqrt(exp(-2.0d0 * kk / 3.0d2)) * exp(cmplx(0.0d0, phase)) 
        end do; end do; end do
        
        call rhsDealias()
        
        ! Apply projection operator to ensure divergence-free field:
        call rhsProject()    
        
        ! Set initial kinetic energy to 1.0
        call io_Ekinetic()
        
        factor = sqrt(1.0d0 / Ekin)
        
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
            
            uhat(i, j, k) = uhat(i, j, k) * factor
            vhat(i, j, k) = vhat(i, j, k) * factor
            what(i, j, k) = what(i, j, k) * factor
            
        end do; end do; end do
        
        ! Transform to configuration space:
        call p3dfft_btran_c2r (uhat, u, 'tff')
        call p3dfft_btran_c2r (vhat, v, 'tff')
        call p3dfft_btran_c2r (what, w, 'tff')
                
        ! Rescale modes after inv-Fourier transformation:
        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
            u(i, j, k) = u(i, j, k) * scalemodes
            v(i, j, k) = v(i, j, k) * scalemodes
            w(i, j, k) = w(i, j, k) * scalemodes
        end do; end do; end do            
        
    else
        ! Load initial state:
        call io_loadState('state0000.h5')
    end if

    call p3dfft_ftran_r2c (u, uhat, 'fft')
    call p3dfft_ftran_r2c (v, vhat, 'fft')
    call p3dfft_ftran_r2c (w, what, 'fft')

    call rhsDealias()

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

    call rhsIntFact() ! compute integration factor    
    
    open (99, file='RUNNING')
         write(99,*) 'This file indicates a job running in this directory.'
         write(99,*) 'Delete this file to cleanly terminate the process.'
    close(99)
    inquire(file='RUNNING', exist=running_exist)
    
    do while (running_exist)
    n = n + 1
        
        ! compute the nonlinear term for uhattemp:
        call rhsNonlinear()
        
        ! Predictor
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
            
            ! Predicted next step for corrector calculation:
            uhattemp(i, j, k) = intFact(i, j, k) &
                              * (uhatold(i, j, k) + dt * nonlinuhat(i, j, k))
            
            ! Contribution to the final step from the predictor calculation:
            uhat(i, j, k) = intFact(i, j, k) &
                          * (uhatold(i, j, k) &
                             + dt * nonlinuhat(i, j, k) * 0.5)
            
            ! Predicted next step for corrector calculation:
            vhattemp(i, j, k) = intFact(i, j, k) &
                              * (vhatold(i, j, k) + dt * nonlinvhat(i, j, k))
            
            ! Contribution to the final step from the predictor calculation:
            vhat(i, j, k) = intFact(i, j, k) &
                          * (vhatold(i, j, k) + dt * nonlinvhat(i, j, k) * 0.5)
                        
            ! Predicted next step for corrector calculation:
            whattemp(i, j, k) = intFact(i, j, k) &
                              * (whatold(i, j, k) + dt * nonlinwhat(i, j, k))
            
            ! Contribution to the final step from the predictor calculation:
            what(i, j, k) = intFact(i, j, k) &
                          * (whatold(i, j, k) + dt * nonlinwhat(i, j, k) * 0.5)
                        
        end do; end do; end do    
        
        ! compute the nonlinear term for uhattemp:
        call rhsNonlinear()
        
        ! Corrector
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
        
            uhat(i, j, k) = uhat(i, j, k) + dt * nonlinuhat(i, j, k) * 0.5
            vhat(i, j, k) = vhat(i, j, k) + dt * nonlinvhat(i, j, k) * 0.5
            what(i, j, k) = what(i, j, k) + dt * nonlinwhat(i, j, k) * 0.5
        
            uhattemp(i, j, k) = uhat(i, j, k)
            vhattemp(i, j, k) = vhat(i, j, k)
            whattemp(i, j, k) = what(i, j, k)
            
            uhatold(i, j, k) = uhat(i, j, k) 
            vhatold(i, j, k) = vhat(i, j, k) 
            whatold(i, j, k) = what(i, j, k) 
                
        end do; end do; end do    

        time(n+1) = n * dt
        
        if (proc_id.eq.0) then
            print *, 'time', n * dt
        end if
        
        ! Back to the configuration space:
        call p3dfft_btran_c2r (uhat, u, 'tff')
        call p3dfft_btran_c2r (vhat, v, 'tff')
        call p3dfft_btran_c2r (what, w, 'tff')
        
        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
            
            u(i, j, k) = u(i, j, k) * scalemodes
            v(i, j, k) = v(i, j, k) * scalemodes
            w(i, j, k) = w(i, j, k) * scalemodes
            
        end do; end do; end do
        
        if(modulo(n,iSaveRate1)==0) then
            call io_saveSpectrum()
            call io_saveState()
        end if
        
        if(modulo(n,iSaveRate2)==0) then
            call io_saveStats()
        end if
                
        ! Terminate if maximum time steps are reached or the Courant number 
        ! exceeds the allowed limit:
        if ((n .eq. Nt) .or. (CourantMax .gt. Courant)) then
            if (proc_id.eq.0) then
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

subroutine initMpi()
    
    !initialize mpi:
    call mpi_init (ierr)
    call mpi_comm_size (mpi_comm_world, nproc, ierr)
    call mpi_comm_rank (mpi_comm_world, proc_id, ierr)

end subroutine initMpi

end program nsbox
