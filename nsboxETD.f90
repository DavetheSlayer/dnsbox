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
    
    time(1) = 0.0d0
    factor = sqrt(3.0d0)

    ! initial condition from an analytical solution:
    do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
            
                u(i, j, k) = -0.5 * (factor * cos(x(i)) * sin(y(j)) * sin(z(k)) &
                                     + sin(x(i)) * cos(y(j)) * cos(z(k))) &
                                     * exp(- (factor ** 2) * time(1) * nu)                                     
                
                v(i,j,k) = 0.5 * (factor * sin(x(i)) * cos(y(j)) * sin(z(k)) &
						          - cos(x(i)) * sin(y(j)) * cos(z(k))) &
                                  * exp(-(factor**2)*time(1) * nu)
                
                w(i,j,k) = 1.0 * cos(x(i)) * cos(y(j)) * sin(z(k)) & 
                           * exp(-(factor**2)*time(1) * nu)                 
                
    end do; end do; end do

    ! randomly fill fields:
    do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
            
                call random_number(u(i, j, k))
                u(i, j, k) = u(i, j, k) - 0.5
                myuSum = myuSum + u(i, j, k)
                
                call random_number(u(i, j, k))
                v(i, j, k) = v(i, j, k) - 0.5
                myvSum = myvSum + u(i, j, k)

                call random_number(u(i, j, k))
                w(i, j, k) = w(i, j, k) - 0.5
                mywSum = mywSum + u(i, j, k)
                
    end do; end do; end do
    
    call mpi_allreduce(myuSum, uSum, 1, mpi_double_precision, &
                       mpi_sum, mpi_comm_world, ierr)
    
    call mpi_allreduce(myvSum, vSum, 1, mpi_double_precision, &
                       mpi_sum, mpi_comm_world, ierr)
    
    call mpi_allreduce(mywSum, wSum, 1, mpi_double_precision, &
                       mpi_sum, mpi_comm_world, ierr)
    
    ! ensure zero-mean
    do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)

                u(i, j, k) = (u(i, j, k) - uSum * scalemodes) * 0.1d0
                
                v(i, j, k) = (v(i, j, k) - vSum * scalemodes) * 0.1d0

                w(i, j, k) = (w(i, j, k) - wSum * scalemodes) * 0.1d0
                
    end do; end do; end do
    
    ! Apply projection operator to ensure divergence-free field:
    call rhsProject()    
        
!    call io_loadState('state0000.h5')
    
    call p3dfft_ftran_r2c (u, uhat, 'fft')
    call p3dfft_ftran_r2c (v, vhat, 'fft')
    call p3dfft_ftran_r2c (w, what, 'fft')
    
    if(proc_id .eq. 0) then 
        print *, 'Starting time-stepping'
    endif
    
!    call io_saveState()
    
    ! copy uhat -> uhatold and uhattemp
    do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
        
        uhattemp(i, j, k) = uhat(i, j, k) * scalemodes
        vhattemp(i, j, k) = vhat(i, j, k) * scalemodes
        whattemp(i, j, k) = what(i, j, k) * scalemodes
        
        uhatold(i, j, k) = uhat(i, j, k) * scalemodes
        vhatold(i, j, k) = vhat(i, j, k) * scalemodes
        whatold(i, j, k) = what(i, j, k) * scalemodes
        
        uhat(i, j, k) = uhatold(i, j, k)
        vhat(i, j, k) = vhatold(i, j, k)
        what(i, j, k) = whatold(i, j, k)
        
    end do; end do; end do
    
!    if(proc_id .eq. 0) then 
!        print *, 'check 1'
!    endif    
    
    call rhsIntFact() ! compute integration factor    
    do n = 1, Nt
        
        ! compute the nonlinear term for uhattemp:
        call rhsNonlinear()

!    if(proc_id .eq. 0) then 
!        print *, 'check 2'
!    endif    
        
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

!    if(proc_id .eq. 0) then 
!        print *, 'check 3'
!    endif    
        
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
            call io_saveState()
        end if
        
        if(modulo(n,iSaveRate2)==0) then
            call io_saveStats()
        end if
            
    end do
    
    ! Error in final numerical solution:
    
    do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
                
        utemp(i, j, k) = u(i, j, k) &
                        - (-0.5*( factor*cos(x(i))*sin(y(j))*sin(z(k))&
						+sin(x(i))*cos(y(j))*cos(z(k)) )*exp(-(factor**2)*time(Nt+1) * nu))
                
        vtemp(i, j, k) = v(i, j, k) &
                        - (0.5*(  factor*sin(x(i))*cos(y(j))*sin(z(k))&
						-cos(x(i))*sin(y(j))*cos(z(k)) )*exp(-(factor**2)*time(Nt+1) * nu))
                
        wtemp(i, j, k) = w(i, j, k) &
                        - (cos(x(i))*cos(y(j))*sin(z(k))*exp(-(factor**2)*time(Nt+1) * nu))
                            
    end do; end do; end do        
    
    mychg(1) = maxval(abs(utemp))
    mychg(2) = maxval(abs(vtemp))
    mychg(3) = maxval(abs(wtemp))
    
    call mpi_allreduce(mychg, allchg, 3, MPI_DOUBLE_PRECISION, &
                       MPI_MAX, MPI_COMM_WORLD, ierr)
    chg = allchg(1) + allchg(2) + allchg(3)
    
    if (proc_id.eq.0) then
        print *, 'The error at the final timestep is ', chg
    end if
    
	deallocate(x,y,z,time,mychg,allchg,u,v,w,ux,uy,uz,vx,vy,vz,wx,wy,wz,&
               !uold,vold,wold,&
               utemp,vtemp,wtemp,&
               temp_r,kx,ky,kz,uhat,vhat,what,&
               uhattemp,vhattemp,whattemp,&
               uhatold,vhatold,whatold,&
               phat,nonlinuhat,nonlinvhat,nonlinwhat,temp_c,&
               intFact, stat=AllocateStatus)		
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
