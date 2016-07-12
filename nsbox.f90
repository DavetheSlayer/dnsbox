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
                                     * exp(- (factor ** 2) * time(1) / Re)                                     
                
                v(i,j,k) = 0.5 * (factor * sin(x(i)) * cos(y(j)) * sin(z(k)) &
						          - cos(x(i)) * sin(y(j)) * cos(z(k))) &
                                  * exp(-(factor**2)*time(1)/Re)
                
                w(i,j,k) = 1.0 * cos(x(i)) * cos(y(j)) * sin(z(k)) & 
                           * exp(-(factor**2)*time(1)/Re)                 
                
    end do; end do; end do

    ! randomly fill fields:
!    do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
            
!                call random_number(u(i, j, k))
!                u(i, j, k) = u(i, j, k) - 0.5
!                myuSum = myuSum + u(i, j, k)
                
!                call random_number(u(i, j, k))
!                v(i, j, k) = v(i, j, k) - 0.5
!                myvSum = myvSum + u(i, j, k)

!                call random_number(u(i, j, k))
!                w(i, j, k) = w(i, j, k) - 0.5
!                mywSum = mywSum + u(i, j, k)
                
!    end do; end do; end do
    
!    call mpi_allreduce(myuSum, uSum, 1, mpi_double_precision, &
!                       mpi_sum, mpi_comm_world, ierr)
    
!    call mpi_allreduce(myvSum, vSum, 1, mpi_double_precision, &
!                       mpi_sum, mpi_comm_world, ierr)
    
!    call mpi_allreduce(mywSum, wSum, 1, mpi_double_precision, &
!                       mpi_sum, mpi_comm_world, ierr)
    
!    ! ensure zero-mean
!    do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)

!                u(i, j, k) = (u(i, j, k) - uSum * scalemodes) * 0.1d0
                
!                v(i, j, k) = (v(i, j, k) - vSum * scalemodes) * 0.1d0

!                w(i, j, k) = (w(i, j, k) - wSum * scalemodes) * 0.1d0
                
!    end do; end do; end do
    
!    ! Apply projection operator to ensure divergence-free field:
!    call rhsProject()    
    
!    call io_loadState('state0000.h5')
    
    call p3dfft_ftran_r2c (u, uhat, 'fft')
    call p3dfft_ftran_r2c (v, vhat, 'fft')
    call p3dfft_ftran_r2c (w, what, 'fft')

    
    ! test output:
    if (nproc .eq. 1) then
        ! print *, kz(:)
        print *, fstart(1), fend(1), fstart(2), fend(2), fstart(3), fend(3)
        print *, ux(1, 1, :)
        
    end if
    
    ! Derivative of u with respect to x, y, z:
    do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
        temp_c(i, j, k) = uhat(i, j, k) * kx(k) * scalemodes
    end do; end do; end do
    call p3dfft_btran_c2r (temp_c, ux, 'tff')    
    do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
        temp_c(i, j, k) = uhat(i, j, k) * ky(j) * scalemodes
    end do; end do; end do
    call p3dfft_btran_c2r (temp_c, uy, 'tff')
    do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
       temp_c(i, j, k) = uhat(i, j, k) * kz(i) * scalemodes
    end do; end do; end do
    call p3dfft_btran_c2r (temp_c, uz, 'tff')
    
    ! Derivative of v with respect to x, y, z:
    do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
        temp_c(i, j, k) = vhat(i, j, k) * kx(k) * scalemodes
    end do; end do; end do
    call p3dfft_btran_c2r (temp_c, vx, 'tff')
    do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
        temp_c(i, j, k) = vhat(i, j, k) * ky(j) * scalemodes
    end do; end do; end do
    call p3dfft_btran_c2r (temp_c, vy, 'tff')
    do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
        temp_c(i, j, k) = vhat(i, j, k) * kz(i) * scalemodes
    end do; end do; end do
    call p3dfft_btran_c2r (temp_c, vz, 'tff')
    
    ! Derivative of w with respect to x, y, z:
    do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
        temp_c(i, j, k) = what(i, j, k) * kx(k) * scalemodes
    end do; end do; end do
    call p3dfft_btran_c2r (temp_c, wx, 'tff')
    do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
        temp_c(i, j, k) = what(i, j, k) * ky(j) * scalemodes
    end do; end do; end do
    call p3dfft_btran_c2r (temp_c, wy, 'tff')
    do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
        temp_c(i, j, k) = what(i, j, k) * kz(i) * scalemodes
    end do; end do; end do
    call p3dfft_btran_c2r (temp_c, wz, 'tff')

    ! print divergence:
    do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)

                print *, ux(i,j,k) + vy(i,j,k) + wz(i,j,k)
                
    end do; end do; end do

    ! omegax:
    do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
                
        realtemp(i, j, k) = wy(i, j, k) - vz(i, j, k)
                            
    end do; end do; end do
        
    ! some save routine    
    ! omegay:
    do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
                
        realtemp(i, j, k) = uz(i, j, k) - wx(i, j, k)
                            
    end do; end do; end do 

    ! some save routine
    ! omegaz:
    do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
                
        realtemp(i, j, k) = vx(i, j, k) - uy(i, j, k)
                            
    end do; end do; end do  
    ! some save routine

    if(proc_id .eq. 0) then 
        print *, 'Starting time-stepping'
    endif
    
    call io_saveState()
    do n = 1, Nt
    
        ! fixed point
        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
                    
            uold(i, j, k) = u(i, j, k)
            uxold(i, j, k) = ux(i, j, k)
            uyold(i, j, k) = uy(i, j, k)
            uzold(i, j, k) = uz(i, j, k)
            
            vold(i, j, k) = v(i, j, k)
            vxold(i, j, k) = vx(i, j, k)
            vyold(i, j, k) = vy(i, j, k)
            vzold(i, j, k) = vz(i, j, k)
            
            wold(i, j, k) = w(i, j, k)
            wxold(i, j, k) = wx(i, j, k)
            wyold(i, j, k) = wy(i, j, k)
            wzold(i, j, k) = wz(i, j, k)
                    
        end do; end do; end do            
        
        call rhsfix()
 
        chg = 1
        ! starting fixed-point iteration:
        do while (chg .gt. tol)
            
            call rhsRest()
            
            do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
                
                uhat(i, j, k) = (rhsuhatfix(i, j, k) - nonlinuhat(i, j, k) &
                                 - kx(k) * phat(i, j, k)) &
                                / (dtInv - Q - (0.5d0 * ReInv) &
                                  * (kx(k) * kx(k) &
                                     + ky(j) * ky(j) &
                                     + kz(i) * kz(i)))
                
                vhat(i, j, k) = (rhsvhatfix(i, j, k) - nonlinvhat(i, j, k) &
                                 - ky(j) * phat(i, j, k)) &
                                / (dtInv - Q - (0.5d0 * ReInv) &
                                  * (kx(k) * kx(k) &
                                     + ky(j) * ky(j) &
                                     + kz(i) * kz(i)))
                                           
                what(i, j, k) = (rhswhatfix(i, j, k) - nonlinwhat(i, j, k) &
                                 - kz(i) * phat(i, j, k)) &
                                / (dtInv - Q - (0.5d0 * ReInv) &
                                  * (kx(k) * kx(k) &
                                     + ky(j) * ky(j) &
                                     + kz(i) * kz(i)))
                                           
            end do; end do; end do     
            
            ! Derivative of u with respect to x, y, z:
            do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
                temp_c(i, j, k) = uhat(i, j, k) * kx(k) * scalemodes
            end do; end do; end do
            call p3dfft_btran_c2r (temp_c, ux, 'tff')
            do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
                temp_c(i, j, k) = uhat(i, j, k) * ky(j) * scalemodes
            end do; end do; end do
            call p3dfft_btran_c2r (temp_c, uy, 'tff')
            do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
               temp_c(i, j, k) = uhat(i, j, k) * kz(i) * scalemodes
            end do; end do; end do
            call p3dfft_btran_c2r (temp_c, uz, 'tff')
            
            ! Derivative of v with respect to x, y, z:
            do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
                temp_c(i, j, k) = vhat(i, j, k) * kx(k) * scalemodes
            end do; end do; end do
            call p3dfft_btran_c2r (temp_c, vx, 'tff')
            do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
                temp_c(i, j, k) = vhat(i, j, k) * ky(j) * scalemodes
            end do; end do; end do
            call p3dfft_btran_c2r (temp_c, vy, 'tff')
            do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
                temp_c(i, j, k) = vhat(i, j, k) * kz(i) * scalemodes
            end do; end do; end do
            call p3dfft_btran_c2r (temp_c, vz, 'tff')
            
            ! Derivative of w with respect to x, y, z:
            do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
                temp_c(i, j, k) = what(i, j, k) * kx(k) * scalemodes
            end do; end do; end do
            call p3dfft_btran_c2r (temp_c, wx, 'tff')
            do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
                temp_c(i, j, k) = what(i, j, k) * ky(j) * scalemodes
            end do; end do; end do
            call p3dfft_btran_c2r (temp_c, wy, 'tff')
            do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
                temp_c(i, j, k) = what(i, j, k) * kz(i) * scalemodes
            end do; end do; end do
            call p3dfft_btran_c2r (temp_c, wz, 'tff')
            
            do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
                
                utemp(i, j, k) = u(i, j, k)
                vtemp(i, j, k) = v(i, j, k)
                wtemp(i, j, k) = w(i, j, k)
                
            end do; end do; end do
            
            call p3dfft_btran_c2r (uhat, u, 'tff')
            call p3dfft_btran_c2r (vhat, v, 'tff')
            call p3dfft_btran_c2r (what, w, 'tff')
            
            do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
                
                u(i, j, k) = u(i, j, k) * scalemodes
                v(i, j, k) = v(i, j, k) * scalemodes
                w(i, j, k) = w(i, j, k) * scalemodes
                
            end do; end do; end do
            
            mychg(1) = maxval(abs(utemp - u))
            mychg(2) = maxval(abs(vtemp - v))
            mychg(3) = maxval(abs(wtemp - w))
            
            call mpi_allreduce (mychg, allchg, 3, mpi_double_precision, &
                                mpi_max, mpi_comm_world, ierr)
            
            chg = allchg(1) + allchg(2) + allchg(3)
            if (proc_id.eq.0) then
                print *, 'chg:', chg
            end if
        end do
        
        time(n+1) = n * dt
        
        if (proc_id.eq.0) then
            print *, 'time', n * dt
        end if
        
        if(modulo(n,iSaveRate1)==0) then
            call io_saveState()
        end if
        
        if(modulo(n,iSaveRate2)==0) then
            call io_saveStats()
        end if
        
    ! omegax:
    do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
                
        realtemp(i, j, k) = wy(i, j, k) - vz(i, j, k)
                            
    end do; end do; end do
    ! some save routine    
    ! omegay:
    do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
                
        realtemp(i, j, k) = uz(i, j, k) - wx(i, j, k)
                            
    end do; end do; end do 
    ! some save routine
    ! omegaz:
    do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
                
        realtemp(i, j, k) = vx(i, j, k) - uy(i, j, k)
                            
    end do; end do; end do  
    ! some save routine
        
    end do
    
    ! Error in final numerical solution:
    
    do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
                
        utemp(i, j, k) = u(i, j, k) &
                        - (-0.5*( factor*cos(x(i))*sin(y(j))*sin(z(k))&
						+sin(x(i))*cos(y(j))*cos(z(k)) )*exp(-(factor**2)*time(Nt+1)/Re))
                
        vtemp(i, j, k) = v(i, j, k) &
                        - (0.5*(  factor*sin(x(i))*cos(y(j))*sin(z(k))&
						-cos(x(i))*sin(y(j))*cos(z(k)) )*exp(-(factor**2)*time(Nt+1)/Re))
                
        wtemp(i, j, k) = w(i, j, k) &
                        - (cos(x(i))*cos(y(j))*sin(z(k))*exp(-(factor**2)*time(Nt+1)/Re))
                            
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
    
	deallocate(x,y,z,time,mychg,allchg,u,v,w,ux,uy,uz,vx,vy,vz,wx,wy,wz,uold,uxold,uyold,uzold,&
               vold,vxold,vyold,vzold,wold,wxold,wyold,wzold,utemp,vtemp,wtemp,&
               temp_r,kx,ky,kz,uhat,vhat,what,rhsuhatfix,rhsvhatfix,&
               rhswhatfix,phat,nonlinuhat,nonlinvhat,nonlinwhat,temp_c,&
               realtemp,stat=AllocateStatus)		
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
