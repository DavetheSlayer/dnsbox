module state
    
    use variables
    implicit none
    
    contains
    
    subroutine stateProject()
        ! Apply projection on (u,v,w)hat to make it divergence-free:
        
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
            
            uhattemp(i, j, k) = uhat(i, j, k)
            vhattemp(i, j, k) = vhat(i, j, k)
            whattemp(i, j, k) = what(i, j, k)
            
        end do; end do; end do
         
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1) 
            ! k -> kx - index, j -> ky - index, i -> kz - index
            if((kx(k) .eq. cmplx(0.0d0, 0.0d0)) .and. &
               (ky(j) .eq. cmplx(0.0d0, 0.0d0)) .and. &
               (kz(i) .eq. cmplx(0.0d0, 0.0d0))) then
                cycle
            else
            
                uhat(i, j, k) = uhattemp(i, j, k) &
                              - kx(k) * kx(k) * uhattemp(i, j, k) / &
                                (kx(k) * kx(k) &
                               + ky(j) * ky(j) & 
                               + kz(i) * kz(i)) &
                              - kx(k) * ky(j) * vhattemp(i, j, k) / &
                                (kx(k) * kx(k) &
                               + ky(j) * ky(j) & 
                               + kz(i) * kz(i)) &
                              - kx(k) * kz(i) * whattemp(i, j, k) / &
                                (kx(k) * kx(k) &
                               + ky(j) * ky(j) & 
                               + kz(i) * kz(i)) 
                
                vhat(i, j, k) = vhattemp(i, j, k) &
                              - ky(j) * kx(k) * uhattemp(i, j, k) / &
                                (kx(k) * kx(k) &
                               + ky(j) * ky(j) & 
                               + kz(i) * kz(i)) &
                              - ky(j) * ky(j) * vhattemp(i, j, k) / &
                                (kx(k) * kx(k) &
                               + ky(j) * ky(j) & 
                               + kz(i) * kz(i)) &
                              - ky(j) * kz(i) * whattemp(i, j, k) / &
                                (kx(k) * kx(k) &
                               + ky(j) * ky(j) & 
                               + kz(i) * kz(i)) 
                
                what(i, j, k) = whattemp(i, j, k) &
                              - kz(i) * kx(k) * uhattemp(i, j, k) / &
                                (kx(k) * kx(k) &
                               + ky(j) * ky(j) & 
                               + kz(i) * kz(i) + eps) &
                              - kz(i) * ky(j) * vhattemp(i, j, k) / &
                                (kx(k) * kx(k) &
                               + ky(j) * ky(j) & 
                               + kz(i) * kz(i) + eps) &
                              - kz(i) * kz(i) * whattemp(i, j, k) / &
                                (kx(k) * kx(k) &
                               + ky(j) * ky(j) & 
                               + kz(i) * kz(i)) 
                
            end if
                                       
        end do; end do; end do
        
    end subroutine stateProject

    subroutine stateDealias()
    ! Set k > 2/3 k_max elements to 0
            
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1) 
            
            ! Truncation:
            if((spherical .and. &
               (abs(kx(k)) ** 2.0d0 + abs(ky(j)) ** 2.0d0 + abs(kz(i)) ** 2.0d0 & 
                .ge.  ((real(Nx, kind=8) / 2.0d0) &
                     * (2.0d0 * alpha_x / 3.0d0)) ** 2)) .or. &
                (abs(kx(k)) .ge. (real(Nx, kind=8) / 2.0d0) &
                               * (2.0d0 * alpha_x / 3.0d0)  &
            .or. abs(ky(j)) .ge. (real(Ny, kind=8) / 2.0d0) &
                               * (2.0d0 * alpha_y / 3.0d0)  &
            .or. abs(kz(i)) .ge. (real(Nz, kind=8) / 2.0d0) &
                               * (2.0d0 * alpha_z / 3.0d0))) then

                  uhat(i, j, k) = cmplx(0.0d0, 0.0d0)
                  vhat(i, j, k) = cmplx(0.0d0, 0.0d0)
                  what(i, j, k) = cmplx(0.0d0, 0.0d0)          

            end if
            
        end do; end do; end do
        
    end subroutine stateDealias
    
    subroutine stateInitRand()
        ! Initiate fields in momentum space:
        ! Initiation method: Rosales & Meneveau (2005), Phys. Fluids. eq.9
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)    
            ! skip k_i = 0
!            if((kx(k) .eq. cmplx(0.0d0, 0.0d0)) .or. &
!               (ky(j) .eq. cmplx(0.0d0, 0.0d0)) .or. &
!               (kz(i) .eq. cmplx(0.0d0, 0.0d0))) cycle 
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
            ! Fourier mode (eq9 / 4pi k^2):
            uhat(i, j, k) = sqrt(4 * (uzero ** 2) * kk * (kzero ** (-5))  &
                                 * (pi ** (-3.0d0/2.0d0)) * (2.0d0 ** (-0.5)) &
                                 * exp(-2.0d0 * kk / (kzero ** 2))) &
                                 * exp(cmplx(0.0d0, phase)) * scalemodes ** (-1)
            call random_number(phase)
            phase = phase * 2 * pi ! a random phase for the Fourier mode
            ! Assign a gaussian distributed amplitude and a random phase to the 
            ! Fourier mode:
            vhat(i, j, k) = sqrt(4 * (uzero ** 2) * kk * (kzero ** (-5))  &
                                 * (pi ** (-3.0d0/2.0d0)) * (2.0d0 ** (-0.5)) &
                                 * exp(-2.0d0 * kk / (kzero ** 2))) &
                                 * exp(cmplx(0.0d0, phase)) * scalemodes ** (-1)
            call random_number(phase)
            phase = phase * 2 * pi ! a random phase for the Fourier mode
            
            ! Assign a gaussian distributed amplitude and a random phase to the 
            ! Fourier mode:
            what(i, j, k) = sqrt(4 * (uzero ** 2) * kk * (kzero ** (-5))  &
                                 * (pi ** (-3.0d0/2.0d0)) * (2.0d0 ** (-0.5)) &
                                 * exp(-2.0d0 * kk / (kzero ** 2))) &
                                 * exp(cmplx(0.0d0, phase)) * scalemodes ** (-1)
        end do; end do; end do
        
        call stateDealias()
        
        ! Apply projection operator to ensure divergence-free field:
        call stateProject()    
               
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
    
    end subroutine stateInitRand
    
    subroutine stateInitRandConf()
        ! Initiate a random field in the configuration space, then transform 
        ! to Fourier space and act on with the projection operator to make 
        ! divergence free. Finally adjust the spectrum to a predefined one.
        ! for description of the procedure see Yoffe 2012, s. 2.2.5
        
        ! constants:                                S7
        real(kind=8) :: Cone =    8.0d-5            ! 8.0d-2
        real(kind=8) :: Ctwo =    2.0d0             ! 2.0d0
        real(kind=8) :: Cthree =  1d-3 !8.2352309d-2      ! 8.2352309d-2
        real(kind=8) :: Cfour =   2.0d0             ! 2.0d0
        
        if(proc_id .eq. 0) then 
            print *, 'Yoffe initial condition'
            print *, 'C_i = ', Cone, Ctwo, Cthree, Cfour
        endif
        
        ! configuration space loop:
        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
                    
                    call random_number(u(i, j, k))  
                    u(i, j, k) = u(i, j, k) - 0.5d0 
                    call random_number(v(i, j, k))  
                    v(i, j, k) = v(i, j, k) - 0.5d0 
                    call random_number(u(i, j, k))  
                    w(i, j, k) = w(i, j, k) - 0.5d0 
                    
        end do; end do; end do                
        
        call p3dfft_ftran_r2c (u, uhat, 'fft')
        call p3dfft_ftran_r2c (v, vhat, 'fft')
        call p3dfft_ftran_r2c (w, what, 'fft')
        
        call stateDealias() ! Dealias
        call stateProject() ! Make state solenoidal
        call stateComputeSpectrum() ! Compute the spectrum
        
        
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
            
            ! Find which window current k belongs to:
            absk = sqrt(real(conjg(kx(k)) * kx(k) &
                           + conjg(ky(j)) * ky(j) &
                           + conjg(kz(i)) * kz(i))) ! |k|
            if (absk < real(int(absk / Deltak)) + 0.5 * Deltak) then
                nk = int(absk / Deltak)
            else
                nk = int(absk / Deltak) + 1
            end if 
            
            if (nk > Nspec .or. nk .eq. 0) then
                ! If |k| = 0 or |k|>k_max:
                uhat(i, j, k) = 0.0d0
                vhat(i, j, k) = 0.0d0
                what(i, j, k) = 0.0d0
            else 
                
                uhat(i, j, k) = uhat(i, j, k) &
                              * sqrt(Cone * (nk * Deltak) ** Ctwo &
                                   * exp(- Cthree * (nk * Deltak) ** Cfour) &
                                     / Espec(nk))
                vhat(i, j, k) = vhat(i, j, k) &
                              * sqrt(Cone * (nk * Deltak) ** Ctwo &
                                   * exp(- Cthree * (nk * Deltak) ** Cfour) &
                                     / Espec(nk))
                what(i, j, k) = what(i, j, k) &
                              * sqrt(Cone * (nk * Deltak) ** Ctwo &
                                   * exp(- Cthree * (nk * Deltak) ** Cfour) &
                                     / Espec(nk))
                
            end if
            
            
            
        end do; end do; end do                       
        
        ! Back to the configuration space:
        call p3dfft_btran_c2r (uhat, u, 'tff')
        call p3dfft_btran_c2r (vhat, v, 'tff')
        call p3dfft_btran_c2r (what, w, 'tff')
        
        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
            
            u(i, j, k) = u(i, j, k) * scalemodes
            v(i, j, k) = v(i, j, k) * scalemodes
            w(i, j, k) = w(i, j, k) * scalemodes
            
        end do; end do; end do
    
    end subroutine stateInitRandConf
    
    subroutine stateInitAnalytic()
        ! initial condition from an analytical solution:
    
        factor = sqrt(3.0d0)

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
        
    end subroutine stateInitAnalytic
    
    subroutine stateInitEig()
        ! initiate state as the eigenvector
        factor = 1.0d-6
        
        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
                    
                    u(i, j, k) = factor * cos(x(i) + y(j))                                     
                    
                    v(i, j, k) = -1.0d0 * factor * cos(x(i) + y(j))                                     
                                        
                    w(i,j,k) = 0.0d0                 
                    
        end do; end do; end do        
        
        
    end subroutine stateInitEig
    
    subroutine stateCheckError()
        
        ! Back to the configuration space:
        call p3dfft_btran_c2r (uhat, u, 'tff')
        call p3dfft_btran_c2r (vhat, v, 'tff')
        call p3dfft_btran_c2r (what, w, 'tff')
        
        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
            
            u(i, j, k) = u(i, j, k) * scalemodes
            v(i, j, k) = v(i, j, k) * scalemodes
            w(i, j, k) = w(i, j, k) * scalemodes
            
        end do; end do; end do
        
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
        
    end subroutine stateCheckError
    
    subroutine stateComputeSpectrum()
        ! Save window-averaged spectrum
        myEspec = 0d0
        myEZero = 0d0
        myEOne  = 0d0
        myEsqrt2 = 0d0
        
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
            
            ! Find which window current k belongs to:
            absk = sqrt(real(conjg(kx(k)) * kx(k) &
                           + conjg(ky(j)) * ky(j) &
                           + conjg(kz(i)) * kz(i))) ! |k|^2
            if (absk < real(int(absk / Deltak)) + 0.5 * Deltak) then
                nk = int(absk / Deltak)
            else
                nk = int(absk / Deltak) + 1
            end if 
            
            if (absk .eq. 0) then
                myEZero = myEZero &
                        + real(conjg(uhat(i, j, k)) * uhat(i, j, k)) &
                        + real(conjg(vhat(i, j, k)) * vhat(i, j, k)) &
                        + real(conjg(what(i, j, k)) * what(i, j, k))
            end if

            if (absk .eq. 1.0) then
                myEOne = myEOne &
                        + real(conjg(uhat(i, j, k)) * uhat(i, j, k)) &
                        + real(conjg(vhat(i, j, k)) * vhat(i, j, k)) &
                        + real(conjg(what(i, j, k)) * what(i, j, k))
            end if

            if (absk - sqrt(2.0d0) .lt. 1.0d-9) then
                myEsqrt2 = myEsqrt2 &
                        + real(conjg(uhat(i, j, k)) * uhat(i, j, k)) &
                        + real(conjg(vhat(i, j, k)) * vhat(i, j, k)) &
                        + real(conjg(what(i, j, k)) * what(i, j, k))
            end if
            
            if (nk > Nspec .or. nk .eq. 0) cycle
            
            myEspec(nk) = myEspec(nk) &
                        + (real(conjg(uhat(i, j, k)) * uhat(i, j, k)) &
                         + real(conjg(vhat(i, j, k)) * vhat(i, j, k)) &
                         + real(conjg(what(i, j, k)) * what(i, j, k))) &
                          * (scalemodes ** 2) / Deltak
        end do; end do; end do               
        
        call mpi_allreduce(myEZero, EZero, 1, mpi_double_precision, &
                           mpi_sum, mpi_comm_world, ierr)
        
        call mpi_allreduce(myEOne, EOne, 1, mpi_double_precision, &
                           mpi_sum, mpi_comm_world, ierr)
        
        call mpi_allreduce(myEsqrt2, Esqrt2, 1, mpi_double_precision, &
                           mpi_sum, mpi_comm_world, ierr)
                           
        call mpi_allreduce(myEspec, Espec, Nspec, mpi_double_precision, &
                           mpi_sum, mpi_comm_world, ierr)
    
    end subroutine stateComputeSpectrum

    subroutine stateDivergence()

        ! Derivatives are already computed in the time-loop
        
        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
            
            temp_r(i, j, k) = ux(i, j, k) + vy(i, j, k) + wz(i, j, k)
                                
        end do; end do; end do                
        
        mydivMax = maxval(temp_r)
        call mpi_allreduce(mydivMax, divMax, 1, MPI_DOUBLE_PRECISION, &
                           MPI_MAX, MPI_COMM_WORLD, ierr)
        
    end subroutine stateDivergence
    
end module state
