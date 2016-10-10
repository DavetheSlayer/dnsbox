module state
    
    use variables
    implicit none
    
    real(kind=8) :: scalekx   !scaling to avoid counting kx=0 modes twice
    
    contains
    
    subroutine state_u2utemp()

        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
                    
            utemp(i, j, k) = u(i, j, k)
            vtemp(i, j, k) = v(i, j, k)
            wtemp(i, j, k) = w(i, j, k)
                    
        end do; end do; end do                        
        
    end subroutine state_u2utemp
    
    subroutine state_utemp2u()

        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
                    
            u(i, j, k) = utemp(i, j, k)
            v(i, j, k) = vtemp(i, j, k)
            w(i, j, k) = wtemp(i, j, k)
                    
        end do; end do; end do                        
        
    end subroutine state_utemp2u
    
    subroutine state_uhat2uhattemp()
        
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
            
            uhattemp(i, j, k) = uhat(i, j, k)
            vhattemp(i, j, k) = vhat(i, j, k)
            whattemp(i, j, k) = what(i, j, k)
            
        end do; end do; end do
        
    end subroutine state_uhat2uhattemp
    
    subroutine state_uhattemp2uhat()
        
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
            
            uhat(i, j, k) = uhattemp(i, j, k)
            vhat(i, j, k) = vhattemp(i, j, k)
            what(i, j, k) = whattemp(i, j, k)
            
        end do; end do; end do
        
    end subroutine state_uhattemp2uhat
    
    subroutine state_utemp2uhattemp()
        
        call p3dfft_ftran_r2c (utemp, uhattemp, 'fft')
        call p3dfft_ftran_r2c (vtemp, vhattemp, 'fft')
        call p3dfft_ftran_r2c (wtemp, whattemp, 'fft')
        
    end subroutine state_utemp2uhattemp
    
    subroutine state_uhattemp2utemp()
        
        call p3dfft_btran_c2r (uhattemp, utemp, 'tff')
        call p3dfft_btran_c2r (vhattemp, vtemp, 'tff')
        call p3dfft_btran_c2r (whattemp, wtemp, 'tff')
        
    end subroutine state_uhattemp2utemp
    
    subroutine state_u2uhat()
        
        call state_u2utemp()
        call state_utemp2uhattemp()
        call state_uhattemp2uhat()
        
    end subroutine state_u2uhat
    
    subroutine state_uhat2u()
        
        call state_uhat2uhattemp()
        call state_uhattemp2utemp()
        call state_utemp2u()
                
        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
            
            u(i, j, k) = u(i, j, k) * scalemodes
            v(i, j, k) = v(i, j, k) * scalemodes
            w(i, j, k) = w(i, j, k) * scalemodes
            
        end do; end do; end do
        
    end subroutine state_uhat2u

    subroutine state_copy_conf(a, b)
        ! Copy configuration space field a to b
        real(p3dfft_type), &
        dimension(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)):: a, b
        
        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
            b(i, j, k) = a(i, j, k)
        end do; end do; end do        
        
    end subroutine state_copy_conf

    subroutine state_copy_fourier(a, b)
        ! Copy fourier space field a to b
        complex(p3dfft_type), &
        dimension(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)):: a, b
        
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
            b(i, j, k) = a(i, j, k)
        end do; end do; end do        
        
    end subroutine state_copy_fourier
    
    subroutine state_project()
        ! Apply projection on (u,v,w)hat to make it divergence-free:
        
        call state_uhat2uhattemp()
        
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
        
    end subroutine state_project

    subroutine state_dealias()
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
                               * (2.0d0 * alpha_z / 3.0d0)) & !)then  
            .or. ((kx(k) .eq. cmplx(0.0d0, 0.0d0)) .and. &
                  (ky(j) .eq. cmplx(0.0d0, 0.0d0)) .and. &
                  (kz(i) .eq. cmplx(0.0d0, 0.0d0)))) then

                  uhat(i, j, k) = cmplx(0.0d0, 0.0d0)
                  vhat(i, j, k) = cmplx(0.0d0, 0.0d0)
                  what(i, j, k) = cmplx(0.0d0, 0.0d0)          

            end if
            
        end do; end do; end do
        
    end subroutine state_dealias
    
    subroutine state_init_rand()
        ! Initiate a random field in the configuration space, then transform 
        ! to Fourier space and act on with the projection operator to make 
        ! divergence free. Finally adjust the spectrum to a predefined one.
        ! for description of the procedure see Yoffe 2012, s. 2.2.5
        
        ! constants:                                S5             S7
        real(kind=8) :: Cone =  8.0d-2            ! 0.001702     ! 8.0d-2
        real(kind=8) :: Ctwo = 2.0d0              ! 4.0d0        ! 2.0d0
        real(kind=8) :: Cthree = 8.2352309d-2     ! 0.08         ! 8.2352309d-2
        real(kind=8) :: Cfour = 2.0d0             ! 2.0d0        ! 2.0d0
        
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
        
        call state_u2uhat()
        
        call state_dealias() ! Dealias
        call state_project() ! Make state solenoidal
        call state_spectrum() ! Compute the spectrum
                
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
            
            if (nk .ge. Nspec .or. nk .eq. 0) then
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
        call state_uhat2u()
    
    end subroutine state_init_rand
    
    subroutine state_init_analytic()
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
        
    end subroutine state_init_analytic
    
    subroutine state_check_error()
        
        ! Back to the configuration space:
        call state_uhat2u()
        factor = sqrt(3.0d0)
        
        if (proc_id.eq.0) then
            print *, 'n = ', n
            print *, 'Nt = ', Nt
        endif
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
        
    end subroutine state_check_error
    
    subroutine state_spectrum()
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
            
            if (nk .ge. Nspec .or. nk .eq. 0) cycle
            
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
    
    end subroutine state_spectrum
    
    subroutine state_divergence()

        ! Derivatives are already computed in the time-loop
        
        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
            
            temp_r(i, j, k) = ux(i, j, k) + vy(i, j, k) + wz(i, j, k)
                                
        end do; end do; end do                
        
        mydivMax = maxval(temp_r)
        call mpi_allreduce(mydivMax, divMax, 1, MPI_DOUBLE_PRECISION, &
                           MPI_MAX, MPI_COMM_WORLD, ierr)
        
    end subroutine state_divergence

    subroutine state_kinetic()
        !*************************!
        ! Average kinetic energy: !
        !*************************!
        myEkin = 0d0
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
        !   print *, k
            if (kx(k).eq.cmplx(0.0d0, 0.0d0)) then              
                scalekx = 0.5d0
            else 
                scalekx = 1.0d0 
            end if
            
            myEkin = myEkin &
                   + real(conjg(uhat(i, j, k)) * uhat(i, j, k)) * scalekx &
                   + real(conjg(vhat(i, j, k)) * vhat(i, j, k)) * scalekx &
                   + real(conjg(what(i, j, k)) * what(i, j, k)) * scalekx
            
            
        end do; end do; end do   
        
        myEkin = myEkin * (scalemodes ** 2)
        call mpi_allreduce(myEkin, Ekin, 1, mpi_double_precision, &
                           mpi_sum, mpi_comm_world, ierr)
                
    end subroutine state_kinetic
        
    subroutine state_derivatives()
        ! Compute space-derivatives of u,v,w from (u,v,w)hattemp
    
        ! Derivative of u with respect to x, y, z:
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
            temp_c(i, j, k) = uhattemp(i, j, k) * kx(k) * scalemodes
        end do; end do; end do
        call p3dfft_btran_c2r (temp_c, ux, 'tff')
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
            temp_c(i, j, k) = uhattemp(i, j, k) * ky(j) * scalemodes
        end do; end do; end do
        call p3dfft_btran_c2r (temp_c, uy, 'tff')
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
            temp_c(i, j, k) = uhattemp(i, j, k) * kz(i) * scalemodes
        end do; end do; end do
        call p3dfft_btran_c2r (temp_c, uz, 'tff')
        
        ! Derivative of v with respect to x, y, z:
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
            temp_c(i, j, k) = vhattemp(i, j, k) * kx(k) * scalemodes
        end do; end do; end do
        call p3dfft_btran_c2r (temp_c, vx, 'tff')
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
            temp_c(i, j, k) = vhattemp(i, j, k) * ky(j) * scalemodes
        end do; end do; end do
        call p3dfft_btran_c2r (temp_c, vy, 'tff')
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
            temp_c(i, j, k) = vhattemp(i, j, k) * kz(i) * scalemodes
        end do; end do; end do
        call p3dfft_btran_c2r (temp_c, vz, 'tff')        
        
        ! Derivative of w with respect to x, y, z:
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
            temp_c(i, j, k) = whattemp(i, j, k) * kx(k) * scalemodes
        end do; end do; end do
        call p3dfft_btran_c2r (temp_c, wx, 'tff')
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
            temp_c(i, j, k) = whattemp(i, j, k) * ky(j) * scalemodes
        end do; end do; end do
        call p3dfft_btran_c2r (temp_c, wy, 'tff')
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
            temp_c(i, j, k) = whattemp(i, j, k) * kz(i) * scalemodes
        end do; end do; end do
        call p3dfft_btran_c2r (temp_c, wz, 'tff')
    end subroutine state_derivatives
    
    subroutine state_vorticity()
        ! computes vorticity
        ! should be called after state_derivatives()
        
        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
            
            ! omegax:        
            omegax(i, j, k) = wy(i, j, k) - vz(i, j, k)
            omegay(i, j, k) = uz(i, j, k) - wx(i, j, k)
            omegaz(i, j, k) = vx(i, j, k) - uy(i, j, k)
                                
        end do; end do; end do        
        
    end subroutine state_vorticity

    subroutine state_dissipation()
        
        !*******************!
        ! Dissipation rate: !
        !*******************!        
        call state_derivatives() ! Compute derivatives for dissipation calculation
        
        myDisp = 0d0
        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
            
            myDisp = myDisp + ux(i, j, k) * ux(i, j, k) &
                            + uy(i, j, k) * uy(i, j, k) &
                            + uz(i, j, k) * uz(i, j, k) &
                            + vx(i, j, k) * vx(i, j, k) &
                            + vy(i, j, k) * vy(i, j, k) &
                            + vz(i, j, k) * vz(i, j, k) &
                            + wx(i, j, k) * wx(i, j, k) &
                            + wy(i, j, k) * wy(i, j, k) &
                            + wz(i, j, k) * wz(i, j, k) 
                            
        end do; end do; end do        
        
        myDisp = myDisp * nu * scalemodes

        call mpi_allreduce(myDisp, Disp, 1, mpi_double_precision, &
                           mpi_sum, mpi_comm_world, ierr)                             
        
    end subroutine state_dissipation
    
    
end module state
