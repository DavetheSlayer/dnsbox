module state
    
    use variables
    implicit none
    
    real(kind=8) :: scalekx   !scaling to avoid counting kx=0 modes twice
    
    contains
    
    subroutine state_utemp2uhattemp()
        
        ! Forward fft
        call dfftw_execute_(plan_forward_u)
        call dfftw_execute_(plan_forward_v)
        call dfftw_execute_(plan_forward_w)
    
    end subroutine state_utemp2uhattemp
    
    subroutine state_uhattemp2utemp()
        
        ! Forward fft
        call dfftw_execute_(plan_backward_u)
        call dfftw_execute_(plan_backward_v)
        call dfftw_execute_(plan_backward_w)
    
    end subroutine state_uhattemp2utemp
    
    subroutine state_u2utemp()
        
        do i=1,Nx; do j=1,Ny; do k=1,Nz
            utemp(i, j, k) = u(i, j, k)
            vtemp(i, j, k) = v(i, j, k)
            wtemp(i, j, k) = w(i, j, k)
        end do; end do; end do
        
    end subroutine state_u2utemp
    
    subroutine state_utemp2u()
        
        do i=1,Nx; do j=1,Ny; do k=1,Nz
            u(i, j, k) = utemp(i, j, k)
            v(i, j, k) = vtemp(i, j, k)
            w(i, j, k) = wtemp(i, j, k)
        end do; end do; end do
        
    end subroutine state_utemp2u
    
    subroutine state_uhat2uhattemp()
        
        do i=1,Nh; do j=1,Ny; do k=1,Nz
            uhattemp(i, j, k) = uhat(i, j, k)
            vhattemp(i, j, k) = vhat(i, j, k)
            whattemp(i, j, k) = what(i, j, k)
        end do; end do; end do
        
    end subroutine state_uhat2uhattemp
    
    subroutine state_uhattemp2uhat()
        
        do i=1,Nh; do j=1,Ny; do k=1,Nz
            uhat(i, j, k) = uhattemp(i, j, k)
            vhat(i, j, k) = vhattemp(i, j, k)
            what(i, j, k) = whattemp(i, j, k)
        end do; end do; end do
        
    end subroutine state_uhattemp2uhat
    
    subroutine state_u2uhat()
        
        call state_u2utemp()
        call state_utemp2uhattemp()
        call state_uhattemp2uhat()
        
    end subroutine state_u2uhat

    subroutine state_uhat2u()
        
        call state_uhat2uhattemp()
        call state_uhattemp2utemp()
        call state_utemp2u()
        
        do k=1,Nz; do j=1,Ny; do i=1,Nx
                    
                    u(i, j, k) = u(i, j, k) * scalemodes
                    v(i, j, k) = v(i, j, k) * scalemodes
                    w(i, j, k) = w(i, j, k) * scalemodes
                    
        end do; end do; end do               
        
    end subroutine state_uhat2u

    subroutine state_copy_conf(a, b)
        ! Copy configuration space field a to b
        real(kind=8), dimension(Nx, Ny, Nz):: a, b
        
        do i=1,Nx; do j=1,Ny; do k=1,Nz
            b(i, j, k) = a(i, j, k)
        end do; end do; end do        
        
    end subroutine state_copy_conf

    subroutine state_copy_fourier(a, b)
        ! Copy fourier space field a to b
        complex(kind=8), dimension(Nh, Ny, Nz):: a, b
        
        do i=1,Nh; do j=1,Ny; do k=1,Nz
            b(i, j, k) = a(i, j, k)
        end do; end do; end do        
        
    end subroutine state_copy_fourier
    
    subroutine state_project()
        ! Apply projection on (u,v,w)hat to make it divergence-free:
        
        call state_uhat2uhattemp()
         
        do i=1,Nh; do j=1,Ny; do k=1,Nz 

            if(kx(i) .eq. 0.0d0 .and. &
               ky(j) .eq. 0.0d0 .and. &
               kz(k) .eq. 0.0d0) then
               cycle
            else
                kk = kx(i) * kx(i) + ky(j) * ky(j) + kz(k) * kz(k)
                
                uhat(i, j, k) = uhattemp(i, j, k) &
                              - kx(k) * kx(k) * uhattemp(i, j, k) / kk &
                              - kx(k) * ky(j) * vhattemp(i, j, k) / kk &
                              - kx(k) * kz(i) * whattemp(i, j, k) / kk 
                
                vhat(i, j, k) = vhattemp(i, j, k) &
                              - ky(j) * kx(k) * uhattemp(i, j, k) / kk &
                              - ky(j) * ky(j) * vhattemp(i, j, k) / kk &
                              - ky(j) * kz(i) * whattemp(i, j, k) / kk 
                
                what(i, j, k) = whattemp(i, j, k) &
                              - kz(i) * kx(k) * uhattemp(i, j, k) / kk &
                              - kz(i) * ky(j) * vhattemp(i, j, k) / kk &
                              - kz(i) * kz(i) * whattemp(i, j, k) / kk 
                
            end if
                                       
        end do; end do; end do
        
    end subroutine state_project

    subroutine state_dealias()
    ! Set k > 2/3 k_max elements to 0
            
        do i=1,Nh; do j=1,Ny; do k=1,Nz 
            
            ! Truncation:
            if((spherical .and. &
               (kx(i) ** 2.0d0 + ky(j) ** 2.0d0 + kz(k) ** 2.0d0 & 
                .ge.  ((real(Nx, kind=8) / 2.0d0) &
                     * (2.0d0 * alpha_x / 3.0d0)) ** 2)) .or. &
                (kx(i) .ge. (real(Nx, kind=8) / 2.0d0) &
                          * (2.0d0 * alpha_x / 3.0d0)  &
            .or. ky(j) .ge. (real(Ny, kind=8) / 2.0d0) &
                          * (2.0d0 * alpha_y / 3.0d0)  &
            .or. kz(k) .ge. (real(Nz, kind=8) / 2.0d0) &
                          * (2.0d0 * alpha_z / 3.0d0)) & 
            .or. ((kx(i) .eq. 0.0d0) .and. &
                  (ky(j) .eq. 0.0d0) .and. &
                  (kz(k) .eq. 0.0d0))) then

                  uhat(i, j, k) = 0.0d0
                  vhat(i, j, k) = 0.0d0
                  what(i, j, k) = 0.0d0          

            end if
            
        end do; end do; end do
        
    end subroutine state_dealias

    subroutine state_spectrum()
        ! Compute window-averaged spectrum
        Espec = 0d0
        Ezero = 0d0
        
        do i=1,Nh; do j=1,Ny; do k=1,Nz
            
            ! Spectrum scale factor:
            if (kx(i).eq. 0.0d0) then              
                scalekx = 0.5d0
            else 
                scalekx = 1.0d0 
            end if
            
            ! Find which window current k belongs to:
            absk = sqrt(kx(i) * kx(i) + ky(j) * ky(j) + kz(k) * kz(k)) ! |k|
            
            if (absk < real(int(absk / Deltak)) + 0.5 * Deltak) then
                nk = int(absk / Deltak)
            else
                nk = int(absk / Deltak) + 1
            end if 
            
            if (absk .eq. 0) then
                Ezero = Ezero &
                      + real(conjg(uhat(i, j, k)) * uhat(i, j, k)) * scalekx &
                      + real(conjg(vhat(i, j, k)) * vhat(i, j, k)) * scalekx &
                      + real(conjg(what(i, j, k)) * what(i, j, k)) * scalekx 
            end if
            
            if (nk .gt. Nspec .or. nk .eq. 0) cycle
            
            Espec(nk) = Espec(nk) &
                      + (real(conjg(uhat(i, j, k)) * uhat(i, j, k)) &
                       + real(conjg(vhat(i, j, k)) * vhat(i, j, k)) &
                       + real(conjg(what(i, j, k)) * what(i, j, k))) &
                      * (scalemodes ** 2)  * scalekx / Deltak
        end do; end do; end do               
            
    end subroutine state_spectrum
        
    subroutine state_init_rand()
        ! Initiate a random field in the configuration space, then transform 
        ! to Fourier space and act on with the projection operator to make 
        ! divergence free. Finally adjust the spectrum to a predefined one.
        ! for description of the procedure see Yoffe 2012, s. 2.2.5
        
        ! constants:                                S7
        real(kind=8) :: Cone =    0.001702          ! 8.0d-2
        real(kind=8) :: Ctwo =    4.0d0             ! 2.0d0
        real(kind=8) :: Cthree =  0.08              ! 8.2352309d-2
        real(kind=8) :: Cfour =   2.0d0             ! 2.0d0
        
        print *, 'Yoffe initial condition'
        print *, 'C_i = ', Cone, Ctwo, Cthree, Cfour
        
        ! configuration space loop:
        do k=1,Nz; do j=1,Ny; do i=1,Nx
                    
                    call random_number(u(i, j, k))  
                    u(i, j, k) = u(i, j, k) - 0.5d0 
                    call random_number(v(i, j, k))  
                    v(i, j, k) = v(i, j, k) - 0.5d0 
                    call random_number(u(i, j, k))  
                    w(i, j, k) = w(i, j, k) - 0.5d0 
                    
        end do; end do; end do                
        
        call state_u2uhat()  ! Transform into Fourier space
        call state_dealias() ! Dealias
        call state_project()  ! Make state solenoidal
        call state_spectrum() ! Compute the spectrum
        
        do i=1,Nh; do j=1,Ny; do k=1,Nz
            
            ! Spectrum scale factor:
            if (kx(i).eq. 0.0d0) then              
                scalekx = 0.5d0
            else 
                scalekx = 1.0d0 
            end if
            
            ! Find which window current k belongs to:
            absk = sqrt(kx(k) * kx(k) + ky(j) * ky(j) + kz(i) * kz(i)) ! |k|
            
            if (absk < real(int(absk / Deltak)) + 0.5 * Deltak) then
                nk = int(absk / Deltak)
            else
                nk = int(absk / Deltak) + 1
            end if
            
            if (nk .gt. Nspec .or. nk .eq. 0) then
                ! If |k| = 0 or |k|>k_max:
                uhat(i, j, k) = 0.0d0
                vhat(i, j, k) = 0.0d0
                what(i, j, k) = 0.0d0
            else 
                
                uhat(i, j, k) = uhat(i, j, k) * scalekx &
                              * sqrt(Cone * (nk * Deltak) ** Ctwo &
                                   * exp(- Cthree * (nk * Deltak) ** Cfour) &
                                     / Espec(nk))
                vhat(i, j, k) = vhat(i, j, k) * scalekx  &
                              * sqrt(Cone * (nk * Deltak) ** Ctwo &
                                   * exp(- Cthree * (nk * Deltak) ** Cfour) &
                                     / Espec(nk))
                what(i, j, k) = what(i, j, k) * scalekx  &
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

        do k=1,Nz; do j=1,Ny; do i=1,Nx
                
            u(i,j,k) = -0.5 * (factor * cos(x(i)) * sin(y(j)) * sin(z(k)) &
                             + sin(x(i)) * cos(y(j)) * cos(z(k))) &
                             * exp(- (factor ** 2) * time(1) * nu)                                     
            
            v(i,j,k) = 0.5 * (factor * sin(x(i)) * cos(y(j)) * sin(z(k)) &
                            - cos(x(i)) * sin(y(j)) * cos(z(k))) &
                            * exp(-(factor**2)*time(1) * nu)
            
            w(i,j,k) = 1.0 * cos(x(i)) * cos(y(j)) * sin(z(k)) & 
                           * exp(-(factor**2)*time(1) * nu)                 
                    
        end do; end do; end do
        
    end subroutine state_init_analytic
        
    subroutine state_kinetic()
        !*************************!
        ! Mean kinetic energy:    !
        !*************************!
        Ekin = 0d0
        do i=1,Nh; do j=1,Ny; do k=1,Nz
            
            ! Spectrum scale factor:
            if (kx(i) .eq. 0.0d0) then              
                scalekx = 0.5d0
            else 
                scalekx = 1.0d0 
            end if
            
            Ekin = Ekin &
                 + (real(conjg(uhat(i, j, k)) * uhat(i, j, k)) &
                  + real(conjg(vhat(i, j, k)) * vhat(i, j, k)) &
                  + real(conjg(what(i, j, k)) * what(i, j, k))) &
                  * scalekx * (scalemodes ** 2)
            
            
        end do; end do; end do   
                
    end subroutine state_kinetic
    
    subroutine state_check_error()
        
        factor = sqrt(3.0d0)
        ! Back to the configuration space:
        call state_uhat2u()
                
        ! Error in final numerical solution:
                
        do k=1,Nz; do j=1,Ny; do i=1,Nx
                    
            utemp(i, j, k) = u(i, j, k) &
                            - (-0.5*( factor*cos(x(i))*sin(y(j))*sin(z(k))&
                            +sin(x(i))*cos(y(j))*cos(z(k)) )*exp(-(factor**2)*time(Nt+1) * nu))
                    
            vtemp(i, j, k) = v(i, j, k) &
                            - (0.5*(  factor*sin(x(i))*cos(y(j))*sin(z(k))&
                            -cos(x(i))*sin(y(j))*cos(z(k)) )*exp(-(factor**2)*time(Nt+1) * nu))
                    
            wtemp(i, j, k) = w(i, j, k) &
                            - (cos(x(i))*cos(y(j))*sin(z(k))*exp(-(factor**2)*time(Nt+1) * nu))
                                
        end do; end do; end do        

        chg = maxval(abs(utemp)) + maxval(abs(vtemp)) + maxval(abs(wtemp))
        
        print *, 'The error at the final timestep is ', chg
        
    end subroutine state_check_error    
    
    
end module state
