module rhs
    
    use state
    implicit none
    
    contains
        
    subroutine rhs_derivatives()
        ! Compute space-derivatives of u,v,w from (u,v,w)hattemp
    
        ! Derivative of u with respect to x, y, z:
        do i=1,Nh; do j=1,Ny; do k=1,Nz
            temp_c(i, j, k) = uhattemp(i, j, k) * cmplx(0.0d0, 1.0d0) * kx(i) &
                            * scalemodes
        end do; end do; end do
        call dfftw_execute_(plan_backward_temp)
        call state_copy_conf(temp_r, ux)
        
        do i=1,Nh; do j=1,Ny; do k=1,Nz
            temp_c(i, j, k) = uhattemp(i, j, k) * cmplx(0.0d0, 1.0d0) * ky(j) &
                            * scalemodes
        end do; end do; end do
        call dfftw_execute_(plan_backward_temp)
        call state_copy_conf(temp_r, uy)
                
        do i=1,Nh; do j=1,Ny; do k=1,Nz
            temp_c(i, j, k) = uhattemp(i, j, k) * cmplx(0.0d0, 1.0d0) * kz(k) &
                            * scalemodes
        end do; end do; end do
        call dfftw_execute_(plan_backward_temp)
        call state_copy_conf(temp_r, uz)

        ! Derivative of v with respect to x, y, z:
        do i=1,Nh; do j=1,Ny; do k=1,Nz
            temp_c(i, j, k) = vhattemp(i, j, k) * cmplx(0.0d0, 1.0d0) * kx(i) &
                            * scalemodes
        end do; end do; end do
        call dfftw_execute_(plan_backward_temp)
        call state_copy_conf(temp_r, vx)
        
        do i=1,Nh; do j=1,Ny; do k=1,Nz
            temp_c(i, j, k) = vhattemp(i, j, k) * cmplx(0.0d0, 1.0d0) * ky(j) &
                            * scalemodes
        end do; end do; end do
        call dfftw_execute_(plan_backward_temp)
        call state_copy_conf(temp_r, vy)
                
        do i=1,Nh; do j=1,Ny; do k=1,Nz
            temp_c(i, j, k) = vhattemp(i, j, k) * cmplx(0.0d0, 1.0d0) * kz(k) &
                            * scalemodes
        end do; end do; end do
        call dfftw_execute_(plan_backward_temp)
        call state_copy_conf(temp_r, vz)
           
        ! Derivative of w with respect to x, y, z:
        do i=1,Nh; do j=1,Ny; do k=1,Nz
            temp_c(i, j, k) = whattemp(i, j, k) * cmplx(0.0d0, 1.0d0) * kx(i) &
                            * scalemodes
        end do; end do; end do
        call dfftw_execute_(plan_backward_temp)
        call state_copy_conf(temp_r, wx)
        
        do i=1,Nh; do j=1,Ny; do k=1,Nz
            temp_c(i, j, k) = whattemp(i, j, k) * cmplx(0.0d0, 1.0d0) * ky(j) &
                            * scalemodes
        end do; end do; end do
        call dfftw_execute_(plan_backward_temp)
        call state_copy_conf(temp_r, wy)
                
        do i=1,Nh; do j=1,Ny; do k=1,Nz
            temp_c(i, j, k) = whattemp(i, j, k) * cmplx(0.0d0, 1.0d0) * kz(k) &
                            * scalemodes
        end do; end do; end do
        call dfftw_execute_(plan_backward_temp)
        call state_copy_conf(temp_r, wz)
        
    end subroutine rhs_derivatives
    
    subroutine rhs_vorticity()
        ! computes vorticity
        ! should be called after rhs_derivatives()
        
        do i=1,Nx; do j=1,Ny; do k=1,Nz
            
            ! omegax:        
            omegax(i, j, k) = wy(i, j, k) - vz(i, j, k)
            omegay(i, j, k) = uz(i, j, k) - wx(i, j, k)
            omegaz(i, j, k) = vx(i, j, k) - uy(i, j, k)
                                
        end do; end do; end do        
        
    end subroutine rhs_vorticity

    subroutine rhs_nonlinear()
        ! Compute nonlinear term of the Navier-Stokes equation 
        ! in rotation form for (u,v,w)hattemp:
        
        call rhs_derivatives() 
        call rhs_vorticity()
        
        if (bandlim) then
            call rhs_Eband()
        end if
            
        ! Configuration space velocity fields fields:
        call state_uhattemp2utemp() ! Now has a factor N^3
        
        do i=1,Nx; do j=1,Ny; do k=1,Nz
            
            temp_r(i, j, k) = (omegay(i, j, k) * wtemp(i, j, k) &
                             - omegaz(i, j, k) * vtemp(i, j, k)) &
                             * scalemodes ! divide by N^3
            
        end do; end do; end do        
        
        call dfftw_execute_(plan_forward_temp)
        call state_copy_fourier(temp_c, nonlinuhat)
                     
        do i=1,Nx; do j=1,Ny; do k=1,Nz
        
            temp_r(i, j, k) = (omegaz(i, j, k) * utemp(i, j, k) &
                             - omegax(i, j, k) * wtemp(i, j, k)) &
                             * scalemodes ! divide by N^3
            
        end do; end do; end do        

        call dfftw_execute_(plan_forward_temp)
        call state_copy_fourier(temp_c, nonlinvhat)
                             
        do i=1,Nx; do j=1,Ny; do k=1,Nz
        
            temp_r(i, j, k) = (omegax(i, j, k) * vtemp(i, j, k) &
                             - omegay(i, j, k) * utemp(i, j, k)) &
                             * scalemodes ! divide by N^3
            
        end do; end do; end do        

        call dfftw_execute_(plan_forward_temp)
        call state_copy_fourier(temp_c, nonlinwhat)
                
        do i=1,Nh; do j=1,Ny; do k=1,Nz            
        
            absk = sqrt(kx(i) * kx(i) + ky(j) * ky(j) + kz(k) * kz(k)) ! |k|
            ! Dealiasing:
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

                nonlinuhat(i, j, k) = cmplx(0.0d0, 0.0d0)
                nonlinvhat(i, j, k) = cmplx(0.0d0, 0.0d0)
                nonlinwhat(i, j, k) = cmplx(0.0d0, 0.0d0)          
                phat(i, j, k) = cmplx(0.0d0, 0.0d0)
                    
            else
                phat(i, j, k) = cmplx(0.0d0, 1.0d0) &
                              * (kx(i) * nonlinuhat(i, j, k) &
                               + ky(j) * nonlinvhat(i, j, k) & 
                               + kz(k) * nonlinwhat(i, j, k)) / absk ** 2
            end if

            
            ! from n_k(u) to N_k(u)
            if (bandlim .and. absk .lt. kCutOff) then
                nonlinuhat(i, j, k) = - nonlinuhat(i, j, k) &
                               - cmplx(0.0d0, 1.0d0) * kx(i) * phat(i, j, k) &
                               + (Pin / (2.0d0 * Eband)) * uhattemp(i, j, k)
                nonlinvhat(i, j, k) = - nonlinvhat(i, j, k) &
                               - cmplx(0.0d0, 1.0d0) * ky(j) * phat(i, j, k) &
                               + (Pin / (2.0d0 * Eband)) * vhattemp(i, j, k)  
                nonlinwhat(i, j, k) = - nonlinwhat(i, j, k) &
                               - cmplx(0.0d0, 1.0d0) * kz(k) * phat(i, j, k) &
                               + (Pin / (2.0d0 * Eband)) * whattemp(i, j, k)  
            else
                nonlinuhat(i, j, k) = - nonlinuhat(i, j, k) &
                                      - kx(i) * phat(i, j, k) ! &
                                      ! + Q * uhattemp(i, j, k)
                nonlinvhat(i, j, k) = - nonlinvhat(i, j, k) &
                                      - ky(j) * phat(i, j, k) ! &
                                      ! + Q * vhattemp(i, j, k)
                nonlinwhat(i, j, k) = - nonlinwhat(i, j, k) &
                                      - kz(k) * phat(i, j, k) ! &
                                      ! + Q * whattemp(i, j, k)
            end if
            
        end do; end do; end do
        
    end subroutine rhs_nonlinear
            
    subroutine rhs_Eband()
        
        ! Compute total energy in the input band
        real(kind=8) :: scalex   !scaling to avoid counting kx=0 modes twice
        Eband = 0d0
        
        do i=1,Nh; do j=1,Ny; do k=1,Nz
!            print *, k
            absk = sqrt(kx(i) * kx(i) + ky(j) * ky(j) + kz(k) * kz(k)) ! |k|
               
            if (absk .ge. kCutOff) cycle

            if (kx(i).eq. 0.0d0) then              
                scalex = 0.5d0
            else 
                scalex = 1.0d0 
            end if
            
            Eband = Eband &
                  + (real(conjg(uhattemp(i, j, k)) * uhattemp(i, j, k))  &
                   + real(conjg(vhattemp(i, j, k)) * vhattemp(i, j, k))  &
                   + real(conjg(whattemp(i, j, k)) * whattemp(i, j, k))) &
                  * scalex * (scalemodes ** 2)
            
        end do; end do; end do   
        
    end subroutine rhs_Eband 
    
    
    
end module rhs
