module rhs
    
    use variables
    implicit none
    
    contains
        
    subroutine rhsDerivatives()
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
    end subroutine rhsDerivatives
    
    subroutine rhsVorticity()
        ! computes vorticity
        ! should be called after rhsDerivatives()
        
        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
            
            ! omegax:        
            omegax(i, j, k) = wy(i, j, k) - vz(i, j, k)
            omegay(i, j, k) = uz(i, j, k) - wx(i, j, k)
            omegaz(i, j, k) = vx(i, j, k) - uy(i, j, k)
                                
        end do; end do; end do        
        
    end subroutine rhsVorticity

    subroutine rhsIntFact()
        ! Compute integrating factor for given time-step
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
            intFact(i, j, k) = exp((nu * (kx(k) * kx(k) &
                                        + ky(j) * ky(j) &
                                        + kz(i) * kz(i)) + Q) * dt)
        end do; end do; end do
        
    end subroutine rhsIntFact
      
    subroutine rhsNonlinear()
        ! Compute nonlinear term of the Navier-Stokes equation 
        ! in rotation form for (u,v,w)hattemp:
        
        call rhsDerivatives() 
        call rhsVorticity()
            
        ! Configuration space velocity fields fields:
        call p3dfft_btran_c2r (uhattemp, utemp, 'tff') ! Now has a factor N^3
        call p3dfft_btran_c2r (vhattemp, vtemp, 'tff') 
        call p3dfft_btran_c2r (whattemp, wtemp, 'tff')        

        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
            
            temp_r(i, j, k) = (omegay(i, j, k) * wtemp(i, j, k) &
                             - omegaz(i, j, k) * vtemp(i, j, k)) &
                             * scalemodes ! divide by N^3
            
        end do; end do; end do        
                     
        call p3dfft_ftran_r2c (temp_r, nonlinuhat, 'fft')
                     
        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
        
            temp_r(i, j, k) = (omegaz(i, j, k) * utemp(i, j, k) &
                             - omegax(i, j, k) * wtemp(i, j, k)) &
                             * scalemodes ! divide by N^3
            
        end do; end do; end do        
        call p3dfft_ftran_r2c (temp_r, nonlinvhat, 'fft')
        
        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
        
            temp_r(i, j, k) = (omegax(i, j, k) * vtemp(i, j, k) &
                             - omegay(i, j, k) * utemp(i, j, k)) &
                             * scalemodes ! divide by N^3
            
        end do; end do; end do        
        call p3dfft_ftran_r2c (temp_r, nonlinwhat, 'fft')
                
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1) 
            ! k -> kx - index, j -> ky - index, i -> kz - index

            ! 2/3 dealiasing:
            
            ! Cubic truncation:
!            if  ((abs(kx(k)) .gt.  (real(Nx, kind=8) / 2.0d0) &
!                                 * (2.0d0 * alpha_x / 3.0d0))  &
!            .or. (abs(ky(j)) .gt.  (real(Ny, kind=8) / 2.0d0) &
!                                 * (2.0d0 * alpha_y / 3.0d0))  &
!            .or. (abs(kz(i)) .gt.  (real(Nz, kind=8) / 2.0d0) &
!                                 * (2.0d0 * alpha_z / 3.0d0))) then
            
            ! Isotropic truncation:
            if(abs(kx(k)) ** 2.0d0 + abs(ky(j)) ** 2.0d0 + abs(kz(i)) ** 2.0d0 & 
               .gt.  ((real(Nx, kind=8) / 2.0d0) &
                    * (2.0d0 * alpha_x / 3.0d0)) ** 2) then

                      nonlinuhat(i, j, k) = cmplx(0.0d0, 0.0d0)
                      nonlinvhat(i, j, k) = cmplx(0.0d0, 0.0d0)
                      nonlinwhat(i, j, k) = cmplx(0.0d0, 0.0d0)          
                
            else 
                if((kx(k) .eq. cmplx(0.0d0, 0.0d0)) .and. &
                   (ky(j) .eq. cmplx(0.0d0, 0.0d0)) .and. &
                   (kz(i) .eq. cmplx(0.0d0, 0.0d0))) then
                    eps = 0.1d0 ** 13
                else
                    eps = 0.0d0
                end if                
                
                
                phat(i, j, k) = -1.0d0 * (kx(k) * nonlinuhat(i, j, k) &
                                        + ky(j) * nonlinvhat(i, j, k) & 
                                        + kz(i) * nonlinwhat(i, j, k)) & 
                                         /  (kx(k) * kx(k) &
                                           + ky(j) * ky(j) & 
                                           + kz(i) * kz(i) + eps)
                
                ! from n_k(u) to N_k(u)
                nonlinuhat(i, j, k) = - nonlinuhat(i, j, k) &
                                      - kx(k) * phat(i, j, k)
                nonlinvhat(i, j, k) = - nonlinvhat(i, j, k) &
                                      - ky(j) * phat(i, j, k)
                nonlinwhat(i, j, k) = - nonlinwhat(i, j, k) &
                                      - kz(i) * phat(i, j, k)
            
            end if
            
        end do; end do; end do
        
    end subroutine rhsNonlinear
    
!    subroutine rhsAll()
!        ! Compute RHS for uhattemp
!        ! utilde in Adjoint descent
!        call rhsNonlinear()
!        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1) 
!            rhsuhat(i, j, k) = (nu *(kx(k) * kx(k) &
!                                   + ky(j) * ky(j) & 
!                                   + kz(i) * kz(i)) + Q) * uhattemp(i, j, k)&
!                               + nonlinuhat(i, j, k)
!            rhsvhat(i, j, k) = (nu *(kx(k) * kx(k) &
!                                   + ky(j) * ky(j) & 
!                                   + kz(i) * kz(i)) + Q) * vhattemp(i, j, k)&
!                               + nonlinvhat(i, j, k)
!            rhswhat(i, j, k) = (nu *(kx(k) * kx(k) &
!                                   + ky(j) * ky(j) & 
!                                   + kz(i) * kz(i)) + Q) * whattemp(i, j, k)&
!                               + nonlinwhat(i, j, k)
            
!        end do; end do; end do
        
        
!    end subroutine rhsAll
    
end module rhs
