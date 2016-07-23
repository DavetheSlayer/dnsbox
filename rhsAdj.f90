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
        ! in Fourier space for (u,v,w)hattemp:
        
        call rhsDerivatives() 
        
!        if(proc_id .eq. 0) then 
!            print *, 'check nonlin 1'
!        endif    
            
        ! Configuration space velocity fields fields:
        call p3dfft_btran_c2r (uhattemp, utemp, 'tff') ! Now has a factor N^3
        call p3dfft_btran_c2r (vhattemp, vtemp, 'tff') 
        call p3dfft_btran_c2r (whattemp, wtemp, 'tff')        
        
!        if(proc_id .eq. 0) then 
!            print *, 'check nonlin 2'
!        endif    
        
        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
            
            temp_r(i, j, k) = (utemp(i, j, k) * ux(i, j, k) &
                             + vtemp(i, j, k) * uy(i, j, k) &
                             + wtemp(i, j, k) * uz(i, j, k)) &
                             * scalemodes ! divide by N^3
            
        end do; end do; end do        
                   
!        if(proc_id .eq. 0) then 
!            print *, 'check nonlin 3.1.1'
!        endif               
                     
        call p3dfft_ftran_r2c (temp_r, nonlinuhat, 'fft')
           
!        if(proc_id .eq. 0) then 
!            print *, 'check nonlin 3.1.2'
!        endif               
                     
        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
        
            temp_r(i, j, k) = (utemp(i, j, k) * vx(i, j, k) &
                    + vtemp(i, j, k) * vy(i, j, k) &
                    + wtemp(i, j, k) * vz(i, j, k)) * scalemodes
            
        end do; end do; end do        
        call p3dfft_ftran_r2c (temp_r, nonlinvhat, 'fft')
           
!        if(proc_id .eq. 0) then 
!            print *, 'check nonlin 3.2'
!        endif               
        
        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
        
            temp_r(i, j, k) = (utemp(i, j, k) * wx(i, j, k) &
                    + vtemp(i, j, k) * wy(i, j, k) &
                    + wtemp(i, j, k) * wz(i, j, k)) * scalemodes
            
        end do; end do; end do        
        call p3dfft_ftran_r2c (temp_r, nonlinwhat, 'fft')
        
!        if(proc_id .eq. 0) then 
!            print *, 'check nonlin 3.3'
!        endif               
                
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1) 
            ! k -> kx - index, j -> ky - index, i -> kz - index
            
            if  ((abs(kx(k)) .gt.  (real(Nx, kind=8) / 2.0d0) &
                                 * (2.0d0 * alpha_x / 3.0d0))  &
            .or. (abs(ky(j)) .gt.  (real(Ny, kind=8) / 2.0d0) &
                                 * (2.0d0 * alpha_y / 3.0d0))  &
            .or. (abs(kz(i)) .gt.  (real(Nz, kind=8) / 2.0d0) &
                                 * (2.0d0 * alpha_z / 3.0d0))) then
            ! 2/3 dealiasing:
            
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

!        if(proc_id .eq. 0) then 
!            print *, 'check nonlin 4'
!        endif               
        
    end subroutine rhsNonlinear
    
    subroutine rhsupptilde()
        ! Compute RHS of NSe (utilde) for uhattemp
        call rhsNonlinear()
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1) 
            
            if (hminone) then
                factor = 1.0d0 / (1.0d0 - real((kx(k) * kx(k) &
                                              + ky(j) * ky(j) & 
                                              + kz(i) * kz(i)), kind=8))
            else
                factor = 1.0d0
            end if
        
            upptildehat(i, j, k) = factor &
                                   * ((nu *(kx(k) * kx(k) &
                                          + ky(j) * ky(j) & 
                                          + kz(i) * kz(i)) + Q) &
                                       * uhattemp(i, j, k) &
                                      + nonlinuhat(i, j, k))

            vpptildehat(i, j, k) = factor &
                                   * ((nu *(kx(k) * kx(k) &
                                          + ky(j) * ky(j) & 
                                          + kz(i) * kz(i)) + Q) &
                                       * vhattemp(i, j, k) &
                                      + nonlinvhat(i, j, k))

            wpptildehat(i, j, k) = factor &
                                   * ((nu *(kx(k) * kx(k) &
                                          + ky(j) * ky(j) & 
                                          + kz(i) * kz(i)) + Q) &
                                       * whattemp(i, j, k) &
                                      + nonlinwhat(i, j, k))
            
        end do; end do; end do
        
    end subroutine rhsupptilde

    subroutine rhsDerivativespp()
        ! Compute space-derivatives of u,v,w from (u,v,w)hattemp
    
        ! Derivative of u with respect to x, y, z:
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
            temp_c(i, j, k) = upptildehat(i, j, k) * kx(k) * scalemodes
        end do; end do; end do
        call p3dfft_btran_c2r (temp_c, upptildex, 'tff')
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
            temp_c(i, j, k) = upptildehat(i, j, k) * ky(j) * scalemodes
        end do; end do; end do
        call p3dfft_btran_c2r (temp_c, upptildey, 'tff')
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
            temp_c(i, j, k) = upptildehat(i, j, k) * kz(i) * scalemodes
        end do; end do; end do
        call p3dfft_btran_c2r (temp_c, upptildez, 'tff')
        
        ! Derivative of v with respect to x, y, z:
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
            temp_c(i, j, k) = vpptildehat(i, j, k) * kx(k) * scalemodes
        end do; end do; end do
        call p3dfft_btran_c2r (temp_c, vpptildex, 'tff')
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
            temp_c(i, j, k) = vpptildehat(i, j, k) * ky(j) * scalemodes
        end do; end do; end do
        call p3dfft_btran_c2r (temp_c, vpptildey, 'tff')
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
            temp_c(i, j, k) = vpptildehat(i, j, k) * kz(i) * scalemodes
        end do; end do; end do
        call p3dfft_btran_c2r (temp_c, vpptildez, 'tff')        
        
        ! Derivative of w with respect to x, y, z:
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
            temp_c(i, j, k) = wpptildehat(i, j, k) * kx(k) * scalemodes
        end do; end do; end do
        call p3dfft_btran_c2r (temp_c, wpptildex, 'tff')
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
            temp_c(i, j, k) = wpptildehat(i, j, k) * ky(j) * scalemodes
        end do; end do; end do
        call p3dfft_btran_c2r (temp_c, wpptildey, 'tff')
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
            temp_c(i, j, k) = wpptildehat(i, j, k) * kz(i) * scalemodes
        end do; end do; end do
        call p3dfft_btran_c2r (temp_c, wpptildez, 'tff')
    end subroutine rhsDerivativespp
    
    subroutine rhsNonlinearAdj()
        ! Pseudospectral computation of nonlinear part of the adjoint equation
        call rhsupptilde()
        call rhsDerivativespp() 
            
        ! Configuration space velocity fields fields:
        call p3dfft_btran_c2r (uhattemp, utemp, 'tff') ! Now has a factor N^3
        call p3dfft_btran_c2r (vhattemp, vtemp, 'tff') 
        call p3dfft_btran_c2r (whattemp, wtemp, 'tff')        
        
        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
            
            temp_r(i, j, k) = (upptildex(i,j,k) * utemp(i,j,k) &
                             + vpptildex(i,j,k) * vtemp(i,j,k) &
                             + wpptildex(i,j,k) * wtemp(i,j,k) &
                             + upptildex(i,j,k) * utemp(i,j,k) &
                             + upptildey(i,j,k) * vtemp(i,j,k) &
                             + upptildez(i,j,k) * wtemp(i,j,k)) &
                            * scalemodes ! divide by N^3
            
        end do; end do; end do        
                     
        call p3dfft_ftran_r2c (temp_r, nonlinuhatadj, 'fft')
                                
        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
        
            temp_r(i, j, k) = (upptildey(i,j,k) * utemp(i,j,k) &
                             + vpptildey(i,j,k) * vtemp(i,j,k) &
                             + wpptildey(i,j,k) * wtemp(i,j,k) &
                             + vpptildex(i,j,k) * utemp(i,j,k) &
                             + vpptildey(i,j,k) * vtemp(i,j,k) &
                             + vpptildez(i,j,k) * wtemp(i,j,k)) &
                            * scalemodes ! divide by N^3
            
        end do; end do; end do        
        call p3dfft_ftran_r2c (temp_r, nonlinvhatadj, 'fft')
        
        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
        
            temp_r(i, j, k) = (upptildez(i,j,k) * utemp(i,j,k) &
                             + vpptildez(i,j,k) * vtemp(i,j,k) &
                             + wpptildez(i,j,k) * wtemp(i,j,k) &
                             + wpptildex(i,j,k) * utemp(i,j,k) &
                             + wpptildey(i,j,k) * vtemp(i,j,k) &
                             + wpptildez(i,j,k) * wtemp(i,j,k)) &
                            * scalemodes ! divide by N^3
            
        end do; end do; end do        
        call p3dfft_ftran_r2c (temp_r, nonlinwhatadj, 'fft')
                
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1) 
            ! k -> kx - index, j -> ky - index, i -> kz - index
            
            if  ((abs(kx(k)) .gt.  (real(Nx, kind=8) / 2.0d0) &
                                 * (2.0d0 * alpha_x / 3.0d0))  &
            .or. (abs(ky(j)) .gt.  (real(Ny, kind=8) / 2.0d0) &
                                 * (2.0d0 * alpha_y / 3.0d0))  &
            .or. (abs(kz(i)) .gt.  (real(Nz, kind=8) / 2.0d0) &
                                 * (2.0d0 * alpha_z / 3.0d0))) then
            ! 2/3 dealiasing:
            
                      nonlinuhatadj(i, j, k) = cmplx(0.0d0, 0.0d0)
                      nonlinvhatadj(i, j, k) = cmplx(0.0d0, 0.0d0)
                      nonlinwhatadj(i, j, k) = cmplx(0.0d0, 0.0d0)          
                
            else 
                if((kx(k) .eq. cmplx(0.0d0, 0.0d0)) .and. &
                   (ky(j) .eq. cmplx(0.0d0, 0.0d0)) .and. &
                   (kz(i) .eq. cmplx(0.0d0, 0.0d0))) then
                    eps = 0.1d0 ** 13
                else
                    eps = 0.0d0
                end if                
                
                phat(i, j, k) = -1.0d0 * (kx(k) * nonlinuhatadj(i, j, k) &
                                        + ky(j) * nonlinvhatadj(i, j, k) & 
                                        + kz(i) * nonlinwhatadj(i, j, k)) & 
                                         /  (kx(k) * kx(k) &
                                           + ky(j) * ky(j) & 
                                           + kz(i) * kz(i) + eps)
                
                ! from n_k(u) to N_k(u)
                nonlinuhatadj(i, j, k) = - nonlinuhatadj(i, j, k) &
                                      - kx(k) * phatadj(i, j, k)
                nonlinvhatadj(i, j, k) = - nonlinvhatadj(i, j, k) &
                                      - ky(j) * phatadj(i, j, k)
                nonlinwhatadj(i, j, k) = - nonlinwhatadj(i, j, k) &
                                      - kz(i) * phatadj(i, j, k)
            
            end if
            
        end do; end do; end do
        
    end subroutine rhsNonlinearAdj

    subroutine rhsAdj()
        ! Compute RHS of NSe (utilde) for uhattemp
        call rhsNonlinearAdj()
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1) 

            rhsadjuhat(i, j, k) = (nu *(kx(k) * kx(k) &
                                      + ky(j) * ky(j) & 
                                      + kz(i) * kz(i)) + Q) &
                                  * upptildehat(i, j, k) &
                                + nonlinuhatadj(i, j, k)

            rhsadjvhat(i, j, k) = (nu *(kx(k) * kx(k) &
                                      + ky(j) * ky(j) & 
                                      + kz(i) * kz(i)) + Q) &
                                  * vpptildehat(i, j, k) &
                                + nonlinvhatadj(i, j, k)

            rhsadjwhat(i, j, k) = (nu *(kx(k) * kx(k) &
                                      + ky(j) * ky(j) & 
                                      + kz(i) * kz(i)) + Q) &
                                  * wpptildehat(i, j, k) &
                                + nonlinwhatadj(i, j, k)
            
        end do; end do; end do
        
    end subroutine rhsAdj
    
end module rhs
