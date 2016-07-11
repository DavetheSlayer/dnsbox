module rhs
    
    use variables
    implicit none
    
    contains
    
    subroutine rhsProject()
        ! Apply projection on (u,v,w)hat to make it divergence-free:
        
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
            
            uhattemp(i, j, k) = uhat(i, j, k)
            vhattemp(i, j, k) = vhat(i, j, k)
            whattemp(i, j, k) = what(i, j, k)
            
        end do; end do; end do
         
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1) 
            ! k -> kx - index, j -> ky - index, i -> kz - index
          
            uhat(i, j, k) = uhattemp(i, j, k) &
                          - kx(k) * kx(k) * uhattemp(i, j, k) / &
                            (kx(k) * kx(k) &
                           + ky(j) * ky(j) & 
                           + kz(i) * kz(i) + 0.1d0 ** 13) &
                          - kx(k) * ky(j) * vhattemp(i, j, k) / &
                            (kx(k) * kx(k) &
                           + ky(j) * ky(j) & 
                           + kz(i) * kz(i) + 0.1d0 ** 13) &
                          - kx(k) * kz(i) * whattemp(i, j, k) / &
                            (kx(k) * kx(k) &
                           + ky(j) * ky(j) & 
                           + kz(i) * kz(i) + 0.1d0 ** 13) 
            
            vhat(i, j, k) = vhattemp(i, j, k) &
                          - ky(j) * kx(k) * uhattemp(i, j, k) / &
                            (kx(k) * kx(k) &
                           + ky(j) * ky(j) & 
                           + kz(i) * kz(i) + 0.1d0 ** 13) &
                          - ky(j) * ky(j) * vhattemp(i, j, k) / &
                            (kx(k) * kx(k) &
                           + ky(j) * ky(j) & 
                           + kz(i) * kz(i) + 0.1d0 ** 13) &
                          - ky(j) * kz(i) * whattemp(i, j, k) / &
                            (kx(k) * kx(k) &
                           + ky(j) * ky(j) & 
                           + kz(i) * kz(i) + 0.1d0 ** 13) 
            
            what(i, j, k) = whattemp(i, j, k) &
                          - kz(i) * kx(k) * uhattemp(i, j, k) / &
                            (kx(k) * kx(k) &
                           + ky(j) * ky(j) & 
                           + kz(i) * kz(i) + 0.1d0 ** 13) &
                          - kz(i) * ky(j) * vhattemp(i, j, k) / &
                            (kx(k) * kx(k) &
                           + ky(j) * ky(j) & 
                           + kz(i) * kz(i) + 0.1d0 ** 13) &
                          - kz(i) * kz(i) * whattemp(i, j, k) / &
                            (kx(k) * kx(k) &
                           + ky(j) * ky(j) & 
                           + kz(i) * kz(i) + 0.1d0 ** 13) 
                                       
        end do; end do; end do
        
    end subroutine rhsProject
        
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
        call p3dfft_btran_c2r (uhattemp, utemp, 'tff') ! Now has a factor 
        call p3dfft_btran_c2r (vhattemp, vtemp, 'tff') ! N^(3/2)
        call p3dfft_btran_c2r (whattemp, wtemp, 'tff')        
        
!        if(proc_id .eq. 0) then 
!            print *, 'check nonlin 2'
!        endif    
        
        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
            
            temp_r(i, j, k) = (utemp(i, j, k) * ux(i, j, k) &
                             + vtemp(i, j, k) * uy(i, j, k) &
                             + wtemp(i, j, k) * uz(i, j, k)) &
                             * scalemodessquare ! divide by N^3
            
        end do; end do; end do        
                   
!        if(proc_id .eq. 0) then 
!            print *, 'check nonlin 3.1.1'
!        endif               
                     
        call p3dfft_ftran_r2c (temp_r, nonlinuhat, 'fft') ! another N^(3/2)
           
!        if(proc_id .eq. 0) then 
!            print *, 'check nonlin 3.1.2'
!        endif               
                     
        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
        
            temp_r(i, j, k) = (utemp(i, j, k) * vx(i, j, k) &
                    + vtemp(i, j, k) * vy(i, j, k) &
                    + wtemp(i, j, k) * vz(i, j, k)) * scalemodessquare
            
        end do; end do; end do        
        call p3dfft_ftran_r2c (temp_r, nonlinvhat, 'fft')
           
!        if(proc_id .eq. 0) then 
!            print *, 'check nonlin 3.2'
!        endif               
        
        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
        
            temp_r(i, j, k) = (utemp(i, j, k) * wx(i, j, k) &
                    + vtemp(i, j, k) * wy(i, j, k) &
                    + wtemp(i, j, k) * wz(i, j, k)) * scalemodessquare
            
        end do; end do; end do        
        call p3dfft_ftran_r2c (temp_r, nonlinwhat, 'fft')
        
!        if(proc_id .eq. 0) then 
!            print *, 'check nonlin 3.3'
!        endif               
                
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1) 
            ! k -> kx - index, j -> ky - index, i -> kz - index
            
            if  (abs(kx(k)) >  (real(Nx, kind=8) / 2.0d0) &
                                 * (2.0d0 * alpha_x / 3.0d0)  &
                .or. abs(ky(j)) >  (real(Ny, kind=8) / 2.0d0) &
                                 * (2.0d0 * alpha_y / 3.0d0)  &
                .or. abs(kz(i)) >  (real(Nz, kind=8) / 2.0d0) &
                                 * (2.0d0 * alpha_z / 3.0d0)) then
            ! 2/3 dealiasing:
            
                      nonlinuhat(i, j, k) = cmplx(0.0d0, 0.0d0)
                      nonlinvhat(i, j, k) = cmplx(0.0d0, 0.0d0)
                      nonlinwhat(i, j, k) = cmplx(0.0d0, 0.0d0)          
                
            else 
            
                phat(i, j, k) = -1.0d0 * (kx(k) * nonlinuhat(i, j, k) &
                                        + ky(j) * nonlinvhat(i, j, k) & 
                                        + kz(i) * nonlinwhat(i, j, k)) & 
                                         /  (kx(k) * kx(k) &
                                           + ky(j) * ky(j) & 
                                           + kz(i) * kz(i) + 0.1d0 ** 13)
                
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
    
!    subroutine rhsAll()
!        ! Compute RHS for uhattemp
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

    subroutine rhsDealias()
    ! Set k > 2/3 k_max elements to 0
            
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1) 
            
            if  (abs(kx(k)) >  (real(Nx, kind=8) / 2.0d0) &
                             * (2.0d0 * alpha_x / 3.0d0)  &
            .or. abs(ky(j)) >  (real(Ny, kind=8) / 2.0d0) &
                             * (2.0d0 * alpha_y / 3.0d0)  &
            .or. abs(kz(i)) >  (real(Nz, kind=8) / 2.0d0) &
                             * (2.0d0 * alpha_z / 3.0d0)) then

                  uhat(i, j, k) = cmplx(0.0d0, 0.0d0)
                  vhat(i, j, k) = cmplx(0.0d0, 0.0d0)
                  what(i, j, k) = cmplx(0.0d0, 0.0d0)          

            end if
            
        end do; end do; end do
        
    end subroutine rhsDealias
    
end module rhs
