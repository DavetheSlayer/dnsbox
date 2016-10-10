module rhs
    
    use state
    implicit none
    
    contains
    
    subroutine rhs_int_fact()
        ! Compute integrating factor for given time-step
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
            intFact(i, j, k) = exp((nu * (kx(k) * kx(k) &
                                        + ky(j) * ky(j) &
                                        + kz(i) * kz(i)) + Q) * dt)
        end do; end do; end do
        
    end subroutine rhs_int_fact
      
    subroutine rhs_nonlinear()
        ! Compute nonlinear term of the Navier-Stokes equation 
        ! in rotation form for (u,v,w)hattemp:
        
        if (bandlim) then
            call rhs_Eband()
        end if
        
        call state_derivatives() 
        call state_vorticity()
            
        ! Configuration space velocity fields fields:
        call state_copy_fourier (uhattemp, temp_c)
        call p3dfft_btran_c2r (temp_c, utemp, 'tff') ! Now has a factor N^3
        
        call state_copy_fourier (vhattemp, temp_c)
        call p3dfft_btran_c2r (temp_c, vtemp, 'tff') ! Now has a factor N^3
        
        call state_copy_fourier (whattemp, temp_c)
        call p3dfft_btran_c2r (temp_c, wtemp, 'tff') ! Now has a factor N^3       

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

            absk = sqrt(real(conjg(kx(k)) * kx(k) &
                           + conjg(ky(j)) * ky(j) &
                           + conjg(kz(i)) * kz(i))) ! |k|
            ! Dealiasing:
            if((spherical .and. &
               (abs(kx(k)) ** 2.0d0 &
              + abs(ky(j)) ** 2.0d0 &
              + abs(kz(i)) ** 2.0d0 & 
            .ge.  ((real(Nx, kind=8) / 2.0d0) * 0.5d0 * alpha_x) ** 2 )) .or. &
                !* (2.0d0 * alpha_x / 3.0d0)) ** 2)) .or. &
               (abs(kx(k)) .ge. (real(Nx, kind=8) / 2.0d0) &
                              * (2.0d0 * alpha_x / 3.0d0)  &
           .or. abs(ky(j)) .ge. (real(Ny, kind=8) / 2.0d0) &
                              * (2.0d0 * alpha_y / 3.0d0)  &
           .or. abs(kz(i)) .ge. (real(Nz, kind=8) / 2.0d0) &
                              * (2.0d0 * alpha_z / 3.0d0)) & ! ) then 
            .or. ((kx(k) .eq. cmplx(0.0d0, 0.0d0)) .and. &
                  (ky(j) .eq. cmplx(0.0d0, 0.0d0)) .and. &
                  (kz(i) .eq. cmplx(0.0d0, 0.0d0)))) then

                nonlinuhat(i, j, k) = cmplx(0.0d0, 0.0d0)
                nonlinvhat(i, j, k) = cmplx(0.0d0, 0.0d0)
                nonlinwhat(i, j, k) = cmplx(0.0d0, 0.0d0)          
                    
            else
                if((kx(k) .eq. cmplx(0.0d0, 0.0d0)) .and. &
                   (ky(j) .eq. cmplx(0.0d0, 0.0d0)) .and. &
                   (kz(i) .eq. cmplx(0.0d0, 0.0d0))) then
                    phat(i, j, k) = 0.0d0                
                else                
                            
                    phat(i, j, k) = -1.0d0 * (kx(k) * nonlinuhat(i, j, k) &
                                            + ky(j) * nonlinvhat(i, j, k) & 
                                            + kz(i) * nonlinwhat(i, j, k)) & 
                                             /  (kx(k) * kx(k) &
                                               + ky(j) * ky(j) & 
                                               + kz(i) * kz(i))
                end if
                ! from n_k(u) to N_k(u)
                if (bandlim .and. absk .le. kCutOff) then
                    nonlinuhat(i, j, k) = - nonlinuhat(i, j, k) &
                      - kx(k) * phat(i, j, k) &
                      + (Pin / (2.0d0 * Eband)) * uhattemp(i, j, k)
                    nonlinvhat(i, j, k) = - nonlinvhat(i, j, k) &
                      - ky(j) * phat(i, j, k) &
                      + (Pin / (2.0d0 * Eband)) * vhattemp(i, j, k)  
                    nonlinwhat(i, j, k) = - nonlinwhat(i, j, k) &
                      - kz(i) * phat(i, j, k) &
                      + (Pin / (2.0d0 * Eband)) * whattemp(i, j, k)  
                else
                    nonlinuhat(i, j, k) = - nonlinuhat(i, j, k) &
                                          - kx(k) * phat(i, j, k) ! &
                                          ! + Q * uhattemp(i, j, k)
                    nonlinvhat(i, j, k) = - nonlinvhat(i, j, k) &
                                          - ky(j) * phat(i, j, k) ! &
                                          ! + Q * vhattemp(i, j, k)
                    nonlinwhat(i, j, k) = - nonlinwhat(i, j, k) &
                                          - kz(i) * phat(i, j, k) ! &
                                          ! + Q * whattemp(i, j, k)
                end if
                
            end if                
            
        end do; end do; end do
        
    end subroutine rhs_nonlinear
  
    subroutine rhs_Eband()
        
        ! Compute total energy in the input band
        real(kind=8) :: scalex   !scaling to avoid counting kx=0 modes twice
        myEband = 0d0
        
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
!            print *, k
            absk = sqrt(real(conjg(kx(k)) * kx(k) &
                           + conjg(ky(j)) * ky(j) &
                           + conjg(kz(i)) * kz(i))) 
               
            if (absk .gt. kCutOff) cycle

            if (kx(k).eq.cmplx(0.0d0, 0.0d0, kind=8)) then              
                scalex = 0.5d0
            else 
                scalex = 1.0d0 
            end if
            
            myEband = myEband &
                   + real(conjg(uhattemp(i, j, k)) * uhattemp(i, j, k)) * scalex &
                   + real(conjg(vhattemp(i, j, k)) * vhattemp(i, j, k)) * scalex &
                   + real(conjg(whattemp(i, j, k)) * whattemp(i, j, k)) * scalex
            
        end do; end do; end do   
        
        myEband = myEband * (scalemodes ** 2)
        call mpi_allreduce(myEband, Eband, 1, mpi_double_precision, &
                           mpi_sum, mpi_comm_world, ierr)
        
    end subroutine rhs_Eband 

  
!    subroutine rhsNonlinearConv()
!        ! Compute nonlinear term of the Navier-Stokes equation 
!        ! in Fourier space for (u,v,w)hattemp:
!        ! in convective form
        
!        call rhsDerivatives() 
            
!        ! Configuration space velocity fields fields:
!        call p3dfft_btran_c2r (uhattemp, utemp, 'tff') ! Now has a factor N^3
!        call p3dfft_btran_c2r (vhattemp, vtemp, 'tff') 
!        call p3dfft_btran_c2r (whattemp, wtemp, 'tff')        
        
!        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
            
!            temp_r(i, j, k) = (utemp(i, j, k) * ux(i, j, k) &
!                             + vtemp(i, j, k) * uy(i, j, k) &
!                             + wtemp(i, j, k) * uz(i, j, k)) &
!                             * scalemodes ! divide by N^3
            
!        end do; end do; end do        
                     
!        call p3dfft_ftran_r2c (temp_r, nonlinuhat, 'fft')
                     
!        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
        
!            temp_r(i, j, k) = (utemp(i, j, k) * vx(i, j, k) &
!                    + vtemp(i, j, k) * vy(i, j, k) &
!                    + wtemp(i, j, k) * vz(i, j, k)) * scalemodes
            
!        end do; end do; end do        
!        call p3dfft_ftran_r2c (temp_r, nonlinvhat, 'fft')
        
!        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
        
!            temp_r(i, j, k) = (utemp(i, j, k) * wx(i, j, k) &
!                    + vtemp(i, j, k) * wy(i, j, k) &
!                    + wtemp(i, j, k) * wz(i, j, k)) * scalemodes
            
!        end do; end do; end do        
!        call p3dfft_ftran_r2c (temp_r, nonlinwhat, 'fft')
        
!        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1) 
!            ! k -> kx - index, j -> ky - index, i -> kz - index
            
!            ! Dealiasing:
!            if((spherical .and. &
!               (abs(kx(k)) ** 2.0d0 &
!              + abs(ky(j)) ** 2.0d0 &
!              + abs(kz(i)) ** 2.0d0 & 
!              .ge. ((real(Nx, kind=8) / 2.0d0) &
!                  * (2.0d0 * alpha_x / 3.0d0)) ** 2)) .or. &
!               (abs(kx(k)) .ge. (real(Nx, kind=8) / 2.0d0) &
!                              * (2.0d0 * alpha_x / 3.0d0)  &
!           .or. abs(ky(j)) .ge. (real(Ny, kind=8) / 2.0d0) &
!                              * (2.0d0 * alpha_y / 3.0d0)  &
!           .or. abs(kz(i)) .ge. (real(Nz, kind=8) / 2.0d0) &
!                              * (2.0d0 * alpha_z / 3.0d0)) .or. &
!               ((kx(k) .eq. cmplx(0.0d0, 0.0d0)) &
!          .and. (ky(j) .eq. cmplx(0.0d0, 0.0d0)) &
!          .and. (kz(i) .eq. cmplx(0.0d0, 0.0d0)))) then

!                nonlinuhat(i, j, k) = cmplx(0.0d0, 0.0d0)
!                nonlinvhat(i, j, k) = cmplx(0.0d0, 0.0d0)
!                nonlinwhat(i, j, k) = cmplx(0.0d0, 0.0d0)          
                    
!            else
            
!                phat(i, j, k) = -1.0d0 * (kx(k) * nonlinuhat(i, j, k) &
!                                        + ky(j) * nonlinvhat(i, j, k) & 
!                                        + kz(i) * nonlinwhat(i, j, k)) & 
!                                         /  (kx(k) * kx(k) &
!                                           + ky(j) * ky(j) & 
!                                           + kz(i) * kz(i))
                
!                nonlinuhat(i, j, k) = - nonlinuhat(i, j, k) &
!                                      - kx(k) * phat(i, j, k)
!                nonlinvhat(i, j, k) = - nonlinvhat(i, j, k) &
!                                      - ky(j) * phat(i, j, k)
!                nonlinwhat(i, j, k) = - nonlinwhat(i, j, k) &
!                                      - kz(i) * phat(i, j, k)
            
!            end if
            
!        end do; end do; end do

!    end subroutine rhsNonlinearConv    

    
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
