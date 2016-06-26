module rhs
    
    use variables
    implicit none
    
    contains
    
    
    subroutine rhsFix()
        
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
                            
            rhsuhatfix(i, j, k) = (dtInv + (0.5 * ReInv) &
                                * (kx(k) * kx(k) &
                                 + ky(j) * ky(j) &
                                 + kz(i) * kz(i))) * uhat(i, j, k)
            
            rhsvhatfix(i, j, k) = (dtInv + (0.5 * ReInv) &
                                * (kx(k) * kx(k) &
                                 + ky(j) * ky(j) &
                                 + kz(i) * kz(i))) * vhat(i, j, k)
            
            rhswhatfix(i, j, k) = (dtInv + (0.5 * ReInv) &
                                * (kx(k) * kx(k) &
                                 + ky(j) * ky(j) &
                                 + kz(i) * kz(i))) * what(i, j, k)
                                 
        end do; end do; end do        
            
    end subroutine rhsFix


    subroutine rhsRest()
        
        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
        
            temp_r(i, j, k) = 0.25d0 * ((u(i, j, k) + uold(i, j, k)) &
                                        * (ux(i, j, k) + uxold(i, j, k)) &
                                       + (v(i, j, k) + vold(i, j, k)) &
                                        * (uy(i, j, k) + uyold(i, j, k)) &
                                       + (w(i, j, k) + wold(i, j, k)) &
                                        * (uz(i, j, k) + uzold(i, j, k)))
            
        end do; end do; end do
        
        call p3dfft_ftran_r2c (temp_r, nonlinuhat, 'fft')
        
        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
        
            temp_r(i, j, k) = 0.25d0 * ((u(i, j, k) + uold(i, j, k)) &
                                        * (vx(i, j, k) + vxold(i, j, k)) &
                                       + (v(i, j, k) + vold(i, j, k)) &
                                        * (vy(i, j, k) + vyold(i, j, k)) &
                                       + (w(i, j, k) + wold(i, j, k)) &
                                        * (vz(i, j, k) + vzold(i, j, k)))
            
        end do; end do; end do
        
        call p3dfft_ftran_r2c (temp_r, nonlinvhat, 'fft')
        
        do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
        
            temp_r(i, j, k) = 0.25d0 * ((u(i, j, k) + uold(i, j, k)) &
                                        * (wx(i, j, k) + wxold(i, j, k)) &
                                       + (v(i, j, k) + vold(i, j, k)) &
                                        * (wy(i, j, k) + wyold(i, j, k)) &
                                       + (w(i, j, k) + wold(i, j, k)) &
                                        * (wz(i, j, k) + wzold(i, j, k)))
            
        end do; end do; end do
        
        call p3dfft_ftran_r2c (temp_r, nonlinwhat, 'fft')
        
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
            
            phat(i, j, k) = -1.0d0 * (kx(k) * nonlinuhat(i, j, k) &
                                     + ky(j) * nonlinvhat(i, j, k) & 
                                     + kz(i) * nonlinwhat(i, j, k)) & 
                                     / (kx(k) * kx(k) &
                                       + ky(j) * ky(j) & 
                                       + kz(i) * kz(i) + 0.1d0 ** 13)
                                       
        end do; end do; end do
        
    end subroutine rhsRest
    
end module rhs
