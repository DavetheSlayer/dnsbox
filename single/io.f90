!**************************************************************************
!  IN/OUT 
!
! netcdf reference:
! http://www.unidata.ucar.edu/software/netcdf/netcdf-4/newdocs/netcdf-f90.html
!**************************************************************************
 module io
!**************************************************************************
    use state
    use netcdf  

    integer(kind=8) :: io_save_count   ! save counter
        
    integer,     private  :: io_stats, io_Spec, io_Pars 
    character(4) :: cnum        ! Save count
    
    ! netcdf variables:
    integer :: ncid, fields_id, x_dimid, y_dimid, z_dimid, dimids(3), &
               x_varid, y_varid, z_varid, &
               u_varid, v_varid, w_varid

    ! Variables for reading parameters:
    integer(kind=4) :: io_Nx, io_Ny, io_Nz
    real(kind=8) :: io_alpha_x, io_alpha_y, io_alpha_z
    
contains
    
    subroutine io_init()
        
        ! Initialize input/output
        io_save_count = 0
        io_Ekin=20
        call io_save_info()
        
    end subroutine io_init
    
    subroutine io_finalize()
        ! Close files
    
        if(proc_id.eq.0) then
            close(io_Ekin)
        end if     
        if (proc_id.eq.0) then
            print *,'Saving the final state'
        end if        
        call io_save_state()
        
    end subroutine io_finalize
    
    subroutine io_save_state()
        !
        ! Saves the current state
        !
                   
        write(cnum,'(I4.4)') io_save_count
        print*, ' saving state'//cnum//'.nc  t=', time(n+1)
        
        ! Create a new state file, overwrite the existing one (NF90_ClOBBER)
        call io_check( nf90_create('state'//cnum//'.nc', NF90_CLOBBER, ncid) )
        ! Define dimensions:
        call io_check( nf90_def_dim(ncid, "x", Nx, x_dimid) )
        call io_check( nf90_def_dim(ncid, "y", Ny, y_dimid) )
        call io_check( nf90_def_dim(ncid, "z", Nz, z_dimid) )
        ! Define dimension variables:
        call io_check( nf90_def_var(ncid, "x", NF90_DOUBLE, x_dimid, x_varid) )
        call io_check( nf90_def_var(ncid, "y", NF90_DOUBLE, y_dimid, y_varid) )
        call io_check( nf90_def_var(ncid, "z", NF90_DOUBLE, z_dimid, z_varid) )
        dimids = (/ x_dimid, y_dimid, z_dimid /)
        ! Define fields
        call io_check( nf90_def_var(ncid, "u", NF90_DOUBLE, dimids, u_varid) )
        call io_check( nf90_def_var(ncid, "v", NF90_DOUBLE, dimids, v_varid) )
        call io_check( nf90_def_var(ncid, "w", NF90_DOUBLE, dimids, w_varid) )

        ! End define mode
        call io_check( nf90_enddef(ncid))
        ! Write coordinate variables:
        call io_check( nf90_put_var(ncid, x_varid, x) )
        call io_check( nf90_put_var(ncid, y_varid, y) )
        call io_check( nf90_put_var(ncid, z_varid, z) )
        ! Write fields
        call io_check( nf90_put_var(ncid, u_varid, u) )
        call io_check( nf90_put_var(ncid, v_varid, v) )
        call io_check( nf90_put_var(ncid, w_varid, w) )

        ! Close netcdf file
        call io_check( nf90_close(ncid) )
        
        
        io_save_count = io_save_count + 1
    
    end subroutine io_save_state
    
    
    subroutine io_load_state(stateName)
        !
        ! loads the stateName
                
        character(12), intent(in) :: stateName 
         
        print*, ' loading '//stateName//' '
        
        call io_check( nf90_open(stateName, NF90_NOWRITE, ncid) )
        ! Get the varid s for u, v, w
        call io_check( nf90_inq_varid(ncid, "u", u_varid) )
        call io_check( nf90_inq_varid(ncid, "v", v_varid) )
        call io_check( nf90_inq_varid(ncid, "w", w_varid) )
        ! Read data
        call io_check( nf90_get_var(ncid, u_varid, u) )
        call io_check( nf90_get_var(ncid, v_varid, v) )
        call io_check( nf90_get_var(ncid, w_varid, w) )
                
        call io_check( nf90_close(ncid) )
        
    end subroutine io_load_state

    subroutine io_check(status)
        integer, intent ( in) :: status
    
        if(status /= nf90_noerr) then 
            print *, trim(nf90_strerror(status))
            stop "Stopped"
        end if
    end subroutine io_check  

    subroutine io_save_stats()
        !
        ! Measures some statistics and saves into text files
        !

        call state_kinetic()
        
        call io_courant()
        
        call state_dissipation()
        
        call state_divergence()

        open(io_stats, status='unknown', access='append', file='stats.dat')
        
        print *, ' Saving stats  '
        print *, 't = ', time(n+1)
        print *, 'Ekin = ', Ekin
        print *, 'Dissipation = ', Disp
        print *, 'Courant = ', Courant
        print *, 'Divergence = ', divMax
        print *, 'EZero = ', EZero
        print *, 'Eband = ', Eband
        write(io_Ekin,'(8e20.12)')  time(n+1), &     ! Instance
                                    Ekin, &          ! Total kinetic energ.
                                    Disp, &          ! Dissipation rate 
                                    2.0d0 * Q * Ekin, & ! Energy input rate
                                    Courant, &       ! maximum Courant num.
                                    divMax, &        ! maximum divergence
                                    EZero, &         ! Energy at k=0
                                    Eband            ! Energy of power-band
        
        close(io_stats)
        
    end subroutine io_save_stats

    
    subroutine io_Courant()
        
        !*****************!
        ! Courant number: !
        !*****************!
        Courant = maxval(abs(u)) * dt / (Lx / real(Nx, kind = 8)) &
                + maxval(abs(v)) * dt / (Ly / real(Ny, kind = 8)) &
                + maxval(abs(w)) * dt / (Lz / real(Nz, kind = 8))

    end subroutine io_Courant  

    subroutine io_save_spectrum()
        
        call state_spectrum()
       
        write(cnum,'(I4.4)') io_save_count
             
        print *, ' Saving spectrum'
        open(io_Spec, status='unknown', access='append', &
             file='spec'//cnum//'.dat')
            
        do nk = 1, Nspec
            
            write(io_Spec, '(2e20.12)') kSpec(nk), Espec(nk)
            
        end do   
        close(io_Spec)
              
    end subroutine io_save_spectrum
    
    subroutine io_save_info()
        ! Save parameters
        open(io_Pars, status='unknown', access='append', file='info.dat')
                
            write(io_Pars, '(''Nx,y,z = '', 3I4)') Nx, Ny, Nz
            write(io_Pars, '(''iSaveRate1 = '', 1I6)') iSaveRate1
            write(io_Pars, '(''iSaveRate2 = '', 1I6)' ) iSaveRate2
            write(io_Pars, '(''alpha_x,y,z = '', 3F12.9)') alpha_x, alpha_y, alpha_z
            write(io_Pars, '(''dt = '', F12.9)') dt
            write(io_Pars, '(''nu = '', F12.9)') nu
            if (bandlim) then
                write(io_Pars, '(''Pin = '', F12.9)') Pin
                write(io_Pars, '(''kF = '', F12.9)') kCutOff
            else
                write(io_Pars, '(''Q = '', F12.9)') Q
            end if
            write(io_Pars, '(''Deltak = '', F12.9)') Deltak
            write(io_Pars, '(''Courant = '', F12.9)') Courant
            if (initrand) then
                write(io_Pars, *) 'Random initial condition'
            else
                write(io_Pars, *) 'Initial condition read from file'
            end if
                
        close(io_Pars)
    end subroutine io_save_info
    
!**************************************************************************
 end module io
!**************************************************************************
