!**************************************************************************
!  IN/OUT 
!
!**************************************************************************
 module io
!**************************************************************************
    use state
    use rhs
    use hdf5
    
    integer(kind=8) :: io_iSaveCount1, io_iSaveCount2   ! save counters
    
    integer(HID_T) :: io_file_id       ! File identifier
    integer(HID_T) :: io_dset_id       ! Dataset identifier
    integer(HID_T) :: io_filespace     ! Dataspace identifier in file
    integer(HID_T) :: io_memspace      ! Dataspace identifier in memory
    integer(HID_T) :: io_plist_id      ! Property list identifier
    integer :: io_error                
    integer :: io_error_n         
    
    integer :: io_rank = 3     ! State file dataset rank
    integer(HSIZE_T), dimension(3) :: io_dimsf ! Dataset dimensions in the file
    integer(HSIZE_T), dimension(3) :: io_dimsfi
    integer(HSIZE_T), dimension(3) :: io_chunk_dims ! Chunks dimensions
    
    integer(HSIZE_T), dimension(3)  :: io_count
    integer(HSSIZE_T), dimension(3) :: io_offset
    integer(hsize_t), dimension(3)  :: io_stride
    integer(hsize_t), dimension(3)  :: io_block
    
    integer,     private  :: io_Ekin, io_Spec, io_Pars 
    character(4) :: cnum        ! Save count
    real(kind=8) :: scalekx   !scaling to avoid counting kx=0 modes twice
    
    ! Variables for reading parameters:
    integer(kind=4) :: io_Nx, io_Ny, io_Nz
    real(kind=8) :: io_alpha_x, io_alpha_y, io_alpha_z
    
contains
    
    subroutine io_init()
        character(10), save :: ss = 'unknown', aa = 'sequential'
        ! Initialize input/output
        io_iSaveCount1 = 0
        io_iSaveCount2 = 0        
        !
        ! Initialize FORTRAN predefined datatypes
        !
        call h5open_f(io_error)
        
        ! Set dataset dimensions:
        io_dimsf = (/ Nx, Ny, Nz /)
        io_dimsfi = (/ Nx, Ny, Nz /)
        
        ! Set chunk dimensions:
        io_chunk_dims = (/ iend(1) - istart(1) + 1, &
                           iend(2) - istart(2) + 1, & 
                           iend(3) - istart(3) + 1 /)
                
        ! Set the offsets for the processor:
        io_offset(1) = istart(1) - 1   ! x-offset
        io_offset(2) = istart(2) - 1   ! y-offset
        io_offset(3) = istart(3) - 1   ! z-offset
        
        ! Stride, count, and block for selecting hyperslab:
        io_stride(1) = 1
        io_stride(2) = 1
        io_stride(3) = 1
        io_count(1) = 1
        io_count(2) = 1
        io_count(3) = 1
        io_block(1) = io_chunk_dims(1)
        io_block(2) = io_chunk_dims(2)
        io_block(3) = io_chunk_dims(3)
        
        io_Ekin=20
        !**************************!
        ! Generate .dat files      !
        !**************************!
        if(proc_id.eq.0) then
            open(io_Ekin, status=ss, access=aa, file='Ekin.dat')
            aa = 'append'
        end if 
        
    end subroutine io_init
    
    subroutine io_finalize()
        ! Close files
    
        if(proc_id.eq.0) then
            close(io_Ekin)
        end if     
        if (proc_id.eq.0) then
            print *,'Saving the final state'
        end if        
        call io_saveState()
        
    end subroutine io_finalize
    
    subroutine io_saveState()
        !
        ! Saves the current state
        ! Based on:
        !
        ! https://www.hdfgroup.org/HDF5/Tutor/phypechk.html
        !

        write(cnum,'(I4.4)') io_iSaveCount1
        
        if(proc_id .eq. 0) then 
            print*, ' saving state'//cnum//'  t=', time(n+1)
        end if
        
        !
        ! Setup file access property list with parallel I/O access
        !
        call h5pcreate_f(H5P_FILE_ACCESS_F, io_plist_id, io_error)
        call h5pset_fapl_mpio_f(io_plist_id, mpi_comm_world, mpi_info_null, & 
                                io_error)
        !
        ! Create the file collectively
        !
        call h5fcreate_f('state'//cnum//'.h5', H5F_ACC_TRUNC_F, io_file_id, &
                         io_error, access_prp = io_plist_id)
        call h5pclose_f(io_plist_id, io_error)
        !
        ! Create the data space for the  dataset. 
        !        
        call h5screate_simple_f(io_rank, io_dimsf, io_filespace, io_error)
        call h5screate_simple_f(io_rank, io_chunk_dims, io_memspace, io_error)
        
        !
        ! Create chunked dataset
        ! 
        call h5pcreate_f(H5P_DATASET_CREATE_F, io_plist_id, io_error)
        call h5pset_chunk_f(io_plist_id, io_rank, io_chunk_dims, io_error)
        call h5dcreate_f (io_file_id, 'u', H5T_NATIVE_DOUBLE, io_filespace, &
                          io_dset_id, io_error, io_plist_id)
        call h5sclose_f(io_filespace, io_error)
        !
        ! Select hyperslab in the file  
        !
        call h5dget_space_f(io_dset_id, io_filespace, io_error)
        call h5sselect_hyperslab_f (io_filespace, H5S_SELECT_SET_F, &
                                    io_offset, io_count, io_error, &
                                    io_stride, io_block)
        !
        ! Create property list for collective dataset write
        !
        call h5pcreate_f (H5P_DATASET_XFER_F, io_plist_id, io_error)
        call h5pset_dxpl_mpio_f(io_plist_id, H5FD_MPIO_COLLECTIVE_F, io_error)
        !
        ! Write the dataset collectively
        !
        call h5dwrite_f (io_dset_id, H5T_NATIVE_DOUBLE, &
                  u(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &  
                  io_dimsfi, io_error, file_space_id = io_filespace, &
                  mem_space_id = io_memspace, xfer_prp = io_plist_id)
        !
        ! Close dataspaces
        !
        call h5sclose_f(io_filespace, io_error)
        call h5sclose_f(io_memspace, io_error)
        !
        ! Close the dataset
        !
        call h5dclose_f (io_dset_id, io_error)        

        !
        ! Close property list
        ! 
        call h5pclose_f(io_plist_id, io_error)
        
        ! REPEAT FOR v:
        !
        ! Create the data space for the  dataset. 
        !        
        call h5screate_simple_f(io_rank, io_dimsf, io_filespace, io_error)
        call h5screate_simple_f(io_rank, io_chunk_dims, io_memspace, io_error)
        
        !
        ! Create chunked dataset
        ! 
        call h5pcreate_f(H5P_DATASET_CREATE_F, io_plist_id, io_error)
        call h5pset_chunk_f(io_plist_id, io_rank, io_chunk_dims, io_error)
        call h5dcreate_f (io_file_id, 'v', H5T_NATIVE_DOUBLE, io_filespace, &
                          io_dset_id, io_error, io_plist_id)
        call h5sclose_f(io_filespace, io_error)
        !
        ! Select hyperslab in the file  
        !
        call h5dget_space_f(io_dset_id, io_filespace, io_error)
        call h5sselect_hyperslab_f (io_filespace, H5S_SELECT_SET_F, &
                                    io_offset, io_count, io_error, &
                                    io_stride, io_block)
        !
        ! Create property list for collective dataset write
        !
        call h5pcreate_f (H5P_DATASET_XFER_F, io_plist_id, io_error)
        call h5pset_dxpl_mpio_f(io_plist_id, H5FD_MPIO_COLLECTIVE_F, io_error)
        !
        ! Write the dataset collectively
        !
        call h5dwrite_f (io_dset_id, H5T_NATIVE_DOUBLE, &
                  v(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &  
                  io_dimsfi, io_error, file_space_id = io_filespace, &
                  mem_space_id = io_memspace, xfer_prp = io_plist_id)
        !
        ! Close dataspaces
        !
        call h5sclose_f(io_filespace, io_error)
        call h5sclose_f(io_memspace, io_error)
        !
        ! Close the dataset
        !
        call h5dclose_f (io_dset_id, io_error)        

        !
        ! Close property list
        ! 
        call h5pclose_f(io_plist_id, io_error)
        
        ! REPEAT FOR w:
        !
        ! Create the data space for the  dataset. 
        !        
        call h5screate_simple_f(io_rank, io_dimsf, io_filespace, io_error)
        call h5screate_simple_f(io_rank, io_chunk_dims, io_memspace, io_error)
        
        !
        ! Create chunked dataset
        ! 
        call h5pcreate_f(H5P_DATASET_CREATE_F, io_plist_id, io_error)
        call h5pset_chunk_f(io_plist_id, io_rank, io_chunk_dims, io_error)
        call h5dcreate_f (io_file_id, 'w', H5T_NATIVE_DOUBLE, io_filespace, &
                          io_dset_id, io_error, io_plist_id)
        call h5sclose_f(io_filespace, io_error)
        !
        ! Select hyperslab in the file  
        !
        call h5dget_space_f(io_dset_id, io_filespace, io_error)
        call h5sselect_hyperslab_f (io_filespace, H5S_SELECT_SET_F, &
                                    io_offset, io_count, io_error, &
                                    io_stride, io_block)
        !
        ! Create property list for collective dataset write
        !
        call h5pcreate_f (H5P_DATASET_XFER_F, io_plist_id, io_error)
        call h5pset_dxpl_mpio_f(io_plist_id, H5FD_MPIO_COLLECTIVE_F, io_error)
        !
        ! Write the dataset collectively
        !
        call h5dwrite_f (io_dset_id, H5T_NATIVE_DOUBLE, &
                  w(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &  
                  io_dimsfi, io_error, file_space_id = io_filespace, &
                  mem_space_id = io_memspace, xfer_prp = io_plist_id)
        !
        ! Close dataspaces
        !
        call h5sclose_f(io_filespace, io_error)
        call h5sclose_f(io_memspace, io_error)
        !
        ! Close the dataset
        !
        call h5dclose_f (io_dset_id, io_error)        

        !
        ! Close property list
        ! 
        call h5pclose_f(io_plist_id, io_error)
        !
        ! Close the file
        !
        call h5fclose_f(io_file_id, io_error)
        
        io_iSaveCount1 = io_iSaveCount1 + 1
    
    end subroutine io_saveState
    
    subroutine io_saveVorticity()
        !
        ! Saves the current state
        ! Based on:
        !
        ! https://www.hdfgroup.org/HDF5/Tutor/phypechk.html
        !
        
        if(proc_id .eq. 0) then 
            print*, 'saving vorticity'
        end if     
        !
        ! Setup file access property list with parallel I/O access
        !
        call h5pcreate_f(H5P_FILE_ACCESS_F, io_plist_id, io_error)
        call h5pset_fapl_mpio_f(io_plist_id, mpi_comm_world, mpi_info_null, & 
                                io_error)
        !
        ! Create the file collectively
        !
        call h5fcreate_f('vorticity.h5', H5F_ACC_TRUNC_F, io_file_id, &
                         io_error, access_prp = io_plist_id)
        call h5pclose_f(io_plist_id, io_error)
        !
        ! Create the data space for the  dataset. 
        !        
        call h5screate_simple_f(io_rank, io_dimsf, io_filespace, io_error)
        call h5screate_simple_f(io_rank, io_chunk_dims, io_memspace, io_error)
        
        !
        ! Create chunked dataset
        ! 
        call h5pcreate_f(H5P_DATASET_CREATE_F, io_plist_id, io_error)
        call h5pset_chunk_f(io_plist_id, io_rank, io_chunk_dims, io_error)
        call h5dcreate_f (io_file_id, 'omegax', H5T_NATIVE_DOUBLE, io_filespace, &
                          io_dset_id, io_error, io_plist_id)
        call h5sclose_f(io_filespace, io_error)
        !
        ! Select hyperslab in the file  
        !
        call h5dget_space_f(io_dset_id, io_filespace, io_error)
        call h5sselect_hyperslab_f (io_filespace, H5S_SELECT_SET_F, &
                                    io_offset, io_count, io_error, &
                                    io_stride, io_block)
        !
        ! Create property list for collective dataset write
        !
        call h5pcreate_f (H5P_DATASET_XFER_F, io_plist_id, io_error)
        call h5pset_dxpl_mpio_f(io_plist_id, H5FD_MPIO_COLLECTIVE_F, io_error)
        !
        ! Write the dataset collectively
        !
        call h5dwrite_f (io_dset_id, H5T_NATIVE_DOUBLE, &
                  omegax(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &  
                  io_dimsfi, io_error, file_space_id = io_filespace, &
                  mem_space_id = io_memspace, xfer_prp = io_plist_id)
        !
        ! Close dataspaces
        !
        call h5sclose_f(io_filespace, io_error)
        call h5sclose_f(io_memspace, io_error)
        !
        ! Close the dataset
        !
        call h5dclose_f (io_dset_id, io_error)        

        !
        ! Close property list
        ! 
        call h5pclose_f(io_plist_id, io_error)
        
        ! REPEAT FOR v:
        !
        ! Create the data space for the  dataset. 
        !        
        call h5screate_simple_f(io_rank, io_dimsf, io_filespace, io_error)
        call h5screate_simple_f(io_rank, io_chunk_dims, io_memspace, io_error)
        
        !
        ! Create chunked dataset
        ! 
        call h5pcreate_f(H5P_DATASET_CREATE_F, io_plist_id, io_error)
        call h5pset_chunk_f(io_plist_id, io_rank, io_chunk_dims, io_error)
        call h5dcreate_f (io_file_id, 'omegay', H5T_NATIVE_DOUBLE, io_filespace, &
                          io_dset_id, io_error, io_plist_id)
        call h5sclose_f(io_filespace, io_error)
        !
        ! Select hyperslab in the file  
        !
        call h5dget_space_f(io_dset_id, io_filespace, io_error)
        call h5sselect_hyperslab_f (io_filespace, H5S_SELECT_SET_F, &
                                    io_offset, io_count, io_error, &
                                    io_stride, io_block)
        !
        ! Create property list for collective dataset write
        !
        call h5pcreate_f (H5P_DATASET_XFER_F, io_plist_id, io_error)
        call h5pset_dxpl_mpio_f(io_plist_id, H5FD_MPIO_COLLECTIVE_F, io_error)
        !
        ! Write the dataset collectively
        !
        call h5dwrite_f (io_dset_id, H5T_NATIVE_DOUBLE, &
                  omegay(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &  
                  io_dimsfi, io_error, file_space_id = io_filespace, &
                  mem_space_id = io_memspace, xfer_prp = io_plist_id)
        !
        ! Close dataspaces
        !
        call h5sclose_f(io_filespace, io_error)
        call h5sclose_f(io_memspace, io_error)
        !
        ! Close the dataset
        !
        call h5dclose_f (io_dset_id, io_error)        

        !
        ! Close property list
        ! 
        call h5pclose_f(io_plist_id, io_error)
        
        ! REPEAT FOR w:
        !
        ! Create the data space for the  dataset. 
        !        
        call h5screate_simple_f(io_rank, io_dimsf, io_filespace, io_error)
        call h5screate_simple_f(io_rank, io_chunk_dims, io_memspace, io_error)
        
        !
        ! Create chunked dataset
        ! 
        call h5pcreate_f(H5P_DATASET_CREATE_F, io_plist_id, io_error)
        call h5pset_chunk_f(io_plist_id, io_rank, io_chunk_dims, io_error)
        call h5dcreate_f (io_file_id, 'omegaz', H5T_NATIVE_DOUBLE, io_filespace, &
                          io_dset_id, io_error, io_plist_id)
        call h5sclose_f(io_filespace, io_error)
        !
        ! Select hyperslab in the file  
        !
        call h5dget_space_f(io_dset_id, io_filespace, io_error)
        call h5sselect_hyperslab_f (io_filespace, H5S_SELECT_SET_F, &
                                    io_offset, io_count, io_error, &
                                    io_stride, io_block)
        !
        ! Create property list for collective dataset write
        !
        call h5pcreate_f (H5P_DATASET_XFER_F, io_plist_id, io_error)
        call h5pset_dxpl_mpio_f(io_plist_id, H5FD_MPIO_COLLECTIVE_F, io_error)
        !
        ! Write the dataset collectively
        !
        call h5dwrite_f (io_dset_id, H5T_NATIVE_DOUBLE, &
                  omegaz(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &  
                  io_dimsfi, io_error, file_space_id = io_filespace, &
                  mem_space_id = io_memspace, xfer_prp = io_plist_id)
        !
        ! Close dataspaces
        !
        call h5sclose_f(io_filespace, io_error)
        call h5sclose_f(io_memspace, io_error)
        !
        ! Close the dataset
        !
        call h5dclose_f (io_dset_id, io_error)        

        !
        ! Close property list
        ! 
        call h5pclose_f(io_plist_id, io_error)
        !
        ! Close the file
        !
        call h5fclose_f(io_file_id, io_error)
        
    end subroutine io_saveVorticity    
    
    subroutine io_loadState(stateName)
        !
        ! loads the stateName
        ! Based on:
                
        character(12), intent(in) :: stateName 
         
        if(proc_id .eq. 0) then 
            print*, ' loading '//stateName//' '
        end if
        
        !
        ! Setup file access property list with parallel I/O access
        !
        call h5pcreate_f(H5P_FILE_ACCESS_F, io_plist_id, io_error)
        call h5pset_fapl_mpio_f(io_plist_id, mpi_comm_world, mpi_info_null, & 
                                io_error)
        
        !
        ! Open the file collectively
        !
        call h5fopen_f(stateName, H5F_ACC_RDONLY_F, io_file_id, &
                       io_error, io_plist_id)
        call h5pclose_f(io_plist_id, io_error)
        !
        ! Create the data space for the  dataset. 
        !        
        call h5screate_simple_f(io_rank, io_dimsf, io_filespace, io_error)
        call h5screate_simple_f(io_rank, io_chunk_dims, io_memspace, io_error)
        
        !
        ! Open chunked dataset
        ! 
        call h5pcreate_f(H5P_DATASET_ACCESS_F, io_plist_id, io_error)
!        call h5pget_chunk_f(io_plist_id, io_rank, io_chunk_dims, io_error)
        call h5dopen_f (io_file_id, 'u', io_dset_id, io_error, io_plist_id)
        call h5sclose_f(io_filespace, io_error)
        !
        ! Select hyperslab in the file  
        !
        call h5dget_space_f(io_dset_id, io_filespace, io_error)
        call h5sselect_hyperslab_f (io_filespace, H5S_SELECT_SET_F, &
                                    io_offset, io_count, io_error, &
                                    io_stride, io_block)
        !
        ! Create property list for collective dataset read
        !
        call h5pcreate_f (H5P_DATASET_XFER_F, io_plist_id, io_error)
        call h5pset_dxpl_mpio_f(io_plist_id, H5FD_MPIO_COLLECTIVE_F, io_error)
        !
        ! Read the dataset collectively
        !
        call h5dread_f (io_dset_id, H5T_NATIVE_DOUBLE, &
                  u(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &  
                  io_dimsfi, io_error, file_space_id = io_filespace, &
                  mem_space_id = io_memspace, xfer_prp = io_plist_id)
        !
        ! Close dataspaces
        !
        call h5sclose_f(io_filespace, io_error)
        call h5sclose_f(io_memspace, io_error)
        !
        ! Close the dataset
        !
        call h5dclose_f (io_dset_id, io_error)        

        !
        ! Close property list
        ! 
        call h5pclose_f(io_plist_id, io_error)

        ! REPEAT V
        !
        ! Create the data space for the  dataset. 
        !        
        call h5screate_simple_f(io_rank, io_dimsf, io_filespace, io_error)
        call h5screate_simple_f(io_rank, io_chunk_dims, io_memspace, io_error)
        
        !
        ! Open chunked dataset
        ! 
        call h5pcreate_f(H5P_DATASET_ACCESS_F, io_plist_id, io_error)
!        call h5pget_chunk_f(io_plist_id, io_rank, io_chunk_dims, io_error)
        call h5dopen_f (io_file_id, 'v', io_dset_id, io_error, io_plist_id)
        call h5sclose_f(io_filespace, io_error)
        !
        ! Select hyperslab in the file  
        !
        call h5dget_space_f(io_dset_id, io_filespace, io_error)
        call h5sselect_hyperslab_f (io_filespace, H5S_SELECT_SET_F, &
                                    io_offset, io_count, io_error, &
                                    io_stride, io_block)
        !
        ! Create property list for collective dataset read
        !
        call h5pcreate_f (H5P_DATASET_XFER_F, io_plist_id, io_error)
        call h5pset_dxpl_mpio_f(io_plist_id, H5FD_MPIO_COLLECTIVE_F, io_error)
        !
        ! Read the dataset collectively
        !
        call h5dread_f (io_dset_id, H5T_NATIVE_DOUBLE, &
                  v(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &  
                  io_dimsfi, io_error, file_space_id = io_filespace, &
                  mem_space_id = io_memspace, xfer_prp = io_plist_id)
        !
        ! Close dataspaces
        !
        call h5sclose_f(io_filespace, io_error)
        call h5sclose_f(io_memspace, io_error)
        !
        ! Close the dataset
        !
        call h5dclose_f (io_dset_id, io_error)        

        !
        ! Close property list
        ! 
        call h5pclose_f(io_plist_id, io_error)

        ! REPEAT w
        !
        ! Create the data space for the  dataset. 
        !        
        call h5screate_simple_f(io_rank, io_dimsf, io_filespace, io_error)
        call h5screate_simple_f(io_rank, io_chunk_dims, io_memspace, io_error)
        
        !
        ! Open chunked dataset
        ! 
        call h5pcreate_f(H5P_DATASET_ACCESS_F, io_plist_id, io_error)
!        call h5pget_chunk_f(io_plist_id, io_rank, io_chunk_dims, io_error)
        call h5dopen_f (io_file_id, 'w', io_dset_id, io_error, io_plist_id)
        call h5sclose_f(io_filespace, io_error)
        !
        ! Select hyperslab in the file  
        !
        call h5dget_space_f(io_dset_id, io_filespace, io_error)
        call h5sselect_hyperslab_f (io_filespace, H5S_SELECT_SET_F, &
                                    io_offset, io_count, io_error, &
                                    io_stride, io_block)
        !
        ! Create property list for collective dataset read
        !
        call h5pcreate_f (H5P_DATASET_XFER_F, io_plist_id, io_error)
        call h5pset_dxpl_mpio_f(io_plist_id, H5FD_MPIO_COLLECTIVE_F, io_error)
        !
        ! Read the dataset collectively
        !
        call h5dread_f (io_dset_id, H5T_NATIVE_DOUBLE, &
                  w(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), &  
                  io_dimsfi, io_error, file_space_id = io_filespace, &
                  mem_space_id = io_memspace, xfer_prp = io_plist_id)
        !
        ! Close dataspaces
        !
        call h5sclose_f(io_filespace, io_error)
        call h5sclose_f(io_memspace, io_error)
        !
        ! Close the dataset
        !
        call h5dclose_f (io_dset_id, io_error)        

        !
        ! Close property list
        ! 
        call h5pclose_f(io_plist_id, io_error)
        !
        ! Close file
        ! 
        call h5fclose_f(io_file_id, io_error)
        
    end subroutine io_loadState

    subroutine io_saveStats()
        !
        ! Measures some statistics and saves into text files
        !
        call io_Ekinetic()
        
        call io_Courant()
        
        call io_Dissipation()

        if(proc_id.eq.0) then
            print *, ' Saving stats  '
            print *, 'Ekin = ', Ekin
            print *, 'Dissipation = ', Disp
            print *, 'Input = ', 2.0d0 * Q * Ekin
            print *, 'Courant = ', Courant
            write(io_Ekin,'(5e20.12)')  time(n+1), &     ! Instance
                                        Ekin, &          ! Total kinetic energ.
                                        Disp, &          ! Dissipation rate 
                                        2.0d0 * Q * Ekin, & ! Energy input rate
                                        Courant          ! maximum Courant num.
        end if 
    
    end subroutine io_saveStats
    
    subroutine io_Ekinetic()
        !*************************!
        ! Average kinetic energy: !
        !*************************!
        myEkin = 0d0
        do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
!            print *, k
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
                
    end subroutine io_Ekinetic
    
    subroutine io_Dissipation()
        
        !*******************!
        ! Dissipation rate: !
        !*******************!        
        call rhsDerivatives() ! Compute derivatives for dissipation calculation
        
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
        
    end subroutine io_Dissipation
    
    subroutine io_Courant()
        
        !*****************!
        ! Courant number: !
        !*****************!
        myCourant = maxval(abs(u)) * dt / (Lx / real(Nx, kind(0d0))) &
                  + maxval(abs(v)) * dt / (Ly / real(Ny, kind(0d0))) &
                  + maxval(abs(w)) * dt / (Lz / real(Nz, kind(0d0)))

        call mpi_allreduce(myCourant, Courant, 1, mpi_double_precision, &
                           mpi_max, mpi_comm_world, ierr)                     
        
    end subroutine io_Courant
    
    

    subroutine io_saveSpectrum()
        
        call stateComputeSpectrum()
       
        if(proc_id.eq.0) then
        
            write(cnum,'(I4.4)') io_iSaveCount1
             
            print *, ' Saving spectrum'
            open(io_Spec, status='unknown', access='append', &
                 file='spec'//cnum//'.dat')
            
            do nk = 1, Nspec
                
                write(io_Spec, '(2e20.12)') kSpec(nk), Espec(nk)
                
            end do   

            close(io_Spec)
        
        end if
            
    end subroutine io_saveSpectrum
    
    subroutine io_saveInfo()
        ! Save parameters
        if(proc_id.eq.0) then
            open(io_Pars, status='unknown', access='append', file='info.dat')
                
                write(io_Pars, '(''Nx,y,z = '', 3I4)') Nx, Ny, Nz
                write(io_Pars, '(''iSaveRate1 = '', 1I6)') iSaveRate1
                write(io_Pars, '(''iSaveRate2 = '', 1I6)' ) iSaveRate2
                write(io_Pars, '(''alpha_x,y,z = '', 3F12.9)') alpha_x, alpha_y, alpha_z
                write(io_Pars, '(''dt = '', F12.9)') dt
                write(io_Pars, '(''nu = '', F12.9)') nu
                write(io_Pars, '(''Q = '', F12.9)') Q
                write(io_Pars, '(''Deltak = '', F12.9)') Deltak
                write(io_Pars, '(''Courant = '', F12.9)') Courant
                if (initrand) then
                    write(io_Pars, *) 'Random initial condition'
                else
                    write(io_Pars, *) 'Initial condition read from file'
                end if
                
            close(io_Pars)
        end if
    end subroutine io_saveInfo
    
    subroutine io_loadInfo()
        
        logical :: infoExist
        character (len=200) :: dummyone
        character (len=200) :: dummytwo
        character (len=200) :: dummythree
        character (len=200) :: dummyfour
        character (len=200) :: dummyfive
        character (len=200) :: dummysix
        
        inquire(file='info.dat', exist=infoExist)
        if(.not. infoExist) then
            if(proc_id.eq.0) then
                print *, 'info.dat not found '
            end if
            call abort
        end if
        
        open(unit = 666, file='info.dat', status='old')
            
            read(666, *) dummyone, dummytwo, dummythree, dummyfour, &
                         io_Nx, io_Ny, io_Nz
            
            read(666, *)
            read(666, *)
            read(666, *) dummyone, dummytwo, dummythree, dummyfour, &
                         io_alpha_x, io_alpha_y, io_alpha_z
             
            if(proc_id.eq.0) then
                write(*, *) 'Nx = ', io_Nx  
                write(*, *) 'Ny = ', io_Ny  
                write(*, *) 'Nz = ', io_Nz  
            
                write(*, *) 'alpha_x = ', alpha_x  
                write(*, *) 'alpha_y = ', alpha_y  
                write(*, *) 'alpha_z = ', alpha_z  
            end if
                        
        close(666)
                
    end subroutine io_loadInfo
    
!**************************************************************************
 end module io
!**************************************************************************
