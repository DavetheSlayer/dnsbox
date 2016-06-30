!**************************************************************************
!  IN/OUT 
!
!**************************************************************************
 module io
!**************************************************************************
    use variables
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
    
contains
    
    subroutine io_init()
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
                
        !**************************!
        !**************************!
        ! Generate .dat files here !
        !**************************!
        !**************************!
    end subroutine io_init
    
    subroutine io_saveState()
        !
        ! Saves the current state
        ! Based on:
        !
        ! https://www.hdfgroup.org/HDF5/Tutor/phypechk.html
        !
        character(4) :: cnum        ! Save count
                
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

!**************************************************************************
 end module io
!**************************************************************************
