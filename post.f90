! Post processing 

program post
    
    use rhs
    use state
    use io
    
    character (len=12) :: fileName
    
    ! Initialize:
    call initMpi()
    call var_init()
    call io_init()
    
    ! Read info.dat:
    call io_loadInfo()    
    
    if(proc_id.eq.0) then
        print*, "enter file name:"
        read(*, *) fileName
    end if
    
    call mpi_bcast(fileName, 12, MPI_CHARACTER, 0, mpi_comm_world, ierr)

    if(proc_id.eq.0) then
        print*, "Loading the state file"
    end if    

    call io_loadState(fileName)

    if(proc_id.eq.0) then
        print*, "Done!"
    end if    

    call p3dfft_ftran_r2c (u, uhat, 'fft')
    call p3dfft_ftran_r2c (v, vhat, 'fft')
    call p3dfft_ftran_r2c (w, what, 'fft')

    ! copy uhat -> uhatold and uhattemp
    do k=fstart(3),fend(3); do j=fstart(2),fend(2); do i=fstart(1),fend(1)
        
        uhattemp(i, j, k) = uhat(i, j, k)
        vhattemp(i, j, k) = vhat(i, j, k)
        whattemp(i, j, k) = what(i, j, k)
        
        uhatold(i, j, k) = uhat(i, j, k)
        vhatold(i, j, k) = vhat(i, j, k)
        whatold(i, j, k) = what(i, j, k)
        
    end do; end do; end do
    
    call rhsDerivatives()
    call rhsVorticity()
    ! call io_saveState()
    call io_saveVorticity()
    
	deallocate(x,y,z,time,mychg,allchg,u,v,w,ux,uy,uz,vx,vy,vz,wx,wy,wz,&
               !uold,vold,wold,&
               utemp,vtemp,wtemp,&
               temp_r,kx,ky,kz,uhat,vhat,what,&
               uhattemp,vhattemp,whattemp,&
               uhatold,vhatold,whatold,&
               phat,nonlinuhat,nonlinvhat,nonlinwhat,temp_c,&
               intFact, stat=AllocateStatus)		
	if (AllocateStatus .ne. 0) stop
	if (proc_id.eq.0) then
		print *,'Program execution complete'
	end if
    
    ! Finalize
    call h5close_f(io_error)
    call p3dfft_clean
    call io_finalize()
    call MPI_FINALIZE (ierr)
    
                    
contains

subroutine initMpi()
    
    !initialize mpi:
    call mpi_init (ierr)
    call mpi_comm_size (mpi_comm_world, nproc, ierr)
    call mpi_comm_rank (mpi_comm_world, proc_id, ierr)

end subroutine initMpi

end program post
