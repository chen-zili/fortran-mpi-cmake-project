program fortran_mpi
    use mpi
    implicit none

    integer(4) :: numtasks, taskid, len, ierr
    character(MPI_MAX_PROCESSOR_NAME) :: hostname

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)

    call MPI_GET_PROCESSOR_NAME(hostname, len, ierr)

    write(*, *) taskid, hostname
    if (taskid .eq. 0) then
        write(*, *) numtasks
    end if

    call MPI_FINALIZE(ierr)

end program fortran_mpi