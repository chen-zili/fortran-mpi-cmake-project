program collective_nonblock_comm
    use mpi
    implicit none

    integer(4), parameter :: Nx = 5
    integer(4) :: size, rank, ierr, recv_count
    integer(4) :: array_1d(Nx)
    integer(4) :: orig_group, new_group
    integer(4), allocatable :: comms(:)
    integer(4) :: ranks(2)
    integer(4) :: status(MPI_STATUS_SIZE, 2), reqs(2)
    integer(4) :: i, tmp_send, tmp_recv, comms_index

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

    write(*, *) rank, "/", size

    call MPI_COMM_GROUP(MPI_COMM_WORLD, orig_group, ierr)

    allocate(comms(size-1))
    do i = 1, size-1
        ranks(1) = i-1
        ranks(2) = i

        call MPI_GROUP_INCL(orig_group, 2, ranks, new_group, ierr)
        call MPI_COMM_CREATE(MPI_COMM_WORLD, new_group, comms(i), ierr)
    end do

    array_1d = 1

    write(*, *) "Begin all reduce"

    do i = 1, size-1
        if (rank == i-1 .or. rank == i) then
            if (rank == i-1) then
                tmp_send = array_1d(Nx)
            else
                tmp_send = array_1d(1)
            end if

            call MPI_IALLREDUCE(tmp_send, tmp_recv, 1, MPI_INTEGER, MPI_SUM, comms(i), reqs(1), ierr)
            call MPI_WAITALL(1, reqs(1), status(:, 1), ierr)

            if (rank == i-1) then
                array_1d(Nx) = tmp_recv
            else
                array_1d(1)  = tmp_recv
            end if

        end if
    end do

    write(*, '(*(i4, 1x))') rank, array_1d

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    call MPI_FINALIZE(ierr)
end