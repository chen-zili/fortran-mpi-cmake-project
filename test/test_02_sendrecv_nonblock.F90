program sendrecv_nonblock_comm
    use mpi
    implicit none

    integer(4) :: size, rank, ierr
    integer(4) :: status_left(MPI_STATUS_SIZE, 2), status_right(MPI_STATUS_SIZE, 2), reqs_left(2), reqs_right(2)
    integer(4) :: send_buff(2), recv_buff(2), recv_count
    integer(4) :: rank_left, rank_right

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

    write(*, *) rank, "/", size

    if (0 == rank) then
        rank_left = -1
    else
        rank_left = rank - 1
    end if

    if (size-1 == rank) then
        rank_right = -1
    else
        rank_right = rank + 1
    end if

    send_buff = rank

    if (rank_left /= -1) then
        call MPI_ISEND(send_buff(1), 1, MPI_INTEGER, rank_left, 1, MPI_COMM_WORLD, reqs_left(1), ierr)
        call MPI_IRECV(recv_buff(1), 1, MPI_INTEGER, rank_left, 1, MPI_COMM_WORLD, reqs_left(2), ierr)
    end if

    if (rank_right /= -1) then
        call MPI_ISEND(send_buff(2), 1, MPI_INTEGER, rank_right, 1, MPI_COMM_WORLD, reqs_right(1), ierr)
        call MPI_IRECV(recv_buff(2), 1, MPI_INTEGER, rank_right, 1, MPI_COMM_WORLD, reqs_right(2), ierr)
    end if

    if (rank_left /= -1) then
        call MPI_WAITALL(2, reqs_left, status_left, ierr)
        call MPI_GET_COUNT(status_left(:, 2), MPI_INTEGER, recv_count, ierr)
        write(*, *) "Recv count: ", recv_count
        write(*, *) "Recv ", recv_buff(1), " in ", rank
    end if

    if (rank_right /= -1) then
        call MPI_WAITALL(2, reqs_right, status_right, ierr)
        call MPI_GET_COUNT(status_right(:, 2), MPI_INTEGER, recv_count, ierr)
        write(*, *) "Recv count: ", recv_count
        write(*, *) "Recv ", recv_buff(2), " in ", rank
    end if

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    call MPI_FINALIZE(ierr)
end