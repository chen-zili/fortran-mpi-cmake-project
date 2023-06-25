program sendrecv_comm
    use mpi
    implicit none

    integer :: size, rank, ierr, status(MPI_STATUS_SIZE)
    integer :: send_buff, recv_buff, recv_count

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

    write(*, *) rank, "/", size

    if (0 == rank) then
        send_buff = 0
        call MPI_SEND(send_buff, 1, MPI_INTEGER, 1, 1, MPI_COMM_WORLD, ierr)
        write(*, *) "Send ", send_buff, " in ", rank

    else if (1 == rank) then
        call MPI_PROBE(0, 1, MPI_COMM_WORLD, status, ierr)
        call MPI_GET_COUNT(status, MPI_INTEGER, recv_count, ierr)
        write(*, *) "Recv count: ", recv_count

        call MPI_RECV(recv_buff, recv_count, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, status, ierr)

        write(*, *) "Recv ", recv_buff, " in ", rank
    end if

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    call MPI_FINALIZE(ierr)
end