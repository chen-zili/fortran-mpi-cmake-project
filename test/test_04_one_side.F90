program sendrecv_nonblock_comm
    use mpi
    implicit none

    integer(4), parameter :: Nx = 5
    integer(4) :: size, rank, ierr
    integer(4) :: array_1d(Nx)
    integer(4) :: send_buffer(2)
    integer(4) :: recv_buffer(2)
    integer(4) :: win
    integer(kind=MPI_ADDRESS_KIND) :: win_size, disp_aint

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

    write(*, *) rank, "/", size

    array_1d = 1
    send_buffer(1) = array_1d(1)
    send_buffer(2) = array_1d(Nx)

    win_size = 4 * 2
    call MPI_WIN_CREATE(send_buffer, win_size, 4, MPI_INFO_NULL, MPI_COMM_WORLD, win, ierr)

    if (0 /= rank) then
        disp_aint = 0
        call MPI_WIN_LOCK(MPI_LOCK_SHARED, rank-1, 0, win, ierr)
        call MPI_GET(recv_buffer(1), 1, MPI_INTEGER, rank-1, disp_aint, 1, MPI_INTEGER, win, ierr)
        call MPI_WIN_UNLOCK(rank-1, win, ierr)

        array_1d(1) = array_1d(1) + recv_buffer(1)
    end if

    if (size-1 /= rank) then
        disp_aint = 1
        call MPI_WIN_LOCK(MPI_LOCK_SHARED, rank+1, 0, win, ierr)
        call MPI_GET(recv_buffer(2), 1, MPI_INTEGER, rank+1, disp_aint, 1, MPI_INTEGER, win, ierr)
        call MPI_WIN_UNLOCK(rank+1, win, ierr)

        array_1d(Nx) = array_1d(Nx) + recv_buffer(2)
    end if

    write(*, '(*(i4, 1x))') rank, array_1d

    call MPI_WIN_FREE(win, ierr)
    call MPI_FINALIZE(ierr)
end