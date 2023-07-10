program field_nonblock__comm
    use mpi
    implicit none

    integer(4), parameter :: deps_col = 3
    integer(4), parameter :: Lx = 5
    integer(4), parameter :: Ly = 5
    integer(4) :: field(Lx, Ly)
    integer(4), allocatable :: all_field(:, :), field_out(:, :)

    integer(4), parameter :: negb_type_boundary = 0
    integer(4), parameter :: negb_type_domain = 1

    integer(4) :: negb_type(8)
    integer(4) :: negb_rank(8)

    integer(4) :: send_buff_edge_left(1:Ly)
    integer(4) :: send_buff_edge_right(1:Ly)
    integer(4) :: send_buff_edge_top(1:Lx)
    integer(4) :: send_buff_edge_bottom(1:Lx)
    integer(4) :: send_buff_corner(4)

    integer(4) :: recv_buff_edge_left(1:Ly)
    integer(4) :: recv_buff_edge_right(1:Ly)
    integer(4) :: recv_buff_edge_top(1:Lx)
    integer(4) :: recv_buff_edge_bottom(1:Lx)
    integer(4) :: recv_buff_corner(4)

    integer(4) :: reqs_left(2), reqs_right(2), reqs_top(2), reqs_bottom(2)
    integer(4) :: reqs_corner(8)

    integer(4) :: status_left(MPI_STATUS_SIZE, 2), status_right(MPI_STATUS_SIZE, 2)
    integer(4) :: status_top(MPI_STATUS_SIZE, 2), status_bottom(MPI_STATUS_SIZE, 2)
    integer(4) :: status_corner(MPI_STATUS_SIZE, 8)

    integer(4) :: recv_count, recv_count_corner(4)
    integer(4) :: size, rank, ierr
    integer(4) :: row, col, i, j

    ! 通信初始化
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

    ! 二维数组沿行和列分解，必须恰好整数
    if (mod(size, deps_col) /= 0) then
        write(*, *) "The domian decomposition is error."
        stop -1
    end if

    row = rank / deps_col
    col = mod(rank, deps_col)

    do i = 0, size-1
        if (rank == i) then
            write(*, '(a, i2, a, i2, a, i2, i2)') "Start: ", rank, "/", size, ", ", row, col
        end if

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    end do
    
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    ! domain邻近domain的rank和type
    negb_rank = -1
    negb_type = negb_type_boundary

    if (row - 1 >= 0 .and. col - 1 >= 0) then                       ! pos 1
        negb_rank(1) = (row-1) * deps_col + col - 1
    end if

    if (row - 1 >= 0) then                                          ! pos 2
        negb_rank(2) = (row-1) * deps_col + col
    end if

    if (row - 1 >= 0 .and. col + 1 < deps_col) then                 ! pos 3
        negb_rank(3) = (row-1) * deps_col + col + 1
    end if

    if (col + 1 < deps_col) then                                    ! pos 4
        negb_rank(4) = row * deps_col + col + 1
    end if

    if (row + 1 < size / deps_col .and. col + 1 < deps_col) then    ! pos 5
        negb_rank(5) = (row+1) * deps_col + col + 1
    end if

    if (row + 1 < size / deps_col) then                             ! pos 6
        negb_rank(6) = (row+1) * deps_col + col
    end if

    if (row + 1 < size / deps_col .and. col - 1 >= 0) then          ! pos 7
        negb_rank(7) = (row+1) * deps_col + col - 1
    end if

    if (col - 1 >= 0) then                                          ! pos 8
        negb_rank(8) = row * deps_col + col - 1
    end if

    do i = 1, 8
        if (negb_rank(i) >= 0) then
            negb_type(i) = negb_type_domain
        end if
    end do

    do i = 0, size-1
        if (rank == i) then
            write(*, '(i2, a, 8(i3), a, 8(i3))') rank, ": ", negb_rank, " | ", negb_type
        end if

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    end do

    ! 设置field和recv buffer
    field = 1
    send_buff_edge_left = field(1, 1:Ly)
    send_buff_edge_right = field(Lx, 1:Ly)
    send_buff_edge_top = field(1:Lx, 1)
    send_buff_edge_bottom = field(1:Lx, Ly)
    send_buff_corner(1) = field(1, 1)
    send_buff_corner(2) = field(Lx, 1)
    send_buff_corner(3) = field(Lx, Ly)
    send_buff_corner(4) = field(1, Ly)

    recv_buff_edge_left = 0
    recv_buff_edge_right = 0
    recv_buff_edge_top = 0
    recv_buff_edge_bottom = 0
    recv_buff_corner = 0

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    ! send and recv edge
    if (negb_type(8) == negb_type_domain) then  ! left
        call MPI_ISEND(send_buff_edge_left, Ly, MPI_INTEGER, negb_rank(8), 1, MPI_COMM_WORLD, reqs_left(1), ierr)
        call MPI_IRECV(recv_buff_edge_left, Ly, MPI_INTEGER, negb_rank(8), 1, MPI_COMM_WORLD, reqs_left(2), ierr)
    end if

    if (negb_type(4) == negb_type_domain) then  ! right
        call MPI_ISEND(send_buff_edge_right, Ly, MPI_INTEGER, negb_rank(4), 1, MPI_COMM_WORLD, reqs_right(1), ierr)
        call MPI_IRECV(recv_buff_edge_right, Ly, MPI_INTEGER, negb_rank(4), 1, MPI_COMM_WORLD, reqs_right(2), ierr)
    end if

    if (negb_type(2) == negb_type_domain) then  ! top
        call MPI_ISEND(send_buff_edge_top, Lx, MPI_INTEGER, negb_rank(2), 1, MPI_COMM_WORLD, reqs_top(1), ierr)
        call MPI_IRECV(recv_buff_edge_top, Lx, MPI_INTEGER, negb_rank(2), 1, MPI_COMM_WORLD, reqs_top(2), ierr)
    end if

    if (negb_type(6) == negb_type_domain) then  ! bottom
        call MPI_ISEND(send_buff_edge_bottom, Lx, MPI_INTEGER, negb_rank(6), 1, MPI_COMM_WORLD, reqs_bottom(1), ierr)
        call MPI_IRECV(recv_buff_edge_bottom, Lx, MPI_INTEGER, negb_rank(6), 1, MPI_COMM_WORLD, reqs_bottom(2), ierr)
    end if

    ! send and recv corner
    if (negb_type(1) == negb_type_domain) then ! corner left top
        call MPI_ISEND(send_buff_corner(1), 1, MPI_INTEGER, negb_rank(1), 1, MPI_COMM_WORLD, reqs_corner(1), ierr)
        call MPI_IRECV(recv_buff_corner(1), 1, MPI_INTEGER, negb_rank(1), 1, MPI_COMM_WORLD, reqs_corner(2), ierr)
    end if

    if (negb_type(3) == negb_type_domain) then ! corner right top
        call MPI_ISEND(send_buff_corner(2), 1, MPI_INTEGER, negb_rank(3), 1, MPI_COMM_WORLD, reqs_corner(3), ierr)
        call MPI_IRECV(recv_buff_corner(2), 1, MPI_INTEGER, negb_rank(3), 1, MPI_COMM_WORLD, reqs_corner(4), ierr)
    end if

    if (negb_type(5) == negb_type_domain) then ! corner right bottom
        call MPI_ISEND(send_buff_corner(3), 1, MPI_INTEGER, negb_rank(5), 1, MPI_COMM_WORLD, reqs_corner(5), ierr)
        call MPI_IRECV(recv_buff_corner(3), 1, MPI_INTEGER, negb_rank(5), 1, MPI_COMM_WORLD, reqs_corner(6), ierr)
    end if

    if (negb_type(7) == negb_type_domain) then ! corner left bottom
        call MPI_ISEND(send_buff_corner(4), 1, MPI_INTEGER, negb_rank(7), 1, MPI_COMM_WORLD, reqs_corner(7), ierr)
        call MPI_IRECV(recv_buff_corner(4), 1, MPI_INTEGER, negb_rank(7), 1, MPI_COMM_WORLD, reqs_corner(8), ierr)
    end if

    ! wait and accumulation
    recv_count_corner = 0
    if (negb_type(8) == negb_type_domain) then  ! left
        call MPI_WAITALL(2, reqs_left, status_left, ierr)
        field(1, 1:Ly) = field(1, 1:Ly) + recv_buff_edge_left
        recv_count_corner(1) = recv_count_corner(1) + 1
        recv_count_corner(4) = recv_count_corner(4) + 1
    end if

    if (negb_type(4) == negb_type_domain) then  ! right
        call MPI_WAITALL(2, reqs_right, status_right, ierr)
        field(Lx, 1:Ly) = field(Lx, 1:Ly) + recv_buff_edge_right
        recv_count_corner(2) = recv_count_corner(2) + 1
        recv_count_corner(3) = recv_count_corner(3) + 1
    end if

    if (negb_type(2) == negb_type_domain) then  ! top
        call MPI_WAITALL(2, reqs_top, status_top, ierr)
        field(1:Lx, 1) = field(1:Lx, 1) + recv_buff_edge_top
        recv_count_corner(1) = recv_count_corner(1) + 1
        recv_count_corner(2) = recv_count_corner(2) + 1
    end if

    if (negb_type(6) == negb_type_domain) then  ! bottom
        call MPI_WAITALL(2, reqs_bottom, status_bottom, ierr)
        field(1:Lx, Ly) = field(1:Lx, Ly) + recv_buff_edge_bottom
        recv_count_corner(3) = recv_count_corner(3) + 1
        recv_count_corner(4) = recv_count_corner(4) + 1
    end if

    if (negb_type(1) == negb_type_domain) then  ! left top
        call MPI_WAITALL(2, reqs_corner(1:2), status_corner(:, 1:2), ierr)
        field(1, 1) = field(1, 1) + recv_buff_corner(1)
        recv_count_corner(1) = recv_count_corner(1) + 1
    end if
    
    if (negb_type(3) == negb_type_domain) then  ! right top
        call MPI_WAITALL(2, reqs_corner(3:4), status_corner(:, 3:4), ierr)
        field(Lx, 1) = field(Lx, 1) + recv_buff_corner(2)
        recv_count_corner(2) = recv_count_corner(2) + 1
    end if

    if (negb_type(5) == negb_type_domain) then  ! right bottom
        call MPI_WAITALL(2, reqs_corner(5:6), status_corner(:, 5:6), ierr)
        field(Lx, Ly) = field(Lx, Ly) + recv_buff_corner(3)
        recv_count_corner(3) = recv_count_corner(3) + 1
    end if

    if (negb_type(7) == negb_type_domain) then  ! left bottom
        call MPI_WAITALL(2, reqs_corner(7:8), status_corner(:, 7:8), ierr)
        field(1, Ly) = field(1, Ly) + recv_buff_corner(4)
        recv_count_corner(4) = recv_count_corner(4) + 1
    end if

    ! out result
    allocate(all_field(Lx*deps_col, Ly*size/deps_col))
    allocate(field_out(Lx*deps_col, Ly*size/deps_col))

    all_field = 0
    all_field(col*Lx+1:(col+1)*Lx, row*Ly+1:(row+1)*Ly) = field
    
    call MPI_ALLREDUCE(all_field, field_out, Lx*Ly*size, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    if (0 == rank) then
        open(10, file="./result_field_nonblock.txt")
            do i = 1, Ly*size/deps_col
                write(10, '(*(i10, 1x))') field_out(:, i)
            end do
        close(10)
    end if

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call MPI_FINALIZE(ierr)
end