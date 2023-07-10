program field_one_side_comm
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

    integer(4) :: recv_buff_edge_left(1:Ly)
    integer(4) :: recv_buff_edge_right(1:Ly)
    integer(4) :: recv_buff_edge_top(1:Lx)
    integer(4) :: recv_buff_edge_bottom(1:Lx)
    integer(4) :: recv_buff_corner(4)

    integer(4) :: win_left, win_right, win_top, win_bottom
    integer(kind=MPI_ADDRESS_KIND) :: win_size, disp_aint

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

    ! 创建远程访问窗口
    win_size = 4 * Ly
    call MPI_WIN_CREATE(send_buff_edge_left, win_size, 4, MPI_INFO_NULL, MPI_COMM_WORLD, win_left, ierr)

    win_size = 4 * Ly
    call MPI_WIN_CREATE(send_buff_edge_right, win_size, 4, MPI_INFO_NULL, MPI_COMM_WORLD, win_right, ierr)

    win_size = 4 * Lx
    call MPI_WIN_CREATE(send_buff_edge_top, win_size, 4, MPI_INFO_NULL, MPI_COMM_WORLD, win_top, ierr)

    win_size = 4 * Lx
    call MPI_WIN_CREATE(send_buff_edge_bottom, win_size, 4, MPI_INFO_NULL, MPI_COMM_WORLD, win_bottom, ierr)

    ! 设置field和recv buffer
    field = 1
    send_buff_edge_left = field(1, 1:Ly)
    send_buff_edge_right = field(Lx, 1:Ly)
    send_buff_edge_top = field(1:Lx, 1)
    send_buff_edge_bottom = field(1:Lx, Ly)

    recv_buff_edge_left = 0
    recv_buff_edge_right = 0
    recv_buff_edge_top = 0
    recv_buff_edge_bottom = 0
    recv_buff_corner = 0

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    ! send and recv edge
    if (negb_type(8) == negb_type_domain) then  ! left
        disp_aint = 0
        call MPI_WIN_LOCK(MPI_LOCK_SHARED, negb_type(8), 0, win_left, ierr)
        call MPI_GET(recv_buff_edge_left, Ly, MPI_INTEGER, negb_type(8), disp_aint, Ly, MPI_INTEGER, win_left, ierr)
        call MPI_WIN_UNLOCK(negb_type(8), win_left, ierr)

        field(1, 1:Ly) = field(1, 1:Ly) + recv_buff_edge_left
    end if

    if (negb_type(4) == negb_type_domain) then  ! right
        disp_aint = 0
        call MPI_WIN_LOCK(MPI_LOCK_SHARED, negb_type(4), 0, win_right, ierr)
        call MPI_GET(recv_buff_edge_right, Ly, MPI_INTEGER, negb_type(4), disp_aint, Ly, MPI_INTEGER, win_right, ierr)
        call MPI_WIN_UNLOCK(negb_type(4), win_right, ierr)

        field(Lx, 1:Ly) = field(Lx, 1:Ly) + recv_buff_edge_right
    end if

    if (negb_type(2) == negb_type_domain) then  ! top
        disp_aint = 0
        call MPI_WIN_LOCK(MPI_LOCK_SHARED, negb_type(2), 0, win_top, ierr)
        call MPI_GET(recv_buff_edge_top, Lx, MPI_INTEGER, negb_type(2), disp_aint, Lx, MPI_INTEGER, win_top, ierr)
        call MPI_WIN_UNLOCK(negb_type(2), win_top, ierr)

        field(1:Lx, 1) = field(1:Lx, 1) + recv_buff_edge_top
    end if

    if (negb_type(6) == negb_type_domain) then  ! bottom
        disp_aint = 0
        call MPI_WIN_LOCK(MPI_LOCK_SHARED, negb_type(6), 0, win_bottom, ierr)
        call MPI_GET(recv_buff_edge_bottom, Lx, MPI_INTEGER, negb_type(6), disp_aint, Lx, MPI_INTEGER, win_bottom, ierr)
        call MPI_WIN_UNLOCK(negb_type(6), win_bottom, ierr)

        field(1:Lx, Ly) = field(1:Lx, Ly) + recv_buff_edge_bottom
    end if

    ! send and recv corner
    if (negb_type(1) == negb_type_domain) then ! corner left top
        disp_aint = Lx-1
        call MPI_WIN_LOCK(MPI_LOCK_SHARED, negb_type(1), 0, win_bottom, ierr)
        call MPI_GET(recv_buff_corner(1), 1, MPI_INTEGER, negb_type(1), disp_aint, 1, MPI_INTEGER, win_bottom, ierr)
        call MPI_WIN_UNLOCK(negb_type(1), win_bottom, ierr)

        field(1, 1) = field(1, 1) + recv_buff_corner(1)
    end if

    if (negb_type(3) == negb_type_domain) then ! corner right top
        disp_aint = 0
        call MPI_WIN_LOCK(MPI_LOCK_SHARED, negb_type(3), 0, win_bottom, ierr)
        call MPI_GET(recv_buff_corner(2), 1, MPI_INTEGER, negb_type(3), disp_aint, 1, MPI_INTEGER, win_bottom, ierr)
        call MPI_WIN_UNLOCK(negb_type(3), win_bottom, ierr)

        field(Lx, 1) = field(Lx, 1) + recv_buff_corner(2)
    end if

    if (negb_type(5) == negb_type_domain) then ! corner right bottom
        disp_aint = 0
        call MPI_WIN_LOCK(MPI_LOCK_SHARED, negb_type(5), 0, win_top, ierr)
        call MPI_GET(recv_buff_corner(3), 1, MPI_INTEGER, negb_type(5), disp_aint, 1, MPI_INTEGER, win_top, ierr)
        call MPI_WIN_UNLOCK(negb_type(5), win_top, ierr)

        field(Lx, Ly) = field(Lx, Ly) + recv_buff_corner(3)
    end if

    if (negb_type(7) == negb_type_domain) then ! corner left bottom
        disp_aint = Lx-1
        call MPI_WIN_LOCK(MPI_LOCK_SHARED, negb_type(7), 0, win_top, ierr)
        call MPI_GET(recv_buff_corner(4), 1, MPI_INTEGER, negb_type(7), disp_aint, 1, MPI_INTEGER, win_top, ierr)
        call MPI_WIN_UNLOCK(negb_type(7), win_top, ierr)

        field(1, Ly) = field(1, Ly) + recv_buff_corner(4)
    end if

    call MPI_WIN_FREE(win_left, ierr)
    call MPI_WIN_FREE(win_right, ierr)
    call MPI_WIN_FREE(win_top, ierr)
    call MPI_WIN_FREE(win_bottom, ierr)

    ! out result
    allocate(all_field(Lx*deps_col, Ly*size/deps_col))
    allocate(field_out(Lx*deps_col, Ly*size/deps_col))

    all_field = 0
    all_field(col*Lx+1:(col+1)*Lx, row*Ly+1:(row+1)*Ly) = field
    
    call MPI_ALLREDUCE(all_field, field_out, Lx*Ly*size, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    if (0 == rank) then
        open(10, file="./result_field_one_side.txt")
            do i = 1, Ly*size/deps_col
                write(10, '(*(i10, 1x))') field_out(:, i)
            end do
        close(10)
    end if

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call MPI_FINALIZE(ierr)
end