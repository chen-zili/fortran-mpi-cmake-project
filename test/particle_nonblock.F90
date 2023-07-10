program particle_nonblock_comm
    use mpi
    implicit none

    type ParticleOne
        real(8) :: X, Y, Z, R, Vx, Vy, Vz, Ax, Ay, Az, WQ
    end type

    integer(4), parameter :: deps_col = 3
    real(8) :: xstart, xend, ystart, yend
    real(8) :: x_lb, x_ub, y_lb, y_ub
    real(8) :: R

    integer(4), parameter :: negb_type_boundary = 0
    integer(4), parameter :: negb_type_domain = 1

    integer(4) :: negb_type(8)
    integer(4) :: negb_rank(8)

    integer(4) :: reqs_send(8), reqs_recv(8)
    integer(4) :: status_send(MPI_STATUS_SIZE, 8), status_recv(MPI_STATUS_SIZE, 8)

    integer(4) :: recv_count
    integer(4) :: size, rank, ierr
    integer(4) :: row, col, i, j

    integer(4), parameter :: particle_number_max = 10000
    integer(4), parameter :: particle_number_max_bd = 80000
    integer(4) :: particle_number = 1000
    integer(4) :: particle_number_send(8)
    integer(4) :: particle_number_recv(8)
    type(ParticleOne) :: particle_bundle(particle_number_max)
    type(ParticleOne) :: particle_send_buff(particle_number_max_bd)     ! 不能用二维数组
    type(ParticleOne) :: particle_recv_buff(particle_number_max_bd)
    integer(4) :: index

    integer(4) :: mpi_type_particle_one

    character(len=99) :: filename

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

    ! 创建新的数据类型
    call MPI_TYPE_CONTIGUOUS(11, MPI_DOUBLE, mpi_type_particle_one, ierr)
    call MPI_TYPE_COMMIT(mpi_type_particle_one, ierr)

    ! 区域范围设置
    xstart = dble(col)
    xend = dble(col+1)
    ystart = dble(row)
    yend = dble(row+1)

    x_lb = 0.d0
    x_ub = dble(deps_col)
    y_lb = 0.d0
    y_ub = dble(size/deps_col)

    do i = 0, size-1
        if (rank == i) then
            write(*, '(i2, a, *(f10.4))') rank, ": ", xstart, xend, ystart, yend, x_lb, x_ub, y_lb, y_ub
        end if

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    end do

    ! 生成粒子
    do i = 1, particle_number
        call random_number(R)
        particle_bundle(i)%X = xstart - 1.d0 + (xend - xstart + 2.d0) * R

        call random_number(R)
        particle_bundle(i)%Y = ystart - 1.d0 + (yend - ystart + 2.d0) * R
    end do

    ! 边界处理
    do i = particle_number, 1, -1
        if (particle_bundle(i)%X <= x_lb .or. &
            particle_bundle(i)%X >= x_ub .or. &
            particle_bundle(i)%Y <= y_lb .or. &
            particle_bundle(i)%Y >= y_ub) then

            particle_bundle(i) = particle_bundle(particle_number)
            particle_number = particle_number - 1

        end if
    end do

    ! dump
    write(filename, '(i1)') rank
    open(10, file="raw_par_"//trim(filename)//".txt")
        do i = 1, particle_number
            write(10, '(*(f10.4, 1x))') particle_bundle(i)%X, particle_bundle(i)%Y
        end do
    close(10)

    ! 处理超过domain的粒子
    particle_number_send = 0
    do i = particle_number, 1, -1
        index = 0
        if (particle_bundle(i)%X < xstart) then
            if (particle_bundle(i)%Y < ystart) then
                index = 1
            else if (particle_bundle(i)%Y > yend) then
                index = 7
            else
                index = 8
            end if
        else if (particle_bundle(i)%X > xend) then
            if (particle_bundle(i)%Y < ystart) then
                index = 3
            else if (particle_bundle(i)%Y > yend) then
                index = 5
            else
                index = 4
            end if
        else
            if (particle_bundle(i)%Y < ystart) then
                index = 2
            else if (particle_bundle(i)%Y > yend) then
                index = 6
            else
                index = 0
            end if
        end if

        if (index > 0) then
            particle_number_send(index) = particle_number_send(index) + 1
            particle_send_buff((index-1)*particle_number_max + particle_number_send(index)) = particle_bundle(i)
            particle_bundle(i) = particle_bundle(particle_number)
            particle_number = particle_number - 1
        end if
    end do

    do i = 0, size-1
        if (rank == i) then
            write(*, '(i2, a, *(i6))') rank, ": ", particle_number_send
        end if

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    end do

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    ! dump
    write(filename, '(i1)') rank
    open(10, file="abort_par_"//trim(filename)//".txt")
        do i = 1, particle_number
            write(10, '(*(f10.4, 1x))') particle_bundle(i)%X, particle_bundle(i)%Y
        end do
    close(10)

    ! send and recv
    do i = 1, 8
        if (negb_type(i) == negb_type_domain) then
            call MPI_ISEND(particle_send_buff((i-1)*particle_number_max + 1: &
                            (i-1)*particle_number_max + particle_number_send(i)), &
                            particle_number_send(i), mpi_type_particle_one, negb_rank(i), &
                            1, MPI_COMM_WORLD, reqs_send(i), ierr)
            call MPI_IRECV(particle_recv_buff((i-1)*particle_number_max + 1: &
                            (i-1)*particle_number_max + particle_number_max), particle_number_max, &
                            mpi_type_particle_one, negb_rank(i), 1, MPI_COMM_WORLD, reqs_recv(i), ierr)
        end if
    end do

    ! wait
    do i = 1, 8
        if (negb_type(i) == negb_type_domain) then
            call MPI_WAIT(reqs_recv(i), status_recv(:, i), ierr)
        end if
    end do

    particle_number_recv = 0
    do i = 1, 8
        if (negb_type(i) == negb_type_domain) then
            call MPI_GET_COUNT(status_recv(:, i), mpi_type_particle_one, recv_count, ierr)
            if (recv_count > 0) particle_number_recv(i) = recv_count
        end if
    end do

    do i = 0, size-1
        if (rank == i) then
            write(*, '(i2, a, *(i6))') rank, ": ", particle_number_recv
        end if

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    end do

    do i = 1, 8
        if (negb_type(i) == negb_type_domain .and. particle_number_recv(i) > 0) then
            particle_bundle(particle_number+1:particle_number+particle_number_recv(i)) = &
            particle_recv_buff((i-1)*particle_number_max + 1:(i-1)*particle_number_max + particle_number_recv(i))
            particle_number = particle_number + particle_number_recv(i)
        end if
    end do

    ! dump
    write(filename, '(i1)') rank
    open(10, file="final_par_"//trim(filename)//".txt")
        do i = 1, particle_number
            write(10, '(*(f10.4, 1x))') particle_bundle(i)%X, particle_bundle(i)%Y
        end do
    close(10)

    call MPI_TYPE_FREE(mpi_type_particle_one, ierr)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call MPI_FINALIZE(ierr)
end