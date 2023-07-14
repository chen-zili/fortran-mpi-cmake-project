module ModuleFieldCommunication
    use mpi
    implicit none


    integer(4), parameter :: negb_type_boundary = 0
    integer(4), parameter :: negb_type_domain   = 1

    integer(4), parameter :: negb_left_bottom   = 1
    integer(4), parameter :: negb_bottom        = 2
    integer(4), parameter :: negb_right_bottom  = 3
    integer(4), parameter :: negb_right         = 4
    integer(4), parameter :: negb_right_top     = 5
    integer(4), parameter :: negb_top           = 6
    integer(4), parameter :: negb_left_top      = 7
    integer(4), parameter :: negb_left          = 8

    integer(4), parameter :: negb_num           = 8
    integer(4), parameter :: corner_num         = 4

    integer(4), parameter :: com_type_nonblock  = 0
    integer(4), parameter :: com_type_onesided  = 1

    integer(4), parameter :: com_opt_sum        = 0
    integer(4), parameter :: com_opt_ext        = 1


    type :: FieldCom2D
        integer(4) :: lx, ly
        integer(4) :: dcol, drow
        integer(4) :: size, rank
        integer(4) :: com_type

        logical :: is_init = .False.
        integer(4) :: col, row

        ! for non-block com
        integer(4) :: reqs_left(2), reqs_right(2), reqs_top(2), reqs_bottom(2)
        integer(4) :: reqs_corner(8)

        integer(4) :: status_left(MPI_STATUS_SIZE, 2), status_right(MPI_STATUS_SIZE, 2)
        integer(4) :: status_top(MPI_STATUS_SIZE, 2), status_bottom(MPI_STATUS_SIZE, 2)
        integer(4) :: status_corner(MPI_STATUS_SIZE, 8)
        integer(4) :: recv_count, recv_count_corner(4)

        ! for one-sided com
        integer(4) :: win_left, win_right, win_top, win_bottom
        integer(kind=MPI_ADDRESS_KIND) :: win_size, disp_aint

        ! negb info
        integer(4) :: negb_type(negb_num)
        integer(4) :: negb_rank(negb_num)

        ! buffer
        real(8), allocatable :: send_buff_edge_left(:)
        real(8), allocatable :: send_buff_edge_right(:)
        real(8), allocatable :: send_buff_edge_bottom(:)
        real(8), allocatable :: send_buff_edge_top(:)
        real(8) :: send_buff_corner(corner_num)

        real(8), allocatable :: recv_buff_edge_left(:)
        real(8), allocatable :: recv_buff_edge_right(:)
        real(8), allocatable :: recv_buff_edge_bottom(:)
        real(8), allocatable :: recv_buff_edge_top(:)
        real(8) :: recv_buff_corner(corner_num)

    contains

        procedure :: init       => initFieldCom2D
        procedure :: destroy    => destroyFieldCom2D
        procedure :: com        => communicationFieldCom2D
        procedure :: gather     => gatherFieldCom2D

    end type FieldCom2D


    contains

        subroutine initFieldCom2D(this, x_local_size, y_local_size, deps_col, deps_row, com_type)
            class(FieldCom2D), intent(inout) :: this
            integer(4), intent(in) :: x_local_size, y_local_size, deps_col, deps_row
            integer(4), optional, intent(in) :: com_type
            integer(4) :: ierr, i

            call MPI_COMM_SIZE(MPI_COMM_WORLD, this%size, ierr)
            call MPI_COMM_RANK(MPI_COMM_WORLD, this%rank, ierr)

            if (x_local_size > 0 .and. y_local_size > 0 .and. deps_col > 0 .and. deps_row > 0) then
                this%lx = x_local_size
                this%ly = y_local_size
                this%dcol = deps_col
                this%drow = deps_row
                this%com_type = com_type_nonblock

                if (present(com_type)) then
                    if (com_type == com_type_nonblock) this%com_type = com_type_nonblock
                    if (com_type == com_type_onesided) this%com_type = com_type_onesided
                end if

                ! 二维数组沿行和列分解，必须恰好整数
                if (this%dcol * this%drow /= this%size) then
                    write(*, *) "The domian decomposition is error."
                    stop
                end if

                this%row = this%rank / this%dcol
                this%col = mod(this%rank, this%dcol)

                do i = 0, this%size-1
                    if (this%rank == i) then
                        write(*, '(a, i2, a, i2, a, i2, i2)') "Start: ", this%rank, "/", this%size, ", ", this%row, this%col
                    end if

                    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                end do
                
                ! 邻近domain的rank和type
                this%negb_rank = -1
                this%negb_type = negb_type_boundary

                if (this%row - 1 >= 0 .and. this%col - 1 >= 0) then                 ! pos 1
                    this%negb_rank(negb_left_bottom) = (this%row-1) * this%dcol + this%col - 1
                end if

                if (this%row - 1 >= 0) then                                         ! pos 2
                    this%negb_rank(negb_bottom) = (this%row-1) * this%dcol + this%col
                end if

                if (this%row - 1 >= 0 .and. this%col + 1 < this%dcol) then          ! pos 3
                    this%negb_rank(negb_right_bottom) = (this%row-1) * this%dcol + this%col + 1
                end if

                if (this%col + 1 < this%dcol) then                                  ! pos 4
                    this%negb_rank(negb_right) = this%row * this%dcol + this%col + 1
                end if

                if (this%row + 1 < this%drow .and. this%col + 1 < this%dcol) then   ! pos 5
                    this%negb_rank(negb_right_top) = (this%row+1) * this%dcol + this%col + 1
                end if

                if (this%row + 1 < this%drow) then                                  ! pos 6
                    this%negb_rank(negb_top) = (this%row+1) * this%dcol + this%col
                end if

                if (this%row + 1 < this%drow .and. this%col - 1 >= 0) then          ! pos 7
                    this%negb_rank(negb_left_top) = (this%row+1) * this%dcol + this%col - 1
                end if

                if (this%col - 1 >= 0) then                                         ! pos 8
                    this%negb_rank(negb_left) = this%row * this%dcol + this%col - 1
                end if

                do i = 1, negb_num
                    if (this%negb_rank(i) >= 0) then
                        this%negb_type(i) = negb_type_domain
                    end if
                end do

                do i = 0, this%size-1
                    if (this%rank == i) then
                        write(*, '(i2, a, 8(i3), a, 8(i3))') this%rank, ": ", this%negb_rank, " | ", this%negb_type
                    end if

                    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                end do

                call this%destroy()
                allocate(this%send_buff_edge_left(this%ly))
                allocate(this%send_buff_edge_right(this%ly))
                allocate(this%send_buff_edge_bottom(this%lx))
                allocate(this%send_buff_edge_top(this%lx))

                allocate(this%recv_buff_edge_left(this%ly))
                allocate(this%recv_buff_edge_right(this%ly))
                allocate(this%recv_buff_edge_bottom(this%lx))
                allocate(this%recv_buff_edge_top(this%lx))

                this%is_init = .True.

            else
                if (0 == this%rank) then
                    write(*, *) "The input parameters of the FieldCom2D are invalid."
                end if
                stop
            end if

        end subroutine initFieldCom2D


        subroutine destroyFieldCom2D(this)
            class(FieldCom2D), intent(inout) :: this

            if (allocated(this%send_buff_edge_left))    deallocate(this%send_buff_edge_left)
            if (allocated(this%send_buff_edge_right))   deallocate(this%send_buff_edge_right)
            if (allocated(this%send_buff_edge_bottom))  deallocate(this%send_buff_edge_bottom)
            if (allocated(this%send_buff_edge_top))     deallocate(this%send_buff_edge_top)

            if (allocated(this%recv_buff_edge_left))    deallocate(this%recv_buff_edge_left)
            if (allocated(this%recv_buff_edge_right))   deallocate(this%recv_buff_edge_right)
            if (allocated(this%recv_buff_edge_bottom))  deallocate(this%recv_buff_edge_bottom)
            if (allocated(this%recv_buff_edge_top))     deallocate(this%recv_buff_edge_top)

            this%is_init = .False.

        end subroutine destroyFieldCom2D


        subroutine communicationFieldCom2D(this, array2d, opt_type, xstart, xend, ystart, yend)
            class(FieldCom2D), intent(inout) :: this
            real(8), allocatable, intent(inout) :: array2d(:, :)
            integer(4), intent(in) :: opt_type
            integer(4), intent(in) :: xstart, xend, ystart, yend
            integer(4) :: ierr

            if (this%is_init) then
                if (xend-xstart+1 /= this%lx .or. yend-ystart+1 /= this%ly) then
                    if (0 == this%rank) then
                        write(*, *) "The input parameters (xstart, xend, ystart, yend) are invalid."
                    end if
                    stop
                end if

                if (com_opt_sum == opt_type) then
                    if (com_type_nonblock == this%com_type) then
                        this%send_buff_edge_left    = array2d(xstart, ystart:yend)
                        this%send_buff_edge_right   = array2d(xstart, ystart:yend)
                        this%send_buff_edge_bottom  = array2d(xstart:xend, ystart)
                        this%send_buff_edge_top     = array2d(xstart:xend, yend)

                        this%send_buff_corner(1)    = array2d(xstart, ystart)
                        this%send_buff_corner(2)    = array2d(xend, ystart)
                        this%send_buff_corner(3)    = array2d(xend, yend)
                        this%send_buff_corner(4)    = array2d(xstart, yend)

                        this%recv_buff_edge_left   = 0
                        this%recv_buff_edge_right  = 0
                        this%recv_buff_edge_bottom = 0
                        this%recv_buff_edge_top    = 0
                        this%recv_buff_corner      = 0

                        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

                        ! send and recv edge
                        if (this%negb_type(negb_left) == negb_type_domain) then         ! left
                            call MPI_ISEND(this%send_buff_edge_left, this%ly, MPI_DOUBLE, this%negb_rank(negb_left), 1, MPI_COMM_WORLD, this%reqs_left(1), ierr)
                            call MPI_IRECV(this%recv_buff_edge_left, this%ly, MPI_DOUBLE, this%negb_rank(negb_left), 1, MPI_COMM_WORLD, this%reqs_left(2), ierr)
                        end if

                        if (this%negb_type(negb_right) == negb_type_domain) then        ! right
                            call MPI_ISEND(this%send_buff_edge_right, this%ly, MPI_DOUBLE, this%negb_rank(negb_right), 1, MPI_COMM_WORLD, this%reqs_right(1), ierr)
                            call MPI_IRECV(this%recv_buff_edge_right, this%ly, MPI_DOUBLE, this%negb_rank(negb_right), 1, MPI_COMM_WORLD, this%reqs_right(2), ierr)
                        end if

                        if (this%negb_type(negb_bottom) == negb_type_domain) then       ! bottom
                            call MPI_ISEND(this%send_buff_edge_bottom, this%lx, MPI_DOUBLE, this%negb_rank(negb_bottom), 1, MPI_COMM_WORLD, this%reqs_bottom(1), ierr)
                            call MPI_IRECV(this%recv_buff_edge_bottom, this%lx, MPI_DOUBLE, this%negb_rank(negb_bottom), 1, MPI_COMM_WORLD, this%reqs_bottom(2), ierr)
                        end if

                        if (this%negb_type(negb_top) == negb_type_domain) then          ! top
                            call MPI_ISEND(this%send_buff_edge_top, this%lx, MPI_DOUBLE, this%negb_rank(negb_top), 1, MPI_COMM_WORLD, this%reqs_top(1), ierr)
                            call MPI_IRECV(this%recv_buff_edge_top, this%lx, MPI_DOUBLE, this%negb_rank(negb_top), 1, MPI_COMM_WORLD, this%reqs_top(2), ierr)
                        end if

                        ! send and recv corner
                        if (this%negb_type(negb_left_bottom) == negb_type_domain) then  ! corner left bottom
                            call MPI_ISEND(this%send_buff_corner(1), 1, MPI_DOUBLE, this%negb_rank(negb_left_bottom), 1, MPI_COMM_WORLD, this%reqs_corner(1), ierr)
                            call MPI_IRECV(this%recv_buff_corner(1), 1, MPI_DOUBLE, this%negb_rank(negb_left_bottom), 1, MPI_COMM_WORLD, this%reqs_corner(2), ierr)
                        end if

                        if (this%negb_type(negb_right_bottom) == negb_type_domain) then ! corner right bottom
                            call MPI_ISEND(this%send_buff_corner(2), 1, MPI_DOUBLE, this%negb_rank(negb_right_bottom), 1, MPI_COMM_WORLD, this%reqs_corner(3), ierr)
                            call MPI_IRECV(this%recv_buff_corner(2), 1, MPI_DOUBLE, this%negb_rank(negb_right_bottom), 1, MPI_COMM_WORLD, this%reqs_corner(4), ierr)
                        end if

                        if (this%negb_type(negb_right_top) == negb_type_domain) then    ! corner right top
                            call MPI_ISEND(this%send_buff_corner(3), 1, MPI_DOUBLE, this%negb_rank(negb_right_top), 1, MPI_COMM_WORLD, this%reqs_corner(5), ierr)
                            call MPI_IRECV(this%recv_buff_corner(3), 1, MPI_DOUBLE, this%negb_rank(negb_right_top), 1, MPI_COMM_WORLD, this%reqs_corner(6), ierr)
                        end if

                        if (this%negb_type(negb_left_top) == negb_type_domain) then     ! corner left top
                            call MPI_ISEND(this%send_buff_corner(4), 1, MPI_DOUBLE, this%negb_rank(negb_left_top), 1, MPI_COMM_WORLD, this%reqs_corner(7), ierr)
                            call MPI_IRECV(this%recv_buff_corner(4), 1, MPI_DOUBLE, this%negb_rank(negb_left_top), 1, MPI_COMM_WORLD, this%reqs_corner(8), ierr)
                        end if

                        ! wait and accumulation
                        this%recv_count_corner = 0
                        if (this%negb_type(negb_left) == negb_type_domain) then         ! left
                            call MPI_WAITALL(2, this%reqs_left, this%status_left, ierr)
                            array2d(xstart, ystart:yend) = array2d(xstart, ystart:yend) + this%recv_buff_edge_left
                            this%recv_count_corner(1) = this%recv_count_corner(1) + 1
                            this%recv_count_corner(4) = this%recv_count_corner(4) + 1
                        end if

                        if (this%negb_type(negb_right) == negb_type_domain) then        ! right
                            call MPI_WAITALL(2, this%reqs_right, this%status_right, ierr)
                            array2d(xend, ystart:yend) = array2d(xend, ystart:yend) + this%recv_buff_edge_right
                            this%recv_count_corner(2) = this%recv_count_corner(2) + 1
                            this%recv_count_corner(3) = this%recv_count_corner(3) + 1
                        end if

                        if (this%negb_type(negb_bottom) == negb_type_domain) then       ! bottom
                            call MPI_WAITALL(2, this%reqs_bottom, this%status_bottom, ierr)
                            array2d(xstart:xend, ystart) = array2d(xstart:xend, ystart) + this%recv_buff_edge_bottom
                            this%recv_count_corner(1) = this%recv_count_corner(1) + 1
                            this%recv_count_corner(2) = this%recv_count_corner(2) + 1
                        end if

                        if (this%negb_type(negb_top) == negb_type_domain) then          ! top
                            call MPI_WAITALL(2, this%reqs_top, this%status_top, ierr)
                            array2d(xstart:xend, yend) = array2d(xstart:xend, yend) + this%recv_buff_edge_top
                            this%recv_count_corner(4) = this%recv_count_corner(4) + 1
                            this%recv_count_corner(3) = this%recv_count_corner(3) + 1
                        end if

                        if (this%negb_type(negb_left_bottom) == negb_type_domain) then  ! left bottom
                            call MPI_WAITALL(2, this%reqs_corner(1:2), this%status_corner(:, 1:2), ierr)
                            array2d(xstart, ystart) = array2d(xstart, ystart) + this%recv_buff_corner(1)
                            this%recv_count_corner(1) = this%recv_count_corner(1) + 1
                        end if
                        
                        if (this%negb_type(negb_right_bottom) == negb_type_domain) then  ! right bottom
                            call MPI_WAITALL(2, this%reqs_corner(3:4), this%status_corner(:, 3:4), ierr)
                            array2d(xend, ystart) = array2d(xend, ystart) + this%recv_buff_corner(2)
                            this%recv_count_corner(2) = this%recv_count_corner(2) + 1
                        end if

                        if (this%negb_type(negb_right_top) == negb_type_domain) then    ! right top
                            call MPI_WAITALL(2, this%reqs_corner(5:6), this%status_corner(:, 5:6), ierr)
                            array2d(xend, yend) = array2d(xend, yend) + this%recv_buff_corner(3)
                            this%recv_count_corner(3) = this%recv_count_corner(3) + 1
                        end if

                        if (this%negb_type(negb_left_top) == negb_type_domain) then     ! left top
                            call MPI_WAITALL(2, this%reqs_corner(7:8), this%status_corner(:, 7:8), ierr)
                            array2d(xstart, yend) = array2d(xstart, yend) + this%recv_buff_corner(4)
                            this%recv_count_corner(4) = this%recv_count_corner(4) + 1
                        end if

                    else if (com_type_onesided == this%com_type) then

                    end if

                else if (com_opt_ext == opt_type) then
                    if (com_type_nonblock == this%com_type) then

                    else if (com_type_onesided == this%com_type) then

                    end if

                else
                    if (0 == this%rank) then
                        write(*, *) "The input parameter of the opt_type is invalid."
                    end if
                    stop
                end if

            else
                write(*, *) "The FieldCom2D is not initialized."
                stop
            end if

        end subroutine communicationFieldCom2D


        subroutine gatherFieldCom2D(this, array2d_local, array2d_global, xstart, xend, ystart, yend, x_global_size, y_global_size)
            class(FieldCom2D), intent(inout) :: this
            real(8), allocatable, intent(inout) :: array2d_local(:, :)
            real(8), allocatable, intent(inout) :: array2d_global(:, :)
            integer(4), intent(in) :: xstart, xend, ystart, yend, x_global_size, y_global_size
            real(8) :: array2d_tmp(x_global_size, y_global_size)
            integer(4) :: ierr

            if (this%is_init) then
                array2d_tmp = 0
                array2d_tmp(this%col*this%lx+1:(this%col+1)*this%lx, this%row*this%ly+1:(this%row+1)*this%ly) = array2d_local(xstart:xend, ystart:yend)
                
                call MPI_ALLREDUCE(array2d_tmp, array2d_global, this%lx*this%ly*this%size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)

            else
                write(*, *) "The FieldCom2D is not initialized."
                stop
            end if

        end subroutine gatherFieldCom2D

end module ModuleFieldCommunication