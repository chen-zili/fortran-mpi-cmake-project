program fortran_mpi
    use mpi
    use ModulePICCommunication
    implicit none

    integer(4) :: size, rank, ierr, i, j
    type(PICCom2D) :: mycom
    real(8), allocatable :: array(:, :), array_ext(:, :)
    real(8), allocatable :: array_out(:, :)
    character(len=99) :: filename

    integer(4) :: lx = 4, ly = 5
    integer(4) :: px = 3, py = 3

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

    call mycom%init(lx, ly, px, py, com_type_nonblock)

    allocate(array(lx, ly))
    allocate(array_ext(0:lx+1, 0:ly+1))
    allocate(array_out(lx*px, ly*py))

    do i = 1, lx
        do j = 1, ly
            array(i, j) = rank + i + j
        end do
    end do

    array_ext = -1.d0
    array_ext(1:lx, 1:ly) = array(1:lx, 1:ly)

    mycom%left_ext_type = com_field_ext_type_symmetry

    call mycom%comf(array, com_field_opt_sum, 1, lx, 1, ly)
    call mycom%comf(array_ext, com_field_opt_ext, 1, lx, 1, ly)

    call mycom%gather(array, array_out, 1, lx, 1, ly, lx*px, ly*py)

    if (0 == rank) then
        open(10, file="./result_field_nonblock.txt")
            do i = 1, ly*py
                write(10, '(*(f10.4, 1x))') array_out(:, i)
            end do
        close(10)
    end if

    write(filename, '(i2)') rank
    open(10, file="./result_field_nonblock_"//trim(filename)//".txt")
        do i = 0, ly+1
            write(10, '(*(f10.4, 1x))') array_ext(:, i)
        end do
    close(10)

    deallocate(array)
    deallocate(array_ext)
    deallocate(array_out)

    call mycom%destroy()
    call MPI_FINALIZE(ierr)

end program fortran_mpi