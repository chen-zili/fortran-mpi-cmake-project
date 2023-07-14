program fortran_mpi
    use mpi
    use ModuleFieldCommunication
    implicit none

    integer(4) :: size, rank, ierr, i
    type(FieldCom2D) :: mycom
    real(8), allocatable :: array(:, :)
    real(8), allocatable :: array_out(:, :)

    integer(4) :: lx = 4, ly = 5
    integer(4) :: px = 3, py = 3

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

    call mycom%init(lx, ly, px, py, com_type_nonblock)

    allocate(array(lx, ly))
    allocate(array_out(lx*px, ly*py))

    array = 1.d0
    call mycom%com(array, com_opt_sum, 1, lx, 1, ly)

    call mycom%gather(array, array_out, 1, lx, 1, ly, lx*px, ly*py)

    if (0 == rank) then
        open(10, file="./result_field_nonblock.txt")
            do i = 1, ly*py
                write(10, '(*(f10.4, 1x))') array_out(:, i)
            end do
        close(10)
    end if

    deallocate(array)
    deallocate(array_out)

    call mycom%destroy()
    call MPI_FINALIZE(ierr)

end program fortran_mpi