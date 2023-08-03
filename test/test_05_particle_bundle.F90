program fortran_mpi
    use mpi
    use ModulePICCommunication
    implicit none

    type(ParticleBundle) :: pb, b2
    type(ParticleOne) :: one
    integer(4) :: i

    ! add bun
    ! call pb%init(100, 10, 1000, 0.25d0, 0.5d0, 1.5d0)
    ! call b2%init(990, 10, 10000)
    ! pb%npar = 50
    ! b2%npar = 950
    ! call pb%addbun(b2)

    ! add one
    ! call pb%init(100, 10, 1000, 0.25d0, 0.5d0, 1.5d0)
    ! do i = 1, 10001
    !     call pb%addone(one)
    !     write(*, *) pb%npar, pb%size
    ! end do

    ! del one
    call pb%init(1000, 10, 1000, 0.25d0, 0.5d0, 1.5d0)
    pb%npar = 1000
    do i = 1000, 1, -1
        call pb%delone(i)
        write(*, *) pb%npar, pb%size
    end do

    call pb%destroy()
    call b2%destroy()
end