# CTest
enable_testing()

set(nprocess 4)

# test 01 测试阻塞点对点通信
set(test_01 "test_01")
add_executable(${test_01} test_01_sendrecv.F90)
target_link_libraries(${test_01} PUBLIC MPI::MPI_Fortran)

add_test(NAME ${test_01} COMMAND "mpirun" "-np" "2" "./${test_01}")


# test 02 测试非阻塞点对点通信
set(test_02 "test_02")
add_executable(${test_02} test_02_sendrecv_nonblock.F90)
target_link_libraries(${test_02} PUBLIC MPI::MPI_Fortran)

add_test(NAME ${test_02} COMMAND "mpirun" "-np" "4" "./${test_02}")

# test 03 测试非阻塞集体通信
set(test_03 "test_03")
add_executable(${test_03} test_03_collective_nonblock.F90)
target_link_libraries(${test_03} PUBLIC MPI::MPI_Fortran)

add_test(NAME ${test_03} COMMAND "mpirun" "-np" "4" "./${test_03}")

# test 04 测试单边通信
set(test_04 "test_04")
add_executable(${test_04} test_04_one_side.F90)
target_link_libraries(${test_04} PUBLIC MPI::MPI_Fortran)

add_test(NAME ${test_04} COMMAND "mpirun" "-np" "4" "./${test_04}")