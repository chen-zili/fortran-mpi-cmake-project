cmake_minimum_required(VERSION 3.1)

project(fortran-mpi LANGUAGES Fortran)

find_package(MPI)

set(CMAKE_BUILD_TYPE Release CACHE STRING "set build type to Release")

add_executable(${PROJECT_NAME} ${PROJECT_SOURCE_DIR}/src/main.F90)

target_link_libraries(${PROJECT_NAME} PUBLIC MPI::MPI_Fortran)

add_subdirectory(test)
