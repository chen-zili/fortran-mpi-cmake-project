rm -rf build
rm -rf run

export FC=ifort
export F9X=ifort

cmake -S . -B build
cmake --build build -j 8

# ctest --test-dir build/test
