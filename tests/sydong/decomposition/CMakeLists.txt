cmake_minimum_required(VERSION 3.1...3.29)
project(
  ReadTensor
  VERSION 1.0
  LANGUAGES CXX C)
set(CMAKE_BUILD_TYPE Debug)
set(Eigen3_DIR "/home/sydong/lib/eigen_3.3.9/share/eigen3/cmake/")
set(ambit_DIR /home/sydong/bin/ambit-Release/share/cmake/ambit/)

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(HDF5 COMPONENTS CXX)
find_package(ambit 0.6 REQUIRED)

add_executable(decomposite decomp_final.cpp)
target_link_libraries(decomposite PRIVATE ambit::ambit)
target_link_libraries(decomposite PRIVATE HDF5::HDF5)


# add_executable(svd svd.cpp)
# target_link_libraries(svd PRIVATE Eigen3::Eigen)