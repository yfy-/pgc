cmake_minimum_required(VERSION 3.0)

project(PGC)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(MPI REQUIRED)
find_package(METIS REQUIRED)
include_directories(${MPI_CXX_INCLUDE_DIRS} ${METIS_INCLUDE_DIRS})

add_executable(pgc main.cpp)
target_link_libraries(pgc ${MPI_CXX_LIBRARIES} ${METIS_LIBRARIES})
