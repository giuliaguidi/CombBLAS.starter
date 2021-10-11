cmake_minimum_required(VERSION 3.11)
project(bridget)
set(CMAKE_CXX_STANDARD 14)

include_directories(
	  $ENV{COMBBLAS_HOME}/CombBLAS/install/include)
  
link_directories(
    $ENV{COMBBLAS_HOME}/CombBLAS/install/lib)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(DISCOVER_LIBC)
  find_library(libc c)
else()
  set(libc "")
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-std=c++11 -fpermissive -w -O3 -mavx2")

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

find_package(ZLIB)

add_executable(bridget
        src/main.cpp
        )
 
target_link_libraries(bridget PUBLIC ZLIB::ZLIB CombBLAS GraphGenlib Usortlib)

if (MPI_FOUND)
    include_directories(${MPI_INCLUDE_PATH})
    target_link_libraries(bridget PUBLIC ${MPI_LIBRARIES})
endif()

if(OPENMP_FOUND OR OpenMP_CXX_FOUND)
    # Set THREADED if OpenMP is found
    target_compile_definitions(bridget PUBLIC THREADED)
    if(TARGET OpenMP::OpenMP_CXX)
        target_link_libraries(bridget PUBLIC OpenMP::OpenMP_CXX)
    else()
        target_compile_options(bridget PUBLIC "${OpenMP_CXX_FLAGS}")
        target_link_libraries(bridget PUBLIC "${OpenMP_CXX_FLAGS}")
    endif()
endif()
