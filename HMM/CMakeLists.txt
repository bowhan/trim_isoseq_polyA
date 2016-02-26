cmake_minimum_required(VERSION 3.3)
project(HMM)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    char_traits.hpp
    fasta.hpp
    fastq.hpp
    format.hpp
    hmm_model.hpp
    hmm_utilities.h
    kernel_color.h
    main.cpp
    matrix.hpp
    polyA_hmm_model.hpp
    quality.hpp
    sequence.hpp
    type_policy.h)

find_package (Threads)
find_package(Boost COMPONENTS system program_options iostreams REQUIRED)
find_package(Zlib REQUIRED)
find_package(BZip2 REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})
add_executable(trim_isoseq_poly ${SOURCE_FILES})
target_link_libraries(trim_isoseq_poly ${ZLIB_LIBRARIES} ${BZIP2_LIBRARIES} ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})