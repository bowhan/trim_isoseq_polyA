file(MAKE_DIRECTORY ${TrimIsoseqPolyA_TestsDir}/bin)
file(MAKE_DIRECTORY ${TrimIsoseqPolyA_TestsDir}/out)

# Generate path for test data
configure_file(
    ${TrimIsoseqPolyA_TestsDir}/src/TestData.h.in
    ${TrimIsoseqPolyA_TestsDir}/src/TestData.h
)

include(files.cmake)
include_directories(
    ${TrimIsoseqPolyA_SourceDir}
    ${GTEST_SourceDir}
    ${Boost_INCLUDE_DIRS}
    ${ZLIB_INCLUDE_DIRS}
)

set(SOURCES
    ${TrimIsoseqPolyA_Test_H}
    ${TrimIsoseqPolyA_Test_CPP}
)

# quash gtest warnings
set(TrimIsoseqPolyA_CXX_FLAGS "${TrimIsoseqPolyA_CXX_FLAGS} -Wno-sign-compare")

# CXX and LINKER FLAGS
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TrimIsoseqPolyA_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${TrimIsoseqPolyA_EXE_LINKER_FLAGS}")

add_executable(test_trim_isoseq_polyA ${SOURCES})

set_target_properties(test_trim_isoseq_polyA PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${TrimIsoseqPolyA_TestsDir}/bin
)

if (NOT APPLE)
    set(MY_LIBRT -lrt)
else()
endif()

target_link_libraries(test_trim_isoseq_polyA
    trima
    libgtest
    libgmock
    ${Boost_LIBRARIES}
    ${ZLIB_LIBRARIES}
    ${BZIP2_LIBRARIES} 
    ${CMAKE_THREAD_LIBS_INIT} # quirky pthreads
    #${MY_LIBRT}
)

add_test(
    NAME UnitTests
    WORKING_DIRECTORY ${TrimIsoseqPolyA_TestsDir}/bin
    COMMAND test_trim_isoseq_polyA
)
