include_directories(
        .
        ${Boost_INCLUDE_DIRS}
        ${ZLIB_INCLUDE_DIRS}
)

# source files
set(LIB_SOURCE_FILES
        char_traits.hpp
        fasta.hpp
        fastq.hpp
        format.hpp
        hmm_model.cpp
        hmm_model.hpp
        hmm_utilities.hpp
        kernel_color.hpp
        matrix.hpp
        polyA_hmm_model.cpp
        polyA_hmm_model.hpp
        quality.hpp
        sequence.hpp
        type_policy.hpp
        thread.hpp
        common.hpp
        common.cpp
        flnc.cpp
        flnc.hpp
        )

set(EXE_SOURCE_FILES
        ${LIB_SOURCE_FILES}
        main.cpp
        )

# CXX and LINKER FLAGS
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TrimIsoseqPolyA_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${TrimIsoseqPolyA_EXE_LINKER_FLAGS}")

add_library(trima ${LIB_SOURCE_FILES})
set_target_properties(trima PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY ${TrimIsoseqPolyA_LibDir}
        RUNTIME_OUTPUT_DIRECTORY ${TrimIsoseqPolyA_LibDir}
        LIBRARY_OUTPUT_DIRECTORY ${TrimIsoseqPolyA_LibDir}
        )

add_executable(trim_isoseq_polyA ${EXE_SOURCE_FILES})
set_target_properties(trim_isoseq_polyA PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${TrimIsoseqPolyA_BinDir}
        )

message(status "Boost_LIBRARIES: " ${Boost_LIBRARIES})

target_link_libraries(trim_isoseq_polyA
        ${Boost_LIBRARIES}
        ${ZLIB_LIBRARIES}
        ${BZIP2_LIBRARIES}
        Threads::Threads
        )
