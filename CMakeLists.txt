cmake_minimum_required (VERSION 3.4)
project (seqan3_tutorial CXX)

add_subdirectory(Kseq_Cpp)
# add seqan3 to search path
list (APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/seqan3/build_system")
list(APPEND CMAKE_PREFIX_PATH "C:/Users/cuixinyuan/Documents/develope/IEContig/zlib")

# require seqan3 with a version between >=3.0.0 and <4.0.0
find_package (seqan3 3.0 REQUIRED)
find_package(ZLIB REQUIRED)
find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ZLIB_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

# build app with seqan3
add_executable (hello_world hello_world.cpp)
target_link_libraries (hello_world seqan3::seqan3)


add_executable (IEContig main.cpp)
target_include_directories(IEContig PUBLIC zlib/include)
target_include_directories(IEContig PUBLIC Kseq_Cpp)
target_link_libraries(IEContig  Kseq_Cpp_LIB)
target_link_libraries (IEContig seqan3::seqan3)
#target_link_libraries(IEContig  ZLIB::ZLIB)
