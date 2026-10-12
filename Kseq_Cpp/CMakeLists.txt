add_library(Kseq_Cpp_LIB Kseq_Cpp.cpp)
target_include_directories(Kseq_Cpp_LIB PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(Kseq_Cpp_LIB PUBLIC ZLIB::ZLIB)
