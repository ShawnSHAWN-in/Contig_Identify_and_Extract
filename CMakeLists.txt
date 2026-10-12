cmake_minimum_required(VERSION 3.24)
project(IEContig VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(ZLIB REQUIRED)
find_package(OpenMP REQUIRED COMPONENTS CXX)
find_package(Threads REQUIRED)

option(IECONTIG_USE_SYSTEM_SEQAN3 "Use an installed SeqAn3 3.3.0 instead of downloading it" OFF)
if(IECONTIG_USE_SYSTEM_SEQAN3)
    find_package(seqan3 3.3.0 EXACT REQUIRED CONFIG)
else()
    include(FetchContent)
    FetchContent_Declare(seqan3
        GIT_REPOSITORY https://github.com/seqan/seqan3.git
        # SeqAn3 3.3.0; the argument_parser API is used by this application.
        GIT_TAG 6b681fb2eae5ab2997d293e99fc6a7f869a20316
        GIT_SUBMODULES submodules/cereal submodules/sdsl-lite
        GIT_SUBMODULES_RECURSE FALSE
        # Fetch sources, then use the upstream package configuration below.
        # This directory has no CMakeLists.txt, so upstream installation is not added.
        SOURCE_SUBDIR build_system
    )
    FetchContent_MakeAvailable(seqan3)
    find_package(seqan3 3.3.0 EXACT REQUIRED CONFIG
        PATHS "${seqan3_SOURCE_DIR}/build_system" NO_DEFAULT_PATH)
endif()

add_subdirectory(Kseq_Cpp)
add_executable(IEContig main.cpp)
target_link_libraries(IEContig PRIVATE
    Kseq_Cpp_LIB seqan3::seqan3 OpenMP::OpenMP_CXX Threads::Threads)

option(BUILD_TESTING "Enable CLI regression tests (requires Python 3)" OFF)
include(CTest)
if(BUILD_TESTING)
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    add_test(NAME cli_regression
        COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_cli.py"
                "$<TARGET_FILE:IEContig>" -v)
    set_tests_properties(cli_regression PROPERTIES TIMEOUT 120)
endif()

include(GNUInstallDirs)
install(TARGETS IEContig RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
if(NOT IECONTIG_USE_SYSTEM_SEQAN3)
    set(license_dir "${CMAKE_INSTALL_DATADIR}/doc/IEContig")
    install(FILES "${seqan3_SOURCE_DIR}/LICENSE.md"
        DESTINATION "${license_dir}" RENAME seqan3-LICENSE.md)
    install(FILES "${seqan3_SOURCE_DIR}/submodules/sdsl-lite/LICENSE"
        DESTINATION "${license_dir}" RENAME sdsl-LICENSE)
    install(FILES "${seqan3_SOURCE_DIR}/submodules/cereal/LICENSE"
        DESTINATION "${license_dir}" RENAME cereal-LICENSE)
    install(FILES "${seqan3_SOURCE_DIR}/include/seqan3/contrib/std/LICENSE"
        DESTINATION "${license_dir}" RENAME seqan3-std-LICENSE)
endif()
