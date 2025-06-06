#########################
# Important note: Minisat does not properly annotate its API with
# dllimport/dllexport, so the library build is forced to be static.
#########################
cmake_minimum_required(VERSION 3.5)
project(MiniSat VERSION 2.2 LANGUAGES CXX)

# Check if minisat is being used directly or via add_subdirectory, but allow overriding
if (NOT DEFINED MINISAT_MASTER_PROJECT)
    if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
        set(MINISAT_MASTER_PROJECT ON)
    else()
        set(MINISAT_MASTER_PROJECT OFF)
    endif()
endif ()

# options
option(MINISAT_BUILD_BINARIES "minisat: build binaries" ${MINISAT_MASTER_PROJECT})
option(MINISAT_BUILD_TESTING "Build and run MiniSat's tests" ON)
option(MINISAT_TEST_BENCHMARKS "Register benchmarks with CTest" OFF)
option(MINISAT_INSTALL "minisat: generate the install target" ON)


include(GNUInstallDirs)

add_library(libminisat STATIC
    # Impl files
    minisat/core/Solver.cc
    minisat/core/SolverTypes.cc
    minisat/utils/Options.cc
    minisat/utils/System.cc
    minisat/simp/SimpSolver.cc
    # Header files for IDEs
    minisat/core/Dimacs.h
    minisat/core/Solver.h
    minisat/core/SolverTypes.h
    minisat/mtl/Alg.h
    minisat/mtl/Alloc.h
    minisat/mtl/Heap.h
    minisat/mtl/IntTypes.h
    minisat/mtl/Map.h
    minisat/mtl/Queue.h
    minisat/mtl/Sort.h
    minisat/mtl/Vec.h
    minisat/mtl/XAlloc.h
    minisat/utils/Options.h
    minisat/utils/ParseUtils.h
    minisat/utils/System.h
    minisat/simp/SimpSolver.h
)

# Keep the library named as either libminisat.a or minisat.lib
# While having the target's logical name be distinct from minisat (the binary)
set_target_properties(libminisat
    PROPERTIES
      OUTPUT_NAME "minisat"
)


## Add a namespace alias.
# This is useful to abstract over use of the library as installed vs subdirectory build
add_library(MiniSat::libminisat ALIAS libminisat)

target_compile_features(libminisat
    PUBLIC
      cxx_attributes
      cxx_defaulted_functions
      cxx_deleted_functions
      cxx_final
)

target_include_directories(libminisat
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

set(targets libminisat)
if (MINISAT_BUILD_BINARIES OR (MINISAT_BUILD_TESTING AND BUILD_TESTING) OR (MINISAT_TEST_BENCHMARKS AND BUILD_TESTING))
    # Also build two MiniSat executables
    # NOTE: `minisat` is used in tests
    add_executable(minisat
        minisat/core/Main.cc
    )
    target_link_libraries(minisat libminisat)


    add_executable(minisat-simp
        minisat/simp/Main.cc
    )
    target_link_libraries(minisat-simp libminisat)
    list(APPEND targets minisat minisat-simp)
endif()

# Workaround for libstdc++ + Clang + -std=gnu++11 bug.
set_target_properties(${targets}
    PROPERTIES
      CXX_EXTENSIONS OFF
)


foreach (target ${targets})
    if ( CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang|GNU" )
        target_compile_options( ${target} PRIVATE -Wall -Wextra )
    endif()
    if ( CMAKE_CXX_COMPILER_ID MATCHES "MSVC" )
        target_compile_options( ${target} PRIVATE /W4 /wd4267 )
        target_compile_definitions( ${target} PRIVATE _CRT_SECURE_NO_WARNINGS )
    endif()
endforeach()

###############
# Testing
##
include(CTest)

if (MINISAT_BUILD_TESTING AND BUILD_TESTING)
    message(STATUS "Registering integration tests")
    # Read all easy instances from a file
    file(READ "${PROJECT_SOURCE_DIR}/tests/inputs/easy.txt" MINISAT_INTEGRATION_TESTS)
    string(REGEX REPLACE ";" "\\\\;" MINISAT_INTEGRATION_TESTS "${MINISAT_INTEGRATION_TESTS}")
    string(REGEX REPLACE "\n" ";" MINISAT_INTEGRATION_TESTS "${MINISAT_INTEGRATION_TESTS}")

    # Add 1 test for each easy instance
    foreach(INTEGRATION_TEST ${MINISAT_INTEGRATION_TESTS})
        add_test(NAME "integration:${INTEGRATION_TEST}"
            COMMAND minisat -verb=0 "tests/inputs/${INTEGRATION_TEST}"
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        )
        if ("${INTEGRATION_TEST}" MATCHES "^SAT")
            # The output can contain multiple lines, so we cannot do a full-line match
            # This means that we have to positive match SAT, but negative match UNSAT
            # because SAT is a substring of UNSAT.
            set_tests_properties("integration:${INTEGRATION_TEST}" PROPERTIES PASS_REGULAR_EXPRESSION "SAT\n"
                                                                              FAIL_REGULAR_EXPRESSION "UNSAT\n")
        else()
            set_tests_properties("integration:${INTEGRATION_TEST}" PROPERTIES PASS_REGULAR_EXPRESSION "SAT\n")
        endif()
        set_tests_properties("integration:${INTEGRATION_TEST}" PROPERTIES
            TIMEOUT 30
        ) # 30s timeout
    endforeach(INTEGRATION_TEST)
endif() # TESTING


if (MINISAT_TEST_BENCHMARKS AND BUILD_TESTING)
    message(STATUS "Registering benchmarks")

    file(READ "${PROJECT_SOURCE_DIR}/tests/inputs/benchmarks.txt" MINISAT_BENCHMARKS)
    string(REGEX REPLACE ";" "\\\\;" MINISAT_BENCHMARKS "${MINISAT_BENCHMARKS}")
    string(REGEX REPLACE "\n" ";" MINISAT_BENCHMARKS "${MINISAT_BENCHMARKS}")

    # Add 1 test for each easy instance
    foreach(BENCHMARK ${MINISAT_BENCHMARKS})
        add_test(NAME "benchmark:${BENCHMARK}"
            COMMAND minisat -verb=0 "tests/inputs/${BENCHMARK}"
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        )
        if ("${BENCHMARK}" MATCHES "^SAT")
            # The output can contain multiple lines, so we cannot do a full-line match
            # This means that we have to positive match SAT, but negative match UNSAT
            # because SAT is a substring of UNSAT.
            set_tests_properties("benchmark:${BENCHMARK}" PROPERTIES PASS_REGULAR_EXPRESSION "SAT\n"
                                                                     FAIL_REGULAR_EXPRESSION "UNSAT\n")
        else()
            set_tests_properties("benchmark:${BENCHMARK}" PROPERTIES PASS_REGULAR_EXPRESSION "SAT\n")
        endif()
        set_tests_properties("benchmark:${BENCHMARK}" PROPERTIES
            TIMEOUT 86400
        ) # 1 day timeout
    endforeach(BENCHMARK)
endif()


###############
# Installation
##
if (MINISAT_INSTALL)
    set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/MiniSat)

    install(
        TARGETS
          libminisat
        EXPORT MiniSatTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )

    if (MINISAT_BUILD_BINARIES)
        install(
            TARGETS
              minisat
              minisat-simp
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        )
    endif()

    install(DIRECTORY minisat/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/minisat FILES_MATCHING PATTERN "*.h*")

    install(EXPORT MiniSatTargets
        FILE MiniSatTargets.cmake
        NAMESPACE MiniSat::
        DESTINATION ${INSTALL_CONFIGDIR}
    )


    #####################
    # ConfigVersion file
    ##
    include(CMakePackageConfigHelpers)
    write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/MiniSatConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion
    )

    configure_package_config_file(
        ${CMAKE_CURRENT_LIST_DIR}/CMake/MiniSatConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/MiniSatConfig.cmake
        INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
    )

    ## Install all the helper files
    install(
        FILES
          ${CMAKE_CURRENT_BINARY_DIR}/MiniSatConfig.cmake
          ${CMAKE_CURRENT_BINARY_DIR}/MiniSatConfigVersion.cmake
        DESTINATION ${INSTALL_CONFIGDIR}
    )
endif()
