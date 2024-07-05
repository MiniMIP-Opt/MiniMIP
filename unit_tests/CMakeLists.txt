# unit_tests/CMakeLists.txt

# Add utils library
add_library(utils utils.cc utils.h)
if(APPLE)
  set_target_properties(utils PROPERTIES
    INSTALL_RPATH "@loader_path")
elseif(UNIX)
  set_target_properties(utils PROPERTIES
    INSTALL_RPATH "$ORIGIN")
endif()
target_link_libraries(utils PUBLIC
        data_structures::strong_sparse_vector
        lp_interface::lpi
        absl::base
        absl::status
        absl::flags
        absl::log
        GTest::gtest
        GTest::gmock
)

# Add gtest_main library
add_library(custom_gtest_main gtest_main.cc)
if(APPLE)
  set_target_properties(custom_gtest_main PROPERTIES
    INSTALL_RPATH "@loader_path")
elseif(UNIX)
  set_target_properties(custom_gtest_main PROPERTIES
    INSTALL_RPATH "$ORIGIN")
endif()
target_link_libraries(custom_gtest_main PRIVATE
        data_structures::strong_sparse_vector
        GTest::gtest
        ortools::ortools
)

include(GNUInstallDirs)
# Apply different linker flags based on platform
if(APPLE)
    set_target_properties(custom_gtest_main PROPERTIES
            LINK_FLAGS "-Wl,-all_load"
    )
endif()

# Add solve_loop_test executable
add_executable(solve_loop_test solve_loop_test.cc)
if(APPLE)
  set_target_properties(solve_loop_test PROPERTIES
    INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
elseif(UNIX)
  set_target_properties(solve_loop_test PROPERTIES
    INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()
target_link_libraries(solve_loop_test PRIVATE
        solver
        custom_gtest_main
        utils
        absl::status
        GTest::gtest
        GTest::gmock
)
add_test(NAME solve_loop_test COMMAND solve_loop_test)

# Add solver_test executable
add_executable(solver_test solver_test.cc)
if(APPLE)
  set_target_properties(solver_test PROPERTIES
    INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
elseif(UNIX)
  set_target_properties(solver_test PROPERTIES
    INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()
target_link_libraries(solver_test PRIVATE
        solver
        custom_gtest_main
        utils
        absl::status
        GTest::gtest
        GTest::gmock
)
add_test(NAME solver_test COMMAND solver_test)

# Include subdirectories
add_subdirectory(branching_interface)
add_subdirectory(cutting_interface)
add_subdirectory(data_structures)
add_subdirectory(lp_interface)
