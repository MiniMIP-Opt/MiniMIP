# unit_tests/cutting_interface/CMakeLists.txt

# Add runner_test executable
add_executable(runner_test runner_test.cc)
if(APPLE)
  set_target_properties(runner_test PROPERTIES
    INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
elseif(UNIX)
  set_target_properties(runner_test PROPERTIES
    INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()
target_link_libraries(runner_test PRIVATE
        cutting_interface::cuts_runner
        solver
        custom_gtest_main
        utils
        absl::status
        GTest::gtest
        GTest::gmock
)
add_test(NAME runner_test COMMAND runner_test)

# Add aggregating_generator_test executable
add_executable(aggregating_generator_test aggregating_generator_test.cc)
if(APPLE)
  set_target_properties(aggregating_generator_test PROPERTIES
    INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
elseif(UNIX)
  set_target_properties(aggregating_generator_test PROPERTIES
    INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()
target_link_libraries(aggregating_generator_test PRIVATE
        cutting_interface::aggregating_generator
        custom_gtest_main
        utils
        solver
        absl::status
        GTest::gtest
        GTest::gmock
)
add_test(NAME aggregating_generator_test COMMAND aggregating_generator_test)

# Add hybrid_selector_test executable
add_executable(hybrid_selector_test hybrid_selector_test.cc)
if(APPLE)
  set_target_properties(hybrid_selector_test PROPERTIES
    INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
elseif(UNIX)
  set_target_properties(hybrid_selector_test PROPERTIES
    INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()
target_link_libraries(hybrid_selector_test PRIVATE
        cutting_interface::hybrid_selector
        solver
        custom_gtest_main
        utils
        absl::status
        GTest::gtest
        GTest::gmock
)
add_test(NAME hybrid_selector_test COMMAND hybrid_selector_test)
