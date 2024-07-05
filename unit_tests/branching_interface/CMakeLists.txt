# unit_tests/branching_interface/CMakeLists.txt

find_package(GTest REQUIRED)
find_package(Protobuf REQUIRED)
find_package(ortools REQUIRED)

# Add lpi_test executable
add_executable(branching_test branching_test.cc)
if(APPLE)
  set_target_properties(branching_test PROPERTIES
    INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
elseif(UNIX)
  set_target_properties(branching_test PROPERTIES
    INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()
target_link_libraries(branching_test PRIVATE
        branching_interface::branching_interface
        branching_interface::branching_factory
        solver
        custom_gtest_main
        utils
        absl::status
        GTest::gtest
        GTest::gmock
)
add_test(NAME branching_test COMMAND branching_test)
