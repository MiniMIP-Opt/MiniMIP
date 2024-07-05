# unit_tests/lp_interface/CMakeLists.txt

# Add lpi_test executable
add_executable(lpi_test lpi_test.cc)
if(APPLE)
  set_target_properties(lpi_test PROPERTIES
    INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
elseif(UNIX)
  set_target_properties(lpi_test PROPERTIES
    INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()
target_link_libraries(lpi_test PRIVATE
        data_structures::strong_sparse_vector
        lp_interface::lpi
        lp_interface::lpi_factory
        custom_gtest_main
        utils
        absl::status
        GTest::gtest
        GTest::gmock
)
add_test(NAME lpi_test COMMAND lpi_test)
