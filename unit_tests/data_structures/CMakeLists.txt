# unit_tests/data_structures/CMakeLists.txt

# Add problem_test executable
add_executable(problem_test problem_test.cc)
if(APPLE)
  set_target_properties(problem_test PROPERTIES
    INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
elseif(UNIX)
  set_target_properties(problem_test PROPERTIES
    INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()
target_link_libraries(problem_test PRIVATE
        data_structures::problem
        custom_gtest_main
        utils
        absl::status
        GTest::gtest
        GTest::gmock
)
add_test(NAME problem_test COMMAND problem_test)

# Add mip_data_test executable
add_executable(mip_data_test mip_data_test.cc)
if(APPLE)
  set_target_properties(mip_data_test PROPERTIES
    INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
elseif(UNIX)
  set_target_properties(mip_data_test PROPERTIES
    INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()
target_link_libraries(mip_data_test PRIVATE
        data_structures::mip_data
        custom_gtest_main
        utils
        absl::status
        GTest::gtest
        GTest::gmock
)
add_test(NAME mip_data_test COMMAND mip_data_test)

# Add mip_tree_test executable
add_executable(mip_tree_test mip_tree_test.cc)
if(APPLE)
  set_target_properties(mip_tree_test PROPERTIES
    INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
elseif(UNIX)
  set_target_properties(mip_tree_test PROPERTIES
    INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()
target_link_libraries(mip_tree_test PRIVATE
        data_structures::mip_tree
        data_structures::strong_sparse_vector
        custom_gtest_main
        utils
        absl::status
        GTest::gtest
        GTest::gmock
)
add_test(NAME mip_tree_test COMMAND mip_tree_test)

# Add cuts_data_test executable
add_executable(cuts_data_test cuts_data_test.cc)
if(APPLE)
  set_target_properties(cuts_data_test PROPERTIES
    INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
elseif(UNIX)
  set_target_properties(cuts_data_test PROPERTIES
    INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()
target_link_libraries(cuts_data_test PRIVATE
        data_structures::cuts_data
        custom_gtest_main
        utils
        absl::status
        GTest::gtest
        GTest::gmock
)
add_test(NAME cuts_data_test COMMAND cuts_data_test)

# Add strong_sparse_vector_test executable
add_executable(strong_sparse_vector_test strong_sparse_vector_test.cc)
if(APPLE)
  set_target_properties(strong_sparse_vector_test PROPERTIES
    INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
elseif(UNIX)
  set_target_properties(strong_sparse_vector_test PROPERTIES
    INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()
target_link_libraries(strong_sparse_vector_test PRIVATE
        data_structures::strong_sparse_vector
        custom_gtest_main
        utils
        GTest::gtest
        GTest::gmock
        ortools::ortools
        # ortools::base
        # ortools::strong_vector
        # ortools::strong_integers
)
add_test(NAME strong_sparse_vector_test COMMAND strong_sparse_vector_test)

# Add strong_sparse_matrix_test executable
add_executable(strong_sparse_matrix_test strong_sparse_matrix_test.cc)
if(APPLE)
  set_target_properties(strong_sparse_matrix_test PROPERTIES
    INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
elseif(UNIX)
  set_target_properties(strong_sparse_matrix_test PROPERTIES
    INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()
target_link_libraries(strong_sparse_matrix_test PRIVATE
        data_structures::strong_sparse_matrix
        custom_gtest_main
        utils
        GTest::gtest
        GTest::gmock
        ortools::ortools
        # ortools::base
        # ortools::strong_vector
        # ortools::lp_data::base
        # ortools::strong_integers
)
add_test(NAME strong_sparse_matrix_test COMMAND strong_sparse_matrix_test)
