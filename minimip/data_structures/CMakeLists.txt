# minimip/data_structures/CMakeLists.txt
cmake_minimum_required(VERSION 3.15)

# Define the problem library
add_library(problem problem.cc problem.h)
if(APPLE)
  set_target_properties(problem PROPERTIES
    INSTALL_RPATH "@loader_path")
elseif(UNIX)
  set_target_properties(problem PROPERTIES
    INSTALL_RPATH "$ORIGIN")
endif()
target_link_libraries(problem PUBLIC
  absl::flat_hash_set
  absl::strings
  ortools::ortools
  #ortools::mps_reader
)
add_library(data_structures::problem ALIAS problem)

# Define the cuts_data library
add_library(cuts_data cuts_data.cc cuts_data.h)
if(APPLE)
  set_target_properties(cuts_data PROPERTIES
    INSTALL_RPATH "@loader_path")
elseif(UNIX)
  set_target_properties(cuts_data PROPERTIES
    INSTALL_RPATH "$ORIGIN")
endif()
target_link_libraries(cuts_data PUBLIC
  mip_data
  strong_sparse_vector
  lp_interface::lpi
)
add_library(data_structures::cuts_data ALIAS cuts_data)

# Define the mip_data library
add_library(mip_data mip_data.cc mip_data.h)
if(APPLE)
  set_target_properties(mip_data PROPERTIES
    INSTALL_RPATH "@loader_path")
elseif(UNIX)
  set_target_properties(mip_data PROPERTIES
    INSTALL_RPATH "$ORIGIN")
endif()
target_link_libraries(mip_data PUBLIC
  problem
  strong_sparse_matrix
)
add_library(data_structures::mip_data ALIAS mip_data)

# Define the mip_types library
add_library(mip_types mip_types.cc mip_types.h)
if(APPLE)
  set_target_properties(mip_types PROPERTIES
    INSTALL_RPATH "@loader_path")
elseif(UNIX)
  set_target_properties(mip_types PROPERTIES
    INSTALL_RPATH "$ORIGIN")
endif()
target_link_libraries(mip_types PUBLIC
  ortools::ortools
  #ortools::strong_vector
  #ortools::strong_integers
)
add_library(data_structures::mip_types ALIAS mip_types)

# Define the mip_tree library
add_library(mip_tree mip_tree.cc mip_tree.h)
if(APPLE)
  set_target_properties(mip_tree PROPERTIES
    INSTALL_RPATH "@loader_path")
elseif(UNIX)
  set_target_properties(mip_tree PROPERTIES
    INSTALL_RPATH "$ORIGIN")
endif()
target_link_libraries(mip_tree PUBLIC
  mip_types
  lp_interface::lpi
  ortools::ortools
  #ortools::base
  #ortools::strong_vector
  #ortools::strong_integers
)
add_library(data_structures::mip_tree ALIAS mip_tree)

# Define the strong_sparse_matrix library
add_library(strong_sparse_matrix strong_sparse_matrix.cc strong_sparse_matrix.h)
target_include_directories(strong_sparse_matrix PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>)
if(APPLE)
  set_target_properties(strong_sparse_matrix PROPERTIES
    INSTALL_RPATH "@loader_path")
elseif(UNIX)
  set_target_properties(strong_sparse_matrix PROPERTIES
    INSTALL_RPATH "$ORIGIN")
endif()
target_link_libraries(strong_sparse_matrix PUBLIC
  strong_sparse_vector
  ortools::ortools
  #ortools::base
  #ortools::strong_vector
  #ortools::lp_data_base
  #ortools::strong_integers
)
add_library(data_structures::strong_sparse_matrix ALIAS strong_sparse_matrix)

# Define the strong_sparse_vector library (assuming it's header-only)
add_library(strong_sparse_vector INTERFACE)
target_include_directories(strong_sparse_vector INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>)
if(APPLE)
  set_target_properties(strong_sparse_vector PROPERTIES
    INSTALL_RPATH "@loader_path")
elseif(UNIX)
  set_target_properties(strong_sparse_vector PROPERTIES
    INSTALL_RPATH "$ORIGIN")
endif()
target_link_libraries(strong_sparse_vector INTERFACE
  mip_types
  ortools::ortools
  #ortools::base
  #ortools::strong_vector
  #ortools::strong_integers
)
add_library(data_structures::strong_sparse_vector ALIAS strong_sparse_vector)
