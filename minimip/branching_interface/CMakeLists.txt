# minimip/branching_interface/CMakeLists.txt

# Header-only library targets
add_library(branching_interface INTERFACE branching_interface.h)
target_include_directories(branching_interface INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(branching_interface INTERFACE
  data_structures::mip_types
  parameters_cc_proto
  absl::status
  absl::statusor
)
add_library(branching_interface::branching_interface ALIAS branching_interface)

add_library(branching_factory INTERFACE branching_factory.h)
target_include_directories(branching_factory INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(branching_factory INTERFACE
  random_branching
  maxfractional_branching
  parameters_cc_proto
  absl::status
  absl::statusor
)
add_library(branching_interface::branching_factory ALIAS branching_factory)

# Source file library targets
add_library(maxfractional_branching maxfractional_branching.cc)
target_include_directories(maxfractional_branching PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
if(APPLE)
  set_target_properties(maxfractional_branching PROPERTIES
    INSTALL_RPATH "@loader_path")
elseif(UNIX)
  set_target_properties(maxfractional_branching PROPERTIES
    INSTALL_RPATH "$ORIGIN")
endif()
target_link_libraries(maxfractional_branching PUBLIC
  solver_context_interface
  parameters_cc_proto
  absl::status
  absl::statusor
)

add_library(random_branching random_branching.cc)
target_include_directories(random_branching PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
if(APPLE)
  set_target_properties(random_branching PROPERTIES
    INSTALL_RPATH "@loader_path")
elseif(UNIX)
  set_target_properties(random_branching PROPERTIES
    INSTALL_RPATH "$ORIGIN")
endif()
target_link_libraries(random_branching PUBLIC
  branching_interface
  solver_context_interface
  parameters_cc_proto
  absl::status
  absl::statusor
)
