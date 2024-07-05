# minimip/lp_interface/CMakeLists.txt

add_library(lpi lpi.cc lpi.h)
target_include_directories(lpi PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
if(APPLE)
  set_target_properties(lpi PROPERTIES
    INSTALL_RPATH "@loader_path")
elseif(UNIX)
  set_target_properties(lpi PROPERTIES
    INSTALL_RPATH "$ORIGIN")
endif()
target_link_libraries(lpi PUBLIC
  parameters_cc_proto
  data_structures::mip_data
  data_structures::mip_types
  data_structures::strong_sparse_vector
  absl::status
  absl::statusor
  ortools::ortools
  #ortools::base
  #ortools::strong_vector
  #ortools::strong_integers
)
add_library(lp_interface::lpi ALIAS lpi)

add_library(lpi_glop lpi_glop.cc lpi_glop.h)
target_include_directories(lpi_glop PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
if(APPLE)
  set_target_properties(lpi_glop PROPERTIES
    INSTALL_RPATH "@loader_path")
elseif(UNIX)
  set_target_properties(lpi_glop PROPERTIES
    INSTALL_RPATH "$ORIGIN")
endif()
target_link_libraries(lpi_glop PUBLIC
  lpi
  parameters_cc_proto
  ortools::ortools
  #ortools::base
  #ortools::glop_lp_solver
  #ortools::glop_revised_simplex
  #ortools::lp_data_lp_data_utils
  #ortools::lp_data_lp_print_utils
  #ortools::lp_data_proto_utils
  #ortools::lp_data_scattered_vector
  #ortools::util_file_util
  #ortools::util_stats
  #ortools::util_time_limit
)
add_library(lp_interface::lpi_glop ALIAS lpi_glop)

add_library(lpi_soplex lpi_soplex.cc lpi_soplex.h)
target_include_directories(lpi_soplex PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
if(APPLE)
  set_target_properties(lpi_soplex PROPERTIES
    INSTALL_RPATH "@loader_path")
elseif(UNIX)
  set_target_properties(lpi_soplex PROPERTIES
    INSTALL_RPATH "$ORIGIN")
endif()
target_link_libraries(lpi_soplex PUBLIC
  lpi
  libsoplex-pic  # Link against the soplex target defined via FetchContent
  parameters_cc_proto
  ortools::ortools
  #ortools::util_file_util
)
add_library(lp_interface::lpi_soplex ALIAS lpi_soplex)

add_library(lpi_factory lpi_factory.cc lpi_factory.h)
target_include_directories(lpi_factory PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
if(APPLE)
  set_target_properties(lpi_factory PROPERTIES
    INSTALL_RPATH "@loader_path")
elseif(UNIX)
  set_target_properties(lpi_factory PROPERTIES
    INSTALL_RPATH "$ORIGIN")
endif()
target_link_libraries(lpi_factory PUBLIC
  lpi
  lpi_glop
  lpi_soplex
  parameters_cc_proto
  ortools::ortools
  #ortools::base
)
add_library(lp_interface::lpi_factory ALIAS lpi_factory)
