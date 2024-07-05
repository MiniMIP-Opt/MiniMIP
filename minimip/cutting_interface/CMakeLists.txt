# minimip/cutting_interface/CMakeLists.txt

# Header-only library targets
add_library(cuts_generator INTERFACE cuts_generator.h)
target_include_directories(cuts_generator INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(cuts_generator INTERFACE
  absl::status
  absl::statusor
  parameters_cc_proto
  data_structures::cuts_data)

add_library(generator_factory INTERFACE generator_factory.h)
target_include_directories(generator_factory INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(generator_factory INTERFACE
  absl::status
  absl::statusor
  parameters_cc_proto
  aggregating_generator)

add_library(cuts_runner INTERFACE cuts_runner.h)
target_include_directories(cuts_runner INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(cuts_runner INTERFACE
  cuts_selector
  cuts_generator
  data_structures::cuts_data
  absl::status
  absl::statusor)
add_library(cutting_interface::cuts_runner ALIAS cuts_runner)

add_library(runner_factory INTERFACE runner_factory.h)
target_include_directories(runner_factory INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(runner_factory INTERFACE
  default_runner
  parameters_cc_proto
  absl::status
  absl::statusor)
add_library(cutting_interface::runner_factory ALIAS runner_factory)

add_library(cuts_selector INTERFACE cuts_selector.h)
target_include_directories(cuts_selector INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(cuts_selector INTERFACE
  data_structures::cuts_data
  absl::status
  absl::statusor
  ortools::ortools
  #ortools::status_macros
)

add_library(selector_factory INTERFACE selector_factory.h)
target_include_directories(selector_factory INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(selector_factory INTERFACE
  hybrid_selector
  parameters_cc_proto
  absl::status
  absl::statusor
)

# Source file library targets
add_library(aggregating_generator aggregating_generator.cc)
target_include_directories(aggregating_generator PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
if(APPLE)
  set_target_properties(aggregating_generator PROPERTIES
    INSTALL_RPATH "@loader_path")
elseif(UNIX)
  set_target_properties(aggregating_generator PROPERTIES
    INSTALL_RPATH "$ORIGIN")
endif()
target_link_libraries(aggregating_generator PUBLIC
  solver_context_interface
  cuts_generator
  lp_interface::lpi
  ortools::ortools
  #ortools::status_macros
)
add_library(cutting_interface::aggregating_generator ALIAS aggregating_generator)


add_library(default_runner default_runner.cc)
target_include_directories(default_runner PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
if(APPLE)
  set_target_properties(default_runner PROPERTIES
    INSTALL_RPATH "@loader_path")
elseif(UNIX)
  set_target_properties(default_runner PROPERTIES
    INSTALL_RPATH "$ORIGIN")
endif()
target_link_libraries(default_runner PUBLIC
  selector_factory
  generator_factory
  cuts_runner
  absl::status
  absl::statusor)
add_library(cutting_interface::default_runner ALIAS default_runner)


add_library(hybrid_selector hybrid_selector.cc)
target_include_directories(hybrid_selector PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
if(APPLE)
  set_target_properties(hybrid_selector PROPERTIES
    INSTALL_RPATH "@loader_path")
elseif(UNIX)
  set_target_properties(hybrid_selector PROPERTIES
    INSTALL_RPATH "$ORIGIN")
endif()
target_link_libraries(hybrid_selector PUBLIC
  solver_context_interface
  cuts_selector
  lp_interface::lpi
  ortools::ortools
  #ortools::status_macros
)
add_library(cutting_interface::hybrid_selector ALIAS hybrid_selector)
