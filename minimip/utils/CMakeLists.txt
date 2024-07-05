# minimip/utils/CMakeLists.txt

# Source file library targets
add_library(reader reader.cc)
target_include_directories(reader PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(reader PUBLIC
  data_structures::problem
  absl::status
  absl::statusor
  parameters_cc_proto
  ortools::ortools
)