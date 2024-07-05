# minimip/CMakeLists.txt

# Include subdirectories first and in order as libraries defined there might be needed
include_directories(${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR})
add_subdirectory(data_structures)
add_subdirectory(lp_interface)
add_subdirectory(cutting_interface)
add_subdirectory(branching_interface)
add_subdirectory(utils)

# Handle the Protobuf file

# When fetching Protobuf the macro protobuf_generate_cpp is not available...

if(NOT CMAKE_CROSSCOMPILING)
  set(PROTOC_PRG protobuf::protoc)
else()
  message(FATAL_ERROR "Not supported here")
endif()
# get_cpp_proto()
# CMake macro to generate Protobuf cpp sources
# Parameters:
#  the proto c++ headers list
#  the proto c++ sources list
# e.g.:
# get_cpp_proto(PROTO_HDRS PROTO_SRCS)
macro(get_cpp_proto PROTO_HDRS PROTO_SRCS)
  file(GLOB_RECURSE PROTO_FILES RELATIVE ${PROJECT_SOURCE_DIR} "*.proto")
  ## Get Protobuf include dir
  get_target_property(protobuf_dirs protobuf::libprotobuf INTERFACE_INCLUDE_DIRECTORIES)
  foreach(dir IN LISTS protobuf_dirs)
    if (NOT "${dir}" MATCHES "INSTALL_INTERFACE|-NOTFOUND")
      message(STATUS "protoc(cc) Adding proto path: ${dir}")
      list(APPEND PROTO_DIRS "--proto_path=${dir}")
    endif()
  endforeach()

  foreach(PROTO_FILE IN LISTS PROTO_FILES)
    message(STATUS "protoc(cc) .proto: ${PROTO_FILE}")
    get_filename_component(PROTO_DIR ${PROTO_FILE} DIRECTORY)
    get_filename_component(PROTO_NAME ${PROTO_FILE} NAME_WE)
    set(PROTO_HDR ${PROJECT_BINARY_DIR}/${PROTO_DIR}/${PROTO_NAME}.pb.h)
    set(PROTO_SRC ${PROJECT_BINARY_DIR}/${PROTO_DIR}/${PROTO_NAME}.pb.cc)
    message(STATUS "protoc(cc) hdr: ${PROTO_HDR}")
    message(STATUS "protoc(cc) src: ${PROTO_SRC}")
    add_custom_command(
      OUTPUT ${PROTO_SRC} ${PROTO_HDR}
      COMMAND ${PROTOC_PRG}
        "--proto_path=${PROJECT_SOURCE_DIR}"
        "--cpp_out=${PROJECT_BINARY_DIR}"
        ${PROTO_DIRS} ${PROTO_FILE}
      DEPENDS ${PROTO_NAME}.proto ${PROTOC_PRG}
      COMMENT "Generate C++ protocol buffer for ${PROTO_FILE}"
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      VERBATIM)
    list(APPEND ${PROTO_HDRS} ${PROTO_HDR})
    list(APPEND ${PROTO_SRCS} ${PROTO_SRC})
  endforeach()
endmacro()
get_cpp_proto(PROTO_HDRS PROTO_SRCS)

#find_package(Protobuf REQUIRED)
#protobuf_generate_cpp(PROTO_SRCS PROTO_HDRS parameters.proto)

add_library(parameters_cc_proto ${PROTO_SRCS} ${PROTO_HDRS})
target_include_directories(parameters_cc_proto PUBLIC
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>)
if(APPLE)
  set_target_properties(parameters_cc_proto PROPERTIES
    INSTALL_RPATH "@loader_path")
elseif(UNIX)
  set_target_properties(parameters_cc_proto PROPERTIES
    INSTALL_RPATH "$ORIGIN")
endif()
target_link_libraries(parameters_cc_proto PUBLIC protobuf::libprotobuf)

# Parameters Factory Header-Only Library
add_library(parameters_factory INTERFACE parameters_factory.h)
target_include_directories(parameters_factory INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(parameters_factory INTERFACE parameters_cc_proto)

# Solver Context Interface Header-Only Library
add_library(solver_context_interface INTERFACE solver_context_interface.h)
target_include_directories(solver_context_interface INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(solver_context_interface INTERFACE
  branching_interface::branching_interface
  cutting_interface::cuts_runner
  data_structures::cuts_data
  data_structures::mip_data
  data_structures::mip_tree
  data_structures::problem
  lp_interface::lpi
  # Ensure these dependencies are defined or findable by CMake
)

# Define the `solver` target
add_library(solver solver.cc solver.h)
target_include_directories(solver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
if(APPLE)
  set_target_properties(solver PROPERTIES
    INSTALL_RPATH "@loader_path")
elseif(UNIX)
  set_target_properties(solver PROPERTIES
    INSTALL_RPATH "$ORIGIN")
endif()
target_link_libraries(solver PUBLIC
  solver_context_interface
  absl::status
  branching_interface::branching_factory
  cutting_interface::runner_factory
  lp_interface::lpi_factory
  parameters_cc_proto
  parameters_factory
  # Ensure these dependencies are defined or findable by CMake
)
