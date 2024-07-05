# unit_tests/utils/CMakeLists.txt

# Add reader_test executable
add_executable(reader_test reader_test.cc)
target_link_libraries(reader_test PRIVATE
        gtest_main
        utils
        utils::reader
        absl::status
        GTest::gtest
)

# Add test data to the executable's runtime environment
add_custom_command(TARGET reader_test POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_SOURCE_DIR}/test_data/50v-10.mps
        $<TARGET_FILE_DIR:reader_test>/test_data/50v-10.mps)

add_test(NAME reader_test COMMAND reader_test)
