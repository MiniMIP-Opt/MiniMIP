#include <gtest/gtest.h>

#include "absl/status/status.h"

#ifndef UNIT_TESTS_UTILS_H_
#define UNIT_TESTS_UTILS_H_

#define ASSERT_OK(x) ASSERT_EQ((x), absl::OkStatus());

#endif  // UNIT_TESTS_UTILS_H_
