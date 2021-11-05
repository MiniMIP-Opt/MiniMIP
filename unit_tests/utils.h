#include <gtest/gtest.h>

#include "absl/status/status.h"

#define ASSERT_OK(x) ASSERT_EQ((x), absl::OkStatus());
