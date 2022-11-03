// Copyright 2022 the MiniMIP Project
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "absl/status/status.h"
#include "ortools/base/status_macros.h"
#include "src/data_structures/strong_sparse_vector.h"

#ifndef UNIT_TESTS_UTILS_H_
#define UNIT_TESTS_UTILS_H_

#define ASSERT_OK(x) ASSERT_EQ((x), absl::OkStatus());
#define EXPECT_OK(x) EXPECT_EQ((x), absl::OkStatus());

// Small variation on ortools' `ASSIGN_OR_RETURN` macro
#define ASSERT_OK_AND_ASSIGN(lhs, rexpr)    \
  STATUS_MACROS_IMPL_ASSERT_OK_AND_ASSIGN_( \
      STATUS_MACROS_IMPL_CONCAT_(_status_or_value, __COUNTER__), lhs, rexpr);
#define STATUS_MACROS_IMPL_ASSERT_OK_AND_ASSIGN_(statusor, lhs, rexpr) \
  auto (statusor) = (rexpr);                                             \
  ASSERT_OK((statusor).status());                                        \
  STATUS_MACROS_IMPL_UNPARENTHESIS(lhs) = std::move(statusor).value()

namespace minimip {

SparseRow CreateSparseRow(const std::vector<std::pair<int, double>>& v) {
  SparseRow res;
  for (auto [index, value] : v) res.AddEntry(ColIndex(index), value);
  res.CleanUpIfNeeded();
  return res;
}

SparseCol CreateSparseCol(const std::vector<std::pair<int, double>>& v) {
  SparseCol res;
  for (auto [index, value] : v) res.AddEntry(RowIndex(index), value);
  res.CleanUpIfNeeded();
  return res;
}

// Simple helper function to allow writing e.g.
// EXPECT_THAT(row.entries(), EntriesAre<ColIndex>({{1, 2.0}, {2, 7.5}}));
// without having to repeatedly specify the index type.
template <typename IndexType>
auto EntriesAre(const std::vector<std::pair<int, double>>& val) {
  std::vector<SparseEntry<IndexType>> vec;
  for (auto [index, value] : val) vec.push_back({IndexType(index), value});
  return testing::ElementsAreArray(vec);
}

}  // namespace minimip

#endif  // UNIT_TESTS_UTILS_H_
