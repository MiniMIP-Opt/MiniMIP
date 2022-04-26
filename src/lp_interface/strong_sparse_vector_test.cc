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

#include "src/lp_interface/strong_sparse_vector.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "ortools/base/commandlineflags.h"
#include "ortools/base/logging.h"

namespace minimip {
namespace {

using ::testing::ElementsAre;

// We only test SparseRow, because SparseCol is just another templated version.

TEST(StrongSparseRow, InitializeAndGetValue) {
  SparseRow row({{4, 40.0}, {3, 0.0}, {1, 10.0}});
  EXPECT_FALSE(row.MayNeedCleaning());
  EXPECT_TRUE(row.IsClean());
  EXPECT_THAT(row.entries(), ElementsAre(RowEntry(1, 10.0), RowEntry(4, 40.0)));
  EXPECT_EQ(row.value(ColIndex(1)), 10.0);
  EXPECT_EQ(row.value(ColIndex(2)), 0.0);
  EXPECT_EQ(row.value(ColIndex(3)), 0.0);
  EXPECT_EQ(row.value(ColIndex(4)), 40.0);
  EXPECT_EQ(row.value(ColIndex(5)), 0.0);
}

TEST(StrongSparseRow, MakeCopy) {
  SparseRow x({{1, 10.0}});
  SparseRow y = x;
  // We clear x to make sure y is independent.
  x.Clear();
  EXPECT_TRUE(x.entries().empty());
  EXPECT_FALSE(y.MayNeedCleaning());
  EXPECT_THAT(y.entries(), ElementsAre(RowEntry(1, 10.0)));
}

TEST(StrongSparseRow, MakeReference) {
  SparseRow x({{1, 10.0}});
  // Note, we don't need a "view" type for sparse vector, just use a reference.
  const SparseRow& y = x;
  x.Clear();
  EXPECT_TRUE(x.entries().empty());
  EXPECT_TRUE(y.entries().empty());
}

TEST(StrongSparseRow, AndEntryAndvalue) {
  SparseRow row;
  row.AddEntry(ColIndex(1), 1.0);
  row.AddEntry(ColIndex(2), 2.0);
  row.AddEntry(ColIndex(3), 3.0);
  EXPECT_THAT(row.entries(), ElementsAre(RowEntry(1, 1.0), RowEntry(2, 2.0),
                                         RowEntry(3, 3.0)));
  row.AddEntry(ColIndex(3), 30.0);
  row.AddEntry(ColIndex(2), 0.0);
  row.AddEntry(ColIndex(1), 10.0);
  row.CleanUpIfNeeded();
  EXPECT_THAT(row.entries(), ElementsAre(RowEntry(1, 10.0), RowEntry(3, 30.0)));
}

TEST(StrongSparseRow, IsNotCleanWhenSortedButContainsDuplicatesIndex) {
  SparseRow row({{1, 1.0}});
  row.AddEntry(ColIndex(1), 2.0);
  EXPECT_TRUE(row.MayNeedCleaning());
}

TEST(StrongSparseRow, IsNotCleanWhenAddsNotSortedNonZero) {
  SparseRow row({{1, 1.0}});
  row.AddEntry(ColIndex(1), 2.0);
  EXPECT_TRUE(row.MayNeedCleaning());
}

TEST(StrongSparseRow, IsNotCleanWhenAddsNotSortedZero) {
  SparseRow row({{1, 1.0}});
  row.AddEntry(ColIndex(1), 0.0);
  EXPECT_TRUE(row.MayNeedCleaning());
}

TEST(StrongSparseRow, IsCleanWhenAddsSortedZero) {
  SparseRow row({{1, 1.0}});
  row.AddEntry(ColIndex(2), 0.0);
  EXPECT_FALSE(row.MayNeedCleaning());
}

TEST(StrongSparseRow, IsCleanWhenAddsSortedNonZero) {
  SparseRow row({{1, 1.0}});
  row.AddEntry(ColIndex(2), 1.0);
  EXPECT_FALSE(row.MayNeedCleaning());
}

TEST(StrongSparseRow, MutateValuesWithCleanUpNeeded) {
  SparseRow row({{1, 10.0}, {2, 20.0}});
  for (auto& e : row.mutable_entries()) {
    e.value /= 10.0;
  }
  EXPECT_TRUE(row.MayNeedCleaning());
  row.CleanUpIfNeeded();
  EXPECT_THAT(row.entries(), ElementsAre(RowEntry(1, 1.0), RowEntry(2, 2.0)));
}

TEST(StrongSparseRow, MutateIndicesWithCleanUpNeeded) {
  SparseRow row({{1, 10.0}, {2, 20.0}});
  for (auto& e : row.mutable_entries()) {
    e.index = ColIndex(10) - e.index;
  }
  EXPECT_TRUE(row.MayNeedCleaning());
  row.CleanUpIfNeeded();
  EXPECT_THAT(row.entries(), ElementsAre(RowEntry(8, 20.0), RowEntry(9, 10.0)));
}

TEST(StrongSparseRowDeathTest, DiesWhenForgotToCleanUpMutatedEntries1) {
  // The test only dies in debug mode.
  if (!DEBUG_MODE) return;
  ASSERT_DEATH(
      {
        SparseRow row({{1, 10.0}, {2, 20.0}});
        // Accessing mutable entries marks the vector for cleaning.
        ASSERT_EQ(row.mutable_entries().size(), 2);
        EXPECT_EQ(row.entries().size(), 2);
      },
      "The vector is not clean!");
}

TEST(StrongSparseRowDeathTest, DiesWhenForgotToCleanUpMutatedEntries2) {
  // The test only dies in debug mode.
  if (!DEBUG_MODE) return;
  ASSERT_DEATH(
      {
        SparseRow row({{1, 10.0}, {2, 20.0}});
        // Accessing mutable entries marks the vector for cleaning.
        ASSERT_EQ(row.mutable_entries().size(), 2);
        EXPECT_EQ(row.mutable_entries().size(), 2);
      },
      "The vector is not clean!");
}

TEST(StrongSparseRowDeathTest, DiesWhenForgotToCleanUpMutatedEntries3) {
  // The test only dies in debug mode.
  if (!DEBUG_MODE) return;
  ASSERT_DEATH(
      {
        SparseRow row({{1, 10.0}, {2, 20.0}});
        // Accessing mutable entries marks the vector for cleaning.
        ASSERT_EQ(row.mutable_entries().size(), 2);
        EXPECT_EQ(row.value(ColIndex(8)), 20.0);
      },
      "The vector is not clean!");
}

}  // namespace
}  // namespace minimip
