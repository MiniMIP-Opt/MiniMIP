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

#include "src/data_structures/strong_sparse_vector.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "ortools/base/commandlineflags.h"
#include "ortools/base/logging.h"
#include "unit_tests/utils.h"

namespace minimip {
namespace {

using ::testing::ElementsAre;
using ::testing::ElementsAreArray;
using ::testing::IsEmpty;

// We only test SparseRow, because SparseCol is just another templated version.

TEST(StrongSparseRow, InitializeAndGetValue) {
  SparseRow row = CreateSparseRow({{4, 40.0}, {3, 0.0}, {1, 10.0}});
  EXPECT_FALSE(row.MayNeedCleaning());
  EXPECT_TRUE(row.IsClean());
  EXPECT_THAT(row.entries(), EntriesAre<ColIndex>({{1, 10.0}, {4, 40.0}}));
  EXPECT_EQ(row.value(ColIndex(1)), 10.0);
  EXPECT_EQ(row.value(ColIndex(2)), 0.0);
  EXPECT_EQ(row.value(ColIndex(3)), 0.0);
  EXPECT_EQ(row.value(ColIndex(4)), 40.0);
  EXPECT_EQ(row.value(ColIndex(5)), 0.0);
  EXPECT_EQ(row[ColIndex(1)], 10.0);
  EXPECT_EQ(row[ColIndex(2)], 0.0);
  EXPECT_EQ(row[ColIndex(3)], 0.0);
  EXPECT_EQ(row[ColIndex(4)], 40.0);
  EXPECT_EQ(row[ColIndex(5)], 0.0);
  EXPECT_THAT(row.indices(), ElementsAre(ColIndex(1), ColIndex(4)));
  EXPECT_THAT(row.values(), ElementsAre(10.0, 40.0));
}

TEST(StrongSparseRow, MakeCopy) {
  SparseRow x = CreateSparseRow({{1, 10.0}});
  SparseRow y = x;
  // We clear x to make sure y is independent.
  x.Clear();
  EXPECT_TRUE(x.entries().empty());
  EXPECT_FALSE(y.MayNeedCleaning());
  EXPECT_THAT(y.entries(), EntriesAre<ColIndex>({{1, 10.0}}));
}

TEST(StrongSparseRow, MakeReference) {
  SparseRow x = CreateSparseRow({{1, 10.0}});
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
  EXPECT_THAT(row.entries(),
              EntriesAre<ColIndex>({{1, 1.0}, {2, 2.0}, {3, 3.0}}));
  row.AddEntry(ColIndex(3), 30.0);
  row.AddEntry(ColIndex(2), 0.0);
  row.AddEntry(ColIndex(1), 10.0);
  row.CleanUpIfNeeded();
  EXPECT_THAT(row.entries(), EntriesAre<ColIndex>({{1, 10.0}, {3, 30.0}}));
}

TEST(StrongSparseRow, IsNotCleanWhenSortedButContainsDuplicatesIndex) {
  SparseRow row = CreateSparseRow({{1, 1.0}});
  row.AddEntry(ColIndex(1), 2.0);
  EXPECT_TRUE(row.MayNeedCleaning());
}

TEST(StrongSparseRow, IsNotCleanWhenAddsNotSortedNonZero) {
  SparseRow row = CreateSparseRow({{1, 1.0}});
  row.AddEntry(ColIndex(1), 2.0);
  EXPECT_TRUE(row.MayNeedCleaning());
}

TEST(StrongSparseRow, IsNotCleanWhenAddsNotSortedZero) {
  SparseRow row = CreateSparseRow({{1, 1.0}});
  row.AddEntry(ColIndex(1), 0.0);
  EXPECT_TRUE(row.MayNeedCleaning());
}

TEST(StrongSparseRow, IsCleanWhenAddsSortedZero) {
  SparseRow row = CreateSparseRow({{1, 1.0}});
  row.AddEntry(ColIndex(2), 0.0);
  EXPECT_FALSE(row.MayNeedCleaning());
}

TEST(StrongSparseRow, IsCleanWhenAddsSortedNonZero) {
  SparseRow row = CreateSparseRow({{1, 1.0}});
  row.AddEntry(ColIndex(2), 1.0);
  EXPECT_FALSE(row.MayNeedCleaning());
}

TEST(StrongSparseRow, IsCleanWhenAddsZeroToEmpty) {
  SparseRow row;
  row.AddEntry(ColIndex(1), 0.0);
  EXPECT_FALSE(row.MayNeedCleaning());
  EXPECT_THAT(row.entries(), IsEmpty());
}

TEST(StrongSparseRow, MutateValuesWithCleanUpNeeded) {
  SparseRow row = CreateSparseRow({{1, 10.0}, {2, 20.0}});
  for (auto& e : row.mutable_entries()) {
    e.value /= 10.0;
  }
  EXPECT_TRUE(row.MayNeedCleaning());
  row.CleanUpIfNeeded();
  EXPECT_THAT(row.entries(), EntriesAre<ColIndex>({{1, 1.0}, {2, 2.0}}));
}

TEST(StrongSparseRow, MutateIndicesWithCleanUpNeeded) {
  SparseRow row = CreateSparseRow({{1, 10.0}, {2, 20.0}});
  for (auto& e : row.mutable_entries()) {
    e.index = ColIndex(10) - e.index;
  }
  EXPECT_TRUE(row.MayNeedCleaning());
  row.CleanUpIfNeeded();
  EXPECT_THAT(row.entries(), EntriesAre<ColIndex>({{8, 20.0}, {9, 10.0}}));
}

TEST(StrongSparseRowDeathTest, DiesWhenForgotToCleanUpMutatedEntries1) {
  // The test only dies in debug mode.
  if (!DEBUG_MODE) return;
  ASSERT_DEATH(
      {
        SparseRow row = CreateSparseRow({{1, 10.0}, {2, 20.0}});
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
        SparseRow row = CreateSparseRow({{1, 10.0}, {2, 20.0}});
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
        SparseRow row = CreateSparseRow({{1, 10.0}, {2, 20.0}});
        // Accessing mutable entries marks the vector for cleaning.
        ASSERT_EQ(row.mutable_entries().size(), 2);
        EXPECT_EQ(row.value(ColIndex(8)), 20.0);
      },
      "The vector is not clean!");
  ASSERT_DEATH(
      {
        SparseRow row = CreateSparseRow({{1, 10.0}, {2, 20.0}});
        // Accessing mutable entries marks the vector for cleaning.
        ASSERT_EQ(row.mutable_entries().size(), 2);
        EXPECT_EQ(row[ColIndex(8)], 20.0);
      },
      "The vector is not clean!");
}

TEST(StrongSparseRow, Negation) {
  const SparseRow row = CreateSparseRow({{2, 3.0}, {3, -2.0}});
  EXPECT_THAT((-row).entries(), EntriesAre<ColIndex>({{2, -3.0}, {3, 2.0}}));
}

TEST(StrongSparseRow, ScalarMultiplication) {
  const SparseRow row = CreateSparseRow({{2, 3.0}, {3, -2.0}});
  EXPECT_THAT((row * 1.5).entries(),
              EntriesAre<ColIndex>({{2, 4.5}, {3, -3.0}}));
  EXPECT_THAT((1.5 * row).entries(),
              EntriesAre<ColIndex>({{2, 4.5}, {3, -3.0}}));
  EXPECT_THAT(((-1.5) * row).entries(),
              EntriesAre<ColIndex>({{2, -4.5}, {3, 3.0}}));
  EXPECT_THAT((row * (-1.5)).entries(),
              EntriesAre<ColIndex>({{2, -4.5}, {3, 3.0}}));
  EXPECT_THAT((row * 0).entries(), IsEmpty());
  EXPECT_THAT((0 * row).entries(), IsEmpty());

  {
    SparseRow modifiable_row = row;
    EXPECT_THAT((modifiable_row *= 1.5).entries(),
                EntriesAre<ColIndex>({{2, 4.5}, {3, -3.0}}));
    EXPECT_THAT(modifiable_row.entries(),
                EntriesAre<ColIndex>({{2, 4.5}, {3, -3.0}}));
  }

  {
    SparseRow modifiable_row = row;
    EXPECT_THAT((modifiable_row *= -1.5).entries(),
                EntriesAre<ColIndex>({{2, -4.5}, {3, 3.0}}));
    EXPECT_THAT(modifiable_row.entries(),
                EntriesAre<ColIndex>({{2, -4.5}, {3, 3.0}}));
  }

  {
    SparseRow modifiable_row = row;
    EXPECT_THAT((modifiable_row *= 0.0).entries(), IsEmpty());
    EXPECT_THAT(modifiable_row.entries(), IsEmpty());
  }
}

TEST(StrongSparseRow, ScalarDivision) {
  const SparseRow row = CreateSparseRow({{2, 6.25}, {3, -5.0}});
  EXPECT_THAT((row / 2.5).entries(),
              EntriesAre<ColIndex>({{2, 2.5}, {3, -2.0}}));
  EXPECT_THAT((row / (-2.5)).entries(),
              EntriesAre<ColIndex>({{2, -2.5}, {3, 2.0}}));

  {
    SparseRow modifiable_row = row;
    EXPECT_THAT((modifiable_row /= 2.5).entries(),
                EntriesAre<ColIndex>({{2, 2.5}, {3, -2.0}}));
    EXPECT_THAT(modifiable_row.entries(),
                EntriesAre<ColIndex>({{2, 2.5}, {3, -2.0}}));
  }

  {
    SparseRow modifiable_row = row;
    EXPECT_THAT((modifiable_row /= (-2.5)).entries(),
                EntriesAre<ColIndex>({{2, -2.5}, {3, 2.0}}));
    EXPECT_THAT(modifiable_row.entries(),
                EntriesAre<ColIndex>({{2, -2.5}, {3, 2.0}}));
  }
}

TEST(StrongSparseRow, VectorAddition) {
  const SparseRow row1 = CreateSparseRow({{2, 6.25}, {3, -5.0}});
  const SparseRow row2 = CreateSparseRow({{1, 2.3}, {3, 1.0}, {5, 1.23}});
  EXPECT_THAT(
      (row1 + row2).entries(),
      EntriesAre<ColIndex>({{1, 2.3}, {2, 6.25}, {3, -4.0}, {5, 1.23}}));
  EXPECT_THAT(
      (row2 + row1).entries(),
      EntriesAre<ColIndex>({{1, 2.3}, {2, 6.25}, {3, -4.0}, {5, 1.23}}));

  {
    SparseRow modifiable_row = row1;
    EXPECT_THAT(
        (modifiable_row += row2).entries(),
        EntriesAre<ColIndex>({{1, 2.3}, {2, 6.25}, {3, -4.0}, {5, 1.23}}));
    EXPECT_THAT(
        modifiable_row.entries(),
        EntriesAre<ColIndex>({{1, 2.3}, {2, 6.25}, {3, -4.0}, {5, 1.23}}));
  }

  {
    SparseRow modifiable_row = row2;
    EXPECT_THAT(
        (modifiable_row += row1).entries(),
        EntriesAre<ColIndex>({{1, 2.3}, {2, 6.25}, {3, -4.0}, {5, 1.23}}));
    EXPECT_THAT(
        modifiable_row.entries(),
        EntriesAre<ColIndex>({{1, 2.3}, {2, 6.25}, {3, -4.0}, {5, 1.23}}));
  }
}

TEST(StrongSparseRow, VectorSubtraction) {
  const SparseRow row1 = CreateSparseRow({{2, 6.25}, {3, -5.0}});
  const SparseRow row2 = CreateSparseRow({{1, 2.3}, {3, 1.0}, {5, 1.23}});
  EXPECT_THAT(
      (row1 - row2).entries(),
      EntriesAre<ColIndex>({{1, -2.3}, {2, 6.25}, {3, -6.0}, {5, -1.23}}));
  EXPECT_THAT(
      (row2 - row1).entries(),
      EntriesAre<ColIndex>({{1, 2.3}, {2, -6.25}, {3, 6.0}, {5, 1.23}}));

  {
    SparseRow modifiable_row = row1;
    EXPECT_THAT(
        (modifiable_row -= row2).entries(),
        EntriesAre<ColIndex>({{1, -2.3}, {2, 6.25}, {3, -6.0}, {5, -1.23}}));
    EXPECT_THAT(
        modifiable_row.entries(),
        EntriesAre<ColIndex>({{1, -2.3}, {2, 6.25}, {3, -6.0}, {5, -1.23}}));
  }

  {
    SparseRow modifiable_row = row2;
    EXPECT_THAT(
        (modifiable_row -= row1).entries(),
        EntriesAre<ColIndex>({{1, 2.3}, {2, -6.25}, {3, 6.0}, {5, 1.23}}));
    EXPECT_THAT(
        modifiable_row.entries(),
        EntriesAre<ColIndex>({{1, 2.3}, {2, -6.25}, {3, 6.0}, {5, 1.23}}));
  }
}

TEST(StrongSparseRow, AddMultipleOfVector) {
  const SparseRow row1 = CreateSparseRow({{1, 2.4}, {4, -3.2}});
  const SparseRow row2 = CreateSparseRow({{0, -1.2}, {4, 3.2}, {5, 1.0}});

  {
    // Note: Element 4 cancels and should be removed in the result vector.
    SparseRow modifiable_row = row1;
    EXPECT_THAT(modifiable_row.AddMultipleOfVector(1.0, row2).entries(),
                EntriesAre<ColIndex>({{0, -1.2}, {1, 2.4}, {5, 1.0}}));
    EXPECT_THAT(modifiable_row.entries(),
                EntriesAre<ColIndex>({{0, -1.2}, {1, 2.4}, {5, 1.0}}));
  }

  {
    SparseRow modifiable_row = row1;
    EXPECT_THAT(modifiable_row.AddMultipleOfVector(0.0, row2).entries(),
                ElementsAreArray(row1.entries()));
    EXPECT_THAT(modifiable_row.entries(), ElementsAreArray(row1.entries()));
  }

  {
    SparseRow modifiable_row = row1;
    EXPECT_THAT(
        modifiable_row.AddMultipleOfVector(-100.0, row2).entries(),
        EntriesAre<ColIndex>({{0, 120}, {1, 2.4}, {4, -323.2}, {5, -100.0}}));
    EXPECT_THAT(
        modifiable_row.entries(),
        EntriesAre<ColIndex>({{0, 120}, {1, 2.4}, {4, -323.2}, {5, -100.0}}));
  }
}

TEST(StrongSparseRow, DotProduct) {
  // DotProduct is templated for different index types, so we simply test that
  // different configurations compile and then run all tests on a single
  // configuration.

  EXPECT_EQ(CreateSparseRow({{0, 1.0}}).DotProduct(CreateSparseRow({{0, 1.0}})),
            1.0);
  EXPECT_EQ(CreateSparseRow({{0, 1.0}}).DotProduct(CreateSparseCol({{0, 1.0}})),
            1.0);
  EXPECT_EQ(CreateSparseCol({{0, 1.0}}).DotProduct(CreateSparseRow({{0, 1.0}})),
            1.0);
  EXPECT_EQ(CreateSparseCol({{0, 1.0}}).DotProduct(CreateSparseCol({{0, 1.0}})),
            1.0);

  EXPECT_EQ(CreateSparseRow({{0, 1.0}, {2, 1.0}})
                .DotProduct(CreateSparseRow({{1, 1.0}, {3, 1.0}})),
            0.0);
  EXPECT_EQ(CreateSparseRow({{0, 1.0}, {2, 1.0}})
                .DotProduct(CreateSparseRow({{0, 1.0}, {2, 1.0}})),
            2.0);
  EXPECT_EQ(CreateSparseRow({{0, 1.0}, {1, -2.0}, {2, 1.0}})
                .DotProduct(CreateSparseRow({{0, 1.0}, {1, 1.0}, {2, 1.0}})),
            0.0);
  EXPECT_EQ(CreateSparseRow({{0, 1.0}, {1, -2.0}, {2, 1.0}})
                .DotProduct(CreateSparseRow({{0, 1.0}, {1, 2.0}, {2, 1.0}})),
            -2.0);
  EXPECT_EQ(CreateSparseRow({{0, 1.0}, {2, 2.5}})
                .DotProduct(CreateSparseRow({{2, 2.5}})),
            6.25);
}

}  // namespace
}  // namespace minimip
