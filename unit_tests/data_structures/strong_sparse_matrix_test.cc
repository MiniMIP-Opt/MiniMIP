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

#include "src/data_structures/strong_sparse_matrix.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "ortools/base/commandlineflags.h"
#include "ortools/base/logging.h"

namespace minimip {
namespace {

using ::testing::ElementsAre;
using DenseRowBasedMatrix =
    absl::StrongVector<RowIndex, absl::StrongVector<ColIndex, double>>;

StrongSparseMatrix DenseToSparseMatrix(const DenseRowBasedMatrix& dense) {
  StrongSparseMatrix sparse;
  sparse.Resize(ColIndex(dense[RowIndex(0)].size()), RowIndex(dense.size()));
  for (RowIndex row(0); row < sparse.num_rows(); ++row) {
    DCHECK_EQ(ColIndex(dense[row].size()), sparse.num_cols());
    for (ColIndex col(0); col < sparse.num_cols(); ++col) {
      if (dense[row][col] != 0.0) {
        sparse.SetCoefficient(col, row, dense[row][col]);
      }
    }
  }
  QCHECK(sparse.AllRowsAreClean());
  QCHECK(sparse.AllColsAreClean());
  return sparse;
}

TEST(StrongSparseMatrix, InitializeAndGetValues) {
  const DenseRowBasedMatrix dense = {{1.0, 0.0, 3.0}, {0.0, 2.0, 0.0}};
  const StrongSparseMatrix sparse = DenseToSparseMatrix(dense);

  // CopyDenseIntoSparseMatrix(dense, &sparse);
  EXPECT_EQ(sparse.num_rows(), RowIndex(2));
  EXPECT_EQ(sparse.num_cols(), ColIndex(3));
  EXPECT_THAT(sparse.row(RowIndex(0)).entries(),
              ElementsAre(RowEntry(ColIndex(0), 1.0), RowEntry(ColIndex(2), 3.0)));
  EXPECT_THAT(sparse.row(RowIndex(1)).entries(), ElementsAre(RowEntry(ColIndex(1), 2.0)));
  EXPECT_THAT(sparse.col(ColIndex(0)).entries(), ElementsAre(ColEntry(RowIndex(0), 1.0)));
  EXPECT_THAT(sparse.col(ColIndex(1)).entries(), ElementsAre(ColEntry(RowIndex(1), 2.0)));
  EXPECT_THAT(sparse.col(ColIndex(2)).entries(), ElementsAre(ColEntry(RowIndex(0), 3.0)));
  for (RowIndex row(0); row < sparse.num_rows(); ++row) {
    for (ColIndex col(0); col < sparse.num_cols(); ++col) {
      EXPECT_EQ(sparse.GetCoefficient(col, row), dense[row][col]);
    }
  }
}

TEST(StrongSparseMatrix, PopulatesColsFromRows) {
  StrongSparseMatrix sparse(ColIndex(3), RowIndex(2));
  SparseRow row1({RowEntry(ColIndex(1), 2.0)});
  sparse.PopulateRow(RowIndex(0), std::move(row1));
  SparseRow row2({RowEntry(ColIndex(0), 1.0), RowEntry(ColIndex(2), 3.0)});
  sparse.PopulateRow(RowIndex(1), std::move(row2));
  EXPECT_THAT(sparse.col(ColIndex(0)).entries(), ElementsAre(ColEntry(RowIndex(1), 1.0)));
  EXPECT_THAT(sparse.col(ColIndex(1)).entries(), ElementsAre(ColEntry(RowIndex(0), 2.0)));
  EXPECT_THAT(sparse.col(ColIndex(2)).entries(), ElementsAre(ColEntry(RowIndex(1), 3.0)));
}

TEST(StrongSparseMatrix, PopulatesRowsFromCols) {
  StrongSparseMatrix sparse(ColIndex(2), RowIndex(3));
  SparseCol col1({ColEntry(RowIndex(1), 2.0)});
  sparse.PopulateCol(ColIndex(0), std::move(col1));
  SparseCol col2({ColEntry(RowIndex(0), 1.0), ColEntry(RowIndex(2), 3.0)});
  sparse.PopulateCol(ColIndex(1), std::move(col2));
  EXPECT_THAT(sparse.row(RowIndex(0)).entries(), ElementsAre(RowEntry(ColIndex(1), 1.0)));
  EXPECT_THAT(sparse.row(RowIndex(1)).entries(), ElementsAre(RowEntry(ColIndex(0), 2.0)));
  EXPECT_THAT(sparse.row(RowIndex(2)).entries(), ElementsAre(RowEntry(ColIndex(1), 3.0)));
}

// TODO(lpawel): Implement more tests after initial review.

}  // namespace
}  // namespace minimip
