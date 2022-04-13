#include "src/lp_interface/strong_sparse_matrix.h"

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
  sparse.Resize(dense[RowIndex(0)].size(), dense.size());
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
              ElementsAre(RowEntry(0, 1.0), RowEntry(2, 3.0)));
  EXPECT_THAT(sparse.row(RowIndex(1)).entries(), ElementsAre(RowEntry(1, 2.0)));
  EXPECT_THAT(sparse.col(ColIndex(0)).entries(), ElementsAre(ColEntry(0, 1.0)));
  EXPECT_THAT(sparse.col(ColIndex(1)).entries(), ElementsAre(ColEntry(1, 2.0)));
  EXPECT_THAT(sparse.col(ColIndex(2)).entries(), ElementsAre(ColEntry(0, 3.0)));
  for (RowIndex row(0); row < sparse.num_rows(); ++row) {
    for (ColIndex col(0); col < sparse.num_cols(); ++col) {
      EXPECT_EQ(sparse.GetCoefficient(col, row), dense[row][col]);
    }
  }
}

// TODO(lpawel): Implement more tests after initial review.

}  // namespace
}  // namespace minimip
