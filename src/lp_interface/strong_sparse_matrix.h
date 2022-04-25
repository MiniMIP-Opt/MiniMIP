// Implementation of `StrongSparseMatrix` that maintains both row-major and
// col-major views, and lazily keeps them in sync recomputing one from the other
// (and back) when needed. Under-the-hood the matrix uses `SparseCol` and
// `SparseRow`, which are strongly typed by `RowIndex` and `ColIndex`,
// respectively.

#ifndef SRC_LP_INTERFACE_STRONG_SPARSE_MATRIX_H_
#define SRC_LP_INTERFACE_STRONG_SPARSE_MATRIX_H_

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include "ortools/base/logging.h"
#include "ortools/base/strong_vector.h"
#include "src/lp_interface/strong_sparse_vector.h"

namespace minimip {

class StrongSparseMatrix {
 public:
  StrongSparseMatrix() {}

  StrongSparseMatrix(ColIndex num_cols, RowIndex num_rows)
      : cols_(num_cols.value()), rows_(num_rows.value()) {}

  // StrongSparseMatrix is not copyable to make sure a copy will not be
  // triggered by accident (copy constructor and assign operator are private).
  // StrongSparseMatrix is (no-throw) moveable.
  StrongSparseMatrix(StrongSparseMatrix&&) noexcept            = default;
  StrongSparseMatrix& operator=(StrongSparseMatrix&&) noexcept = default;

  // Use this to initialize by deep copy from another matrix `m`. Under-the-hood
  // we just use a private copy / move constructor and assignment operator.
  void PopulateFromSparseMatrix(StrongSparseMatrix m);

  // Current size of this matrix. Note, because we're paranoid and want the
  // compiler to catch any bug related to accidently swapping rows and cols
  // in some code, we return appropriate strong types (i.e., `ColIndex` and
  // `RowIndex`) and not `int` or `size_t`.
  ColIndex num_cols() const { return ColIndex(cols_.size()); }
  RowIndex num_rows() const { return RowIndex(rows_.size()); }

  // Getters to rows and cols represented as sparse vectors.
  const SparseRow& row(RowIndex row) const;
  const SparseCol& col(ColIndex col) const;

  // This usually runs in O(log(k)) where k is the number of entries in the
  // underlying container (row or col). We prefer to avoid cleaning up the
  // underlying container, but if that's impossible it will run in O(k log(k)).
  // This never triggers a recomputation of the row / col major view, because
  // one of them is guaranteed to be up to date.
  double GetCoefficient(ColIndex col, RowIndex row) const;

  // Runs in O(1), may cause (lazy) row and/or col clean up on later access.
  void SetCoefficient(ColIndex col, RowIndex row, double value);

  // Overrides the row at `row` with `row_data`.
  void PopulateRow(RowIndex row, SparseRow row_data);

  // Overrides the column at `col` with `col_data`.
  void PopulateCol(ColIndex col, SparseCol col_data);

  // Recomputes `rows_` or `cols_` only if really needed. Prefers to keep
  // `rows_` (over `cols_`) in a consistent state, if possible. We use strong
  // types to prevent bugs related to accidently swapping rows and cols.
  void Resize(ColIndex new_num_cols, RowIndex new_num_rows);

  // Clears both the row and col major view of the matrix.
  void Clear() {
    rows_.clear();
    cols_.clear();
    consistency_ = MatrixConsistency::kRowsAndColsAreUpToDate;
  }

  bool rows_need_recomputation() const {
    return consistency_ == MatrixConsistency::kRowsNeedRecomputation;
  }

  bool cols_need_recomputation() const {
    return consistency_ == MatrixConsistency::kColsNeedRecomputation;
  }

  bool AllRowsAreClean() const {
    if (rows_need_recomputation()) return false;
    return std::all_of(rows_.begin(), rows_.end(), [](const SparseRow& row) {
      return !row.MayNeedCleaning();
    });
  }

  bool AllColsAreClean() const {
    if (cols_need_recomputation()) return false;
    return std::all_of(cols_.begin(), cols_.end(), [](const SparseCol& col) {
      return !col.MayNeedCleaning();
    });
  }

 private:
  // This class maintains both row and col major view and keeps them in sync.
  enum class MatrixConsistency {
    // Both `rows_` and `cols_` are up-to-date. The rows and cols stored there
    // could still require cleaning (which will happen lazily on access).
    kRowsAndColsAreUpToDate = 0,

    // Only `cols_` are up-to-date, but not `rows_`. I.e., `cols_` can be
    // mutated, but `rows_` needs to be recomputed before mutating it (which
    // will happen lazily when needed).
    kRowsNeedRecomputation = 1,

    // Only `rows_` are up-to-date, but not `cols_`. I.e., `rows_` can be
    // mutated, but `cols_` needs to be recomputed before mutating it (which
    // will happen lazily when needed).
    kColsNeedRecomputation = 2,
  };

  // We keep the copy constructor and copy assing operator around to use in
  // `PopulateFromSparseMatrix()`.
  StrongSparseMatrix(const StrongSparseMatrix&)            = default;
  StrongSparseMatrix& operator=(const StrongSparseMatrix&) = default;

  // Helper function to compute the transpose in both directions.
  template <typename TargetIndex, typename TransposeIndex>
  void PopulateFromTranspose(
      absl::StrongVector<TargetIndex,
                         StrongSparseVectorOfDoubles<TransposeIndex>>& target,
      absl::StrongVector<TransposeIndex,
                         StrongSparseVectorOfDoubles<TargetIndex>>& transpose)
      const {
    for (auto& target_entry : target) {
      target_entry.Clear();
    }
    for (TransposeIndex i(0); i < transpose.size(); ++i) {
      transpose[i].CleanUpIfNeeded();
      for (const auto& e : transpose[i].entries()) {
        DCHECK_GE(e.index, TargetIndex(0));
        DCHECK_LT(e.index, TargetIndex(target.size()));
        target[e.index].AddEntry(i, e.value);
      }
    }
    consistency_ = MatrixConsistency::kRowsAndColsAreUpToDate;
    DCHECK(AllRowsAreClean());
    DCHECK(AllColsAreClean());
  }

  // Helper functions to recompute row / col view. They are marked const so that
  // getters can be marked const too.
  void RecomputeRowsFromColsIfNeeded() const;
  void RecomputeColsFromRowsIfNeeded() const;

  // Rows and cols are mutable so that we can automatically clean them (if
  // needed) on access.
  mutable absl::StrongVector<ColIndex, SparseCol> cols_;
  mutable absl::StrongVector<RowIndex, SparseRow> rows_;

  // Indicates which of the `cols_` and `rows_` needs recomputation (if any).
  mutable MatrixConsistency consistency_;
};

}  // namespace minimip

#endif  // SRC_LP_INTERFACE_STRONG_SPARSE_MATRIX_H_
