// Copyright 2024 the MiniMIP Project
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

#include "minimip/data_structures/strong_sparse_matrix.h"

#include <algorithm>
#include <utility>
#include <vector>

#include "ortools/base/logging.h"
#include "ortools/base/strong_vector.h"

namespace minimip {

void StrongSparseMatrix::PopulateFromSparseMatrix(StrongSparseMatrix m) {
  VLOG(10) << "calling PopulateFromSparseMatrix().";
  *this = std::move(m);
}

const SparseRow& StrongSparseMatrix::row(RowIndex row) const {
  DCHECK_GE(row, RowIndex(0));
  DCHECK_LT(row, num_rows());
  RecomputeRowsFromColsIfNeeded();
  rows_[row].CleanUpIfNeeded();
  return rows_[row];
}

const SparseCol& StrongSparseMatrix::col(ColIndex col) const {
  DCHECK_GE(col, ColIndex(0));
  DCHECK_LT(col, num_cols());
  RecomputeColsFromRowsIfNeeded();
  cols_[col].CleanUpIfNeeded();
  return cols_[col];
}

double StrongSparseMatrix::GetCoefficient(ColIndex col, RowIndex row) const {
  VLOG(10) << "calling GetCoefficient().";
  DCHECK_GE(col, ColIndex(0));
  DCHECK_LT(col, num_cols());
  DCHECK_GE(row, RowIndex(0));
  DCHECK_LT(row, num_rows());

  if (!rows_need_recomputation() and !rows_[row].MayNeedCleaning()) {
    return rows_[row].value(col);
  }

  if (!cols_need_recomputation() and !cols_[col].MayNeedCleaning()) {
    return cols_[col].value(row);
  }

  // We try reading the data from rows first (if possible) -- we assume rows
  // will be in cleaned state more often (in contrast to an LP solver, a MIP
  // solver usually accesses data row-by-row).
  if (!rows_need_recomputation()) {
    rows_[row].CleanUpIfNeeded();
    return rows_[row].value(col);
  }

  cols_[col].CleanUpIfNeeded();
  return cols_[col].value(row);
}

void StrongSparseMatrix::SetCoefficient(ColIndex col, RowIndex row,
                                        double value) {
  VLOG(10) << "calling SetCoefficient().";
  DCHECK_GE(col, ColIndex(0));
  DCHECK_LT(col, num_cols());
  DCHECK_GE(row, RowIndex(0));
  DCHECK_LT(row, num_rows());
  rows_[row].AddEntry(col, value);
  cols_[col].AddEntry(row, value);
}

void StrongSparseMatrix::PopulateCol(ColIndex col, SparseCol col_data) {
  VLOG(10) << "calling PopulateCol().";
  RecomputeColsFromRowsIfNeeded();
  DCHECK_GE(col, ColIndex(0));
  DCHECK_LT(col, num_cols());
  cols_[col] = std::move(col_data);
  cols_[col].CleanUpIfNeeded();
  DCHECK(cols_[col].entries().empty() ||
         (cols_[col].first_index() >= RowIndex(0) &&
          cols_[col].last_index() < num_rows()));
  consistency_ = MatrixConsistency::kRowsNeedRecomputation;
}

void StrongSparseMatrix::PopulateRow(RowIndex row, SparseRow row_data) {
  VLOG(10) << "calling PopulateRow().";

  RecomputeRowsFromColsIfNeeded();
  DCHECK_GE(row, RowIndex(0));
  DCHECK_LT(row, num_rows());
  rows_[row] = std::move(row_data);
  rows_[row].CleanUpIfNeeded();
  DCHECK(rows_[row].entries().empty() ||
         (rows_[row].first_index() >= ColIndex(0) &&
          rows_[row].last_index() < num_cols()));
  consistency_ = MatrixConsistency::kColsNeedRecomputation;
}

void StrongSparseMatrix::Resize(ColIndex new_num_cols, RowIndex new_num_rows) {
  VLOG(10) << "calling Resize().";
  // Both row and col dimensions drop to 0.
  if (new_num_cols == ColIndex(0) and new_num_rows == RowIndex(0)) {
    Clear();
    return;
  }

  // Only row dimension drops to 0.
  if (new_num_rows == RowIndex(0)) {
    rows_.clear();
    cols_.resize(new_num_cols.value());
    for (auto& col : cols_) {
      col.Clear();
    }
    consistency_ = MatrixConsistency::kRowsAndColsAreUpToDate;
    return;
  }

  // Only col dimension drops to 0.
  if (new_num_cols == ColIndex(0)) {
    cols_.clear();
    rows_.resize(new_num_rows.value());
    for (auto& row : rows_) {
      row.Clear();
    }
    consistency_ = MatrixConsistency::kRowsAndColsAreUpToDate;
    return;
  }

  // We just add empty rows.
  if (new_num_rows >= num_rows()) {
    rows_.resize(new_num_rows.value());
  }

  // We just add empty cols.
  if (new_num_cols >= num_cols()) {
    cols_.resize(new_num_cols.value());
  }

  // Nothing else to do.
  if (new_num_rows == num_rows() and new_num_cols == num_cols()) {
    return;
  }

  // We are left only with deleting rows.
  if (new_num_rows < num_rows() and new_num_cols == num_cols()) {
    RecomputeRowsFromColsIfNeeded();
    rows_.resize(new_num_rows.value());
    consistency_ = MatrixConsistency::kColsNeedRecomputation;
    return;
  }

  // We are left only with deleting cols.
  if (new_num_cols < num_cols() and new_num_rows == num_rows()) {
    RecomputeColsFromRowsIfNeeded();
    cols_.resize(new_num_cols.value());
    consistency_ = MatrixConsistency::kRowsNeedRecomputation;
    return;
  }

  // We are left with deleting both rows and cols.
  if (!cols_need_recomputation()) {
    // Note, this `if` also applies to the case when both rows and cols do not
    // need recomputation. In this case, we prefer to end up with updated rows
    // (and cols marked for recomputation).
    cols_.resize(new_num_cols.value());
    consistency_ = MatrixConsistency::kRowsNeedRecomputation;
    RecomputeRowsFromColsIfNeeded();
    rows_.resize(new_num_rows.value());
    consistency_ = MatrixConsistency::kColsNeedRecomputation;
  } else {
    DCHECK(!rows_need_recomputation());
    rows_.resize(new_num_rows.value());
    consistency_ = MatrixConsistency::kColsNeedRecomputation;
    RecomputeColsFromRowsIfNeeded();
    cols_.resize(new_num_cols.value());
    consistency_ = MatrixConsistency::kRowsNeedRecomputation;
  }
}

void StrongSparseMatrix::RecomputeRowsFromColsIfNeeded() const {
  VLOG(10) << "calling RecomputeRowsFromColsIfNeeded().";
  if (!rows_need_recomputation()) return;
  DCHECK(!cols_need_recomputation());
  PopulateFromTranspose(rows_, cols_);
}

void StrongSparseMatrix::RecomputeColsFromRowsIfNeeded() const {
  VLOG(10) << "calling RecomputeColsFromRowsIfNeeded().";
  if (!cols_need_recomputation()) return;
  DCHECK(!rows_need_recomputation());
  PopulateFromTranspose(cols_, rows_);
}

}  // namespace minimip
