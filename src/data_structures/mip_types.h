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

#ifndef SRC_DATA_STRUCTURES_MIP_TYPES_H_
#define SRC_DATA_STRUCTURES_MIP_TYPES_H_

#include "ortools/util/strong_integers.h"

namespace minimip {

// Constant for infinity used everywhere in MiniMip.
constexpr double kInfinity = std::numeric_limits<double>::infinity();

// These are strong types for columns (aka variables) and rows (aka
// constraints) used everywhere in MiniMip.
DEFINE_STRONG_INDEX_TYPE(ColIndex);
DEFINE_STRONG_INDEX_TYPE(RowIndex);

// Strongly types dense row and col (see "strong_sparse_vector.h" for
// sparse counterparts).
// TODO(lpawel): Use these inside lpi.h (and everywhere).
using DenseRow = absl::StrongVector<ColIndex, double>;
using DenseCol = absl::StrongVector<RowIndex, double>;

template <typename SparseIndex>
constexpr SparseIndex kInvalidSparseIndex(-1);

constexpr RowIndex kInvalidRow = kInvalidSparseIndex<RowIndex>;
constexpr ColIndex kInvalidCol = kInvalidSparseIndex<ColIndex>;

// This representes an indexs which may be either a column or a row of the
// original problem (useful to, e.g., specify basis).
class ColOrRowIndex {
 public:
  explicit ColOrRowIndex(ColIndex only_col)
      : col_(only_col), row_(kInvalidRow) {}
  explicit ColOrRowIndex(RowIndex only_row)
      : col_(kInvalidCol), row_(only_row) {}

  ColIndex col() const { return col_; }
  RowIndex row() const { return row_; }

 private:
  ColIndex col_;
  RowIndex row_;
};

// A struct to contain a single double value for a variable (used to, e.g.,
// specify sparse implied bounds).
struct ColAndValue {
  ColIndex col = kInvalidCol;
  double value = 0.0;
};

}  // namespace minimip

#endif  // SRC_DATA_STRUCTURES_MIP_TYPES_H_
