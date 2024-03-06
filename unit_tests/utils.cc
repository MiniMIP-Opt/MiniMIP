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

#include "unit_tests/utils.h"

namespace minimip {
namespace {

std::string FormatDouble(double val) { return absl::StrFormat("% 7.2lf", val); }

std::string FormatRow(const SparseRow& row, int num_columns) {
  std::string res;
  for (ColIndex col(0); col < num_columns; ++col) {
    absl::StrAppend(&res, FormatDouble(row[col]));
  }
  return res;
}

std::string FormatLPBasisStatus(LpBasisStatus status) {
  switch (status) {
    case minimip::LpBasisStatus::kAtLowerBound:
      return "LOWER_BOUND";
    case minimip::LpBasisStatus::kAtUpperBound:
      return "UPPER_BOUND";
    case minimip::LpBasisStatus::kBasic:
      return "BASIC";
    case minimip::LpBasisStatus::kFixed:
      return "FIXED";
    case minimip::LpBasisStatus::kFree:
      return "FREE";
  }
}

}  // namespace

std::string LpModelDebugString(const LpInterface* lpi) {
  const ColIndex num_columns = lpi->GetNumberOfColumns();
  const RowIndex num_rows = lpi->GetNumberOfRows();
  std::string s;

  absl::StrAppend(&s, lpi->IsMaximization() ? "max" : "min", ": ");
  SparseRow obj;
  for (ColIndex col(0); col < num_columns; ++col) {
    obj.AddEntry(col, lpi->GetObjectiveCoefficient(col));
  }
  absl::StrAppend(&s, FormatRow(obj, num_columns.value()), "\n");

  for (RowIndex row(0); row < num_rows; ++row) {
    absl::StrAppend(
        &s,
        absl::StrFormat(
            "%s <= %s <= %s\n", FormatDouble(lpi->GetLeftHandSide(row)),
            FormatRow(lpi->GetSparseRowCoefficients(row), num_columns.value()),
            FormatDouble(lpi->GetRightHandSide(row))));
  }

  for (ColIndex col(0); col < num_columns; ++col) {
    absl::StrAppend(
        &s, absl::StrFormat("%s <= x%d <= %s\n",
                            FormatDouble(lpi->GetLowerBound(col)), col.value(),
                            FormatDouble(lpi->GetUpperBound(col))));
  }
  return s;
}

std::string LpStatusDebugString(const LpInterface* lpi) {
  const ColIndex num_columns = lpi->GetNumberOfColumns();
  const RowIndex num_rows = lpi->GetNumberOfRows();
  std::string s;

  if (!lpi->IsSolved()) return "LP isn't solved.";
  if (!lpi->IsOptimal()) return "LP isn't solved to optimality.";

  absl::StrAppend(
      &s, "objective value: ", FormatDouble(lpi->GetObjectiveValue()), "\n");
  absl::StrAppend(
      &s, "primal values: ",
      FormatRow(SparseRow(lpi->GetPrimalValues().value()), num_columns.value()),
      "\n");
  const absl::StrongVector<ColIndex, LpBasisStatus> column_status =
      lpi->GetBasisStatusForColumns().value();
  for (ColIndex col(0); col < num_columns; ++col) {
    absl::StrAppend(&s,
                    absl::StrFormat("x%d: %s\n", col.value(),
                                    FormatLPBasisStatus(column_status[col])));
  }
  const absl::StrongVector<RowIndex, LpBasisStatus> row_status =
      lpi->GetBasisStatusForRows().value();
  for (RowIndex row(0); row < num_rows; ++row) {
    absl::StrAppend(&s, absl::StrFormat("constraint %d: %s\n", row.value(),
                                        FormatLPBasisStatus(row_status[row])));
  }

  std::vector<ColOrRowIndex> basis = lpi->GetColumnsAndRowsInBasis();
  for (RowIndex row(0); row < num_rows; ++row) {
    const std::string name =
        basis[row.value()].col() == kInvalidCol
            ? absl::StrFormat("constraint %d", basis[row.value()].row().value())
            : absl::StrFormat("x%d", basis[row.value()].col().value());
    absl::StrAppend(
        &s, absl::StrFormat("basis inverse for %s: %s\n", name,
                            FormatRow(lpi->GetSparseRowOfBInverted(row).value(),
                                      num_rows.value())));
  }

  return s;
}

}  // namespace minimip
