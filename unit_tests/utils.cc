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

std::string FormatLPBasisStatus(LPBasisStatus status) {
  switch (status) {
    case minimip::LPBasisStatus::kAtLowerBound:
      return "LOWER_BOUND";
    case minimip::LPBasisStatus::kAtUpperBound:
      return "UPPER_BOUND";
    case minimip::LPBasisStatus::kBasic:
      return "BASIC";
    case minimip::LPBasisStatus::kFixed:
      return "FIXED";
    case minimip::LPBasisStatus::kFree:
      return "FREE";
  }
}

}  // namespace

std::string LPModelDebugString(const LPInterface* lpi) {
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

std::string LPStatusDebugString(const LPInterface* lpi) {
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
  const absl::StrongVector<ColIndex, LPBasisStatus> column_status =
      lpi->GetBasisStatusForColumns().value();
  for (ColIndex col(0); col < num_columns; ++col) {
    absl::StrAppend(&s,
                    absl::StrFormat("x%d: %s\n", col.value(),
                                    FormatLPBasisStatus(column_status[col])));
  }
  const absl::StrongVector<RowIndex, LPBasisStatus> row_status =
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