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

#include "src/lp_interface/lpi_glop.h"

#include <limits>
#include <memory>

using operations_research::MPModelProto;
using operations_research::TimeLimit;
using operations_research::glop::BasisState;
using operations_research::glop::ColIndexVector;
using operations_research::glop::ConstraintStatus;
using operations_research::glop::ConstraintStatusColumn;
using operations_research::glop::DenseBooleanColumn;
using operations_research::glop::DenseBooleanRow;
using operations_research::glop::DenseColumn;
using operations_research::glop::DenseRow;
using operations_research::glop::Fractional;
using operations_research::glop::GetProblemStatusString;
using operations_research::glop::GlopParameters;
using operations_research::glop::ProblemStatus;
using operations_research::glop::RowToColIndex;
using operations_research::glop::ScatteredColumn;
using operations_research::glop::ScatteredColumnIterator;
using operations_research::glop::ScatteredRow;
using operations_research::glop::ScatteredRowIterator;
using operations_research::glop::SparseColumn;
using operations_research::glop::SparseMatrix;
using operations_research::glop::VariableStatus;
using operations_research::glop::VariableStatusRow;

using GlopColIndex = operations_research::glop::ColIndex;
using GlopRowIndex = operations_research::glop::RowIndex;

namespace minimip {

namespace {

LpBasisStatus ConvertGlopVariableStatus(VariableStatus status,
                                        double reduced_cost) {
  VLOG(2) << "calling ConvertGlopVariableStatus().";
  switch (status) {
    case VariableStatus::BASIC:
      return LpBasisStatus::kBasic;
    case VariableStatus::AT_UPPER_BOUND:
      return LpBasisStatus::kAtUpperBound;
    case VariableStatus::AT_LOWER_BOUND:
      return LpBasisStatus::kAtLowerBound;
    case VariableStatus::FREE:
      return LpBasisStatus::kFree;
    case VariableStatus::FIXED_VALUE:
      return reduced_cost > 0.0 ? LpBasisStatus::kAtLowerBound
                                : LpBasisStatus::kAtUpperBound;
    default:
      LOG(FATAL) << "Unknown Glop column basis status.";
  }
}

LpBasisStatus ConvertGlopConstraintStatus(ConstraintStatus status,
                                          double dual_value) {
  VLOG(2) << "calling ConvertGlopConstraintStatus().";
  switch (status) {
    case ConstraintStatus::BASIC:
      return LpBasisStatus::kBasic;
    case ConstraintStatus::AT_UPPER_BOUND:
      return LpBasisStatus::kAtUpperBound;
    case ConstraintStatus::AT_LOWER_BOUND:
      return LpBasisStatus::kAtLowerBound;
    case ConstraintStatus::FREE:
      return LpBasisStatus::kFree;
    case ConstraintStatus::FIXED_VALUE:
      return dual_value > 0.0 ? LpBasisStatus::kAtLowerBound
                              : LpBasisStatus::kAtUpperBound;
    default:
      LOG(FATAL) << "Unknown Glop row basis status.";
  }
}

VariableStatus ConvertMiniMIPVariableStatus(LpBasisStatus status) {
  VLOG(2) << "calling ConvertMiniMIPVariableStatus().";
  switch (status) {
    case LpBasisStatus::kBasic:
      return VariableStatus::BASIC;
    case LpBasisStatus::kAtUpperBound:
      return VariableStatus::AT_UPPER_BOUND;
    case LpBasisStatus::kAtLowerBound:
      return VariableStatus::AT_LOWER_BOUND;
    case LpBasisStatus::kFixed:
      return VariableStatus::FIXED_VALUE;
    case LpBasisStatus::kFree:
      return VariableStatus::FREE;
    default:
      LOG(FATAL) << "Unknown MiniMip col basis status.";
  }
}

VariableStatus ConvertMiniMIPConstraintStatusToSlackStatus(
    LpBasisStatus status) {
  VLOG(2) << "calling ConvertMiniMIPConstraintStatusToSlackStatus().";
  // We swap lower and upper bound, because Glop adds slacks with -1.0
  // coefficient, whereas LP interface assumes slacks with +1.0 coefficient.
  switch (status) {
    case LpBasisStatus::kBasic:
      return VariableStatus::BASIC;
    case LpBasisStatus::kAtUpperBound:
      return VariableStatus::AT_LOWER_BOUND;
    case LpBasisStatus::kAtLowerBound:
      return VariableStatus::AT_UPPER_BOUND;
    case LpBasisStatus::kFixed:
      return VariableStatus::FIXED_VALUE;
    case LpBasisStatus::kFree:
      return VariableStatus::FREE;
    default:
      LOG(FATAL) << "Unknown MiniMip row basis status.";
  }
}

bool IsDualBoundValid(ProblemStatus status) {
  VLOG(2) << "calling IsDualBoundValid().";
  return status == ProblemStatus::OPTIMAL or
         status == ProblemStatus::DUAL_FEASIBLE or
         status == ProblemStatus::DUAL_UNBOUNDED;
}

}  // namespace

LpGlopInterface::LpGlopInterface()
    : lp_modified_since_last_solve_(true),
      lp_time_limit_was_reached_(false),
      num_iterations_of_last_solve_(0) {
  VLOG(2) << "calling LpGlopInterface().";
  tmp_row_ = std::make_unique<ScatteredRow>();
  tmp_column_ = std::make_unique<ScatteredColumn>();
  // We set the parameters explicitly to the default values, because
  // `SetLpParameters` decides some default values (e.g., default tolerances).
  // This way all parameter values are consistent from the start and do not
  // depend on Glop implementation.
  LpParameters params;
  params.set_lp_solver_type(LpParameters::LP_GLOP);
  CHECK_OK(SetLpParameters(params));
}

// ==========================================================================
// LP model setters.
// ==========================================================================

absl::Status LpGlopInterface::PopulateFromMipData(const MipData& mip_data) {
  VLOG(2) << "calling PopulateFromMipData().";

  RETURN_IF_ERROR(Clear());
  DCHECK_EQ(mip_data.constraint_names().size(),
            mip_data.left_hand_sides().size());
  DCHECK_EQ(mip_data.left_hand_sides().size(),
            mip_data.right_hand_sides().size());
  DCHECK_EQ(mip_data.variable_names().size(), mip_data.lower_bounds().size());
  DCHECK_EQ(mip_data.lower_bounds().size(), mip_data.upper_bounds().size());

  const RowIndex num_rows(mip_data.left_hand_sides().size());
  for (RowIndex row(0); row < num_rows; ++row) {
    RETURN_IF_ERROR(AddRow({}, mip_data.left_hand_sides()[row],
                           mip_data.right_hand_sides()[row],
                           mip_data.constraint_names()[row]));
  }
  for (ColIndex col(0); col < mip_data.matrix().num_cols(); ++col) {
    RETURN_IF_ERROR(
        AddColumn(mip_data.matrix().col(col), mip_data.lower_bounds()[col],
                  mip_data.upper_bounds()[col], mip_data.objective().value(col),
                  mip_data.variable_names()[col]));
  }
  RETURN_IF_ERROR(SetObjectiveSense(false));
  return absl::OkStatus();
}

absl::Status LpGlopInterface::AddColumn(const SparseCol& col_data,
                                        double lower_bound, double upper_bound,
                                        double objective_coefficient,
                                        const std::string& name) {
  VLOG(2) << "calling AddColumn().";
  DCHECK(!col_data.MayNeedCleaning());
  DCHECK(std::all_of(col_data.entries().begin(), col_data.entries().end(),
                     [num_rows = GetNumberOfRows()](const ColEntry& e) {
                       return RowIndex(0) <= e.index and e.index < num_rows;
                     }));

  const GlopColIndex col = lp_.CreateNewVariable();
  VLOG(3) << "Adding column at index = " << col.value()
          << " with num_non_zeros = " << col_data.entries().size() << ".";

  lp_.SetVariableBounds(col, lower_bound, upper_bound);
  lp_.SetObjectiveCoefficient(col, objective_coefficient);
  lp_.SetVariableName(col, name);
  for (const ColEntry& e : col_data.entries()) {
    DCHECK_GE(e.index, 0);
    DCHECK_LT(e.index, GetNumberOfRows());
    DCHECK(!IsInfinity(std::abs(e.value)));
    const GlopRowIndex row(e.index.value());
    lp_.SetCoefficient(row, col, e.value);
  }

  lp_modified_since_last_solve_ = true;
  return absl::OkStatus();
}

absl::Status LpGlopInterface::AddColumns(
    const StrongSparseMatrix& matrix,
    const absl::StrongVector<ColIndex, double>& lower_bounds,
    const absl::StrongVector<ColIndex, double>& upper_bounds,
    const absl::StrongVector<ColIndex, double>& objective_coefficients,
    const absl::StrongVector<ColIndex, std::string>& names) {
  VLOG(2) << "calling AddColumns().";
  DCHECK_EQ(names.size(), lower_bounds.size());
  DCHECK_EQ(lower_bounds.size(), upper_bounds.size());
  DCHECK_EQ(upper_bounds.size(), matrix.num_cols());
  DCHECK_EQ(matrix.num_rows(), GetNumberOfRows());
  for (ColIndex col(0); col < matrix.num_cols(); ++col) {
    RETURN_IF_ERROR(AddColumn(matrix.col(col), lower_bounds[col],
                              upper_bounds[col], objective_coefficients[col],
                              names[col]));
  }
  return absl::OkStatus();
}

absl::Status LpGlopInterface::DeleteColumns(ColIndex first_col,
                                            ColIndex last_col) {
  VLOG(2) << "calling DeleteColumns().";
  DCHECK_GE(first_col, ColIndex(0));
  DCHECK_GE(last_col, first_col);
  DCHECK_LT(last_col, GetNumberOfColumns());

  VLOG(3) << "Deleting columns from " << first_col << " to " << last_col << ".";

  const GlopColIndex num_cols = lp_.num_variables();
  DenseBooleanRow columns_to_delete(num_cols, false);
  for (ColIndex col = first_col; col <= last_col; ++col) {
    columns_to_delete[GlopColIndex(col.value())] = true;
  }
  lp_.DeleteColumns(columns_to_delete);

  lp_modified_since_last_solve_ = true;
  return absl::OkStatus();
}

absl::Status LpGlopInterface::AddRow(const SparseRow& row_data,
                                     double left_hand_side,
                                     double right_hand_side,
                                     const std::string& name) {
  VLOG(2) << "calling AddRow().";
  DCHECK(!row_data.MayNeedCleaning());
  DCHECK(std::all_of(row_data.entries().begin(), row_data.entries().end(),
                     [num_cols = GetNumberOfColumns()](const RowEntry& e) {
                       return ColIndex(0) <= e.index and e.index < num_cols;
                     }));

  const GlopRowIndex row = lp_.CreateNewConstraint();
  VLOG(3) << "Adding row at index = " << row.value()
          << " with num_non_zeros = " << row_data.entries().size() << ".";

  lp_.SetConstraintBounds(row, left_hand_side, right_hand_side);
  lp_.SetConstraintName(row, name);
  for (const RowEntry& e : row_data.entries()) {
    DCHECK_GE(e.index, 0);
    DCHECK_LT(e.index, GetNumberOfColumns());
    DCHECK(!IsInfinity(std::abs(e.value)));
    const GlopColIndex col(e.index.value());
    lp_.SetCoefficient(row, col, e.value);
  }

  lp_modified_since_last_solve_ = true;
  return absl::OkStatus();
}

absl::Status LpGlopInterface::AddRows(
    const absl::StrongVector<RowIndex, SparseRow>& rows,
    const absl::StrongVector<RowIndex, double>& left_hand_sides,
    const absl::StrongVector<RowIndex, double>& right_hand_sides,
    const absl::StrongVector<RowIndex, std::string>& names) {
  VLOG(2) << "calling AddRows().";
  DCHECK_EQ(names.size(), left_hand_sides.size());
  DCHECK_EQ(left_hand_sides.size(), right_hand_sides.size());
  DCHECK_EQ(right_hand_sides.size(), rows.size());
  for (RowIndex row(0); row < rows.size(); ++row) {
    RETURN_IF_ERROR(AddRow(rows[row], left_hand_sides[row],
                           right_hand_sides[row], names[row]));
  }
  return absl::OkStatus();
}

absl::Status LpGlopInterface::DeleteRows(RowIndex first_row,
                                         RowIndex last_row) {
  VLOG(2) << "calling DeleteRows().";
  DCHECK_GE(first_row, RowIndex(0));
  DCHECK_GE(last_row, first_row);
  DCHECK_LT(last_row, GetNumberOfRows());

  VLOG(3) << "Deleting rows from " << first_row << " to " << last_row << ".";

  DenseBooleanColumn glop_rows_to_delete(lp_.num_constraints(), false);
  for (RowIndex row = first_row; row <= last_row; ++row) {
    glop_rows_to_delete[GlopRowIndex(row.value())] = true;
  }
  DeleteRowsAndUpdateCurrentBasis(glop_rows_to_delete);
  return absl::OkStatus();
}

absl::StatusOr<absl::StrongVector<RowIndex, RowIndex>>
LpGlopInterface::DeleteRowSet(
    const absl::StrongVector<RowIndex, bool>& rows_to_delete) {
  VLOG(2) << "calling DeleteRowSet().";
  DCHECK_EQ(rows_to_delete.size(), GetNumberOfRows());
  DenseBooleanColumn glop_rows_to_delete(rows_to_delete.begin(),
                                         rows_to_delete.end());
  absl::StrongVector<RowIndex, RowIndex> row_mapping(rows_to_delete.size());
  RowIndex next_index(0);
  int num_deleted_rows = 0;
  for (RowIndex row(0); row < GetNumberOfRows(); ++row) {
    if (rows_to_delete[row]) {
      row_mapping[row] = kInvalidRow;
      ++num_deleted_rows;
    } else {
      row_mapping[row] = next_index++;
    }
  }
  VLOG(3) << "Deleting row set of size = " << num_deleted_rows << ".";
  DeleteRowsAndUpdateCurrentBasis(glop_rows_to_delete);
  return row_mapping;
}

void LpGlopInterface::DeleteRowsAndUpdateCurrentBasis(
    const DenseBooleanColumn& rows_to_delete) {
  VLOG(2) << "calling DeleteRowsAndUpdateCurrentBasis().";
  const GlopRowIndex num_rows = lp_.num_constraints();
  const GlopColIndex num_cols = lp_.num_variables();

  // Try to repair basis status if problem size has not changed.
  BasisState state = solver_.GetState();
  if (state.statuses.size() == num_cols.value() + num_rows.value()) {
    // Shift the status of the non-deleted rows. Note that if a deleted row
    // was in the basis (i.e., the constraint was not tight in the solution)
    // then we should still be left with a correct basis after deletion.
    // This should be the most common use case in MiniMIP.
    GlopColIndex new_size = num_cols;
    for (GlopRowIndex row(0); row < num_rows; ++row) {
      if (rows_to_delete[row]) continue;
      state.statuses[new_size++] =
          state.statuses[num_cols + RowToColIndex(row)];
    }
    state.statuses.resize(new_size);
    solver_.LoadStateForNextSolve(state);
  }

  lp_.DeleteRows(rows_to_delete);
  lp_modified_since_last_solve_ = true;
}

absl::Status LpGlopInterface::Clear() {
  VLOG(2) << "calling LpGlopInterface::Clear().";
  lp_.Clear();
  lp_modified_since_last_solve_ = true;
  return absl::OkStatus();
}

absl::Status LpGlopInterface::ClearState() {
  VLOG(2) << "calling LpGlopInterface::ClearState().";
  solver_.ClearStateForNextSolve();
  return absl::OkStatus();
}

absl::Status LpGlopInterface::SetColumnBounds(ColIndex col, double lower_bound,
                                              double upper_bound) {
  VLOG(2) << "calling SetColumnBounds().";
  DCHECK_GE(col, 0);
  DCHECK_LT(col, GetNumberOfColumns());
  DCHECK(!IsInfinity(lower_bound));
  DCHECK(!IsInfinity(-upper_bound));

  VLOG(3) << "Set column bounds: col=" << col << ", lower_bound=" << lower_bound
          << ", upper_bound=" << upper_bound << ".";
  lp_.SetVariableBounds(GlopColIndex(col.value()), lower_bound, upper_bound);
  lp_modified_since_last_solve_ = true;
  return absl::OkStatus();
}

absl::Status LpGlopInterface::SetRowSides(RowIndex row, double left_hand_side,
                                          double right_hand_side) {
  VLOG(2) << "calling SetRowSides().";
  DCHECK_GE(row, 0);
  DCHECK_LT(row, GetNumberOfRows());
  DCHECK(!IsInfinity(left_hand_side));
  DCHECK(!IsInfinity(-right_hand_side));
  DCHECK_LE(left_hand_side, right_hand_side);
  VLOG(3) << "Set row sides: row=" << row
          << ", left_hand_side=" << left_hand_side
          << ", right_hand_side=" << right_hand_side << ".";
  lp_.SetConstraintBounds(GlopRowIndex(row.value()), left_hand_side,
                          right_hand_side);
  lp_modified_since_last_solve_ = true;
  return absl::OkStatus();
}

absl::Status LpGlopInterface::SetObjectiveSense(bool is_maximization) {
  VLOG(2) << "calling SetObjectiveSense().";
  VLOG(3) << "Setting maximize=" << is_maximization << ".";
  lp_.SetMaximizationProblem(is_maximization);
  lp_modified_since_last_solve_ = true;
  return absl::OkStatus();
}

absl::Status LpGlopInterface::SetObjectiveCoefficient(
    ColIndex col, double objective_coefficient) {
  VLOG(2) << "calling SetObjectiveCoefficient().";
  DCHECK_GE(col, 0);
  DCHECK_LT(col, GetNumberOfColumns());
  DCHECK(!IsInfinity(std::abs(objective_coefficient)));
  VLOG(3) << "Setting objective coefficient, col=" << col.value()
          << ", objective_coefficient=" << objective_coefficient << ".";
  lp_.SetObjectiveCoefficient(GlopColIndex(col.value()), objective_coefficient);
  lp_modified_since_last_solve_ = true;
  return absl::OkStatus();
}

// ==========================================================================
// LP model getters.
// ==========================================================================

RowIndex LpGlopInterface::GetNumberOfRows() const {
  VLOG(2) << "calling GetNumberOfRows().";
  return RowIndex(lp_.num_constraints().value());
}

ColIndex LpGlopInterface::GetNumberOfColumns() const {
  VLOG(2) << "calling GetNumberOfColumns().";
  return ColIndex(lp_.num_variables().value());
}

bool LpGlopInterface::IsMaximization() const {
  VLOG(2) << "calling IsMaximization().";
  return lp_.IsMaximizationProblem();
}

int64_t LpGlopInterface::GetNumberOfNonZeros() const {
  VLOG(2) << "calling GetNumberOfNonZeros().";
  return lp_.num_entries().value();
}

SparseCol LpGlopInterface::GetSparseColumnCoefficients(ColIndex col) const {
  VLOG(2) << "calling GetSparseColumnCoefficients().";
  DCHECK_GE(col, ColIndex(0));
  DCHECK_LT(col, GetNumberOfColumns());
  const SparseColumn& column_in_glop =
      lp_.GetSparseColumn(GlopColIndex(col.value()));
  SparseCol col_data;
  col_data.mutable_entries().resize(column_in_glop.num_entries().value());
  // Accessing `mutable_entries()` forces us to clean up -- even if the sparse
  // col is still empty.
  col_data.CleanUpIfNeeded();
  for (const auto& entry : column_in_glop) {
    col_data.AddEntry(RowIndex(entry.row().value()), entry.coefficient());
  }
  DCHECK(!col_data.MayNeedCleaning());
  return col_data;
}

SparseRow LpGlopInterface::GetSparseRowCoefficients(RowIndex row) const {
  VLOG(2) << "calling GetSparseRowCoefficients().";
  DCHECK_GE(row, RowIndex(0));
  DCHECK_LT(row, GetNumberOfRows());
  // Note, there is no casting from col to row in Glop, hence we keep the row
  // (grabbed from the transposed matrix) in the column type.
  const SparseColumn& row_in_glop =
      lp_.GetTransposeSparseMatrix().column(GlopColIndex(row.value()));
  SparseRow row_data;
  row_data.mutable_entries().resize(row_in_glop.num_entries().value());
  // Accessing `mutable_entries()` forces us to clean up -- even if the sparse
  // row is still empty.
  row_data.CleanUpIfNeeded();
  for (const auto& entry : row_in_glop) {
    row_data.AddEntry(ColIndex(entry.row().value()), entry.coefficient());
  }
  DCHECK(!row_data.MayNeedCleaning());
  return row_data;
}

double LpGlopInterface::GetObjectiveCoefficient(ColIndex col) const {
  VLOG(2) << "calling GetObjectiveCoefficient().";
  CHECK_GE(col, 0);
  CHECK_LT(col, GetNumberOfColumns());
  return lp_.objective_coefficients()[GlopColIndex(col.value())];
}

double LpGlopInterface::GetLowerBound(ColIndex col) const {
  VLOG(2) << "calling GetLowerBound().";
  CHECK_GE(col, 0);
  CHECK_LT(col, GetNumberOfColumns());
  return lp_.variable_lower_bounds()[GlopColIndex(col.value())];
}

double LpGlopInterface::GetUpperBound(ColIndex col) const {
  VLOG(2) << "calling GetUpperBound().";
  CHECK_GE(col, 0);
  CHECK_LT(col, GetNumberOfColumns());
  return lp_.variable_upper_bounds()[GlopColIndex(col.value())];
}

double LpGlopInterface::GetLeftHandSide(RowIndex row) const {
  VLOG(2) << "calling GetLeftHandSide().";
  CHECK_GE(row, 0);
  CHECK_LT(row, GetNumberOfRows());
  return lp_.constraint_lower_bounds()[GlopRowIndex(row.value())];
}

double LpGlopInterface::GetRightHandSide(RowIndex row) const {
  VLOG(2) << "calling GetRightHandSide().";
  CHECK_GE(row, 0);
  CHECK_LT(row, GetNumberOfRows());
  return lp_.constraint_upper_bounds()[GlopRowIndex(row.value())];
}

double LpGlopInterface::GetMatrixCoefficient(ColIndex col, RowIndex row) const {
  VLOG(2) << "calling GetMatrixCoefficient().";
  CHECK_GE(col, 0);
  CHECK_LT(col, GetNumberOfColumns());
  CHECK_GE(row, 0);
  CHECK_LT(row, GetNumberOfRows());
  // Caution: Runs in O(num_entries_in_col).
  return lp_.GetSparseMatrix().LookUpValue(GlopRowIndex(row.value()),
                                           GlopColIndex(col.value()));
}

// ============================================================================
// Internal solving methods.
// ============================================================================

// NOLINTNEXTLINE(misc-no-recursion)
absl::Status LpGlopInterface::SolveInternal(bool recursive,
                                            TimeLimit* time_limit) {
  VLOG(2) << "calling SolveInternal().";
  // Recompute `scaled_lp_`.
  if (lp_modified_since_last_solve_) {
    // TODO(lpawel): Avoid doing a copy if there is no scaling.
    scaled_lp_.PopulateFromLinearProgram(lp_);
    scaled_lp_.AddSlackVariablesWhereNecessary(false);

    if (solver_.GetParameters().use_scaling()) {
      // TODO(lpawel): Avoid rescaling if nothing changed.
      scaler_.Scale(&scaled_lp_);
    } else {
      scaler_.Clear();
    }
  }

  if (!recursive) {
    num_iterations_of_last_solve_ = 0;
    lp_time_limit_was_reached_ = false;
  }

  if (solve_from_scratch_) {
    solver_.ClearStateForNextSolve();
  }

  const operations_research::glop::Status glop_solve_status =
      solver_.Solve(scaled_lp_, time_limit);
  if (!glop_solve_status.ok()) {
    return absl::Status(absl::StatusCode::kInternal,
                        absl::StrFormat("LP Error: code=%d, msg=%s",
                                        glop_solve_status.error_code(),
                                        glop_solve_status.error_message()));
  }

  lp_time_limit_was_reached_ = time_limit->LimitReached();
  num_iterations_of_last_solve_ += solver_.GetNumberOfIterations();

  const ProblemStatus status = solver_.GetProblemStatus();
  VLOG(3) << "Glop Solve() finished, "
          << "status=" << GetProblemStatusString(status)
          << ", obj=" << solver_.GetObjectiveValue()
          << ", iterations=" << solver_.GetNumberOfIterations();

  // In case the solution is not feasible wrt original problem, we will attempt
  // to solve the unscaled version from scratch.
  if ((status == ProblemStatus::PRIMAL_FEASIBLE or
       status == ProblemStatus::OPTIMAL) &&
      solver_.GetParameters().use_scaling()) {
    const auto primal_values = GetPrimalValues();
    CHECK_OK(primal_values);
    const operations_research::glop::DenseRow solution(primal_values->begin(),
                                                       primal_values->end());
    if (!lp_.SolutionIsLPFeasible(
            solution, solver_.GetParameters().primal_feasibility_tolerance())) {
      VLOG(1) << "Solution not feasible w.r.t. absolute tolerance "
              << solver_.GetParameters().primal_feasibility_tolerance()
              << ". Will re-optimize.";
      GlopParameters glop_params = solver_.GetParameters();
      glop_params.set_use_scaling(false);
      solver_.SetParameters(glop_params);

      // This is needed to force setting `scaled_lp_ = lp_` in the recursive
      // call.
      // TODO(lpawel): Fix this logic and make it more explicit.
      lp_modified_since_last_solve_ = true;

      RETURN_IF_ERROR(SolveInternal(/*recursive=*/true, time_limit));
      glop_params.set_use_scaling(true);
      solver_.SetParameters(glop_params);
    }
  }

  lp_modified_since_last_solve_ = false;
  return absl::OkStatus();
}

// ==========================================================================
// Solving methods.
// ==========================================================================

absl::Status LpGlopInterface::SolveLpWithPrimalSimplex() {
  VLOG(2) << "calling SolveLpWithPrimalSimplex().";
  VLOG(3) << "Solving with primal simplex: "
          << "num_cols=" << lp_.num_variables().value()
          << ", num_rows=" << lp_.num_constraints().value();

  GlopParameters glop_params = solver_.GetParameters();
  std::unique_ptr<TimeLimit> time_limit =
      TimeLimit::FromParameters(glop_params);
  glop_params.set_use_dual_simplex(false);
  solver_.SetParameters(glop_params);
  return SolveInternal(/*recursive=*/false, time_limit.get());
}

absl::Status LpGlopInterface::SolveLpWithDualSimplex() {
  VLOG(2) << "calling SolveLpWithDualSimplex().";
  VLOG(3) << "Solving with dual simplex: "
          << "num_cols=" << lp_.num_variables().value()
          << ", num_rows=" << lp_.num_constraints().value();
  GlopParameters glop_params = solver_.GetParameters();
  std::unique_ptr<TimeLimit> time_limit =
      TimeLimit::FromParameters(glop_params);
  glop_params.set_use_dual_simplex(true);
  solver_.SetParameters(glop_params);
  return SolveInternal(/*recursive=*/false, time_limit.get());
}

absl::Status LpGlopInterface::StartStrongBranching() {
  VLOG(2) << "calling StartStrongBranching().";
  // TODO(lpawel): Save Glop state and tune Glop towards strong branching (and
  // avoid rescaling completely when in strong branching).
  return absl::OkStatus();
}

absl::Status LpGlopInterface::EndStrongBranching() {
  VLOG(10) << "calling EndStrongBranching().";
  // TODO(lpawel): Restore the saved Glop state.
  return absl::OkStatus();
}

absl::StatusOr<LpGlopInterface::StrongBranchResult>
LpGlopInterface::SolveDownAndUpStrongBranch(ColIndex col, double primal_value,
                                            int iteration_limit) {
  VLOG(10) << "calling SolveDownAndUpStrongBranch().";
  DCHECK_GE(col, ColIndex(0));
  DCHECK_LT(col, GetNumberOfColumns());

  VLOG(3) << "Solving down and up strong branches: col=" << col.value()
          << ", iteration_limit=" << iteration_limit;

  // We work on the scaled problem.
  const GlopColIndex glop_col(col.value());
  const Fractional lower_bound = scaled_lp_.variable_lower_bounds()[glop_col];
  const Fractional upper_bound = scaled_lp_.variable_upper_bounds()[glop_col];
  const double scaled_primal_value =
      primal_value * scaler_.VariableScalingFactor(glop_col);

  // Configure solver.
  // TODO(lpawel): Use the iteration limit and incrementality.
  GlopParameters glop_params = solver_.GetParameters();
  glop_params.set_use_dual_simplex(true);
  solver_.SetParameters(glop_params);
  const Fractional eps = glop_params.primal_feasibility_tolerance();

  std::unique_ptr<TimeLimit> time_limit =
      TimeLimit::FromParameters(glop_params);

  StrongBranchResult result;

  // Down branch.
  const double new_upper_bound = std::ceil(scaled_primal_value - 1.0 - eps);
  if (new_upper_bound >= lower_bound - 0.5) {
    scaled_lp_.SetVariableBounds(glop_col, lower_bound, new_upper_bound);
    if (solver_.Solve(scaled_lp_, time_limit.get()).ok()) {
      result.iterations += solver_.GetNumberOfIterations();
      result.dual_bound_down_branch = solver_.GetObjectiveValue();
      result.down_valid = IsDualBoundValid(solver_.GetProblemStatus());
    } else {
      result.dual_bound_down_branch = 0.0;
      result.down_valid = false;
    }
  } else {
    result.down_valid = true;
    result.dual_bound_down_branch = lp_.IsMaximizationProblem()
                                        ? glop_params.objective_lower_limit()
                                        : glop_params.objective_upper_limit();
  }

  // Up branch.
  const double new_lower_bound = std::floor(scaled_primal_value + 1.0 + eps);
  if (new_lower_bound <= upper_bound + 0.5) {
    scaled_lp_.SetVariableBounds(glop_col, new_lower_bound, upper_bound);
    if (solver_.Solve(scaled_lp_, time_limit.get()).ok()) {
      result.iterations += solver_.GetNumberOfIterations();
      result.dual_bound_up_branch = solver_.GetObjectiveValue();
      result.up_valid = IsDualBoundValid(solver_.GetProblemStatus());
    } else {
      result.dual_bound_up_branch = 0.0;
      result.down_valid = false;
    }
  } else {
    result.up_valid = true;
    result.dual_bound_up_branch = lp_.IsMaximizationProblem()
                                      ? glop_params.objective_lower_limit()
                                      : glop_params.objective_upper_limit();
  }

  //  Restore the bounds.
  scaled_lp_.SetVariableBounds(glop_col, lower_bound, upper_bound);

  return result;
}

// ==========================================================================
// Solution information getters.
// ==========================================================================

bool LpGlopInterface::IsSolved() const {
  VLOG(10) << "calling IsSolved().";
  // TODO(lpawel): Track this to avoid unneeded resolving.
  return (!lp_modified_since_last_solve_);
}

bool LpGlopInterface::ExistsPrimalRay() const {
  VLOG(10) << "calling ExistsPrimalRay().";
  return solver_.GetProblemStatus() == ProblemStatus::PRIMAL_UNBOUNDED;
}

bool LpGlopInterface::HasPrimalRay() const {
  VLOG(10) << "calling HasPrimalRay().";
  return solver_.GetProblemStatus() == ProblemStatus::PRIMAL_UNBOUNDED;
}

bool LpGlopInterface::IsPrimalUnbounded() const {
  VLOG(10) << "calling IsPrimalUnbounded().";
  return solver_.GetProblemStatus() == ProblemStatus::PRIMAL_UNBOUNDED;
}

bool LpGlopInterface::IsPrimalInfeasible() const {
  VLOG(10) << "calling IsPrimalInfeasible().";
  const ProblemStatus status = solver_.GetProblemStatus();
  return status == ProblemStatus::DUAL_UNBOUNDED or
         status == ProblemStatus::PRIMAL_INFEASIBLE;
}

bool LpGlopInterface::IsPrimalFeasible() const {
  VLOG(10) << "calling IsPrimalFeasible().";
  const ProblemStatus status = solver_.GetProblemStatus();
  return status == ProblemStatus::PRIMAL_FEASIBLE or
         status == ProblemStatus::OPTIMAL;
}

bool LpGlopInterface::ExistsDualRay() const {
  VLOG(10) << "calling ExistsDualRay().";
  return solver_.GetProblemStatus() == ProblemStatus::DUAL_UNBOUNDED;
}

bool LpGlopInterface::HasDualRay() const {
  VLOG(10) << "calling HasDualRay().";
  return solver_.GetProblemStatus() == ProblemStatus::DUAL_UNBOUNDED;
}

bool LpGlopInterface::IsDualUnbounded() const {
  VLOG(10) << "calling IsDualUnbounded().";
  return solver_.GetProblemStatus() == ProblemStatus::DUAL_UNBOUNDED;
}

bool LpGlopInterface::IsDualInfeasible() const {
  VLOG(10) << "calling IsDualInfeasible().";
  const ProblemStatus status = solver_.GetProblemStatus();
  return status == ProblemStatus::PRIMAL_UNBOUNDED or
         status == ProblemStatus::DUAL_INFEASIBLE;
}

bool LpGlopInterface::IsDualFeasible() const {
  VLOG(10) << "calling IsDualFeasible().";
  const ProblemStatus status = solver_.GetProblemStatus();
  return status == ProblemStatus::DUAL_FEASIBLE or
         status == ProblemStatus::OPTIMAL;
}

bool LpGlopInterface::IsOptimal() const {
  VLOG(10) << "calling IsOptimal().";
  return solver_.GetProblemStatus() == ProblemStatus::OPTIMAL;
}

bool LpGlopInterface::IsStable() const {
  VLOG(10) << "calling IsStable().";
  // For correctness, we need to report "unstable" if Glop was not able to
  // prove optimality because of numerical issues. Currently, Glop still
  // reports primal/dual feasible if at the end, one status is within the
  // tolerance but not the other.
  const ProblemStatus status = solver_.GetProblemStatus();
  if ((status == ProblemStatus::PRIMAL_FEASIBLE or
       status == ProblemStatus::DUAL_FEASIBLE) &&
      !ObjectiveLimitIsExceeded() and !IterationLimitIsExceeded() &&
      !TimeLimitIsExceeded()) {
    VLOG(3) << "OPTIMAL not reached and no limit: unstable";
    return false;
  }

  if (status == ProblemStatus::ABNORMAL or
      status == ProblemStatus::INVALID_PROBLEM or
      status == ProblemStatus::IMPRECISE) {
    VLOG(3) << "Errors while solving: unstable";
    return false;
  }
  return true;
}

bool LpGlopInterface::ObjectiveLimitIsExceeded() const {
  VLOG(10) << "calling ObjectiveLimitIsExceeded().";
  return solver_.objective_limit_reached();
}

bool LpGlopInterface::TimeLimitIsExceeded() const {
  VLOG(10) << "calling TimeLimitIsExceeded().";
  return lp_time_limit_was_reached_;
}

bool LpGlopInterface::IterationLimitIsExceeded() const {
  VLOG(10) << "calling IterationLimitIsExceeded().";
  // We might have accumulated iterations across 2 recursive solves,
  // hence _GE, and not _EQ.
  DCHECK_GE(num_iterations_of_last_solve_, solver_.GetNumberOfIterations());
  return solver_.GetParameters().max_number_of_iterations() != -1 &&
         num_iterations_of_last_solve_ >=
             solver_.GetParameters().max_number_of_iterations();
}

int64_t LpGlopInterface::GetNumIterations() const {
  VLOG(10) << "calling GetNumIterations().";
  return num_iterations_of_last_solve_;
}

double LpGlopInterface::GetObjectiveValue() const {
  VLOG(10) << "calling GetObjectiveValue().";
  DCHECK(IsOptimal());
  return solver_.GetObjectiveValue();
}

absl::StatusOr<absl::StrongVector<ColIndex, double>>
LpGlopInterface::GetPrimalValues() const {
  VLOG(10) << "calling GetPrimalValues().";
  DCHECK(IsOptimal());
  absl::StrongVector<ColIndex, double> primal_values;
  primal_values.reserve(lp_.num_variables().value());
  for (GlopColIndex col(0); col < lp_.num_variables(); ++col) {
    primal_values.push_back(
        scaler_.UnscaleVariableValue(col, solver_.GetVariableValue(col)));
  }
  return primal_values;
}

absl::StatusOr<absl::StrongVector<RowIndex, double>>
LpGlopInterface::GetDualValues() const {
  VLOG(10) << "calling GetDualValues().";
  DCHECK(IsOptimal());
  absl::StrongVector<RowIndex, double> dual_values;
  dual_values.reserve(lp_.num_constraints().value());
  for (GlopRowIndex row(0); row < lp_.num_constraints(); ++row) {
    dual_values.push_back(
        scaler_.UnscaleDualValue(row, solver_.GetDualValue(row)));
  }
  return dual_values;
}

absl::StatusOr<absl::StrongVector<ColIndex, double>>
LpGlopInterface::GetReducedCosts() const {
  VLOG(10) << "calling GetReducedCosts().";
  DCHECK(IsOptimal());
  absl::StrongVector<ColIndex, double> reduced_costs;
  reduced_costs.reserve(lp_.num_variables().value());
  for (GlopColIndex col(0); col < lp_.num_variables(); ++col) {
    reduced_costs.push_back(
        scaler_.UnscaleReducedCost(col, solver_.GetReducedCost(col)));
  }
  return reduced_costs;
}

absl::StatusOr<absl::StrongVector<RowIndex, double>>
LpGlopInterface::GetRowActivities() const {
  VLOG(10) << "calling GetRowActivities().";
  DCHECK(IsOptimal());
  absl::StrongVector<RowIndex, double> row_activities;
  row_activities.reserve(lp_.num_constraints().value());
  for (GlopRowIndex row(0); row < lp_.num_constraints(); ++row) {
    row_activities.push_back(scaler_.UnscaleConstraintActivity(
        row, solver_.GetConstraintActivity(row)));
  }
  return row_activities;
}

absl::StatusOr<absl::StrongVector<ColIndex, double>>
LpGlopInterface::GetPrimalRay() const {
  VLOG(10) << "calling GetPrimalRay().";
  DCHECK(HasPrimalRay());
  absl::StrongVector<ColIndex, double> primal_ray;
  primal_ray.reserve(lp_.num_variables().value());
  for (GlopColIndex col(0); col < lp_.num_variables(); ++col) {
    // Note, Glop adds slacks with -1.0 coefficient. LP interface assumes slacks
    // added with +1.0 coefficient. Hence, we need to reverse the sign.
    primal_ray.push_back(
        -scaler_.UnscaleVariableValue(col, solver_.GetPrimalRay()[col]));
  }
  return primal_ray;
}

absl::StatusOr<absl::StrongVector<RowIndex, double>>
LpGlopInterface::GetDualRay() const {
  VLOG(10) << "calling GetDualRay().";
  DCHECK(HasDualRay());
  absl::StrongVector<RowIndex, double> dual_ray;
  dual_ray.reserve(lp_.num_constraints().value());
  for (GlopRowIndex row(0); row < lp_.num_constraints(); ++row) {
    // Note, Glop adds slacks with -1.0 coefficient. LP interface assumes slacks
    // added with +1.0 coefficient. Hence, we need to reverse the sign.
    dual_ray.push_back(
        -scaler_.UnscaleDualValue(row, solver_.GetDualRay()[row]));
  }
  return dual_ray;
}

// ==========================================================================
// Getters and setters of the basis.
// ==========================================================================

absl::StatusOr<absl::StrongVector<ColIndex, LpBasisStatus>>
LpGlopInterface::GetBasisStatusForColumns() const {
  VLOG(10) << "calling GetBasisStatusForColumns().";
  DCHECK(IsOptimal());
  absl::StrongVector<ColIndex, LpBasisStatus> statuses;
  statuses.reserve(lp_.num_variables().value());
  for (GlopColIndex col(0); col < lp_.num_variables(); ++col) {
    statuses.push_back(ConvertGlopVariableStatus(solver_.GetVariableStatus(col),
                                                 solver_.GetReducedCost(col)));
  }
  return statuses;
}

absl::StatusOr<absl::StrongVector<RowIndex, LpBasisStatus>>
LpGlopInterface::GetBasisStatusForRows() const {
  VLOG(10) << "calling GetBasisStatusForRows().";
  DCHECK(IsOptimal());
  absl::StrongVector<RowIndex, LpBasisStatus> statuses;
  statuses.reserve(lp_.num_constraints().value());
  for (GlopRowIndex row(0); row < lp_.num_constraints(); ++row) {
    statuses.push_back(ConvertGlopConstraintStatus(
        solver_.GetConstraintStatus(row), solver_.GetDualValue(row)));
  }
  return statuses;
}

absl::Status LpGlopInterface::SetBasisStatusForColumnsAndRows(
    const absl::StrongVector<ColIndex, LpBasisStatus>& column_basis_statuses,
    const absl::StrongVector<RowIndex, LpBasisStatus>& row_basis_statuses) {
  VLOG(10) << "calling SetBasisStatusForColumnsAndRows().";
  BasisState state;
  state.statuses.reserve(lp_.num_variables() +
                         RowToColIndex(lp_.num_constraints()));

  for (ColIndex col(0); col < GetNumberOfColumns(); ++col) {
    state.statuses.push_back(
        ConvertMiniMIPVariableStatus(column_basis_statuses[col]));
  }
  for (RowIndex row(0); row < GetNumberOfRows(); ++row) {
    state.statuses.push_back(
        ConvertMiniMIPConstraintStatusToSlackStatus(row_basis_statuses[row]));
  }

  solver_.LoadStateForNextSolve(state);
  return absl::OkStatus();
}

std::vector<ColOrRowIndex> LpGlopInterface::GetColumnsAndRowsInBasis() const {
  VLOG(10) << "calling GetColumnsAndRowsInBasis().";
  std::vector<ColOrRowIndex> basis;
  basis.reserve(GetNumberOfRows().value());
  // The order in which we populate the `basis` is important!
  for (GlopRowIndex row(0); row < lp_.num_constraints(); ++row) {
    const GlopColIndex col = solver_.GetBasis(row);
    VLOG(3) << "solver basis index: " << col;
    basis.push_back(col < lp_.num_variables()
                        ? ColOrRowIndex(ColIndex(col.value()))
                        : ColOrRowIndex(RowIndex(col.value() -
                                                 lp_.num_variables().value())));
    VLOG(3) << "Basis index: " << basis.back() << " from GlopRowIndex: " << row;
  }
  return basis;
}

// ==========================================================================
// Getters of vectors in the inverted basis matrix.
// ==========================================================================

absl::StatusOr<SparseRow> LpGlopInterface::GetSparseRowOfBInverted(
    RowIndex row_in_basis) const {
  VLOG(10) << "calling GetSparseRowOfBInverted().";
  SparseRow sparse_row;

  solver_.GetBasisFactorization().LeftSolveForUnitRow(
      GlopColIndex(row_in_basis.value()), tmp_row_.get());
  scaler_.UnscaleUnitRowLeftSolve(
      solver_.GetBasis(GlopRowIndex(row_in_basis.value())), tmp_row_.get());

  // DCHECK_EQ(tmp_row_->values.size(), lp_.num_constraints().value());
  absl::StrongVector<RowIndex, double> row_activities = GetRowActivities().value();
  LOG(INFO) << "Row activities: ";
  for (RowIndex row_index(0); row_index < row_activities.size(); ++row_index) {
    LOG(INFO) << "Row " << row_index << ": " << row_activities[row_index];
  }
  // Vectors in Glop might be stored in dense or sparse format depending on
  // the values. If non_zeros are given, we can directly loop over the
  // non_zeros, otherwise we have to collect the nonzeros.
  if (!tmp_row_->non_zeros.empty()) {
    ScatteredRowIterator end = tmp_row_->end();
    for (ScatteredRowIterator iter = tmp_row_->begin(); iter != end; ++iter) {
      const int idx = (*iter).column().value();
      DCHECK_GE(idx, 0);
      DCHECK_LT(idx, lp_.num_constraints().value());
      sparse_row.AddEntry(ColIndex(idx), (*iter).coefficient());
    }
  } else {
    const Fractional eps =
        solver_.GetParameters().primal_feasibility_tolerance();
    for (GlopColIndex col(0); col < RowToColIndex(lp_.num_constraints());
         ++col) {
      const double value = (*tmp_row_)[col];
      if (std::abs(value) >= eps) {
        sparse_row.AddEntry(ColIndex(col.value()), value);
      }
    }
  }
  VLOG(3) << "Sparse row: " << sparse_row;
  WriteLpToFile("lp_g.txt");
  return sparse_row;
}

absl::StatusOr<SparseCol> LpGlopInterface::GetSparseColumnOfBInverted(
    ColIndex col_in_basis) const {
  VLOG(10) << "calling GetSparseColumnOfBInverted().";
  SparseCol sparse_column;
  // We need to loop through the rows to extract the values for `col_in_basis`.
  for (GlopRowIndex row(0); row < lp_.num_constraints(); ++row) {
    solver_.GetBasisFactorization().LeftSolveForUnitRow(RowToColIndex(row),
                                                        tmp_row_.get());
    scaler_.UnscaleUnitRowLeftSolve(solver_.GetBasis(row), tmp_row_.get());

    double value = (*tmp_row_)[GlopColIndex(col_in_basis.value())];
    if (std::abs(value) >=
        solver_.GetParameters().primal_feasibility_tolerance()) {
      sparse_column.AddEntry(RowIndex(row.value()), value);
    }
  }
  return sparse_column;
}

absl::StatusOr<SparseRow> LpGlopInterface::GetSparseRowOfBInvertedTimesA(
    RowIndex row_in_basis) const {
  VLOG(10) << "calling GetSparseRowOfBInvertedTimesA().";
  SparseRow sparse_row;
  solver_.GetBasisFactorization().LeftSolveForUnitRow(
      GlopColIndex(row_in_basis.value()), tmp_row_.get());
  scaler_.UnscaleUnitRowLeftSolve(
      solver_.GetBasis(GlopRowIndex(row_in_basis.value())), tmp_row_.get());
  for (GlopColIndex col(0); col < lp_.num_variables(); ++col) {
    double value = operations_research::glop::ScalarProduct(
        tmp_row_->values, lp_.GetSparseColumn(col));
    if (std::abs(value) >=
        solver_.GetParameters().primal_feasibility_tolerance()) {
      sparse_row.AddEntry(ColIndex(col.value()), value);
    }
  }
  return sparse_row;
}

absl::StatusOr<SparseCol> LpGlopInterface::GetSparseColumnOfBInvertedTimesA(
    ColIndex col_in_basis) const {
  VLOG(10) << "calling GetSparseColumnOfBInvertedTimesA().";
  SparseCol sparse_column;
  solver_.GetBasisFactorization().RightSolveForProblemColumn(
      GlopColIndex(col_in_basis.value()), tmp_column_.get());
  scaler_.UnscaleColumnRightSolve(solver_.GetBasisVector(),
                                  GlopColIndex(col_in_basis.value()),
                                  tmp_column_.get());

  // DCHECK_EQ(tmp_row_->values.size(), lp_.num_constraints());

  // If `non_zeros` are present, we have a sparse representation and can just
  // iterate across entries.
  if (!tmp_column_->non_zeros.empty()) {
    ScatteredColumnIterator end = tmp_column_->end();
    for (ScatteredColumnIterator iter = tmp_column_->begin(); iter != end;
         ++iter) {
      const RowIndex row_in_basis((*iter).row().value());
      DCHECK_GE(row_in_basis, 0);
      DCHECK_LT(row_in_basis, GetNumberOfRows());
      sparse_column.AddEntry(row_in_basis, (*iter).coefficient());
    }
  } else {
    // Otherwise, we need to iterate through all rows.
    for (GlopRowIndex row(0); row < lp_.num_constraints(); ++row) {
      const double value = (*tmp_column_)[row];
      if (std::abs(value) >
          solver_.GetParameters().primal_feasibility_tolerance()) {
        const RowIndex row_in_basis(row.value());
        sparse_column.AddEntry(row_in_basis, value);
      }
    }
  }
  return sparse_column;
}

// ==========================================================================
// Getters and setters of the parameters.
// ==========================================================================

LpParameters LpGlopInterface::GetLpParameters() const {
  VLOG(10) << "calling GetLpParameters().";
  LpParameters params;
  params.set_lp_solver_type(LpParameters::LP_GLOP);

  params.set_solve_from_scratch(solve_from_scratch_);

  const GlopParameters& glop_params = solver_.GetParameters();
  if (glop_params.use_scaling()) {
    switch (glop_params.scaling_method()) {
      case GlopParameters::DEFAULT:
        params.set_scaling_strategy(LpParameters::SCALING_DEFAULT);
        break;
      case GlopParameters::EQUILIBRATION:
        params.set_scaling_strategy(LpParameters::SCALING_EQUILIBRATION);
        break;
      case GlopParameters::LINEAR_PROGRAM:
        params.set_scaling_strategy(LpParameters::SCALING_LINEAR_PROGRAM);
        break;
    }
  } else {
    params.set_scaling_strategy(LpParameters::SCALING_OFF);
  }

  params.set_use_presolve(glop_params.use_preprocessing());

  DCHECK_EQ(glop_params.feasibility_rule(), glop_params.optimization_rule());
  switch (glop_params.feasibility_rule()) {
    case GlopParameters::DANTZIG:
      params.set_pricing_strategy(LpParameters::PRICING_DANTZIG);
      break;
    case GlopParameters::STEEPEST_EDGE:
      params.set_pricing_strategy(LpParameters::PRICING_STEEPEST_EDGE);
      break;
    case GlopParameters::DEVEX:
      params.set_pricing_strategy(LpParameters::PRICING_DEVEX);
      break;
  }

  params.set_feasibility_tolerance(glop_params.primal_feasibility_tolerance());
  params.set_optimality_tolerance(glop_params.dual_feasibility_tolerance());

  // 0 means that Glop determines the interval automatically (and doesn't
  // support user provided values).
  params.set_refactorization_interval(0);

  params.set_objective_lower_limit(glop_params.objective_lower_limit());

  params.set_objective_upper_limit(glop_params.objective_upper_limit());

  params.set_timing_mode(LpParameters::TIMING_OFF);
  if (glop_params.max_deterministic_time() <
      std::numeric_limits<double>::infinity()) {
    params.set_timing_mode(LpParameters::TIMING_DETERMINISTIC);
    params.set_time_limit(glop_params.max_deterministic_time());
    DCHECK_EQ(glop_params.max_time_in_seconds(),
              std::numeric_limits<double>::infinity());
  }
  if (glop_params.max_time_in_seconds() <
      std::numeric_limits<double>::infinity()) {
    params.set_timing_mode(absl::GetFlag(FLAGS_time_limit_use_usertime)
                               ? LpParameters::TIMING_WALLCLOCK
                               : LpParameters::TIMING_CPU);
    params.set_time_limit(glop_params.max_time_in_seconds());
    DCHECK_EQ(glop_params.max_deterministic_time(),
              std::numeric_limits<double>::infinity());
  }

  // LpInterface assumes 0 means no limit, but Glop uses -1 to this end.
  // Hence, we shouldn't ever see 0 == max_number_of_iterations().
  DCHECK_NE(glop_params.max_number_of_iterations(), 0);
  params.set_iteration_limit(glop_params.max_number_of_iterations() == -1
                                 ? 0
                                 : glop_params.max_number_of_iterations());

  params.set_num_threads(glop_params.num_omp_threads());
  params.set_random_seed(glop_params.random_seed());
  params.set_enable_internal_solver_output(glop_params.log_search_progress());

  return params;
}

namespace {

absl::Status LpParametersAreSupportedByGlop(const LpParameters& params) {
  VLOG(10) << "calling LpParametersAreSupportedByGlop().";
  RETURN_IF_ERROR(LpParametersAreValid(params));

  if (params.lp_solver_type() != LpParameters::LP_GLOP) {
    return absl::InvalidArgumentError("Solver type must be LP_GLOP.");
  }
  if (params.scaling_strategy() == LpParameters::SCALING_LEAST_SQUARES) {
    return absl::InvalidArgumentError("Unsupported scaling strategy.");
  }
  if (params.pricing_strategy() ==
          LpParameters::PRICING_STEEPEST_EDGE_QUICK_START or
      params.pricing_strategy() == LpParameters::PRICING_PARTIAL_DANTZIG) {
    return absl::InvalidArgumentError("Unsupported pricing strategy.");
  }
  if (params.min_markowitz_threshold() != -1) {
    return absl::InvalidArgumentError(
        "Glop only supports default Markowitz threshold.");
  }
  if (params.refactorization_interval() != 0) {
    return absl::InvalidArgumentError(
        "Glop only supports automatic refactorization.");
  }
  return absl::OkStatus();
}

}  // namespace

absl::Status LpGlopInterface::SetLpParameters(const LpParameters& params) {
  VLOG(10) << "calling SetLpParameters().";
  RETURN_IF_ERROR(LpParametersAreSupportedByGlop(params));

  solve_from_scratch_ = params.solve_from_scratch();

  GlopParameters glop_params;
  switch (params.scaling_strategy()) {
    case LpParameters::SCALING_OFF:
      glop_params.set_use_scaling(false);
      break;
    case LpParameters::SCALING_DEFAULT:
    case LpParameters::SCALING_EQUILIBRATION:
      glop_params.set_use_scaling(true);
      glop_params.set_scaling_method(GlopParameters::EQUILIBRATION);
      break;
    case LpParameters::SCALING_LINEAR_PROGRAM:
      glop_params.set_use_scaling(true);
      glop_params.set_scaling_method(GlopParameters::LINEAR_PROGRAM);
      break;
    case LpParameters::SCALING_LEAST_SQUARES:
    default:
      LOG(DFATAL) << "[BUG] Unexpected scaling strategy "
                  << params.scaling_strategy();
      return absl::InvalidArgumentError("[BUG] Unexpected scaling strategy.");
  }

  glop_params.set_use_preprocessing(params.use_presolve());

  switch (params.pricing_strategy()) {
    case LpParameters::PRICING_DEFAULT:
    case LpParameters::PRICING_STEEPEST_EDGE:
      glop_params.set_feasibility_rule(
          operations_research::glop::GlopParameters::STEEPEST_EDGE);
      glop_params.set_optimization_rule(
          operations_research::glop::GlopParameters::STEEPEST_EDGE);
      break;
    case LpParameters::PRICING_DANTZIG:
      glop_params.set_feasibility_rule(
          operations_research::glop::GlopParameters::DANTZIG);
      glop_params.set_optimization_rule(
          operations_research::glop::GlopParameters::DANTZIG);
      break;
    case LpParameters::PRICING_DEVEX:
      glop_params.set_feasibility_rule(
          operations_research::glop::GlopParameters::DEVEX);
      glop_params.set_optimization_rule(
          operations_research::glop::GlopParameters::DEVEX);
      break;
    case LpParameters::PRICING_STEEPEST_EDGE_QUICK_START:
    case LpParameters::PRICING_PARTIAL_DANTZIG:
    default:
      LOG(DFATAL) << "[BUG] Unexpected pricing strategy "
                  << params.pricing_strategy();
      return absl::InvalidArgumentError("[BUG] Unexpected pricing strategy.");
  }

  glop_params.set_primal_feasibility_tolerance(
      params.feasibility_tolerance() == -1 ? 1e-8
                                           : params.feasibility_tolerance());

  glop_params.set_dual_feasibility_tolerance(
      params.optimality_tolerance() == -1 ? 1e-8
                                          : params.optimality_tolerance());

  glop_params.set_objective_lower_limit(params.objective_lower_limit());
  glop_params.set_objective_upper_limit(params.objective_upper_limit());

  switch (params.timing_mode()) {
    case LpParameters::TIMING_OFF:
      glop_params.set_max_time_in_seconds(
          std::numeric_limits<double>::infinity());
      glop_params.set_max_deterministic_time(
          std::numeric_limits<double>::infinity());
      break;
    case LpParameters::TIMING_DETERMINISTIC:
      glop_params.set_max_time_in_seconds(
          std::numeric_limits<double>::infinity());
      glop_params.set_max_deterministic_time(params.time_limit());
      break;
    case LpParameters::TIMING_WALLCLOCK:
      glop_params.set_max_time_in_seconds(params.time_limit());
      glop_params.set_max_deterministic_time(
          std::numeric_limits<double>::infinity());
      absl::SetFlag(&FLAGS_time_limit_use_usertime, true);
      break;
    case LpParameters::TIMING_CPU:
      glop_params.set_max_time_in_seconds(params.time_limit());
      glop_params.set_max_deterministic_time(
          std::numeric_limits<double>::infinity());
      absl::SetFlag(&FLAGS_time_limit_use_usertime, false);
      break;
    default:
      LOG(DFATAL) << "[BUG] Unexpected timing mode " << params.timing_mode();
      return absl::InvalidArgumentError("[BUG] Unexpected timing mode.");
  }

  glop_params.set_max_number_of_iterations(
      params.iteration_limit() == 0 ? -1 : params.iteration_limit());

  glop_params.set_num_omp_threads(
      params.num_threads() == 0 ? 1 : params.num_threads());

  glop_params.set_random_seed(params.random_seed());

  glop_params.set_log_search_progress(params.enable_internal_solver_output());

  solver_.SetParameters(glop_params);
  return absl::OkStatus();
}

// ==========================================================================
// Numerical methods.
// ==========================================================================

double LpGlopInterface::Infinity() const {
  VLOG(10) << "calling Infinity().";
  return std::numeric_limits<double>::infinity();
}

bool LpGlopInterface::IsInfinity(double value) const {
  VLOG(10) << "calling IsInfinity().";
  return value == Infinity();
}

// ==========================================================================
// File interface methods.
// ==========================================================================

absl::Status LpGlopInterface::ReadLpFromFile(const std::string& file_path) {
  VLOG(10) << "calling ReadLpFromFile().";
  MPModelProto proto;
  if (!ReadFileToProto(file_path, &proto)) {
    return absl::Status(absl::StatusCode::kInternal,
                        absl::StrFormat("Could not read: %s", file_path));
  }
  lp_.Clear();
  MPModelProtoToLinearProgram(proto, &lp_);
  return absl::OkStatus();
}

absl::StatusOr<std::string> LpGlopInterface::WriteLpToFile(
    const std::string& file_path) const {
  VLOG(10) << "calling WriteLpToFile().";
  MPModelProto proto;
  LinearProgramToMPModelProto(lp_, &proto);
  if (!WriteProtoToFile(file_path, proto,
                        operations_research::ProtoWriteFormat::kProtoText,
                        /*gzipped=*/false)) {
    return absl::Status(absl::StatusCode::kInternal,
                        absl::StrFormat("Could not write: %s", file_path));
  }
  return absl::StrCat(file_path);
}

}  // namespace minimip
