// Copyright 2023 the MiniMIP Project
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

LPBasisStatus ConvertGlopVariableStatus(VariableStatus status,
                                        double reduced_cost) {
  switch (status) {
    case VariableStatus::BASIC:
      return LPBasisStatus::kBasic;
    case VariableStatus::AT_UPPER_BOUND:
      return LPBasisStatus::kAtUpperBound;
    case VariableStatus::AT_LOWER_BOUND:
      return LPBasisStatus::kAtLowerBound;
    case VariableStatus::FREE:
      return LPBasisStatus::kFree;
    case VariableStatus::FIXED_VALUE:
      return reduced_cost > 0.0 ? LPBasisStatus::kAtLowerBound
                                : LPBasisStatus::kAtUpperBound;
    default:
      LOG(FATAL) << "Unknown Glop column basis status.";
  }
}

LPBasisStatus ConvertGlopConstraintStatus(ConstraintStatus status,
                                          double dual_value) {
  switch (status) {
    case ConstraintStatus::BASIC:
      return LPBasisStatus::kBasic;
    case ConstraintStatus::AT_UPPER_BOUND:
      return LPBasisStatus::kAtUpperBound;
    case ConstraintStatus::AT_LOWER_BOUND:
      return LPBasisStatus::kAtLowerBound;
    case ConstraintStatus::FREE:
      return LPBasisStatus::kFree;
    case ConstraintStatus::FIXED_VALUE:
      return dual_value > 0.0 ? LPBasisStatus::kAtLowerBound
                              : LPBasisStatus::kAtUpperBound;
    default:
      LOG(FATAL) << "Unknown Glop row basis status.";
  }
}

VariableStatus ConvertMiniMIPVariableStatus(LPBasisStatus status) {
  switch (status) {
    case LPBasisStatus::kBasic:
      return VariableStatus::BASIC;
    case LPBasisStatus::kAtUpperBound:
      return VariableStatus::AT_UPPER_BOUND;
    case LPBasisStatus::kAtLowerBound:
      return VariableStatus::AT_LOWER_BOUND;
    case LPBasisStatus::kFixed:
      return VariableStatus::FIXED_VALUE;
    case LPBasisStatus::kFree:
      return VariableStatus::FREE;
    default:
      LOG(FATAL) << "Unknown MiniMip col basis status.";
  }
}

VariableStatus ConvertMiniMIPConstraintStatusToSlackStatus(
    LPBasisStatus status) {
  // We swap lower and upper bound, because Glop adds slacks with -1.0
  // coefficient, whereas LP interface assumes slacks with +1.0 coefficient.
  switch (status) {
    case LPBasisStatus::kBasic:
      return VariableStatus::BASIC;
    case LPBasisStatus::kAtUpperBound:
      return VariableStatus::AT_LOWER_BOUND;
    case LPBasisStatus::kAtLowerBound:
      return VariableStatus::AT_UPPER_BOUND;
    case LPBasisStatus::kFixed:
      return VariableStatus::FIXED_VALUE;
    case LPBasisStatus::kFree:
      return VariableStatus::FREE;
    default:
      LOG(FATAL) << "Unknown MiniMip row basis status.";
  }
}

bool IsDualBoundValid(ProblemStatus status) {
  return status == ProblemStatus::OPTIMAL ||
         status == ProblemStatus::DUAL_FEASIBLE ||
         status == ProblemStatus::DUAL_UNBOUNDED;
}

}  // namespace

LPGlopInterface::LPGlopInterface()
    : lp_modified_since_last_solve_(true),
      lp_time_limit_was_reached_(false),
      num_iterations_of_last_solve_(0),
      lp_info_(false),
      pricing_(LPPricing::kDefault),
      from_scratch_(false),
      num_threads_(0),
      timing_(0) {
  tmp_row_ = std::make_unique<ScatteredRow>();
  tmp_column_ = std::make_unique<ScatteredColumn>();
}

// ==========================================================================
// LP model setters.
// ==========================================================================

absl::Status LPGlopInterface::PopulateFromMipData(const MipData& mip_data) {
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
  RETURN_IF_ERROR(SetObjectiveSense(mip_data.is_maximization()));
  return absl::OkStatus();
}

absl::Status LPGlopInterface::AddColumn(const SparseCol& col_data,
                                        double lower_bound, double upper_bound,
                                        double objective_coefficient,
                                        const std::string& name) {
  DCHECK(!col_data.MayNeedCleaning());
  DCHECK(std::all_of(col_data.entries().begin(), col_data.entries().end(),
                     [num_rows = GetNumberOfRows()](const ColEntry& e) {
                       return RowIndex(0) <= e.index && e.index < num_rows;
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

absl::Status LPGlopInterface::AddColumns(
    const StrongSparseMatrix& matrix,
    const absl::StrongVector<ColIndex, double>& lower_bounds,
    const absl::StrongVector<ColIndex, double>& upper_bounds,
    const absl::StrongVector<ColIndex, double>& objective_coefficients,
    const absl::StrongVector<ColIndex, std::string>& names) {
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

absl::Status LPGlopInterface::DeleteColumns(ColIndex first_col,
                                            ColIndex last_col) {
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

absl::Status LPGlopInterface::AddRow(const SparseRow& row_data,
                                     double left_hand_side,
                                     double right_hand_side,
                                     const std::string& name) {
  DCHECK(!row_data.MayNeedCleaning());
  DCHECK(std::all_of(row_data.entries().begin(), row_data.entries().end(),
                     [num_cols = GetNumberOfColumns()](const RowEntry& e) {
                       return ColIndex(0) <= e.index && e.index < num_cols;
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

absl::Status LPGlopInterface::AddRows(
    const absl::StrongVector<RowIndex, SparseRow>& rows,
    const absl::StrongVector<RowIndex, double>& left_hand_sides,
    const absl::StrongVector<RowIndex, double>& right_hand_sides,
    const absl::StrongVector<RowIndex, std::string>& names) {
  DCHECK_EQ(names.size(), left_hand_sides.size());
  DCHECK_EQ(left_hand_sides.size(), right_hand_sides.size());
  DCHECK_EQ(right_hand_sides.size(), rows.size());
  for (RowIndex row(0); row < rows.size(); ++row) {
    RETURN_IF_ERROR(AddRow(rows[row], left_hand_sides[row],
                           right_hand_sides[row], names[row]));
  }
  return absl::OkStatus();
}

absl::Status LPGlopInterface::DeleteRows(RowIndex first_row,
                                         RowIndex last_row) {
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
LPGlopInterface::DeleteRowSet(
    const absl::StrongVector<RowIndex, bool>& rows_to_delete) {
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

void LPGlopInterface::DeleteRowsAndUpdateCurrentBasis(
    const DenseBooleanColumn& rows_to_delete) {
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

absl::Status LPGlopInterface::Clear() {
  VLOG(3) << "Clearing Glop.";
  lp_.Clear();
  lp_modified_since_last_solve_ = true;
  return absl::OkStatus();
}

absl::Status LPGlopInterface::ClearState() {
  VLOG(3) << "Clear Glop state.";
  solver_.ClearStateForNextSolve();
  return absl::OkStatus();
}

absl::Status LPGlopInterface::SetColumnBounds(ColIndex col, double lower_bound,
                                              double upper_bound) {
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

absl::Status LPGlopInterface::SetRowSides(RowIndex row, double left_hand_side,
                                          double right_hand_side) {
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

absl::Status LPGlopInterface::SetObjectiveSense(bool is_maximization) {
  VLOG(3) << "Setting maximize=" << is_maximization << ".";
  lp_.SetMaximizationProblem(is_maximization);
  lp_modified_since_last_solve_ = true;
  return absl::OkStatus();
}

absl::Status LPGlopInterface::SetObjectiveCoefficient(
    ColIndex col, double objective_coefficient) {
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

RowIndex LPGlopInterface::GetNumberOfRows() const {
  return RowIndex(lp_.num_constraints().value());
}

ColIndex LPGlopInterface::GetNumberOfColumns() const {
  return ColIndex(lp_.num_variables().value());
}

bool LPGlopInterface::IsMaximization() const {
  return lp_.IsMaximizationProblem();
}

int64_t LPGlopInterface::GetNumberOfNonZeros() const {
  return lp_.num_entries().value();
}

SparseCol LPGlopInterface::GetSparseColumnCoefficients(ColIndex col) const {
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

SparseRow LPGlopInterface::GetSparseRowCoefficients(RowIndex row) const {
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

double LPGlopInterface::GetObjectiveCoefficient(ColIndex col) const {
  CHECK_GE(col, 0);
  CHECK_LT(col, GetNumberOfColumns());
  return lp_.objective_coefficients()[GlopColIndex(col.value())];
}

double LPGlopInterface::GetLowerBound(ColIndex col) const {
  CHECK_GE(col, 0);
  CHECK_LT(col, GetNumberOfColumns());
  return lp_.variable_lower_bounds()[GlopColIndex(col.value())];
}

double LPGlopInterface::GetUpperBound(ColIndex col) const {
  CHECK_GE(col, 0);
  CHECK_LT(col, GetNumberOfColumns());
  return lp_.variable_upper_bounds()[GlopColIndex(col.value())];
}

double LPGlopInterface::GetLeftHandSide(RowIndex row) const {
  CHECK_GE(row, 0);
  CHECK_LT(row, GetNumberOfRows());
  return lp_.constraint_lower_bounds()[GlopRowIndex(row.value())];
}

double LPGlopInterface::GetRightHandSide(RowIndex row) const {
  CHECK_GE(row, 0);
  CHECK_LT(row, GetNumberOfRows());
  return lp_.constraint_upper_bounds()[GlopRowIndex(row.value())];
}

double LPGlopInterface::GetMatrixCoefficient(ColIndex col, RowIndex row) const {
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
absl::Status LPGlopInterface::SolveInternal(bool recursive,
                                            TimeLimit* time_limit) {
  // Recompute `scaled_lp_`.
  if (lp_modified_since_last_solve_) {
    // TODO(lpawel): Avoid doing a copy if there is no scaling.
    scaled_lp_.PopulateFromLinearProgram(lp_);
    scaled_lp_.AddSlackVariablesWhereNecessary(false);

    if (parameters_.use_scaling()) {
      // TODO(lpawel): Avoid rescaling if nothing changed.
      scaler_.Scale(&scaled_lp_);
    } else {
      scaler_.Clear();
    }
  }

  solver_.SetParameters(parameters_);

  if (!recursive) {
    num_iterations_of_last_solve_ = 0;
    lp_time_limit_was_reached_ = false;
  }

  if (from_scratch_) {
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
  if ((status == ProblemStatus::PRIMAL_FEASIBLE ||
       status == ProblemStatus::OPTIMAL) &&
      parameters_.use_scaling()) {
    const auto primal_values = GetPrimalValues();
    CHECK_OK(primal_values);
    const operations_research::glop::DenseRow solution(primal_values->begin(),
                                                       primal_values->end());
    if (!lp_.SolutionIsLPFeasible(solution,
                                  parameters_.primal_feasibility_tolerance())) {
      VLOG(1) << "Solution not feasible w.r.t. absolute tolerance "
              << parameters_.primal_feasibility_tolerance()
              << ". Will re-optimize.";
      parameters_.set_use_scaling(false);

      // This is needed to force setting `scaled_lp_ = lp_` in the recursive
      // call.
      // TODO(lpawel): Fix this logic and make it more explicit.
      lp_modified_since_last_solve_ = true;

      RETURN_IF_ERROR(SolveInternal(/*recursive=*/true, time_limit));
      parameters_.set_use_scaling(true);
    }
  }

  lp_modified_since_last_solve_ = false;
  return absl::OkStatus();
}

// ==========================================================================
// Solving methods.
// ==========================================================================

absl::Status LPGlopInterface::SolveLPWithPrimalSimplex() {
  VLOG(3) << "Solving with primal simplex: "
          << "num_cols=" << lp_.num_variables().value()
          << ", num_rows=" << lp_.num_constraints().value();
  std::unique_ptr<TimeLimit> time_limit =
      TimeLimit::FromParameters(parameters_);
  parameters_.set_use_dual_simplex(false);
  return SolveInternal(/*recursive=*/false, time_limit.get());
}

absl::Status LPGlopInterface::SolveLPWithDualSimplex() {
  VLOG(3) << "Solving with dual simplex: "
          << "num_cols=" << lp_.num_variables().value()
          << ", num_rows=" << lp_.num_constraints().value();
  std::unique_ptr<TimeLimit> time_limit =
      TimeLimit::FromParameters(parameters_);
  parameters_.set_use_dual_simplex(true);
  return SolveInternal(/*recursive=*/false, time_limit.get());
}

absl::Status LPGlopInterface::StartStrongBranching() {
  // TODO(lpawel): Save Glop state and tune Glop towards strong branching (and
  // avoid rescaling completely when in strong branching).
  return absl::OkStatus();
}

absl::Status LPGlopInterface::EndStrongBranching() {
  // TODO(lpawel): Restore the saved Glop state.
  return absl::OkStatus();
}

absl::StatusOr<LPGlopInterface::StrongBranchResult>
LPGlopInterface::SolveDownAndUpStrongBranch(ColIndex col, double primal_value,
                                            int iteration_limit) {
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
  parameters_.set_use_dual_simplex(true);
  solver_.SetParameters(parameters_);
  const Fractional eps = parameters_.primal_feasibility_tolerance();

  std::unique_ptr<TimeLimit> time_limit =
      TimeLimit::FromParameters(parameters_);

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
    if (lp_.IsMaximizationProblem()) {
      result.dual_bound_down_branch = parameters_.objective_lower_limit();
    } else {
      result.dual_bound_down_branch = parameters_.objective_upper_limit();
    }
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
    if (lp_.IsMaximizationProblem()) {
      result.dual_bound_up_branch = parameters_.objective_lower_limit();
    } else {
      result.dual_bound_up_branch = parameters_.objective_upper_limit();
    }
  }

  //  Restore the bounds.
  scaled_lp_.SetVariableBounds(glop_col, lower_bound, upper_bound);

  return result;
}

// ==========================================================================
// Solution information getters.
// ==========================================================================

bool LPGlopInterface::IsSolved() const {
  // TODO(lpawel): Track this to avoid unneeded resolving.
  return (!lp_modified_since_last_solve_);
}

bool LPGlopInterface::ExistsPrimalRay() const {
  return solver_.GetProblemStatus() == ProblemStatus::PRIMAL_UNBOUNDED;
}

bool LPGlopInterface::HasPrimalRay() const {
  return solver_.GetProblemStatus() == ProblemStatus::PRIMAL_UNBOUNDED;
}

bool LPGlopInterface::IsPrimalUnbounded() const {
  return solver_.GetProblemStatus() == ProblemStatus::PRIMAL_UNBOUNDED;
}

bool LPGlopInterface::IsPrimalInfeasible() const {
  const ProblemStatus status = solver_.GetProblemStatus();
  return status == ProblemStatus::DUAL_UNBOUNDED ||
         status == ProblemStatus::PRIMAL_INFEASIBLE;
}

bool LPGlopInterface::IsPrimalFeasible() const {
  const ProblemStatus status = solver_.GetProblemStatus();
  return status == ProblemStatus::PRIMAL_FEASIBLE ||
         status == ProblemStatus::OPTIMAL;
}

bool LPGlopInterface::ExistsDualRay() const {
  return solver_.GetProblemStatus() == ProblemStatus::DUAL_UNBOUNDED;
}

bool LPGlopInterface::HasDualRay() const {
  return solver_.GetProblemStatus() == ProblemStatus::DUAL_UNBOUNDED;
}

bool LPGlopInterface::IsDualUnbounded() const {
  return solver_.GetProblemStatus() == ProblemStatus::DUAL_UNBOUNDED;
}

bool LPGlopInterface::IsDualInfeasible() const {
  const ProblemStatus status = solver_.GetProblemStatus();
  return status == ProblemStatus::PRIMAL_UNBOUNDED ||
         status == ProblemStatus::DUAL_INFEASIBLE;
}

bool LPGlopInterface::IsDualFeasible() const {
  const ProblemStatus status = solver_.GetProblemStatus();
  return status == ProblemStatus::DUAL_FEASIBLE ||
         status == ProblemStatus::OPTIMAL;
}

bool LPGlopInterface::IsOptimal() const {
  return solver_.GetProblemStatus() == ProblemStatus::OPTIMAL;
}

bool LPGlopInterface::IsStable() const {
  // For correctness, we need to report "unstable" if Glop was not able to
  // prove optimality because of numerical issues. Currently, Glop still
  // reports primal/dual feasible if at the end, one status is within the
  // tolerance but not the other.
  const ProblemStatus status = solver_.GetProblemStatus();
  if ((status == ProblemStatus::PRIMAL_FEASIBLE ||
       status == ProblemStatus::DUAL_FEASIBLE) &&
      !ObjectiveLimitIsExceeded() && !IterationLimitIsExceeded() &&
      !TimeLimitIsExceeded()) {
    VLOG(3) << "OPTIMAL not reached and no limit: unstable";
    return false;
  }

  if (status == ProblemStatus::ABNORMAL ||
      status == ProblemStatus::INVALID_PROBLEM ||
      status == ProblemStatus::IMPRECISE) {
    VLOG(3) << "Errors while solving: unstable";
    return false;
  }
  return true;
}

bool LPGlopInterface::ObjectiveLimitIsExceeded() const {
  return solver_.objective_limit_reached();
}

bool LPGlopInterface::TimeLimitIsExceeded() const {
  return lp_time_limit_was_reached_;
}

bool LPGlopInterface::IterationLimitIsExceeded() const {
  // We might have accumulated iterations across 2 recursive solves,
  // hence _GE, and not _EQ.
  DCHECK_GE(num_iterations_of_last_solve_, solver_.GetNumberOfIterations());
  // The limit of -1 means there is no limit.
  return parameters_.max_number_of_iterations() != -1 &&
         num_iterations_of_last_solve_ >=
             parameters_.max_number_of_iterations();
}

int64_t LPGlopInterface::GetNumIterations() const {
  return num_iterations_of_last_solve_;
}

double LPGlopInterface::GetObjectiveValue() {
  DCHECK(IsOptimal());
  return solver_.GetObjectiveValue();
}

absl::StatusOr<absl::StrongVector<ColIndex, double>>
LPGlopInterface::GetPrimalValues() const {
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
LPGlopInterface::GetDualValues() const {
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
LPGlopInterface::GetReducedCosts() const {
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
LPGlopInterface::GetRowActivities() const {
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
LPGlopInterface::GetPrimalRay() const {
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
LPGlopInterface::GetDualRay() const {
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

absl::StatusOr<absl::StrongVector<ColIndex, LPBasisStatus>>
LPGlopInterface::GetBasisStatusForColumns() const {
  DCHECK(IsOptimal());
  absl::StrongVector<ColIndex, LPBasisStatus> statuses;
  statuses.reserve(lp_.num_variables().value());
  for (GlopColIndex col(0); col < lp_.num_variables(); ++col) {
    statuses.push_back(ConvertGlopVariableStatus(solver_.GetVariableStatus(col),
                                                 solver_.GetReducedCost(col)));
  }
  return statuses;
}

absl::StatusOr<absl::StrongVector<RowIndex, LPBasisStatus>>
LPGlopInterface::GetBasisStatusForRows() const {
  DCHECK(IsOptimal());
  absl::StrongVector<RowIndex, LPBasisStatus> statuses;
  statuses.reserve(lp_.num_constraints().value());
  for (GlopRowIndex row(0); row < lp_.num_constraints(); ++row) {
    statuses.push_back(ConvertGlopConstraintStatus(
        solver_.GetConstraintStatus(row), solver_.GetDualValue(row)));
  }
  return statuses;
}

absl::Status LPGlopInterface::SetBasisStatusForColumnsAndRows(
    const absl::StrongVector<ColIndex, LPBasisStatus>& column_basis_statuses,
    const absl::StrongVector<RowIndex, LPBasisStatus>& row_basis_statuses) {
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

std::vector<ColOrRowIndex> LPGlopInterface::GetColumnsAndRowsInBasis() const {
  std::vector<ColOrRowIndex> basis;
  basis.reserve(GetNumberOfRows().value());
  // The order in which we populate the `basis` is important!
  for (GlopRowIndex row(0); row < lp_.num_constraints(); ++row) {
    const GlopColIndex col = solver_.GetBasis(row);
    basis.push_back(col < lp_.num_variables()
                        ? ColOrRowIndex(ColIndex(col.value()))
                        : ColOrRowIndex(RowIndex(col.value())));
  }
  return basis;
}

// ==========================================================================
// Getters of vectors in the inverted basis matrix.
// ==========================================================================

absl::StatusOr<SparseRow> LPGlopInterface::GetSparseRowOfBInverted(
    RowIndex row_in_basis) const {
  SparseRow sparse_row;

  solver_.GetBasisFactorization().LeftSolveForUnitRow(
      GlopColIndex(row_in_basis.value()), tmp_row_.get());
  scaler_.UnscaleUnitRowLeftSolve(
      solver_.GetBasis(GlopRowIndex(row_in_basis.value())), tmp_row_.get());

  // DCHECK_EQ(tmp_row_->values.size(), lp_.num_constraints().value());

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
    const Fractional eps = parameters_.primal_feasibility_tolerance();
    for (GlopColIndex col(0); col < RowToColIndex(lp_.num_constraints());
         ++col) {
      const double value = (*tmp_row_)[col];
      if (std::abs(value) >= eps) {
        sparse_row.AddEntry(ColIndex(col.value()), value);
      }
    }
  }
  return sparse_row;
}

absl::StatusOr<SparseCol> LPGlopInterface::GetSparseColumnOfBInverted(
    ColIndex col_in_basis) const {
  SparseCol sparse_column;
  // We need to loop through the rows to extract the values for `col_in_basis`.
  for (GlopRowIndex row(0); row < lp_.num_constraints(); ++row) {
    solver_.GetBasisFactorization().LeftSolveForUnitRow(RowToColIndex(row),
                                                        tmp_row_.get());
    scaler_.UnscaleUnitRowLeftSolve(solver_.GetBasis(row), tmp_row_.get());

    double value = (*tmp_row_)[GlopColIndex(col_in_basis.value())];
    if (std::abs(value) >= parameters_.primal_feasibility_tolerance()) {
      sparse_column.AddEntry(RowIndex(row.value()), value);
    }
  }
  return sparse_column;
}

absl::StatusOr<SparseRow> LPGlopInterface::GetSparseRowOfBInvertedTimesA(
    RowIndex row_in_basis) const {
  SparseRow sparse_row;
  solver_.GetBasisFactorization().LeftSolveForUnitRow(
      GlopColIndex(row_in_basis.value()), tmp_row_.get());
  scaler_.UnscaleUnitRowLeftSolve(
      solver_.GetBasis(GlopRowIndex(row_in_basis.value())), tmp_row_.get());
  for (GlopColIndex col(0); col < lp_.num_variables(); ++col) {
    double value = operations_research::glop::ScalarProduct(
        tmp_row_->values, lp_.GetSparseColumn(col));
    if (std::abs(value) >= parameters_.primal_feasibility_tolerance()) {
      sparse_row.AddEntry(ColIndex(col.value()), value);
    }
  }
  return sparse_row;
}

absl::StatusOr<SparseCol> LPGlopInterface::GetSparseColumnOfBInvertedTimesA(
    ColIndex col_in_basis) const {
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
      if (std::abs(value) > parameters_.primal_feasibility_tolerance()) {
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

absl::StatusOr<int> LPGlopInterface::GetIntegerParameter(
    LPParameter type) const {
  int param_val;
  switch (type) {
    case LPParameter::kFromScratch:
      param_val = static_cast<int>(from_scratch_);
      break;
    case LPParameter::kLPInfo:
      param_val = static_cast<int>(lp_info_);
      break;
    case LPParameter::kLPIterationLimit:
      param_val = static_cast<int>(parameters_.max_number_of_iterations());
      break;
    case LPParameter::kPresolving:
      param_val = static_cast<int>(parameters_.use_preprocessing());
      break;
    case LPParameter::kPolishing:
      param_val = 0;
      break;
    case LPParameter::kPricing:
      param_val = static_cast<int>(pricing_);
      break;
    case LPParameter::kRefactor:
      param_val = 0;
      break;
    case LPParameter::kScaling:
      param_val = static_cast<int>(parameters_.use_scaling());
      break;
    case LPParameter::kThreads:
      param_val = num_threads_;
      break;
    case LPParameter::kTiming:
      param_val = timing_;
      break;
    case LPParameter::kRandomSeed:
      param_val = static_cast<int>(parameters_.random_seed());
      break;
    default:
      return util::InvalidArgumentErrorBuilder()
             << "Unknown integer parameter type " << type;
  }

  return param_val;
}

absl::Status LPGlopInterface::SetIntegerParameter(LPParameter type,
                                                  int param_val) {
  switch (type) {
    case LPParameter::kFromScratch:
      from_scratch_ = static_cast<bool>(param_val);
      break;
    case LPParameter::kLPInfo:
      if (param_val == 0) {
        static_cast<void>(google::SetVLOGLevel("*", google::GLOG_INFO));
        lp_info_ = false;
      } else {
        static_cast<void>(google::SetVLOGLevel("*", google::GLOG_ERROR));
        lp_info_ = true;
      }
      break;
    case LPParameter::kLPIterationLimit:
      parameters_.set_max_number_of_iterations(param_val);
      break;
    case LPParameter::kPresolving:
      parameters_.set_use_preprocessing(static_cast<bool>(param_val));
      break;
    case LPParameter::kPolishing:
      if (param_val != 0) {
        return absl::InvalidArgumentError(
            "Polishing is not supported by glop.");
      }
      break;
    case LPParameter::kPricing:
      pricing_ = static_cast<LPPricing>(param_val);
      switch (pricing_) {
        case LPPricing::kDefault:
        case LPPricing::kAuto:
        case LPPricing::kPartial:
        case LPPricing::kSteep:
        case LPPricing::kSteepQStart:
          parameters_.set_feasibility_rule(
              operations_research::glop::
                  GlopParameters_PricingRule_STEEPEST_EDGE);
          break;
        case LPPricing::kFull:
          // Dantzig does not really fit, but use it anyway
          parameters_.set_feasibility_rule(
              operations_research::glop::GlopParameters_PricingRule_DANTZIG);
          break;
        case LPPricing::kDevex:
          parameters_.set_feasibility_rule(
              operations_research::glop::GlopParameters_PricingRule_DEVEX);
          break;
        default:
          return util::InvalidArgumentErrorBuilder()
                 << "Unknown pricing strategy " << pricing_;
      }
      break;
    case LPParameter::kRefactor:
      if (param_val != 0) {
        return absl::InvalidArgumentError(
            "Glop only supports automatic refactoring.");
      }
      break;
    case LPParameter::kScaling:
      parameters_.set_use_scaling(static_cast<bool>(param_val));
      break;
    case LPParameter::kThreads:
      num_threads_ = param_val;
      parameters_.set_num_omp_threads(num_threads_ == 0 ? 1 : num_threads_);
      break;
    case LPParameter::kTiming:
      assert(param_val <= 2);
      timing_ = param_val;
      absl::SetFlag(&FLAGS_time_limit_use_usertime, timing_ == 1);
      break;
    case LPParameter::kRandomSeed:
      parameters_.set_random_seed(param_val);
      break;
    default:
      return util::InvalidArgumentErrorBuilder()
             << "Unknown integer parameter type " << type << " with value "
             << param_val;
  }

  return absl::OkStatus();
}

absl::StatusOr<double> LPGlopInterface::GetRealParameter(
    LPParameter type) const {
  double param_val;

  switch (type) {
    case LPParameter::kFeasibilityTolerance:
      param_val = parameters_.primal_feasibility_tolerance();
      break;
    case LPParameter::kDualFeasibilityTolerance:
      param_val = parameters_.dual_feasibility_tolerance();
      break;
    case LPParameter::kObjectiveLimit:
      param_val = lp_.IsMaximizationProblem()
                      ? parameters_.objective_lower_limit()
                      : parameters_.objective_upper_limit();
      break;
    case LPParameter::kLPTimeLimit:
      param_val = absl::GetFlag(FLAGS_time_limit_use_usertime)
                      ? parameters_.max_time_in_seconds()
                      : parameters_.max_deterministic_time();
      break;
    default:
      return util::InvalidArgumentErrorBuilder()
             << "Unknown real parameter type " << type;
  }
  return param_val;
}

absl::Status LPGlopInterface::SetRealParameter(LPParameter type,
                                               double param_val) {
  switch (type) {
    case LPParameter::kFeasibilityTolerance:
      parameters_.set_primal_feasibility_tolerance(param_val);
      break;
    case LPParameter::kDualFeasibilityTolerance:
      parameters_.set_dual_feasibility_tolerance(param_val);
      break;
    case LPParameter::kObjectiveLimit:
      if (lp_.IsMaximizationProblem()) {
        parameters_.set_objective_lower_limit(param_val);
      } else {
        parameters_.set_objective_upper_limit(param_val);
      }
      break;
    case LPParameter::kLPTimeLimit:
      if (absl::GetFlag(FLAGS_time_limit_use_usertime)) {
        parameters_.set_max_time_in_seconds(param_val);
      } else {
        parameters_.set_max_deterministic_time(param_val);
      }
      break;
    default:
      return util::InvalidArgumentErrorBuilder()
             << "Unknown real parameter type " << type << " with value "
             << param_val;
  }

  return absl::OkStatus();
}

// ==========================================================================
// Numerical methods.
// ==========================================================================

double LPGlopInterface::Infinity() const {
  return std::numeric_limits<double>::infinity();
}

bool LPGlopInterface::IsInfinity(double value) const {
  return value == Infinity();
}

// ==========================================================================
// File interface methods.
// ==========================================================================

absl::Status LPGlopInterface::ReadLPFromFile(const std::string& file_path) {
  MPModelProto proto;
  if (!ReadFileToProto(file_path, &proto)) {
    return absl::Status(absl::StatusCode::kInternal,
                        absl::StrFormat("Could not read: %s", file_path));
  }
  lp_.Clear();
  MPModelProtoToLinearProgram(proto, &lp_);
  return absl::OkStatus();
}

absl::StatusOr<std::string> LPGlopInterface::WriteLPToFile(
    const std::string& file_path) const {
  MPModelProto proto;
  LinearProgramToMPModelProto(lp_, &proto);
  if (!WriteProtoToFile(file_path, proto,
                        operations_research::ProtoWriteFormat::kProtoText,
                        /*gzipped=*/true)) {
    return absl::Status(absl::StatusCode::kInternal,
                        absl::StrFormat("Could not write: %s", file_path));
  }
  return absl::StrCat(file_path, ".gz");
}

}  // namespace minimip
