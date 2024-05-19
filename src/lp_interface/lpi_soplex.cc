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

#include "src/lp_interface/lpi_soplex.h"

#include <soplex/spxout.h>

namespace minimip {

namespace {

// Macro for a single SoPlex call for which exceptions have to be catched -
// return an LP error. We make no distinction between different exception types,
// e.g., between memory allocation and other exceptions.
#define SOPLEX_TRY(x)                                                         \
  {                                                                           \
    try {                                                                     \
      (x);                                                                    \
    } catch (const soplex::SPxMemoryException& E) {                           \
      LOG(WARNING) << "SoPlex threw a memory exception: " << E.what() << "."; \
      return absl::InternalError(E.what());                                   \
    } catch (const soplex::SPxException& E) {                                 \
      LOG(WARNING) << "SoPlex threw an exception: " << E.what() << ".";       \
      return absl::InternalError(E.what());                                   \
    }                                                                         \
  }

bool FileExists(const std::string& file_path) {
  VLOG(10) << "calling FileExists().";
  FILE* file = fopen(file_path.c_str(), "r");
  if (file == nullptr) return false;

  fclose(file);

  return true;
}
}  // namespace

LpSoplexInterface::LpSoplexInterface()
    : spx_(new soplex::SoPlex),
      solve_from_scratch_(false),
      is_solved_(false),
      col_basis_status_(0),
      row_basis_status_(0) {
  VLOG(10) << "calling LpSoplexInterface().";
  // We set the parameters explicitly to the default values, because
  // `SetLpParameters` decides some default values (e.g., default tolerances).
  // This way all parameter values are consistent from the start and do not
  // depend on SoPlex implementation.
  LpParameters params;
  params.set_lp_solver_type(LpParameters::LP_SOPLEX);
  CHECK_OK(SetLpParameters(params));
}

LpSoplexInterface::~LpSoplexInterface() = default;

soplex::DataArray<soplex::SPxSolver::VarStatus>&
LpSoplexInterface::RowsBasisStatus() {
  VLOG(10) << "calling RowsBasisStatus().";
  return row_basis_status_;
}

soplex::DataArray<soplex::SPxSolver::VarStatus>&
LpSoplexInterface::ColumnsBasisStatus() {
  VLOG(10) << "calling ColumnsBasisStatus().";
  return col_basis_status_;
}

bool LpSoplexInterface::CheckConsistentBounds() const {
  VLOG(10) << "calling CheckConsistentBounds().";
  for (int i = 0; i < spx_->numColsReal(); ++i) {
    if (spx_->lowerReal(i) >
        spx_->upperReal(i) + spx_->realParam(soplex::SoPlex::EPSILON_ZERO)) {
      VLOG(3) << "inconsistent bounds on column " << i
              << ": lower=" << spx_->lowerReal(i)
              << ", upper=" << spx_->upperReal(i) << ".";
      return false;
    }
  }
  return true;
}

bool LpSoplexInterface::CheckConsistentSides() const {
  VLOG(10) << "calling CheckConsistentSides().";
  for (int i = 0; i < spx_->numRowsReal(); ++i) {
    if (spx_->lhsReal(i) >
        spx_->rhsReal(i) + spx_->realParam(soplex::SoPlex::EPSILON_ZERO)) {
      VLOG(3) << "inconsistent sides on row " << i
              << ": lhs=" << spx_->lhsReal(i) << ", rhs=" << spx_->rhsReal(i)
              << ".";
      return false;
    }
  }
  return true;
}

soplex::SPxSolver::Status LpSoplexInterface::LpSolve(
    bool print_warning = true) {
  VLOG(10) << "calling LpSolve().";
  CHECK(CheckConsistentBounds());
  CHECK(CheckConsistentSides());

  try {
    static_cast<void>(spx_->optimize());
  } catch (const soplex::SPxException& x) {
    if (print_warning) {
      LOG(WARNING) << "SoPlex threw an exception: " << x.what() << ".";
    }
    // Since it is not clear if the status in SoPlex are set correctly we want
    // to make sure that if an error is thrown the status is not OPTIMAL
    // anymore.
    CHECK_NE(spx_->status(), soplex::SPxSolver::OPTIMAL);
  }
  // NOLINTBEGIN(readability-simplify-boolean-expr)
  CHECK(spx_->intParam(soplex::SoPlex::ITERLIMIT) < 0 or
        spx_->numIterations() <= spx_->intParam(soplex::SoPlex::ITERLIMIT));
  // NOLINTEND(readability-simplify-boolean-expr)

  // Update the time limit.
  const soplex::Real time_spent = spx_->solveTime();
  const soplex::Real time_limit = spx_->realParam(soplex::SoPlex::TIMELIMIT);

  // Get the current time limit.
  if (time_spent > 0.0) {
    spx_->setRealParam(soplex::SoPlex::TIMELIMIT,
                       std::max(0.0, time_limit - time_spent));
  }
  DCHECK_GE(time_limit, 0);

  // This needs to be a CHECK to trigger the function call and set the new time
  // limit
  CHECK(spx_->setRealParam(soplex::SoPlex::TIMELIMIT, time_limit))
      << "SoPlex: unsupported parameter value.\n";

  const soplex::SPxSolver::Status spx_status = spx_->status();

  return spx_status;
}

absl::Status LpSoplexInterface::SoPlexSolve() {
  VLOG(10) << "calling SoPlexSolve().";
  VLOG(3) << "calling SoPlex solve(): " << spx_->numColsReal() << " cols, "
          << spx_->numRowsReal() << " rows.";

  InvalidateSolution();

  CHECK(PreStrongBranchingBasisFreed());

  // Delete the starting basis if solving from scratch.
  if (solve_from_scratch_) {
    try {
      spx_->clearBasis();
    } catch (const soplex::SPxException& x) {
      LOG(WARNING) << "SoPlex threw an exception: " << x.what() << ".";
      // since it is not clear if the status in SoPlex are set correctly we want
      // to make sure that if an error is thrown the status is not OPTIMAL
      // anymore.
      CHECK_NE(spx_->status(), soplex::SPxSolver::OPTIMAL);
      return {absl::StatusCode::kInternal, "Error"};
    }
  }

  // NOLINTBEGIN(readability-simplify-boolean-expr)
  CHECK(!solve_from_scratch_ or
        spx_->status() == soplex::SPxSolver::NO_PROBLEM);
  // NOLINTEND(readability-simplify-boolean-expr)

  const soplex::SPxSolver::Status status = LpSolve();

  VLOG(3) << " -> SoPlex status: " << status
          << ", basis status: " << spx_->basisStatus() << ".";
  is_solved_ = true;

  switch (status) {
    case soplex::SPxSolver::ABORT_TIME:
    case soplex::SPxSolver::ABORT_ITER:
    case soplex::SPxSolver::ABORT_VALUE:
    case soplex::SPxSolver::SINGULAR:
    case soplex::SPxSolver::REGULAR:
    case soplex::SPxSolver::UNKNOWN:
    case soplex::SPxSolver::OPTIMAL:
    case soplex::SPxSolver::OPTIMAL_UNSCALED_VIOLATIONS:
    case soplex::SPxSolver::UNBOUNDED:
    case soplex::SPxSolver::INFEASIBLE:
      return absl::OkStatus();
    default:
      return {absl::StatusCode::kInternal, "Error"};
  }
}

void LpSoplexInterface::SavePreStrongbranchingBasis() {
  VLOG(10) << "calling SavePreStrongbranchingBasis().";
  row_basis_status_.reSize(spx_->numRowsReal());
  col_basis_status_.reSize(spx_->numColsReal());

  try {
    spx_->getBasis(row_basis_status_.get_ptr(), col_basis_status_.get_ptr());
  } catch (const soplex::SPxException& x) {
    LOG(WARNING) << "SoPlex threw an exception: " << x.what() << ".";
    // since it is not clear if the status in SoPlex are set correctly we want
    // to make sure that if an error is thrown the status is not OPTIMAL
    // anymore.
    CHECK_NE(spx_->status(), soplex::SPxSolver::OPTIMAL);
  }
}

void LpSoplexInterface::RestorePreStrongbranchingBasis() {
  VLOG(10) << "calling RestorePreStrongbranchingBasis().";
  DCHECK_EQ(row_basis_status_.size(), spx_->numRowsReal());
  DCHECK_EQ(col_basis_status_.size(), spx_->numColsReal());

  try {
    spx_->setBasis(row_basis_status_.get_ptr(), col_basis_status_.get_ptr());
  } catch (const soplex::SPxException& x) {
    LOG(WARNING) << "SoPlex threw an exception: " << x.what() << ".";
    // since it is not clear if the status in SoPlex are set correctly we want
    // to make sure that if an error is thrown the status is not OPTIMAL
    // anymore.
    CHECK_NE(spx_->status(), soplex::SPxSolver::OPTIMAL);
  }
}

void LpSoplexInterface::InvalidateSolution() {
  VLOG(10) << "calling InvalidateSolution().";
  is_solved_ = false;
}

bool LpSoplexInterface::PreStrongBranchingBasisFreed() const {
  VLOG(10) << "calling PreStrongBranchingBasisFreed().";
  return ((row_basis_status_.size() == 0) and (col_basis_status_.size() == 0));
}

void LpSoplexInterface::FreePreStrongBranchingBasis() {
  VLOG(10) << "calling FreePreStrongBranchingBasis().";
  row_basis_status_.clear();
  col_basis_status_.clear();
}

// Strongbranching is applied to the given column, with the corresponding
// current primal solution value. The double references are used to store the
// dual bound after branching up and down. Additionally, the validity of both
// bounds is stored, if one bound is not valid it can be used as an
// estimate.
absl::StatusOr<LpInterface::StrongBranchResult>
LpSoplexInterface::SolveDownAndUpStrongBranch(ColIndex col, double primal_value,
                                              int iteration_limit) {
  VLOG(10) << "calling SolveDownAndUpStrongBranch().";
  StrongBranchResult result;
  VLOG(3) << "calling StrongBranch() on variable " << col << "("
          << iteration_limit << " iterations).";

  bool error = false;
  int old_iter_limit = spx_->intParam(soplex::SoPlex::ITERLIMIT);

  // Get the current bounds of the column.
  const double old_lb = spx_->lowerReal(col.value());
  const double old_ub = spx_->upperReal(col.value());
  result.down_valid = false;
  result.up_valid = false;
  result.iterations = 0;

  // Set the algorithm type to use the dual simplex.
  static_cast<void>(spx_->setIntParam(soplex::SoPlex::ALGORITHM,
                                      soplex::SoPlex::ALGORITHM_DUAL));

  // Computation for the down branch.
  double new_ub = std::ceil(primal_value - 1.0 - feasibility_tolerance());

  soplex::SPxSolver::Status status;
  bool from_parent_basis;
  if (new_ub >= old_lb - 0.5) {
    VLOG(3) << "strong branching down on x" << col << " (" << primal_value
            << ") with " << iteration_limit << " iterations.";
    spx_->changeUpperReal(col.value(), new_ub);
    DCHECK_LE(spx_->lowerReal(col.value()), spx_->upperReal(col.value()));

    static_cast<void>(
        spx_->setIntParam(soplex::SoPlex::ITERLIMIT, iteration_limit));
    do {
      status = spx_->optimize();

      VLOG(3) << " --> Terminate with status " << status << ".";
      switch (status) {
        case soplex::SPxSolver::OPTIMAL:
          result.dual_bound_down_branch = spx_->objValueReal();
          result.down_valid = true;
          VLOG(3) << " --> Terminate with value "
                  << result.dual_bound_down_branch << ".";
          break;
        case soplex::SPxSolver::ABORT_TIME:
          // SoPlex does not return a proven dual bound, if it is aborted.
        case soplex::SPxSolver::ABORT_ITER:
        case soplex::SPxSolver::ABORT_CYCLING:
        case soplex::SPxSolver::OPTIMAL_UNSCALED_VIOLATIONS:
          result.dual_bound_down_branch = spx_->objValueReal();
          break;
        case soplex::SPxSolver::ABORT_VALUE:
        case soplex::SPxSolver::INFEASIBLE:
          result.dual_bound_down_branch = objective_limit();
          result.down_valid = true;
          break;
        default:
          error = true;
          break;
      }
      result.iterations += spx_->numIterations();

      // We restore the pre-strong-branching basis by default (and don't solve
      // again).
      CHECK(!PreStrongBranchingBasisFreed());
      RestorePreStrongbranchingBasis();
      from_parent_basis = false;
    } while (from_parent_basis);

    spx_->changeUpperReal(col.value(), old_ub);
    DCHECK_LE(spx_->lowerReal(col.value()), spx_->upperReal(col.value()));
  } else {
    result.dual_bound_down_branch = objective_limit();
    result.down_valid = true;
  }

  // Computation for the up branch.
  if (!error) {
    double new_lb = std::floor(primal_value + 1.0 - feasibility_tolerance());
    if (new_lb <= old_ub + 0.5) {
      VLOG(3) << "strong branching  up  on x" << col << " (" << primal_value
              << ") with " << iteration_limit << "iterations.";
      spx_->changeLowerReal(col.value(), new_lb);
      DCHECK_LE(spx_->lowerReal(col.value()), spx_->upperReal(col.value()));

      static_cast<void>(
          spx_->setIntParam(soplex::SoPlex::ITERLIMIT, iteration_limit));
      do {
        status = spx_->optimize();

        VLOG(3) << " --> Terminate with status " << status << ".";
        switch (status) {
          case soplex::SPxSolver::OPTIMAL:
            result.dual_bound_up_branch = spx_->objValueReal();
            result.up_valid = true;
            VLOG(3) << " --> Terminate with value " << spx_->objValueReal()
                    << ".";
            break;
          case soplex::SPxSolver::ABORT_TIME:
            // SoPlex does not return a proven dual bound, if it is aborted.
          case soplex::SPxSolver::ABORT_ITER:
          case soplex::SPxSolver::ABORT_CYCLING:
          case soplex::SPxSolver::OPTIMAL_UNSCALED_VIOLATIONS:
            result.dual_bound_up_branch = spx_->objValueReal();
            break;
          case soplex::SPxSolver::ABORT_VALUE:
          case soplex::SPxSolver::INFEASIBLE:
            result.dual_bound_up_branch = objective_limit();
            result.up_valid = true;
            break;
          default:
            error = true;
            break;
        }
        result.iterations += spx_->numIterations();

        // We restore the pre-strong-branching basis by default (and don't solve
        // again).
        CHECK(!PreStrongBranchingBasisFreed());
        RestorePreStrongbranchingBasis();
        from_parent_basis = false;
      } while (from_parent_basis);

      spx_->changeLowerReal(col.value(), old_lb);
      DCHECK_LE(spx_->lowerReal(col.value()), spx_->upperReal(col.value()));
    } else {
      result.dual_bound_up_branch = objective_limit();
      result.up_valid = true;
    }
  }
  // Reset the old iteration limit.
  static_cast<void>(
      spx_->setIntParam(soplex::SoPlex::ITERLIMIT, old_iter_limit));

  if (error) {
    VLOG(10) << "StrongBranch() returned SoPlex status "
             << static_cast<int>(status) << ".";
    return absl::InternalError("Error");
  }
  return result;
}

// ==========================================================================
// LP model setters.
// ==========================================================================

absl::Status LpSoplexInterface::PopulateFromMipData(const MipData& mip_data) {
  VLOG(10) << "calling PopulateFromMipData().";

  RETURN_IF_ERROR(Clear());
  DCHECK_EQ(mip_data.constraint_names().size(),
            mip_data.left_hand_sides().size());
  DCHECK_EQ(mip_data.left_hand_sides().size(),
            mip_data.right_hand_sides().size());
  DCHECK_EQ(mip_data.variable_names().size(), mip_data.lower_bounds().size());
  DCHECK_EQ(mip_data.lower_bounds().size(), mip_data.upper_bounds().size());

  InvalidateSolution();

  CHECK(PreStrongBranchingBasisFreed());

  try {
    int num_rows = mip_data.left_hand_sides().size();
    soplex::LPRowSet rowset(num_rows);
    soplex::DSVector empty_vector(0);

    spx_->clearLPReal();

    static_cast<void>(spx_->setIntParam(soplex::SoPlex::OBJSENSE,
                                        soplex::SoPlex::OBJSENSE_MINIMIZE));

    // Create empty rows with the given sides.
    for (RowIndex row(0); row < num_rows; ++row) {
      rowset.add(mip_data.left_hand_sides()[row], empty_vector,
                 mip_data.right_hand_sides()[row]);
    }
    spx_->addRowsReal(rowset);

    // Create the column vectors with the given coefficients and bounds.
    absl::StrongVector<ColIndex, double> objective_coefficients(
        mip_data.lower_bounds().size(), 0.0);
    for (auto [col, coefficient] : mip_data.objective().entries()) {
      objective_coefficients[col] = coefficient;
    }
    RETURN_IF_ERROR(AddColumns(mip_data.matrix(), mip_data.lower_bounds(),
                               mip_data.upper_bounds(), objective_coefficients,
                               mip_data.variable_names()));
  } catch (const soplex::SPxException& x) {
    LOG(WARNING) << "SoPlex threw an exception: " << x.what() << ".";
    return {absl::StatusCode::kInternal, "Error"};
  }
  return absl::OkStatus();
}

absl::Status LpSoplexInterface::AddColumn(const SparseCol& col_data,
                                          double lower_bound,
                                          double upper_bound,
                                          double objective_coefficient,
                                          const std::string& /*unused*/) {
  VLOG(10) << "calling AddColumn().";
  DCHECK(!col_data.MayNeedCleaning());
  DCHECK(std::all_of(col_data.entries().begin(), col_data.entries().end(),
                     [num_rows = GetNumberOfRows()](const ColEntry& e) {
                       return RowIndex(0) <= e.index and e.index < num_rows;
                     }));
  InvalidateSolution();

  CHECK(PreStrongBranchingBasisFreed());

  if (DEBUG_MODE) {
    // Perform a check to ensure that no new rows have been added.
    RowIndex num_rows = GetNumberOfRows();
    for (int j = 0; j < col_data.entries().size(); ++j) {
      DCHECK_LE(0, col_data.indices()[j]);
      DCHECK_LT(col_data.indices()[j], num_rows);
      DCHECK_NE(col_data.values()[j], 0.0);
    }
  }

  try {
    soplex::DSVector col_vector;
    for (SparseEntry entry : col_data.entries()) {
      DCHECK(!IsInfinity(std::abs(entry.value)));
      col_vector.add(entry.index.value(), entry.value);
    }
    VLOG(3) << "Adding column with objective coefficient = "
            << objective_coefficient << ", lower bound = " << lower_bound
            << ", upper bound = " << upper_bound << ".";
    spx_->addColReal(soplex::LPCol(objective_coefficient, col_vector,
                                   upper_bound, lower_bound));
  } catch (const soplex::SPxException& x) {
    LOG(WARNING) << "SoPlex threw an exception: " << x.what() << ".";
    return {absl::StatusCode::kInternal, "Error"};
  }

  return absl::OkStatus();
}

absl::Status LpSoplexInterface::AddColumns(
    const StrongSparseMatrix& matrix,
    const absl::StrongVector<ColIndex, double>& lower_bounds,
    const absl::StrongVector<ColIndex, double>& upper_bounds,
    const absl::StrongVector<ColIndex, double>& objective_coefficients,
    const absl::StrongVector<ColIndex, std::string>& /*unused*/) {
  VLOG(10) << "calling AddColumns().";

  InvalidateSolution();

  CHECK(PreStrongBranchingBasisFreed());

  if (DEBUG_MODE) {
    if (!(matrix.num_cols() == 0) and matrix.AllColsAreClean()) {
      // Perform a check to ensure that no new rows have been added.
      int num_rows = spx_->numRowsReal();
      for (ColIndex col(0); col < matrix.num_cols(); ++col) {
        for (int row = 0; row < matrix.col(col).entries().size(); ++row) {
          DCHECK_LE(0, matrix.col(col).indices()[row]);
          DCHECK_LT(matrix.col(col).indices()[row], num_rows);
          DCHECK_NE(matrix.col(col).values()[row], 0.0);
        }
      }
    }
  }

  try {
    soplex::LPColSet columns(matrix.num_cols().value());
    soplex::DSVector col_vector(matrix.num_cols().value());

    // Create column vectors with coefficients and bounds.
    for (ColIndex col(0); col < matrix.num_cols(); ++col) {
      col_vector.clear();
      if (!matrix.col(col).entries().empty()) {
        for (SparseEntry entry : matrix.col(col).entries()) {
          DCHECK(!IsInfinity(std::abs(entry.value)));
          col_vector.add(entry.index.value(), entry.value);
        }
      }
      columns.add(objective_coefficients[col], lower_bounds[col], col_vector,
                  upper_bounds[col]);
    }
    spx_->addColsReal(columns);
  } catch (const soplex::SPxException& x) {
    LOG(WARNING) << "SoPlex threw an exception: " << x.what() << ".";
    return {absl::StatusCode::kInternal, "Error"};
  }

  return absl::OkStatus();
}

absl::Status LpSoplexInterface::DeleteColumns(ColIndex first_col,
                                              ColIndex last_col) {
  VLOG(10) << "calling DeleteColumns().";

  DCHECK_LE(0, first_col);
  DCHECK_LE(first_col, last_col);
  DCHECK_LT(last_col, spx_->numColsReal());

  InvalidateSolution();

  CHECK(PreStrongBranchingBasisFreed());

  try {
    spx_->removeColRangeReal(first_col.value(), last_col.value());
  } catch (const soplex::SPxMemoryException& e) {
    LOG(ERROR) << "SoPlex threw a memory exception: " << e.what() << ".";
    return {absl::StatusCode::kInternal, "Error"};
  } catch (const soplex::SPxException& e) {
    LOG(ERROR) << "SoPlex threw an exception: " << e.what() << ".";
    return {absl::StatusCode::kInternal, "Error"};
  }

  return absl::OkStatus();
}

absl::Status LpSoplexInterface::AddRow(const SparseRow& row_data,
                                       double left_hand_side,
                                       double right_hand_side,
                                       const std::string& /*unused*/) {
  if (DEBUG_MODE) {
    // Perform a check to ensure that no new rows have been added.
    const int num_cols = spx_->numColsReal();
    for (int j = 0; j < row_data.entries().size(); ++j) {
      DCHECK_NE(row_data.values()[j], 0.0);
      DCHECK_LE(0, row_data.indices()[j]);
      DCHECK_LT(row_data.indices()[j], num_cols);
    }
  }

  InvalidateSolution();

  soplex::DSVector row_vector;
  for (SparseEntry entry : row_data.entries()) {
    DCHECK(!IsInfinity(std::abs(entry.value)));
    row_vector.add(entry.index.value(), entry.value);
  }
  spx_->addRowReal(soplex::LPRow(left_hand_side, row_vector, right_hand_side));

  return absl::OkStatus();
}

absl::Status LpSoplexInterface::AddRows(
    const absl::StrongVector<RowIndex, SparseRow>& rows,
    const absl::StrongVector<RowIndex, double>& left_hand_sides,
    const absl::StrongVector<RowIndex, double>& right_hand_sides,
    const absl::StrongVector<RowIndex, std::string>& /*unused*/) {
  VLOG(10) << "calling AddRows().";

  InvalidateSolution();

  CHECK(PreStrongBranchingBasisFreed());

  if (DEBUG_MODE) {
    for (RowIndex i(0); i < rows.size(); ++i) {
      if (!rows[i].entries().empty()) {
        // Perform a check to ensure that no new rows have been added.
        const int num_cols = spx_->numColsReal();
        for (int j = 0; j < rows[i].entries().size(); ++j) {
          DCHECK_NE(rows[i].values()[j], 0.0);
          DCHECK_LE(0, rows[i].indices()[j]);
          DCHECK_LT(rows[i].indices()[j], num_cols);
        }
      }
    }
  }

  try {
    soplex::LPRowSet row_set(rows.size());
    soplex::DSVector row_vector;

    // Create row vectors from the values of the given sides.
    for (RowIndex i(0); i < rows.size(); ++i) {
      row_vector.clear();
      if (!rows[i].entries().empty()) {
        for (SparseEntry entry : rows[i].entries()) {
          row_vector.add(entry.index.value(), entry.value);
        }
      }
      row_set.add(left_hand_sides[i], row_vector, right_hand_sides[i]);
    }
    spx_->addRowsReal(row_set);
  } catch (const soplex::SPxException& x) {
    LOG(WARNING) << "SoPlex threw an exception: " << x.what() << ".";
    return {absl::StatusCode::kInternal, "Error"};
  }
  return absl::OkStatus();
}

absl::Status LpSoplexInterface::DeleteRows(RowIndex first_row,
                                           RowIndex last_row) {
  VLOG(10) << "calling DeleteRows().";

  DCHECK_LE(0, first_row);
  DCHECK_LE(first_row, last_row);
  DCHECK_LT(last_row, spx_->numRowsReal());

  InvalidateSolution();

  CHECK(PreStrongBranchingBasisFreed());

  SOPLEX_TRY(spx_->removeRowRangeReal(first_row.value(), last_row.value()));

  return absl::OkStatus();
}

// Deletes a set of rows from LP; the new position of a row must not be greater
// than its old position.
absl::StatusOr<absl::StrongVector<RowIndex, RowIndex>>
LpSoplexInterface::DeleteRowSet(
    const absl::StrongVector<RowIndex, bool>& rows_to_delete) {
  VLOG(10) << "calling DeleteRowSet().";

  InvalidateSolution();

  CHECK(PreStrongBranchingBasisFreed());

  std::vector<int> int_deletion_status(rows_to_delete.begin(),
                                       rows_to_delete.end());
  int num_rows = spx_->numRowsReal();

  // SoPlex removeRows() method deletes the rows with deletion_status[i] < 0,
  // so we have to negate the values.
  for (int i = 0; i < num_rows; ++i) int_deletion_status[i] *= -1;

  SOPLEX_TRY(spx_->removeRowsReal(int_deletion_status.data()));

  // Note: We cannot use the same vector for the Soplex call above and the row
  // mapping, even though RowIndex internally consists of an int. In addition to
  // being bad practice to depend on a class' internals, the compiler may
  // introduce additional padding meaning RowIndex and int don't have identical
  // representations.
  absl::StrongVector<RowIndex, RowIndex> row_mapping(rows_to_delete.size());
  for (RowIndex row(0); row < row_mapping.size(); ++row) {
    if (rows_to_delete[row]) {
      row_mapping[row] = kInvalidRow;
    } else {
      row_mapping[row] = int_deletion_status[row.value()];
    }
  }
  return row_mapping;
}

// Clears the whole LP.
absl::Status LpSoplexInterface::Clear() {
  VLOG(10) << "calling LpSoplexInterface::Clear().";

  InvalidateSolution();

  CHECK(PreStrongBranchingBasisFreed());
  SOPLEX_TRY(spx_->clearLPReal());

  return absl::OkStatus();
}

// Clears current LP Interface state (like basis information) of the solver.
absl::Status LpSoplexInterface::ClearState() {
  VLOG(10) << "calling LpSoplexInterface::ClearState().";

  InvalidateSolution();

  try {
    spx_->clearBasis();
  } catch (const soplex::SPxException& x) {
    LOG(WARNING) << "SoPlex threw an exception: " << x.what() << ".";
    CHECK_NE(spx_->status(), soplex::SPxSolver::OPTIMAL);
    return {absl::StatusCode::kInternal, "Error"};
  }
  return absl::OkStatus();
}

absl::Status LpSoplexInterface::SetColumnBounds(ColIndex col,
                                                double lower_bound,
                                                double upper_bound) {
  VLOG(10) << "calling SetColumnBounds().";
  DCHECK(!IsInfinity(lower_bound));
  DCHECK(!IsInfinity(-upper_bound));

  InvalidateSolution();

  CHECK(PreStrongBranchingBasisFreed());
  CHECK(0 <= col);
  CHECK(col < spx_->numColsReal());

  try {
    spx_->changeBoundsReal(col.value(), lower_bound, upper_bound);
    DCHECK_LE(spx_->lowerReal(col.value()),
              spx_->upperReal(col.value()) +
                  spx_->realParam(soplex::SoPlex::EPSILON_ZERO));
  } catch (const soplex::SPxException& x) {
    LOG(WARNING) << "SoPlex threw an exception: " << x.what() << ".";
    return {absl::StatusCode::kInternal, "Error"};
  }
  return absl::OkStatus();
}

absl::Status LpSoplexInterface::SetRowSides(RowIndex row, double left_hand_side,
                                            double right_hand_side) {
  VLOG(10) << "calling SetRowSides().";
  DCHECK_GE(row, 0);
  DCHECK_LT(row, GetNumberOfRows());
  DCHECK(!IsInfinity(left_hand_side));
  DCHECK(!IsInfinity(-right_hand_side));

  InvalidateSolution();

  CHECK(PreStrongBranchingBasisFreed());

  try {
    spx_->changeRangeReal(row.value(), left_hand_side, right_hand_side);
    DCHECK_LE(spx_->lhsReal(row.value()),
              spx_->rhsReal(row.value()) +
                  spx_->realParam(soplex::SoPlex::EPSILON_ZERO));
  } catch (const soplex::SPxException& x) {
    LOG(WARNING) << "SoPlex threw an exception: " << x.what() << ".";
    return {absl::StatusCode::kInternal, "LP Error"};
  }
  return absl::OkStatus();
}

absl::Status LpSoplexInterface::SetObjectiveSense(bool is_maximization) {
  VLOG(10) << "calling SetObjectiveSense().";

  InvalidateSolution();

  CHECK(PreStrongBranchingBasisFreed());

  SOPLEX_TRY(static_cast<void>(spx_->setIntParam(
      soplex::SoPlex::OBJSENSE, is_maximization == 0
                                    ? soplex::SoPlex::OBJSENSE_MINIMIZE
                                    : soplex::SoPlex::OBJSENSE_MAXIMIZE)));

  return absl::OkStatus();
}

absl::Status LpSoplexInterface::SetObjectiveCoefficient(
    ColIndex col, double objective_coefficient) {
  VLOG(10) << "calling SetObjectiveCoefficient().";
  DCHECK(!IsInfinity(std::abs(objective_coefficient)));

  InvalidateSolution();

  CHECK(PreStrongBranchingBasisFreed());

  spx_->changeObjReal(col.value(), objective_coefficient);
  return absl::OkStatus();
}

// ==========================================================================
// LP model getters.
// ==========================================================================

RowIndex LpSoplexInterface::GetNumberOfRows() const {
  VLOG(10) << "calling GetNumberOfRows().";

  return RowIndex(spx_->numRowsReal());
}

ColIndex LpSoplexInterface::GetNumberOfColumns() const {
  VLOG(10) << "calling GetNumberOfColumns().";

  return ColIndex(spx_->numColsReal());
}

int64_t LpSoplexInterface::GetNumberOfNonZeros() const {
  VLOG(10) << "calling GetNumberOfNonZeros().";
  // SoPlex has no direct method to return the number of nonzeros, so we have
  // to count them manually.
  int num_non_zeros = 0;

  if (spx_->numRowsReal() < spx_->numColsReal()) {
    for (int i = 0; i < spx_->numRowsReal(); ++i) {
      num_non_zeros += spx_->rowVectorRealInternal(i).size();
    }
  } else {
    for (int i = 0; i < spx_->numColsReal(); ++i) {
      num_non_zeros += spx_->colVectorRealInternal(i).size();
    }
  }

  return num_non_zeros;
}

bool LpSoplexInterface::IsMaximization() const {
  VLOG(10) << "calling IsMaximization().";

  return spx_->intParam(soplex::SoPlex::OBJSENSE) !=
         soplex::SoPlex::OBJSENSE_MINIMIZE;
}

// Either both, lower_bound and upper_bound, have to be 0, or both have to be
// non-0, either n_non_zeroes, begin_cols, indices, and obj_coeffs have to be
// 0, or all of them have to be non-0.
SparseCol LpSoplexInterface::GetSparseColumnCoefficients(ColIndex col) const {
  VLOG(10) << "calling GetSparseColumnCoefficients().";

  DCHECK_LE(0, col);
  DCHECK_LT(col, spx_->numColsReal());

  SparseCol sparse_column;

  if (spx_->boolParam(soplex::SoPlex::PERSISTENTSCALING)) {
    soplex::DSVector cvec;
    spx_->getColVectorReal(col.value(), cvec);
    for (int j = 0; j < cvec.size(); ++j) {
      sparse_column.AddEntry(RowIndex(cvec.index(j)), cvec.value(j));
    }
  } else {
    const soplex::SVector& cvec = spx_->colVectorRealInternal(col.value());
    for (int j = 0; j < cvec.size(); ++j) {
      sparse_column.AddEntry(RowIndex(cvec.index(j)), cvec.value(j));
    }
  }
  sparse_column.CleanUpIfNeeded();
  return sparse_column;
}

// Either both, left_hand_side and right_hand_side, have to be 0, or both have
// to be non-0, either n_non_zeroes, begin_cols, indices, and obj_coeffs have
// to be 0, or all of them have to be non-0.
SparseRow LpSoplexInterface::GetSparseRowCoefficients(RowIndex row) const {
  VLOG(10) << "calling GetSparseRowCoefficients().";

  DCHECK_LE(0, row);
  DCHECK_LT(row, spx_->numRowsReal());

  SparseRow sparse_row;

  if (spx_->boolParam(soplex::SoPlex::PERSISTENTSCALING)) {
    soplex::DSVector rvec;
    spx_->getRowVectorReal(row.value(), rvec);
    for (int j = 0; j < rvec.size(); ++j) {
      sparse_row.AddEntry(ColIndex(rvec.index(j)), rvec.value(j));
    }
  } else {
    const soplex::SVector& rvec = spx_->rowVectorRealInternal(row.value());
    for (int j = 0; j < rvec.size(); ++j) {
      sparse_row.AddEntry(ColIndex(rvec.index(j)), rvec.value(j));
    }
  }
  return sparse_row;
}

double LpSoplexInterface::GetObjectiveCoefficient(ColIndex col) const {
  VLOG(10) << "calling GetObjectiveCoefficient().";

  DCHECK_LE(0, col);
  DCHECK_LT(col, spx_->numColsReal());

  return spx_->objReal(col.value());
}

double LpSoplexInterface::GetLowerBound(ColIndex col) const {
  VLOG(10) << "calling GetLowerBound().";

  DCHECK_LE(0, col);
  DCHECK_LT(col, spx_->numColsReal());

  return spx_->lowerReal(col.value());
}

double LpSoplexInterface::GetUpperBound(ColIndex col) const {
  VLOG(10) << "calling GetUpperBound().";

  DCHECK_LE(0, col);
  DCHECK_LT(col, spx_->numColsReal());

  return spx_->upperReal(col.value());
}

double LpSoplexInterface::GetLeftHandSide(RowIndex row) const {
  VLOG(10) << "calling GetLeftHandSide().";

  DCHECK_LE(0, row);
  DCHECK_LT(row, spx_->numRowsReal());

  return spx_->lhsReal(row.value());
}

double LpSoplexInterface::GetRightHandSide(RowIndex row) const {
  VLOG(10) << "calling GetRightHandSide().";

  DCHECK_LE(0, row);
  DCHECK_LT(row, spx_->numRowsReal());

  return spx_->rhsReal(row.value());
}

double LpSoplexInterface::GetMatrixCoefficient(ColIndex col,
                                               RowIndex row) const {
  VLOG(10) << "calling GetMatrixCoefficient().";

  DCHECK_LE(0, col);
  DCHECK_LT(col, spx_->numColsReal());
  DCHECK_LE(0, row);
  DCHECK_LT(row, spx_->numRowsReal());

  return spx_->coefReal(row.value(), col.value());
}

// ==========================================================================
// Solving methods.
// ==========================================================================

absl::Status LpSoplexInterface::SolveLpWithPrimalSimplex() {
  VLOG(10) << "calling SolveLpWithPrimalSimplex().";

  static_cast<void>(spx_->setIntParam(soplex::SoPlex::ALGORITHM,
                                      soplex::SoPlex::ALGORITHM_PRIMAL));
  return SoPlexSolve();
}

absl::Status LpSoplexInterface::SolveLpWithDualSimplex() {
  VLOG(10) << "calling SolveLPWithDualSimplex().";

  static_cast<void>(spx_->setIntParam(soplex::SoPlex::ALGORITHM,
                                      soplex::SoPlex::ALGORITHM_DUAL));
  return SoPlexSolve();
}

// This call is needed before any strong branching.
absl::Status LpSoplexInterface::StartStrongBranching() {
  VLOG(10) << "calling StartStrongBranching().";
  CHECK(PreStrongBranchingBasisFreed());
  SavePreStrongbranchingBasis();

  return absl::OkStatus();
}

// This call is needed after any strong branching.
absl::Status LpSoplexInterface::EndStrongBranching() {
  VLOG(10) << "calling EndStrongBranching().";
  CHECK(!PreStrongBranchingBasisFreed());
  RestorePreStrongbranchingBasis();
  FreePreStrongBranchingBasis();

  return absl::OkStatus();
}

// ============================================================================
// Solution information getters.
// ============================================================================

// Returns whether a solve method was called after the last modification of
// the LP.
bool LpSoplexInterface::IsSolved() const {
  VLOG(10) << "calling IsSolved().";
  return is_solved_;
}

// Returns true if the LP is proven to have a primal unbounded ray (but not
// necessary a primal feasible point); this does not necessarily mean that the
// solver knows and can return the primal ray.
bool LpSoplexInterface::ExistsPrimalRay() const {
  VLOG(10) << "calling ExistsPrimalRay().";

  return (spx_->status() == soplex::SPxSolver::UNBOUNDED);
}

// Returns true iff LP is proven to have a primal unbounded ray (but not
// necessary a primal feasible point), and the solver knows and can return the
// primal ray.
bool LpSoplexInterface::HasPrimalRay() const {
  VLOG(10) << "calling HasPrimalRay().";

  return spx_->hasPrimalRay();
}

bool LpSoplexInterface::IsPrimalUnbounded() const {
  VLOG(10) << "calling IsPrimalUnbounded().";

  // If SoPlex returns unbounded, this may only mean that an unbounded ray is
  // available, not necessarily a primal
  //* feasible point; hence we have to check the perturbation.
  return spx_->status() == soplex::SPxSolver::UNBOUNDED;
}

bool LpSoplexInterface::IsPrimalInfeasible() const {
  VLOG(10) << "calling IsPrimalInfeasible().";

  return (spx_->status() == soplex::SPxSolver::INFEASIBLE);
}

bool LpSoplexInterface::IsPrimalFeasible() const {
  VLOG(10) << "calling IsPrimalFeasible().";

  return spx_->basisStatus() == soplex::SPxBasis::OPTIMAL or
         spx_->basisStatus() == soplex::SPxBasis::PRIMAL;
}

// Returns true if the LP is proven to have a dual unbounded ray (but not
// necessary a dual feasible point);
// this does not necessarily mean that the solver knows and can return the
// dual ray.
bool LpSoplexInterface::ExistsDualRay() const {
  VLOG(10) << "calling ExistsDualRay().";

  return (spx_->status() == soplex::SPxSolver::INFEASIBLE);
}

// Returns true if the LP is proven to have a dual unbounded ray (but not
// necessary a dual feasible point),
//*  and the solver knows and can return the dual ray
bool LpSoplexInterface::HasDualRay() const {
  VLOG(10) << "calling HasDualRay().";

  return spx_->hasDualFarkas();
}

bool LpSoplexInterface::IsDualUnbounded() const {
  VLOG(10) << "calling IsDualUnbounded().";

  return spx_->status() == soplex::SPxSolver::INFEASIBLE &&
         spx_->basisStatus() == soplex::SPxBasis::DUAL;
}

bool LpSoplexInterface::IsDualInfeasible() const {
  VLOG(10) << "calling IsDualInfeasible().";

  return (spx_->status() == soplex::SPxSolver::UNBOUNDED);
}

bool LpSoplexInterface::IsDualFeasible() const {
  VLOG(10) << "calling IsDualFeasible().";

  return (spx_->basisStatus() == soplex::SPxBasis::OPTIMAL) or
         spx_->basisStatus() == soplex::SPxBasis::DUAL;
}

bool LpSoplexInterface::IsOptimal() const {
  VLOG(10) << "calling IsOptimal().";

  CHECK((spx_->basisStatus() == soplex::SPxBasis::OPTIMAL) ==
        (IsPrimalFeasible() and IsDualFeasible()));

  return (spx_->status() == soplex::SPxSolver::OPTIMAL);
}

// This function should return true if the solution is reliable, i.e., feasible
// and optimal (or proven infeasible/unbounded) with respect to the original
// problem. The optimality status might be with respect to a scaled version of
// the problem, but the solution might not be feasible to the unscaled original
// problem; in this case, minimip::LpInterface.IsStable() should return false.
bool LpSoplexInterface::IsStable() const {
  VLOG(10) << "calling IsStable().";

  if (spx_->status() == soplex::SPxSolver::ERROR or
      spx_->status() == soplex::SPxSolver::SINGULAR) {
    return false;
  }
  if (spx_->status() == soplex::SPxSolver::OPTIMAL_UNSCALED_VIOLATIONS) {
    return false;
  }
  return true;
}

bool LpSoplexInterface::ObjectiveLimitIsExceeded() const {
  VLOG(10) << "calling ObjectiveLimitIsExceeded().";

  return (spx_->status() == soplex::SPxSolver::ABORT_VALUE);
}

bool LpSoplexInterface::IterationLimitIsExceeded() const {
  VLOG(10) << "calling IterationLimitIsExceeded().";

  return (spx_->status() == soplex::SPxSolver::ABORT_ITER);
}

bool LpSoplexInterface::TimeLimitIsExceeded() const {
  VLOG(10) << "calling TimeLimitIsExceeded().";

  return (spx_->status() == soplex::SPxSolver::ABORT_TIME);
}

double LpSoplexInterface::GetObjectiveValue() const {
  VLOG(10) << "calling GetObjectiveValue().";

  return spx_->objValueReal();
}

// Before calling this function, the caller must ensure that the LP has been
// solved to optimality, i.e., that minimip::LpInterface.IsOptimal() returns
// true.

absl::StatusOr<absl::StrongVector<ColIndex, double>>
LpSoplexInterface::GetPrimalValues() const {
  VLOG(10) << "calling GetPrimalValues().";
  CHECK(IsOptimal());
  absl::StrongVector<ColIndex, double> primal_values(spx_->numColsReal());
  try {
    static_cast<void>(
        spx_->getPrimalReal(primal_values.data(), spx_->numColsReal()));
  } catch (const soplex::SPxException& x) {
    LOG(WARNING) << "SoPlex threw an exception: " << x.what() << ".";
    return absl::Status(absl::StatusCode::kInternal, "Error");
  }
  return primal_values;
}

// Before calling this function, the caller must ensure that the LP has been
// solved to optimality, i.e., that minimip::LpInterface.IsOptimal() returns
// true.
absl::StatusOr<absl::StrongVector<RowIndex, double>>
LpSoplexInterface::GetDualValues() const {
  VLOG(10) << "calling GetDualSolution().";
  CHECK(IsOptimal());
  absl::StrongVector<RowIndex, double> dual_values(spx_->numRowsReal());
  try {
    static_cast<void>(
        spx_->getDualReal(dual_values.data(), spx_->numRowsReal()));
  } catch (const soplex::SPxException& x) {
    LOG(WARNING) << "SoPlex threw an exception: " << x.what() << ".";
    return absl::Status(absl::StatusCode::kInternal, "Error");
  }
  return dual_values;
}

absl::StatusOr<absl::StrongVector<RowIndex, double>>
LpSoplexInterface::GetRowActivities() const {
  VLOG(10) << "calling GetRowActivities().";
  absl::StrongVector<RowIndex, double> row_activities(spx_->numRowsReal());
  try {
    static_cast<void>(
        spx_->getSlacksReal(row_activities.data(),
                            spx_->numRowsReal()));  // in SoPlex, the activities
    // are called "slacks"
  } catch (const soplex::SPxException& x) {
    LOG(WARNING) << "SoPlex threw an exception: " << x.what() << ".";
    return absl::Status(absl::StatusCode::kInternal, "Error");
  }
  return row_activities;
}

absl::StatusOr<absl::StrongVector<ColIndex, double>>
LpSoplexInterface::GetReducedCosts() const {
  VLOG(10) << "calling GetReducedCosts().";
  absl::StrongVector<ColIndex, double> reduced_costs(spx_->numColsReal());
  try {
    static_cast<void>(
        spx_->getRedCostReal(reduced_costs.data(), spx_->numColsReal()));
  } catch (const soplex::SPxException& x) {
    LOG(WARNING) << "SoPlex threw an exception: " << x.what() << ".";
    return absl::Status(absl::StatusCode::kInternal, "Error");
  }
  return reduced_costs;
}

absl::StatusOr<absl::StrongVector<ColIndex, double>>
LpSoplexInterface::GetPrimalRay() const {
  VLOG(10) << "calling GetPrimalRay().";
  absl::StrongVector<ColIndex, double> primal_ray;
  CHECK(spx_->hasPrimalRay());
  try {
    static_cast<void>(
        spx_->getPrimalRayReal(primal_ray.data(), spx_->numColsReal()));
  } catch (const soplex::SPxException& x) {
    LOG(WARNING) << "SoPlex threw an exception: " << x.what() << ".";
    return absl::Status(absl::StatusCode::kInternal, "Error");
  }
  return primal_ray;
}

absl::StatusOr<absl::StrongVector<RowIndex, double>>
LpSoplexInterface::GetDualRay() const {
  VLOG(10) << "calling GetDualRay().";
  CHECK(spx_->hasDualFarkas());
  absl::StrongVector<RowIndex, double> dual_ray;
  try {
    static_cast<void>(
        spx_->getDualFarkasReal(dual_ray.data(), spx_->numRowsReal()));
  } catch (const soplex::SPxException& x) {
    LOG(WARNING) << "SoPlex threw an exception: " << x.what() << ".";
    return absl::Status(absl::StatusCode::kInternal, "Error");
  }
  return dual_ray;
}

// Gets the number of LP iterations of the last solve call.
int64_t LpSoplexInterface::GetNumIterations() const {
  VLOG(10) << "calling GetNumIterations().";

  return spx_->numIterations();
}

// ==========================================================================
// Getters and setters of the basis.
// ==========================================================================

absl::StatusOr<absl::StrongVector<ColIndex, LpBasisStatus>>
LpSoplexInterface::GetBasisStatusForColumns() const {
  absl::StrongVector<ColIndex, LpBasisStatus> statuses(spx_->numColsReal());
  VLOG(10) << "calling GetBasisStatusForColumns().";
  CHECK(PreStrongBranchingBasisFreed());
  CHECK(IsOptimal());

  if (!statuses.empty()) {
    for (ColIndex i(0); i < spx_->numColsReal(); ++i) {
      //         double obj_coeffs = 0.0;
      switch (spx_->basisColStatus(i.value())) {
        case soplex::SPxSolver::BASIC:
          statuses[i] = LpBasisStatus::kBasic;
          break;
        case soplex::SPxSolver::FIXED:
          // Get reduced cost estimation. If the estimation is not correct
          // this should not hurt: If the basis is loaded into SoPlex again,
          // the status is converted to FIXED again; in this case there is no
          // problem at all. If the basis is saved and/or used in some other
          // solver, it usually is very cheap to perform the pivots necessary
          // to get an optimal basis.
          // @todo implement getRedCostEst()
          //
          // RETURN_IF_ERROR( getRedCostEst(spx_, i, &obj_coeffs) );
          //             if( obj_coeffs < 0.0 )  // reduced costs < 0 => UPPER
          //             else => LOWER
          //                column_basis_status.push_back(BASESTAT_UPPER;
          //             else
          [[fallthrough]];
        case soplex::SPxSolver::ON_LOWER:
          statuses[i] = LpBasisStatus::kAtLowerBound;
          break;
        case soplex::SPxSolver::ON_UPPER:
          statuses[i] = LpBasisStatus::kAtUpperBound;
          break;
        case soplex::SPxSolver::ZERO:
          statuses[i] = LpBasisStatus::kFree;
          break;
        case soplex::SPxSolver::UNDEFINED:
        default:
          LOG(ERROR) << "invalid basis status.";
          std::abort();
          return absl::Status(absl::StatusCode::kInvalidArgument,
                              "Invalid Data");
      }
    }
  }
  return statuses;
}

absl::StatusOr<absl::StrongVector<RowIndex, LpBasisStatus>>
LpSoplexInterface::GetBasisStatusForRows() const {
  absl::StrongVector<RowIndex, LpBasisStatus> statuses(spx_->numRowsReal());
  VLOG(10) << "calling GetBasisStatusForRows().";
  CHECK(PreStrongBranchingBasisFreed());
  CHECK(IsOptimal());

  if (!statuses.empty()) {
    for (RowIndex i(0); i < spx_->numRowsReal(); ++i) {
      switch (spx_->basisRowStatus(i.value())) {
        case soplex::SPxSolver::BASIC:
          statuses[i] = LpBasisStatus::kBasic;
          break;
        case soplex::SPxSolver::FIXED:
          statuses[i] = LpBasisStatus::kFixed;
        case soplex::SPxSolver::ON_LOWER:
          statuses[i] = LpBasisStatus::kAtLowerBound;
          break;
        case soplex::SPxSolver::ON_UPPER:
          statuses[i] = LpBasisStatus::kAtUpperBound;
          break;
        case soplex::SPxSolver::ZERO:
          LOG(ERROR)
              << "slack variable has basis status ZERO (should not occur).";
          return absl::Status(absl::StatusCode::kInternal, "Error");
        case soplex::SPxSolver::UNDEFINED:
        default:
          LOG(ERROR) << "invalid basis status.";
          std::abort();
          return absl::Status(absl::StatusCode::kInvalidArgument,
                              "Invalid Data");
      }
    }
  }
  return statuses;
}

absl::Status LpSoplexInterface::SetBasisStatusForColumnsAndRows(
    const absl::StrongVector<ColIndex, LpBasisStatus>& column_basis_statuses,
    const absl::StrongVector<RowIndex, LpBasisStatus>& row_basis_statuses) {
  VLOG(10) << "calling SetBasisStatusForColumnsAndRows().";

  ColIndex num_cols = GetNumberOfColumns();
  RowIndex num_rows = GetNumberOfRows();

  CHECK(PreStrongBranchingBasisFreed());
  InvalidateSolution();

  soplex::DataArray<soplex::SPxSolver::VarStatus>& colstat =
      ColumnsBasisStatus();
  soplex::DataArray<soplex::SPxSolver::VarStatus>& rowstat = RowsBasisStatus();

  colstat.reSize(num_cols.value());
  rowstat.reSize(num_rows.value());

  for (int i = 0; i < num_rows; ++i) {
    switch (row_basis_statuses[RowIndex(i)]) {
      case LpBasisStatus::kBasic:
        rowstat[i] = soplex::SPxSolver::BASIC;
        break;
      case LpBasisStatus::kAtLowerBound:
        rowstat[i] = soplex::SPxSolver::ON_LOWER;
        break;
      case LpBasisStatus::kAtUpperBound:
        rowstat[i] = soplex::SPxSolver::ON_UPPER;
        break;
      case LpBasisStatus::kFixed:
        rowstat[i] = soplex::SPxSolver::FIXED;
        break;
      case LpBasisStatus::kFree:
        LOG(ERROR)
            << "slack variable has basis status ZERO (should not occur).";
        return {absl::StatusCode::kInternal, "Error"};
      default:
        LOG(ERROR) << "invalid basis status.";
        std::abort();
        return {absl::StatusCode::kInvalidArgument, "Invalid Data"};
    }
  }

  for (int i = 0; i < num_cols; ++i) {
    switch (column_basis_statuses[ColIndex(i)]) {
      case LpBasisStatus::kBasic:
        colstat[i] = soplex::SPxSolver::BASIC;
        break;
      case LpBasisStatus::kAtLowerBound:
        colstat[i] = soplex::SPxSolver::ON_LOWER;
        break;
      case LpBasisStatus::kAtUpperBound:
        colstat[i] = soplex::SPxSolver::ON_UPPER;
        break;
      case LpBasisStatus::kFixed:
        colstat[i] = soplex::SPxSolver::FIXED;
        break;
      case LpBasisStatus::kFree:
        colstat[i] = soplex::SPxSolver::ZERO;
        break;
      default:
        LOG(ERROR) << "invalid basis status.";
        std::abort();
        return {absl::StatusCode::kInvalidArgument, "Invalid Data"};
    }
  }

  SOPLEX_TRY(spx_->setBasis(rowstat.get_ptr(), colstat.get_ptr()));
  FreePreStrongBranchingBasis();
  return absl::OkStatus();
}

// Returns the indices of the basic columns and rows; basic column n gives
// value n, basic row m gives value -1-m.
std::vector<ColOrRowIndex> LpSoplexInterface::GetColumnsAndRowsInBasis() const {
  std::vector<int> basis_indices(spx_->numRows());
  std::vector<ColOrRowIndex> basis;
  basis.reserve(spx_->numRows());
  VLOG(10) << "calling GetColumnsAndRowsInBasis().";

  CHECK(PreStrongBranchingBasisFreed());

  spx_->getBasisInd(basis_indices.data());

  for (int index : basis_indices) {
    basis.push_back(index < 0 ? ColOrRowIndex(RowIndex(-1 * index - 1))
                              : ColOrRowIndex(ColIndex(index)));
    VLOG(3) << "Basis index: " << basis.back() << " from (row)Index: " << index;
  }
  return basis;
}

// ==========================================================================
// Getters of vectors in the inverted basis matrix.
// ==========================================================================

// NOTE: The LP interface defines slack variables to have coefficient +1. This
// means that if, internally, the LP solver
//       uses a -1 coefficient, then rows associated with slacks variables
//       whose coefficient is -1, should be negated; see also the explanation
//       in lpi.h.
absl::StatusOr<SparseRow> LpSoplexInterface::GetSparseRowOfBInverted(
    RowIndex row_in_basis) const {
  VLOG(10) << "calling GetSparseRowOfBInverted().";
  CHECK(PreStrongBranchingBasisFreed());
  DCHECK_GE(row_in_basis, 0);
  DCHECK_LT(row_in_basis, spx_->numRowsReal());

  int num_indices;

  std::vector<double> dense_row(spx_->numRowsReal());
  std::vector<int> indices(spx_->numRowsReal());

  if (!spx_->getBasisInverseRowReal(row_in_basis.value(), dense_row.data(),
                                    indices.data(), &num_indices)) {
    return absl::Status(absl::StatusCode::kInternal, "Error");
  }

  absl::StrongVector<RowIndex, double> row_activities =
      GetRowActivities().value();

  LOG(INFO) << "Row activities: ";
  for (RowIndex row_index(0); row_index < row_activities.size(); ++row_index) {
    LOG(INFO) << "Row " << row_index << ": " << row_activities[row_index];
  }

  // TODO (CG): Check if this should be applied automatically.
  int slack = row_activities.at(row_in_basis) < 0 ? -1 : 1;
  if (slack == -1) {
    LOG(INFO) << "SLACK IS NEGATIVE";
  }

  SparseRow sparse_row;
  if (num_indices == -1) {
    // This means that Soplex used a dense representation. If so, the result
    // is found directly in dense_row and must be sparsified.
    for (int i = 0; i < spx_->numRowsReal(); ++i) {
      if (dense_row[i] != 0.0) sparse_row.AddEntry(ColIndex(i), dense_row[i]);
    }
  } else {
    for (int i = 0; i < num_indices; ++i) {
      sparse_row.AddEntry(ColIndex(indices[i]), dense_row[indices[i]]);
    }
  }
  sparse_row.CleanUpIfNeeded();
  LOG(INFO) << "Sparse row: " << sparse_row;
  WriteLpToFile("lp_s.txt");
  return sparse_row;
}

// NOTE: The LP interface defines slack variables to have coefficient +1. This
// means that if, internally, the LP solver
//       uses a -1 coefficient, then rows associated with slacks variables
//       whose coefficient is -1, should be negated; see also the explanation
//       in lpi.h.
//
// Use the column number of B^-1 - this is NOT the number of the column in the
// LP; you have to call minimip::LpInterface.GetBasisIndices() to get the array
// which links the B^-1 column numbers to the row and column numbers of the LP!
// c must be between 0 and num_rows-1, since the basis has the size num_rows *
// num_rows.

absl::StatusOr<SparseCol> LpSoplexInterface::GetSparseColumnOfBInverted(
    ColIndex col_in_basis) const {
  VLOG(10) << "calling GetSparseColumnOfBInverted().";

  // TODO(issues/26): Use getBasisInverseColReal when the SoPlex bug is fixed.
  LOG_FIRST_N(WARNING, 1)
      << "Due to a SoPlex bug, LpSoplexInterface::GetSparseColumnOfBInverted "
         "currently retrieves the matrix by rows, which is inefficient. Prefer "
         "LpSoplexInterface::GetSparseRowOfBInverted if possible.";
  CHECK(PreStrongBranchingBasisFreed());
  DCHECK_GE(col_in_basis, 0);
  DCHECK_LT(col_in_basis, spx_->numRowsReal());

  // Reused in the loop below to avoid reallocation.
  std::vector<double> dense_row(spx_->numRowsReal());
  std::vector<int> indices(spx_->numRowsReal());

  SparseCol sparse_column;
  for (RowIndex row(0); row < spx_->numRowsReal(); ++row) {
    int num_indices;

    if (!spx_->getBasisInverseRowReal(row.value(), dense_row.data(),
                                      indices.data(), &num_indices)) {
      return absl::Status(absl::StatusCode::kInternal, "Error");
    }
    if (num_indices == -1) {
      sparse_column.AddEntry(row, dense_row[col_in_basis.value()]);
    } else {
      // The indices aren't guaranteed to be sorted, so we use linear lookup
      auto it = std::find(indices.begin(), indices.end(), col_in_basis.value());
      if (it != indices.end()) {
        int index = std::distance(indices.begin(), it);
        if (index < num_indices) sparse_column.AddEntry(row, dense_row[index]);
      }
    }
  }
  sparse_column.CleanUpIfNeeded();

  return sparse_column;
}

// NOTE: The LP interface defines slack variables to have coefficient +1. This
// means that if, internally, the LP solver
//       uses a -1 coefficient, then rows associated with slacks variables
//       whose coefficient is -1, should be negated; see also the explanation
//       in lpi.h.
absl::StatusOr<SparseRow> LpSoplexInterface::GetSparseRowOfBInvertedTimesA(
    RowIndex row_in_basis) const {
  VLOG(10) << "calling GetSparseRowOfBInvertedTimesA().";

  int num_rows = spx_->numRowsReal();
  int num_cols = spx_->numColsReal();
  SparseRow sparse_row;

  std::vector<double> dense_binv(num_rows);
  std::vector<int> nonzero_indices(num_rows);
  int num_indices;

  CHECK(PreStrongBranchingBasisFreed());

  // Calculate the row in B^-1.
  if (!spx_->getBasisInverseRowReal(row_in_basis.value(), dense_binv.data(),
                                    nonzero_indices.data(), &num_indices)) {
    return absl::Status(absl::StatusCode::kInternal, "Error");
  }
  soplex::Vector binv_vec(num_rows, dense_binv.data());

  // Temporary unscaled column of A.
  soplex::DSVector acol;
  const double eps = feasibility_tolerance();

  for (int col_idx = 0; col_idx < num_cols; ++col_idx) {
    spx_->getColVectorReal(col_idx, acol);
    double val = binv_vec * acol;
    if (std::abs(val) > eps) {
      sparse_row.AddEntry(ColIndex(col_idx), val);
    }
  }
  return sparse_row;
}

// NOTE: The LP interface defines slack variables to have coefficient +1. This
// means that if, internally, the LP solver
//       uses a -1 coefficient, then rows associated with slacks variables
//       whose coefficient is -1, should be negated; see also the explanation
//       in lpi.h.
absl::StatusOr<SparseCol> LpSoplexInterface::GetSparseColumnOfBInvertedTimesA(
    ColIndex col_in_basis) const {
  VLOG(10) << "calling GetSparseColumnOfBInvertedTimesA().";

  std::vector<double> binv_vec(spx_->numRows());
  SparseCol sparse_col;
  int num_rows = spx_->numRowsReal();

  soplex::DVector column(num_rows);

  // Temporary sparse vector used for unscaling (memory is automatically
  // enlarged).
  soplex::DSVector column_sparse;

  CHECK(PreStrongBranchingBasisFreed());

  DCHECK_GE(col_in_basis, 0);
  DCHECK_LT(col_in_basis, spx_->numColsReal());

  // Col needs to be cleared because copying colVectorReal only regard nonzeros.
  column.clear();

  spx_->getColVectorReal(col_in_basis.value(), column_sparse);
  // The copy is necessary to transform the sparse column into a dense vector.
  column = column_sparse;

  if (!spx_->getBasisInverseTimesVecReal(column.get_ptr(), binv_vec.data())) {
    return absl::Status(absl::StatusCode::kInternal, "Error");
  }
  const double eps = feasibility_tolerance();
  for (int row_idx = 0; row_idx < num_rows; ++row_idx) {
    if (std::abs(binv_vec[row_idx]) > eps) {
      sparse_col.AddEntry(RowIndex(row_idx), binv_vec[row_idx]);
    }
  }
  return sparse_col;
}

// ==========================================================================
// Getters and setters of the parameters.
// ==========================================================================

LpParameters LpSoplexInterface::GetLpParameters() const {
  VLOG(10) << "calling GetLpParameters().";
  LpParameters params;

  params.set_lp_solver_type(LpParameters::LP_SOPLEX);

  params.set_solve_from_scratch(solve_from_scratch_);

  switch (spx_->intParam(soplex::SoPlex::SCALER)) {
    case soplex::SoPlex::SCALER_OFF:
      params.set_scaling_strategy(LpParameters::SCALING_OFF);
      break;
    case soplex::SoPlex::SCALER_BIEQUI:
      params.set_scaling_strategy(LpParameters::SCALING_EQUILIBRATION);
      break;
    case soplex::SoPlex::SCALER_LEASTSQ:
      params.set_scaling_strategy(LpParameters::SCALING_LEAST_SQUARES);
      break;
    default:
      LOG(DFATAL) << "[BUG] Unrecognized scaling strategy.";
  }

  params.set_use_presolve(spx_->intParam(soplex::SoPlex::SIMPLIFIER) !=
                          soplex::SoPlex::SIMPLIFIER_OFF);

  switch (spx_->intParam(soplex::SoPlex::PRICER)) {
    case soplex::SoPlex::PRICER_AUTO:
      params.set_pricing_strategy(LpParameters::PRICING_DEFAULT);
      break;
    case soplex::SoPlex::PRICER_STEEP:
      params.set_pricing_strategy(LpParameters::PRICING_STEEPEST_EDGE);
      break;
    case soplex::SoPlex::PRICER_QUICKSTEEP:
      params.set_pricing_strategy(
          LpParameters::PRICING_STEEPEST_EDGE_QUICK_START);
      break;
    case soplex::SoPlex::PRICER_DANTZIG:
      params.set_pricing_strategy(LpParameters::PRICING_DANTZIG);
      break;
    case soplex::SoPlex::PRICER_PARMULT:
      params.set_pricing_strategy(LpParameters::PRICING_PARTIAL_DANTZIG);
      break;
    case soplex::SoPlex::PRICER_DEVEX:
      params.set_pricing_strategy(LpParameters::PRICING_DEVEX);
      break;
    default:
      LOG(DFATAL) << "[BUG] Unrecognized pricing strategy.";
  }

  params.set_feasibility_tolerance(spx_->realParam(soplex::SoPlex::FEASTOL));
  params.set_optimality_tolerance(spx_->realParam(soplex::SoPlex::OPTTOL));
  params.set_min_markowitz_threshold(
      spx_->realParam(soplex::SoPlex::MIN_MARKOWITZ));
  params.set_refactorization_interval(
      spx_->intParam(soplex::SoPlex::FACTOR_UPDATE_MAX));
  params.set_objective_lower_limit(
      spx_->realParam(soplex::SoPlex::OBJLIMIT_LOWER));
  params.set_objective_upper_limit(
      spx_->realParam(soplex::SoPlex::OBJLIMIT_UPPER));

  switch (spx_->intParam(soplex::SoPlex::TIMER)) {
    case soplex::SoPlex::TIMER_OFF:
      params.set_timing_mode(LpParameters::TIMING_OFF);
      break;
    case soplex::SoPlex::TIMER_WALLCLOCK:
      params.set_timing_mode(LpParameters::TIMING_WALLCLOCK);
      break;
    case soplex::SoPlex::TIMER_CPU:
      params.set_timing_mode(LpParameters::TIMING_CPU);
      break;
    default:
      LOG(DFATAL) << "[BUG] Unrecogniczed timing mode.";
  }

  if (params.timing_mode() != LpParameters::TIMING_OFF) {
    params.set_time_limit(spx_->realParam(soplex::SoPlex::TIMELIMIT));
  }

  const int iteration_limit = spx_->intParam(soplex::SoPlex::ITERLIMIT);
  params.set_iteration_limit(iteration_limit == -1 ? 0 : iteration_limit);

  params.set_random_seed(spx_->randomSeed());

  params.set_enable_internal_solver_output(spx_->spxout.getVerbosity() ==
                                           soplex::SPxOut::DEBUG);

  return params;
}

namespace {

absl::Status LpParametersAreSupportedBySoPlex(const LpParameters& params) {
  VLOG(10) << "calling LpParametersAreSupportedBySoPlex().";
  RETURN_IF_ERROR(LpParametersAreValid(params));

  if (params.lp_solver_type() != LpParameters::LP_SOPLEX) {
    return absl::InvalidArgumentError("Solver type must be LP_SOPLEX.");
  }
  if (params.scaling_strategy() == LpParameters::SCALING_LINEAR_PROGRAM) {
    return absl::InvalidArgumentError("Unsupported scaling strategy.");
  }
  if (params.min_markowitz_threshold() != -1) {
    if (params.min_markowitz_threshold() < 1e-4 or
        params.min_markowitz_threshold() > 1.0 - 1e-4) {
      return absl::InvalidArgumentError(
          "Markowitz threshold must be between 1e-4 and 1-1e-4 (or -1).");
    }
  }
  if (params.timing_mode() == LpParameters::TIMING_DETERMINISTIC) {
    return absl::InvalidArgumentError(
        "SoPlex doesn't support deterministic timing.");
  }
  if (params.num_threads() != 0 and params.num_threads() != 1) {
    return absl::InvalidArgumentError(
        "SoPlex doesn't support using more than one thread.");
  }
  return absl::OkStatus();
}

}  // namespace

absl::Status LpSoplexInterface::SetLpParameters(const LpParameters& params) {
  VLOG(10) << "calling SetLpParameters().";
  RETURN_IF_ERROR(LpParametersAreSupportedBySoPlex(params));

  const LpParameters old_params = GetLpParameters();
  std::string error;

  int scaling_strategy = soplex::SoPlex::SCALER_OFF;
  switch (params.scaling_strategy()) {
    case LpParameters::SCALING_OFF:
      scaling_strategy = soplex::SoPlex::SCALER_OFF;
      break;
    case LpParameters::SCALING_DEFAULT:
    case LpParameters::SCALING_EQUILIBRATION:
      scaling_strategy = soplex::SoPlex::SCALER_BIEQUI;
      break;
    case LpParameters::SCALING_LEAST_SQUARES:
      scaling_strategy = soplex::SoPlex::SCALER_LEASTSQ;
      break;
    case LpParameters::SCALING_LINEAR_PROGRAM:
    default:
      error = "[BUG] Unexpected scaling strategy.";
      LOG(DFATAL) << error;
  }
  if (!spx_->setIntParam(soplex::SoPlex::SCALER, scaling_strategy)) {
    error = "SoPlex error when setting scaling strategy.";
  }

  if (!spx_->setIntParam(soplex::SoPlex::SIMPLIFIER,
                         params.use_presolve()
                             ? soplex::SoPlex::SIMPLIFIER_AUTO
                             : soplex::SoPlex::SIMPLIFIER_OFF)) {
    error = "SoPlex error when setting presolve.";
  }

  int pricing_strategy = soplex::SoPlex::PRICER_AUTO;
  switch (params.pricing_strategy()) {
    case LpParameters::PRICING_DEFAULT:
      pricing_strategy = soplex::SoPlex::PRICER_AUTO;
      break;
    case LpParameters::PRICING_STEEPEST_EDGE:
      pricing_strategy = soplex::SoPlex::PRICER_STEEP;
      break;
    case LpParameters::PRICING_STEEPEST_EDGE_QUICK_START:
      pricing_strategy = soplex::SoPlex::PRICER_QUICKSTEEP;
      break;
    case LpParameters::PRICING_DANTZIG:
      pricing_strategy = soplex::SoPlex::PRICER_DANTZIG;
      break;
    case LpParameters::PRICING_PARTIAL_DANTZIG:
      pricing_strategy = soplex::SoPlex::PRICER_PARMULT;
      break;
    case LpParameters::PRICING_DEVEX:
      pricing_strategy = soplex::SoPlex::PRICER_DEVEX;
      break;
    default:
      error = "[BUG] Unexpected pricing strategy.";
      LOG(DFATAL) << error;
  }
  if (!spx_->setIntParam(soplex::SoPlex::PRICER, pricing_strategy)) {
    error = "SoPlex error when setting pricing strategy.";
  }

  if (!spx_->setRealParam(soplex::SoPlex::FEASTOL,
                          params.feasibility_tolerance() == -1
                              ? 1e-6
                              : params.feasibility_tolerance())) {
    error = "SoPlex error when setting feasibility tolerance.";
  }

  if (!spx_->setRealParam(soplex::SoPlex::OPTTOL,
                          params.optimality_tolerance() == -1
                              ? 1e-6
                              : params.optimality_tolerance())) {
    error = "SoPlex error when setting optimality tolerance.";
  }

  if (!spx_->setRealParam(soplex::SoPlex::MIN_MARKOWITZ,
                          params.min_markowitz_threshold() == -1
                              ? 1e-4
                              : params.min_markowitz_threshold())) {
    error = "SoPlex error when setting min markowitz threshold.";
  }

  if (!spx_->setIntParam(soplex::SoPlex::FACTOR_UPDATE_MAX,
                         params.refactorization_interval())) {
    error = "SoPlex error when setting refactorization interval.";
  }

  if (!spx_->setRealParam(
          soplex::SoPlex::OBJLIMIT_LOWER,
          std::max(-Infinity(), params.objective_lower_limit()))) {
    error = "SoPlex error when setting lower objective limit.";
  }

  if (!spx_->setRealParam(
          soplex::SoPlex::OBJLIMIT_UPPER,
          std::min(Infinity(), params.objective_upper_limit()))) {
    error = "SoPlex error when setting upper objective limit.";
  }

  int timing_mode = soplex::SoPlex::TIMER_OFF;
  switch (params.timing_mode()) {
    case LpParameters::TIMING_OFF:
      timing_mode = soplex::SoPlex::TIMER_OFF;
      break;
    case LpParameters::TIMING_WALLCLOCK:
      timing_mode = soplex::SoPlex::TIMER_WALLCLOCK;
      break;
    case LpParameters::TIMING_CPU:
      timing_mode = soplex::SoPlex::TIMER_CPU;
      break;
    case LpParameters::TIMING_DETERMINISTIC:
    default:
      error = "[BUG] Unexpected timing mode";
      LOG(DFATAL) << error;
  }

  if (!spx_->setIntParam(soplex::SoPlex::TIMER, timing_mode)) {
    error = "SoPlex error when setting timing mode.";
  }

  // SoPlex requires 0 < param_val < DEFAULT_INFINITY (= 1e100). So, if we're
  // out of bounds, we set to Infinity() / 2.0 (still a very large number).
  const double time_limit =
      params.time_limit() == 0.0 or params.time_limit() >= Infinity()
          ? Infinity() / 2.0
          : params.time_limit();
  if (!spx_->setRealParam(soplex::SoPlex::TIMELIMIT, time_limit)) {
    error = "SoPlex error when setting timi limit";
  }

  const int iteration_limit =
      (params.iteration_limit() == 0 or params.iteration_limit() >= INT_MAX)
          ? -1
          : params.iteration_limit();
  if (!spx_->setIntParam(soplex::SoPlex::ITERLIMIT, iteration_limit)) {
    error = "SoPlex error when setting iteration limit.";
  }

  spx_->setRandomSeed(params.random_seed());

  spx_->spxout.setVerbosity(params.enable_internal_solver_output()
                                ? soplex::SPxOut::DEBUG
                                : soplex::SPxOut::ERROR);

  if (error.empty()) return absl::OkStatus();
  CHECK_OK(SetLpParameters(old_params));
  return absl::InvalidArgumentError(error);
}

// ==========================================================================
// Numerical methods.
// ==========================================================================

double LpSoplexInterface::Infinity() const {
  VLOG(10) << "calling Infinity().";

  return spx_->realParam(soplex::SoPlex::INFTY);
}

bool LpSoplexInterface::IsInfinity(double value) const {
  VLOG(10) << "calling IsInfinity().";

  return (value >= spx_->realParam(soplex::SoPlex::INFTY));
}

// ==========================================================================
// File interface methods.
// ==========================================================================

absl::Status LpSoplexInterface::ReadLpFromFile(const std::string& file_path) {
  VLOG(10) << "calling ReadLPFromFile().";

  CHECK(!file_path.empty());

  CHECK(PreStrongBranchingBasisFreed());

  if (!FileExists(file_path)) {
    return {absl::StatusCode::kInternal, "Read Error, file does not exist"};
  }

  DCHECK_EQ(spx_->intParam(soplex::SoPlex::READMODE),
            soplex::SoPlex::READMODE_REAL);

  try {
    if (!spx_->readFile(file_path.c_str())) {
      return {absl::StatusCode::kInternal, "Internal SoPlex Read Error"};
    }
  } catch (const soplex::SPxException& x) {
    LOG(WARNING) << "SoPlex threw an exception: " << x.what() << ".";
    return {absl::StatusCode::kInternal, "Read Error: soplex exception"};
  }

  return absl::OkStatus();
}

absl::StatusOr<std::string> LpSoplexInterface::WriteLpToFile(
    const std::string& file_path) const {
  VLOG(10) << "calling WriteLPToFile().";

  CHECK(!file_path.empty());

  bool success;
  try {
    success = spx_->writeFileReal(file_path.c_str());
  } catch (const soplex::SPxException& x) {
    LOG(WARNING) << "SoPlex threw an exception: " << x.what() << ".";
    return absl::Status(absl::StatusCode::kInternal,
                        "Write Error: soplex exception");
  }
  if (!success) {
    return absl::InternalError("SoPlex failed to write LP to file");
  }
  return file_path;
}

}  // namespace minimip
