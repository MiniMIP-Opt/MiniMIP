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

#include "src/lp_interface/lpi_soplex.h"

const int kSoPlexVerbosityLevel = 0;

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

namespace minimip {

namespace {

bool FileExists(const std::string& file_path) {
  FILE* file = fopen(file_path.c_str(), "r");
  if (file == nullptr) return false;

  fclose(file);

  return true;
}
}  // namespace

LPSoplexInterface::LPSoplexInterface()
    : spx_(new soplex::SoPlex),
      pricing_(LPPricing::kDefault),
      is_solved_(false),
      lp_info_(false),
      from_scratch_(false),
      col_basis_status_(0),
      row_basis_status_(0) {}

LPSoplexInterface::~LPSoplexInterface() = default;

soplex::DataArray<soplex::SPxSolver::VarStatus>&
LPSoplexInterface::RowsBasisStatus() {
  return row_basis_status_;
}

soplex::DataArray<soplex::SPxSolver::VarStatus>&
LPSoplexInterface::ColumnsBasisStatus() {
  return col_basis_status_;
}

void LPSoplexInterface::SetFromScratch(bool from_scratch) {
  from_scratch_ = from_scratch;
}

void LPSoplexInterface::SetLPInfo(bool lp_info) { lp_info_ = lp_info; }

double LPSoplexInterface::objective_limit() const {
  return (spx_->intParam(soplex::SoPlex::OBJSENSE) ==
          soplex::SoPlex::OBJSENSE_MINIMIZE)
             ? spx_->realParam(soplex::SoPlex::OBJLIMIT_UPPER)
             : spx_->realParam(soplex::SoPlex::OBJLIMIT_LOWER);
}

double LPSoplexInterface::feasibility_tolerance() const {
  return spx_->realParam(soplex::SoPlex::FEASTOL);
}

double LPSoplexInterface::optimality_tolerance() const {
  return spx_->realParam(soplex::SoPlex::OPTTOL);
}

bool LPSoplexInterface::CheckConsistentBounds() const {
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

bool LPSoplexInterface::CheckConsistentSides() const {
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

soplex::SPxSolver::Status LPSoplexInterface::LPSolve(
    bool print_warning = true) {
  const soplex::SPxOut::Verbosity verbosity = spx_->spxout.getVerbosity();
  spx_->spxout.setVerbosity(
      (soplex::SPxOut::Verbosity)(GetLPInfo() ? kSoPlexVerbosityLevel : 0));

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
  CHECK(spx_->intParam(soplex::SoPlex::ITERLIMIT) < 0 ||
        spx_->numIterations() <= spx_->intParam(soplex::SoPlex::ITERLIMIT));

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

  // Restore the verbosity after running TrySolve.
  spx_->spxout.setVerbosity(verbosity);

  return spx_status;
}

bool LPSoplexInterface::GetFromScratch() const { return from_scratch_; }

bool LPSoplexInterface::GetLPInfo() const { return lp_info_; }

absl::Status LPSoplexInterface::SoPlexSolve() {
  soplex::SPxOut::Verbosity verbosity = spx_->spxout.getVerbosity();
  spx_->spxout.setVerbosity(
      (soplex::SPxOut::Verbosity)(GetLPInfo() ? kSoPlexVerbosityLevel : 0));

  VLOG(3) << "calling SoPlex solve(): " << spx_->numColsReal() << " cols, "
          << spx_->numRowsReal() << " rows.";

  InvalidateSolution();

  CHECK(PreStrongBranchingBasisFreed());

  // Delete the starting basis if solving from scratch.
  if (GetFromScratch()) {
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
  CHECK(!GetFromScratch() || spx_->status() == soplex::SPxSolver::NO_PROBLEM);

  const soplex::SPxSolver::Status status = LPSolve();

  VLOG(3) << " -> SoPlex status: " << spx_->status()
          << ", basis status: " << spx_->basisStatus() << ".";
  is_solved_ = true;

  // Restore the verbosity after running doSolve.
  spx_->spxout.setVerbosity(verbosity);

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

void LPSoplexInterface::SavePreStrongbranchingBasis() {
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

void LPSoplexInterface::RestorePreStrongbranchingBasis() {
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

void LPSoplexInterface::InvalidateSolution() { is_solved_ = false; }

bool LPSoplexInterface::PreStrongBranchingBasisFreed() const {
  return ((row_basis_status_.size() == 0) && (col_basis_status_.size() == 0));
}

void LPSoplexInterface::FreePreStrongBranchingBasis() {
  row_basis_status_.clear();
  col_basis_status_.clear();
}

// Strongbranching is applied to the given column, with the corresponding
// current primal solution value. The double referenes are used to store the
// dual bound after branching up and down. Additionally the validity of both
// bounds is stored, if one bound is not valid it can be used as an estimate.
absl::Status LPSoplexInterface::StrongBranch(int col, double primal_sol,
                                             int iteration_limit,
                                             StrongBranchResult result) {
  soplex::SPxOut::Verbosity verbosity = spx_->spxout.getVerbosity();
  spx_->spxout.setVerbosity(
      (soplex::SPxOut::Verbosity)(GetLPInfo() ? kSoPlexVerbosityLevel : 0));

  VLOG(3) << "calling StrongBranch() on variable " << col << "("
          << iteration_limit << " iterations).";

  bool error = false;
  int old_iter_limit = spx_->intParam(soplex::SoPlex::ITERLIMIT);

  // Get the current bounds of the column.
  const double old_lb = spx_->lowerReal(col);
  const double old_ub = spx_->upperReal(col);
  result.down_valid = false;
  result.up_valid = false;
  result.iterations = 0;

  // Set the algorithm type to use the dual simplex.
  static_cast<void>(spx_->setIntParam(soplex::SoPlex::ALGORITHM,
                                      soplex::SoPlex::ALGORITHM_DUAL));

  // Computation for the down branch.
  double new_ub = std::ceil(primal_sol - 1.0 - feasibility_tolerance());

  soplex::SPxSolver::Status status;
  bool from_parent_basis;
  if (new_ub >= old_lb - 0.5) {
    VLOG(3) << "strong branching down on x" << col << " (" << primal_sol
            << ") with " << iteration_limit << " iterations.";
    spx_->changeUpperReal(col, new_ub);
    DCHECK_LE(spx_->lowerReal(col), spx_->upperReal(col));

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

    spx_->changeUpperReal(col, old_ub);
    DCHECK_LE(spx_->lowerReal(col), spx_->upperReal(col));
  } else {
    result.dual_bound_down_branch = objective_limit();
    result.down_valid = true;
  }

  // Computation for the up branch.
  if (!error) {
    double new_lb = std::floor(primal_sol + 1.0 - feasibility_tolerance());
    if (new_lb <= old_ub + 0.5) {
      VLOG(3) << "strong branching  up  on x" << col << " (" << primal_sol
              << ") with " << iteration_limit << "iterations.";
      spx_->changeLowerReal(col, new_lb);
      DCHECK_LE(spx_->lowerReal(col), spx_->upperReal(col));

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

      spx_->changeLowerReal(col, old_lb);
      DCHECK_LE(spx_->lowerReal(col), spx_->upperReal(col));
    } else {
      result.dual_bound_up_branch = objective_limit();
      result.up_valid = true;
    }
  }
  // Reset the old iteration limit.
  static_cast<void>(
      spx_->setIntParam(soplex::SoPlex::ITERLIMIT, old_iter_limit));

  // Restore the verbosity after running the branching on the up and down
  // branch.
  spx_->spxout.setVerbosity(verbosity);

  if (error) {
    VLOG(2) << "StrongBranch() returned SoPlex status "
            << static_cast<int>(status) << ".";
    return {absl::StatusCode::kInternal, "Error"};
  }
  return absl::OkStatus();
}

// ==========================================================================
// LP model setters.
// ==========================================================================

absl::Status LPSoplexInterface::PopulateFromMipData(const MipData& mip_data) {
  VLOG(2) << "calling LoadColumnLP().";

  DCHECK_EQ(mip_data.constraint_names().size(),
            mip_data.left_hand_sides().size());
  DCHECK_EQ(mip_data.left_hand_sides().size(),
            mip_data.right_hand_sides().size());
  DCHECK_EQ(mip_data.variable_names().size(), mip_data.lower_bounds().size());
  DCHECK_EQ(mip_data.lower_bounds().size(), mip_data.upper_bounds().size());

  InvalidateSolution();

  CHECK(PreStrongBranchingBasisFreed());

  try {
    int num_rows = mip_data.matrix().num_rows().value();
    soplex::LPRowSet rowset(num_rows);
    soplex::DSVector empty_vector(0);

    spx_->clearLPReal();

    static_cast<void>(spx_->setIntParam(
        soplex::SoPlex::OBJSENSE,
        (mip_data.is_maximization() == 0 ? soplex::SoPlex::OBJSENSE_MINIMIZE
                                         : soplex::SoPlex::OBJSENSE_MAXIMIZE)));

    // Create empty rows with the given sides.
    for (int i = 0; i < num_rows; ++i)
      rowset.add(mip_data.left_hand_sides()[i], empty_vector,
                 mip_data.right_hand_sides()[i]);
    spx_->addRowsReal(rowset);

    // Create the column vectors with the given coefficients and bounds.
    RETURN_IF_ERROR(AddColumns(mip_data.matrix(), mip_data.lower_bounds(),
                               mip_data.upper_bounds(), mip_data.objective(),
                               mip_data.variable_names()));
  } catch (const soplex::SPxException& x) {
    LOG(WARNING) << "SoPlex threw an exception: " << x.what() << ".";
    return {absl::StatusCode::kInternal, "Error"};
  }
  return absl::OkStatus();
}

absl::Status LPSoplexInterface::AddColumn(const SparseCol& col_data,
                                          double lower_bound,
                                          double upper_bound,
                                          double objective_coefficient,
                                          const std::string& name) {
  VLOG(2) << "calling AddColumn().";
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
    spx_->addColReal(soplex::LPCol(objective_coefficient, col_vector,
                                   upper_bound, lower_bound));
  } catch (const soplex::SPxException& x) {
    LOG(WARNING) << "SoPlex threw an exception: " << x.what() << ".";
    return {absl::StatusCode::kInternal, "Error"};
  }

  return absl::OkStatus();
}

absl::Status LPSoplexInterface::AddColumns(
    const StrongSparseMatrix& matrix, const std::vector<double>& lower_bounds,
    const std::vector<double>& upper_bounds,
    const SparseRow& objective_coefficients,
    const std::vector<std::string>& names) {
  VLOG(2) << "calling AddColumns().";

  InvalidateSolution();

  CHECK(PreStrongBranchingBasisFreed());

  if (DEBUG_MODE) {
    if (!(matrix.num_cols() == 0) && matrix.AllColsAreClean()) {
      // Perform a check to ensure that no new rows have been added.
      int num_rows = spx_->numRowsReal();
      for (ColIndex i = ColIndex(0); i < matrix.num_cols(); ++i) {
        for (int j = 0; j < matrix.col(i).entries().size(); ++j) {
          DCHECK_LE(0, matrix.col(i).indices()[j]);
          DCHECK_LT(matrix.col(i).indices()[j], num_rows);
          DCHECK_NE(matrix.col(i).values()[j], 0.0);
        }
      }
    }
  }

  try {
    soplex::LPColSet columns(matrix.num_cols().value());
    soplex::DSVector col_Vector(matrix.num_cols().value());

    // Create column vectors with coefficients and bounds.
    for (ColIndex i = ColIndex(0); i < matrix.num_cols(); ++i) {
      col_Vector.clear();
      if (!matrix.col(i).entries().empty()) {
        for (SparseEntry entry : matrix.col(i).entries()) {
          DCHECK(!IsInfinity(std::abs(entry.value)));
          col_Vector.add(entry.index.value(), entry.value);
        }
      }
      columns.add(objective_coefficients.value(i), lower_bounds[i.value()],
                  col_Vector, upper_bounds[i.value()]);
    }
    spx_->addColsReal(columns);
  } catch (const soplex::SPxException& x) {
    LOG(WARNING) << "SoPlex threw an exception: " << x.what() << ".";
    return {absl::StatusCode::kInternal, "Error"};
  }

  return absl::OkStatus();
}

absl::Status LPSoplexInterface::DeleteColumns(ColIndex first_col,
                                              ColIndex last_col) {
  VLOG(2) << "calling DeleteColumns().";

  DCHECK_LE(0, first_col);
  DCHECK_LE(first_col, last_col);
  DCHECK_LT(last_col, spx_->numColsReal());

  InvalidateSolution();

  CHECK(PreStrongBranchingBasisFreed());

  try {
    spx_->removeColRangeReal(first_col.value(), last_col.value());
  } catch (const soplex::SPxMemoryException& E) {
    LOG(ERROR) << "SoPlex threw a memory exception: " << E.what() << ".";
    return {absl::StatusCode::kInternal, "Error"};
  } catch (const soplex::SPxException& E) {
    LOG(ERROR) << "SoPlex threw an exception: " << E.what() << ".";
    return {absl::StatusCode::kInternal, "Error"};
  }

  return absl::OkStatus();
}

absl::Status LPSoplexInterface::AddRow(const SparseRow& row_data,
                                       double left_hand_side,
                                       double right_hand_side,
                                       const std::string& name) {
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

absl::Status LPSoplexInterface::AddRows(
    const absl::StrongVector<RowIndex, SparseRow>& rows,
    const absl::StrongVector<RowIndex, double>& left_hand_sides,
    const absl::StrongVector<RowIndex, double>& right_hand_sides,
    const absl::StrongVector<RowIndex, std::string>& names) {
  VLOG(2) << "calling AddRows().";

  InvalidateSolution();

  CHECK(PreStrongBranchingBasisFreed());

  if (DEBUG_MODE) {
    for (RowIndex i = RowIndex(0); i < rows.size(); ++i) {
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
    for (RowIndex i = RowIndex(0); i < rows.size(); ++i) {
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

absl::Status LPSoplexInterface::DeleteRows(RowIndex first_row,
                                           RowIndex last_row) {
  VLOG(2) << "calling DeleteRows().";

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
LPSoplexInterface::DeleteRowSet(
    const absl::StrongVector<RowIndex, bool>& rows_to_delete) {
  VLOG(2) << "calling DeleteRowSet().";

  InvalidateSolution();

  CHECK(PreStrongBranchingBasisFreed());

  std::vector<int> int_deletion_status(rows_to_delete.begin(),
                                       rows_to_delete.end());
  int num_rows = spx_->numRowsReal();

  // SoPlex removeRows() method deletes the rows with deletion_status[i] < 0,
  // so we have to negate the values.
  for (int i = 0; i < num_rows; ++i) int_deletion_status[i] *= -1;

  SOPLEX_TRY(spx_->removeRowsReal(int_deletion_status.data()));

  // Note: We cannot use the same vector for the soplex call above and the row
  // mapping, even though RowIndex internally consists of an int. In addition to
  // being bad practice to depend on a class' internals, the compiler may
  // introduce additional padding mening RowIndex and int don't have identical
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
absl::Status LPSoplexInterface::Clear() {
  VLOG(2) << "calling Clear().";

  InvalidateSolution();

  CHECK(PreStrongBranchingBasisFreed());
  SOPLEX_TRY(spx_->clearLPReal());

  return absl::OkStatus();
}

// Clears current LP Interface state (like basis information) of the solver.
absl::Status LPSoplexInterface::ClearState() {
  VLOG(2) << "calling ClearState().";

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

absl::Status LPSoplexInterface::SetColumnBounds(ColIndex col,
                                                double lower_bound,
                                                double upper_bound) {
  DCHECK(!IsInfinity(lower_bound));
  DCHECK(!IsInfinity(-upper_bound));

  VLOG(2) << "calling SetColumnBounds().";

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

absl::Status LPSoplexInterface::SetRowSides(RowIndex row, double left_hand_side,
                                            double right_hand_side) {
  DCHECK_GE(row, 0);
  DCHECK_LT(row, GetNumberOfRows());
  DCHECK(!IsInfinity(left_hand_side));
  DCHECK(!IsInfinity(-right_hand_side));

  VLOG(2) << "calling SetRowSides().";

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

absl::Status LPSoplexInterface::SetObjectiveSense(
    bool is_maximization  // new objective sense
) {
  VLOG(2) << "calling SetObjectiveSense().";

  InvalidateSolution();

  CHECK(PreStrongBranchingBasisFreed());

  SOPLEX_TRY(static_cast<void>(spx_->setIntParam(
      soplex::SoPlex::OBJSENSE, is_maximization == 0
                                    ? soplex::SoPlex::OBJSENSE_MINIMIZE
                                    : soplex::SoPlex::OBJSENSE_MAXIMIZE)));

  return absl::OkStatus();
}

absl::Status LPSoplexInterface::SetObjectiveCoefficient(
    ColIndex col, double objective_coefficient) {
  DCHECK(!IsInfinity(std::abs(objective_coefficient)));
  VLOG(2) << "calling SetObjectiveCoefficient().";

  InvalidateSolution();

  CHECK(PreStrongBranchingBasisFreed());

  spx_->changeObjReal(col.value(), objective_coefficient);
  return absl::OkStatus();
}

// ==========================================================================
// LP model getters.
// ==========================================================================

RowIndex LPSoplexInterface::GetNumberOfRows() const {
  VLOG(2) << "calling GetNumberOfRows().";

  return RowIndex(spx_->numRowsReal());
}

ColIndex LPSoplexInterface::GetNumberOfColumns() const {
  VLOG(2) << "calling GetNumberOfColumns().";

  return ColIndex(spx_->numColsReal());
}

int64_t LPSoplexInterface::GetNumberOfNonZeros() const {
  // SoPlex has no direct method to return the number of nonzeros, so we have
  // to count them manually.
  int num_non_zeros = 0;

  VLOG(2) << "calling GetNumberOfNonZeros().";

  if (spx_->numRowsReal() < spx_->numColsReal()) {
    for (int i = 0; i < spx_->numRowsReal(); ++i)
      num_non_zeros += spx_->rowVectorRealInternal(i).size();
  } else {
    for (int i = 0; i < spx_->numColsReal(); ++i)
      num_non_zeros += spx_->colVectorRealInternal(i).size();
  }

  return num_non_zeros;
}

bool LPSoplexInterface::IsMaximization() const {
  VLOG(2) << "calling GetObjectiveSense().";

  return spx_->intParam(soplex::SoPlex::OBJSENSE) !=
         soplex::SoPlex::OBJSENSE_MINIMIZE;
}

// Either both, lower_bound and upper_bound, have to be 0, or both have to be
// non-0, either n_non_zeroes, begin_cols, indices, and obj_coeffs have to be
// 0, or all of them have to be non-0.
SparseCol LPSoplexInterface::GetSparseColumnCoefficients(ColIndex col) const {
  VLOG(2) << "calling GetSparseColumnCoefficients().";

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
SparseRow LPSoplexInterface::GetSparseRowCoefficients(RowIndex row) const {
  VLOG(2) << "calling GetSparseRowCoefficients().";

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

double LPSoplexInterface::GetObjectiveCoefficient(ColIndex col) const {
  VLOG(2) << "calling GetObjectiveCoefficient().";

  DCHECK_LE(0, col);
  DCHECK_LT(col, spx_->numColsReal());

  return spx_->objReal(col.value());
}

double LPSoplexInterface::GetLowerBound(ColIndex col) const {
  VLOG(2) << "calling GetLowerBound().";

  DCHECK_LE(0, col);
  DCHECK_LT(col, spx_->numColsReal());

  return spx_->lowerReal(col.value());
}

double LPSoplexInterface::GetUpperBound(ColIndex col) const {
  VLOG(2) << "calling GetUpperBound().";

  DCHECK_LE(0, col);
  DCHECK_LT(col, spx_->numColsReal());

  return spx_->upperReal(col.value());
}

double LPSoplexInterface::GetLeftHandSide(RowIndex row) const {
  VLOG(2) << "calling GetLeftHandSide().";

  DCHECK_LE(0, row);
  DCHECK_LT(row, spx_->numRowsReal());

  return spx_->lhsReal(row.value());
}

double LPSoplexInterface::GetRightHandSide(RowIndex row) const {
  VLOG(2) << "calling GetRightHandSide().";

  DCHECK_LE(0, row);
  DCHECK_LT(row, spx_->numRowsReal());

  return spx_->rhsReal(row.value());
}

double LPSoplexInterface::GetMatrixCoefficient(ColIndex col,
                                               RowIndex row) const {
  VLOG(2) << "calling GetMatrixCoefficient().";

  DCHECK_LE(0, col);
  DCHECK_LT(col, spx_->numColsReal());
  DCHECK_LE(0, row);
  DCHECK_LT(row, spx_->numRowsReal());

  return spx_->coefReal(row.value(), col.value());
}

// ==========================================================================
// Solving methods.
// ==========================================================================

absl::Status LPSoplexInterface::SolveLPWithPrimalSimplex() {
  VLOG(2) << "calling SolveLPWithPrimalSimplex().";

  static_cast<void>(spx_->setIntParam(soplex::SoPlex::ALGORITHM,
                                      soplex::SoPlex::ALGORITHM_PRIMAL));
  return SoPlexSolve();
}

absl::Status LPSoplexInterface::SolveLPWithDualSimplex() {
  VLOG(2) << "calling SolveLPWithDualSimplex().";

  static_cast<void>(spx_->setIntParam(soplex::SoPlex::ALGORITHM,
                                      soplex::SoPlex::ALGORITHM_DUAL));
  return SoPlexSolve();
}

// This call is needed before any strong branching.
absl::Status LPSoplexInterface::StartStrongBranching() {
  CHECK(PreStrongBranchingBasisFreed());
  SavePreStrongbranchingBasis();

  return absl::OkStatus();
}

// This call is needed after any strong branching.
absl::Status LPSoplexInterface::EndStrongBranching() {
  CHECK(!PreStrongBranchingBasisFreed());
  RestorePreStrongbranchingBasis();
  FreePreStrongBranchingBasis();

  return absl::OkStatus();
}

// Performs strong branching iterations on one branching candidate.
absl::StatusOr<LPInterface::StrongBranchResult>
LPSoplexInterface::SolveDownAndUpStrongBranch(ColIndex col, double primal_value,
                                              int iteration_limit) {
  StrongBranchResult result;
  absl::Status absl_status_code =
      StrongBranch(col.value(), primal_value, iteration_limit, result);

  if (absl_status_code != absl::OkStatus()) {
    return absl::Status(absl::StatusCode::kInternal, "Error");
  } else {
    return result;
  }
}

// ============================================================================
// Solution information getters.
// ============================================================================

// Returns whether a solve method was called after the last modification of
// the LP.
bool LPSoplexInterface::IsSolved() const { return is_solved_; }

// Returns true if the LP is proven to have a primal unbounded ray (but not
// necessary a primal feasible point); this does not necessarily mean that the
// solver knows and can return the primal ray.
bool LPSoplexInterface::ExistsPrimalRay() const {
  VLOG(2) << "calling ExistsPrimalRay().";

  return (spx_->status() == soplex::SPxSolver::UNBOUNDED);
}

// Returns true iff LP is proven to have a primal unbounded ray (but not
// necessary a primal feasible point), and the solver knows and can return the
// primal ray.
bool LPSoplexInterface::HasPrimalRay() const {
  VLOG(2) << "calling HasPrimalRay().";

  return spx_->hasPrimalRay();
}

bool LPSoplexInterface::IsPrimalUnbounded() const {
  VLOG(2) << "calling IsPrimalUnbounded().";

  // If SoPlex returns unbounded, this may only mean that an unbounded ray is
  // available, not necessarily a primal
  //* feasible point; hence we have to check the perturbation.
  return spx_->status() == soplex::SPxSolver::UNBOUNDED;
}

bool LPSoplexInterface::IsPrimalInfeasible() const {
  VLOG(2) << "calling IsPrimalInfeasible().";

  return (spx_->status() == soplex::SPxSolver::INFEASIBLE);
}

bool LPSoplexInterface::IsPrimalFeasible() const {
  VLOG(2) << "calling IsPrimalFeasible().";

  return spx_->basisStatus() == soplex::SPxBasis::OPTIMAL ||
         spx_->basisStatus() == soplex::SPxBasis::PRIMAL;
}

// Returns true if the LP is proven to have a dual unbounded ray (but not
// necessary a dual feasible point);
// this does not necessarily mean that the solver knows and can return the
// dual ray.
bool LPSoplexInterface::ExistsDualRay() const {
  VLOG(2) << "calling ExistsDualRay().";

  return (spx_->status() == soplex::SPxSolver::INFEASIBLE);
}

// Returns true if the LP is proven to have a dual unbounded ray (but not
// necessary a dual feasible point),
//*  and the solver knows and can return the dual ray
bool LPSoplexInterface::HasDualRay() const {
  VLOG(2) << "calling HasDualRay().";

  return spx_->hasDualFarkas();
}

bool LPSoplexInterface::IsDualUnbounded() const {
  VLOG(2) << "calling IsDualUnbounded().";

  return spx_->status() == soplex::SPxSolver::INFEASIBLE &&
         spx_->basisStatus() == soplex::SPxBasis::DUAL;
}

bool LPSoplexInterface::IsDualInfeasible() const {
  VLOG(2) << "calling IsDualInfeasible().";

  return (spx_->status() == soplex::SPxSolver::UNBOUNDED);
}

bool LPSoplexInterface::IsDualFeasible() const {
  VLOG(2) << "calling IsDualFeasible().";

  return (spx_->basisStatus() == soplex::SPxBasis::OPTIMAL) ||
         spx_->basisStatus() == soplex::SPxBasis::DUAL;
}

bool LPSoplexInterface::IsOptimal() const {
  VLOG(2) << "calling IsOptimal().";

  CHECK((spx_->basisStatus() == soplex::SPxBasis::OPTIMAL) ==
        (IsPrimalFeasible() && IsDualFeasible()));

  return (spx_->status() == soplex::SPxSolver::OPTIMAL);
}

// This function should return true if the solution is reliable, i.e., feasible
// and optimal (or proven infeasible/unbounded) with respect to the original
// problem. The optimality status might be with respect to a scaled version of
// the problem, but the solution might not be feasible to the unscaled original
// problem; in this case, minimip::LPInterface.IsStable() should return false.
bool LPSoplexInterface::IsStable() const {
  VLOG(2) << "calling IsStable().";

  if (spx_->status() == soplex::SPxSolver::ERROR ||
      spx_->status() == soplex::SPxSolver::SINGULAR) {
    return false;
  }
  if (spx_->status() == soplex::SPxSolver::OPTIMAL_UNSCALED_VIOLATIONS) {
    return false;
  }
  return true;
}

bool LPSoplexInterface::ObjectiveLimitIsExceeded() const {
  VLOG(2) << "calling ObjectiveLimitIsExceeded().";

  return (spx_->status() == soplex::SPxSolver::ABORT_VALUE);
}

bool LPSoplexInterface::IterationLimitIsExceeded() const {
  VLOG(2) << "calling IterationLimitIsExceeded().";

  return (spx_->status() == soplex::SPxSolver::ABORT_ITER);
}

bool LPSoplexInterface::TimeLimitIsExceeded() const {
  VLOG(2) << "calling TimeLimitIsExceeded().";

  return (spx_->status() == soplex::SPxSolver::ABORT_TIME);
}

double LPSoplexInterface::GetObjectiveValue() {
  VLOG(2) << "calling GetObjectiveValue().";

  return spx_->objValueReal();
}

// Before calling this function, the caller must ensure that the LP has been
// solved to optimality, i.e., that minimip::LPInterface.IsOptimal() returns
// true.

absl::StatusOr<absl::StrongVector<ColIndex, double>>
LPSoplexInterface::GetPrimalValues() const {
  VLOG(2) << "calling GetPrimalSolution().";
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
// solved to optimality, i.e., that minimip::LPInterface.IsOptimal() returns
// true.
absl::StatusOr<absl::StrongVector<RowIndex, double>>
LPSoplexInterface::GetDualValues() const {
  VLOG(2) << "calling GetDualSolution().";
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
LPSoplexInterface::GetRowActivities() const {
  VLOG(2) << "calling GetRowActivity().";
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
LPSoplexInterface::GetReducedCosts() const {
  VLOG(2) << "calling GetReducedCost().";
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
LPSoplexInterface::GetPrimalRay() const {
  VLOG(2) << "calling GetPrimalRay().";
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
LPSoplexInterface::GetDualRay() const {
  VLOG(2) << "calling GetDualFarkasMultiplier().";
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
int64_t LPSoplexInterface::GetNumIterations() const {
  VLOG(2) << "calling GetNumIterations().";

  return spx_->numIterations();
}

// ==========================================================================
// Getters and setters of the basis.
// ==========================================================================

absl::StatusOr<absl::StrongVector<ColIndex, LPBasisStatus>>
LPSoplexInterface::GetBasisStatusForColumns() const {
  absl::StrongVector<ColIndex, LPBasisStatus> statuses(spx_->numColsReal());
  VLOG(2) << "calling GetColumnBasisStatus().";
  CHECK(PreStrongBranchingBasisFreed());
  CHECK(IsOptimal());

  if (!statuses.empty()) {
    for (ColIndex i = ColIndex(0); i < spx_->numColsReal(); ++i) {
      //         double obj_coeffs = 0.0;
      switch (spx_->basisColStatus(i.value())) {
        case soplex::SPxSolver::BASIC:
          statuses[i] = LPBasisStatus::kBasic;
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
          statuses[i] = LPBasisStatus::kAtLowerBound;
          break;
        case soplex::SPxSolver::ON_LOWER:
          statuses[i] = LPBasisStatus::kAtLowerBound;
          break;
        case soplex::SPxSolver::ON_UPPER:
          statuses[i] = LPBasisStatus::kAtUpperBound;
          break;
        case soplex::SPxSolver::ZERO:
          statuses[i] = LPBasisStatus::kFree;
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

absl::StatusOr<absl::StrongVector<RowIndex, LPBasisStatus>>
LPSoplexInterface::GetBasisStatusForRows() const {
  absl::StrongVector<RowIndex, LPBasisStatus> statuses(spx_->numRowsReal());
  VLOG(2) << "calling GetRowBasisStatus().";
  CHECK(PreStrongBranchingBasisFreed());
  CHECK(IsOptimal());

  if (!statuses.empty()) {
    for (RowIndex i = RowIndex(0); i < spx_->numRowsReal(); ++i) {
      switch (spx_->basisRowStatus(i.value())) {
        case soplex::SPxSolver::BASIC:
          statuses[i] = LPBasisStatus::kBasic;
          break;
        case soplex::SPxSolver::FIXED:
          statuses[i] = LPBasisStatus::kFixed;
        case soplex::SPxSolver::ON_LOWER:
          statuses[i] = LPBasisStatus::kAtLowerBound;
          break;
        case soplex::SPxSolver::ON_UPPER:
          statuses[i] = LPBasisStatus::kAtUpperBound;
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

absl::Status LPSoplexInterface::SetBasisStatusForColumnsAndRows(
    const absl::StrongVector<ColIndex, LPBasisStatus>& column_basis_statuses,
    const absl::StrongVector<RowIndex, LPBasisStatus>& row_basis_statuses) {
  VLOG(2) << "calling SetBasisStatus().";

  ColIndex num_cols = GetNumberOfColumns();
  RowIndex num_rows = GetNumberOfRows();

  CHECK(PreStrongBranchingBasisFreed());
  InvalidateSolution();

  soplex::DataArray<soplex::SPxSolver::VarStatus>& _colstat =
      ColumnsBasisStatus();
  soplex::DataArray<soplex::SPxSolver::VarStatus>& _rowstat = RowsBasisStatus();

  _colstat.reSize(num_cols.value());
  _rowstat.reSize(num_rows.value());

  for (int i = 0; i < num_rows; ++i) {
    switch (row_basis_statuses[RowIndex(i)]) {
      case LPBasisStatus::kBasic:
        _rowstat[i] = soplex::SPxSolver::BASIC;
        break;
      case LPBasisStatus::kAtLowerBound:
        _rowstat[i] = soplex::SPxSolver::ON_LOWER;
        break;
      case LPBasisStatus::kAtUpperBound:
        _rowstat[i] = soplex::SPxSolver::ON_UPPER;
        break;
      case LPBasisStatus::kFixed:
        _rowstat[i] = soplex::SPxSolver::FIXED;
        break;
      case LPBasisStatus::kFree:
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
      case LPBasisStatus::kBasic:
        _colstat[i] = soplex::SPxSolver::BASIC;
        break;
      case LPBasisStatus::kAtLowerBound:
        _colstat[i] = soplex::SPxSolver::ON_LOWER;
        break;
      case LPBasisStatus::kAtUpperBound:
        _colstat[i] = soplex::SPxSolver::ON_UPPER;
        break;
      case LPBasisStatus::kFixed:
        _colstat[i] = soplex::SPxSolver::FIXED;
        break;
      case LPBasisStatus::kFree:
        _colstat[i] = soplex::SPxSolver::ZERO;
        break;
      default:
        LOG(ERROR) << "invalid basis status.";
        std::abort();
        return {absl::StatusCode::kInvalidArgument, "Invalid Data"};
    }
  }

  SOPLEX_TRY(spx_->setBasis(_rowstat.get_ptr(), _colstat.get_ptr()));
  FreePreStrongBranchingBasis();
  return absl::OkStatus();
}

// Returns the indices of the basic columns and rows; basic column n gives
// value n, basic row m gives value -1-m.
std::vector<ColOrRowIndex> LPSoplexInterface::GetColumnsAndRowsInBasis() const {
  std::vector<int> basis_indices(spx_->numRows());
  std::vector<ColOrRowIndex> basis;
  basis.reserve(spx_->numRows());
  VLOG(2) << "calling GetBasisInd().";

  CHECK(PreStrongBranchingBasisFreed());

  spx_->getBasisInd(basis_indices.data());

  for (int index : basis_indices) {
    basis.push_back(index < 0 ? ColOrRowIndex(RowIndex(-1 * index + 1))
                              : ColOrRowIndex(ColIndex(index)));
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
absl::StatusOr<SparseRow> LPSoplexInterface::GetSparseRowOfBInverted(
    RowIndex row_in_basis) const {
  VLOG(2) << "calling GetSparseRowOfBInverted().";
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
  return sparse_row;
}

// NOTE: The LP interface defines slack variables to have coefficient +1. This
// means that if, internally, the LP solver
//       uses a -1 coefficient, then rows associated with slacks variables
//       whose coefficient is -1, should be negated; see also the explanation
//       in lpi.h.
//
// Use the column number of B^-1 - this is NOT the number of the column in the
// LP; you have to call minimip::LPInterface.GetBasisIndices() to get the array
// which links the B^-1 column numbers to the row and column numbers of the LP!
// c must be between 0 and num_rows-1, since the basis has the size num_rows *
// num_rows.

absl::StatusOr<SparseCol> LPSoplexInterface::GetSparseColumnOfBInverted(
    ColIndex col_in_basis) const {
  VLOG(2) << "calling GetColumnOfBInverted().";

  // TODO(issues/26): Use getBasisInverseColReal when the SoPlex bug is fixed.
  LOG_FIRST_N(WARNING, 1)
      << "Due to a SoPlex bug, LPSoplexInterface::GetSparseColumnOfBInverted "
         "currently retrieves the matrix by rows, which is inefficient. Prefer "
         "LPSoplexInterface::GetSparseRowOfBInverted if possible.";
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
absl::StatusOr<SparseRow> LPSoplexInterface::GetSparseRowOfBInvertedTimesA(
    RowIndex row_in_basis) const {
  VLOG(2) << "calling GetSparseRowOfBInvertedTimesA().";

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
absl::StatusOr<SparseCol> LPSoplexInterface::GetSparseColumnOfBInvertedTimesA(
    ColIndex col_in_basis) const {
  VLOG(2) << "calling GetColumnOfBInvertedTimesA().";

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

absl::StatusOr<int> LPSoplexInterface::GetIntegerParameter(
    LPParameter type) const {
  int scale_param;
  int param_val;

  VLOG(2) << "calling GetIntegerParameter().";

  switch (type) {
    case LPParameter::kFromScratch:
      param_val = GetFromScratch();
      break;
    case LPParameter::kLPInfo:
      param_val = GetLPInfo();
      break;
    case LPParameter::kLPIterationLimit:
      param_val = spx_->intParam(soplex::SoPlex::ITERLIMIT);
      if (param_val == -1) {
        param_val = INT_MAX;
      }
      break;
    case LPParameter::kPresolving:
      param_val = spx_->intParam(soplex::SoPlex::SIMPLIFIER) ==
                  soplex::SoPlex::SIMPLIFIER_AUTO;
      break;
    case LPParameter::kPricing:
      param_val = static_cast<int>(pricing_);
      break;
    case LPParameter::kScaling:
      scale_param = spx_->intParam(soplex::SoPlex::SCALER);

      if (scale_param == soplex::SoPlex::SCALER_OFF) {
        param_val = 0;
      } else if (scale_param == soplex::SoPlex::SCALER_BIEQUI) {
        param_val = 1;
      } else {
        DCHECK_EQ(scale_param, soplex::SoPlex::SCALER_LEASTSQ);
        param_val = 2;
      }
      break;
    case LPParameter::kTiming:
      param_val = static_cast<int>(spx_->intParam(soplex::SoPlex::TIMER));
      break;
    case LPParameter::kRandomSeed:
      param_val = static_cast<int>(spx_->randomSeed());
      break;
    case LPParameter::kRefactor:
      param_val =
          static_cast<int>(spx_->intParam(soplex::SoPlex::FACTOR_UPDATE_MAX));
      break;
    default:
      return absl::Status(absl::StatusCode::kInvalidArgument,
                          "Parameter Unknown");
  }

  return param_val;
}

absl::Status LPSoplexInterface::SetIntegerParameter(LPParameter type,
                                                    int param_val) {
  VLOG(2) << "calling SetIntegerParameter().";

  switch (type) {
    case LPParameter::kFromScratch:
      CHECK(param_val == true || param_val == false);
      SetFromScratch(static_cast<bool>(param_val));
      break;
    case LPParameter::kLPInfo:
      CHECK(param_val == true || param_val == false);
      SetLPInfo(static_cast<bool>(param_val));
      break;
    case LPParameter::kLPIterationLimit:
      DCHECK_GE(param_val, 0);
      if (param_val >= INT_MAX) {
        static_cast<void>(spx_->setIntParam(soplex::SoPlex::ITERLIMIT, -1));
      }
      break;
    case LPParameter::kPresolving:
      CHECK(param_val == true || param_val == false);
      static_cast<void>(
          spx_->setIntParam(soplex::SoPlex::SIMPLIFIER,
                            (param_val ? soplex::SoPlex::SIMPLIFIER_AUTO
                                       : soplex::SoPlex::SIMPLIFIER_OFF)));
      break;
    case LPParameter::kPricing:
      pricing_ = (LPPricing)param_val;
      switch (pricing_) {
        case LPPricing::kDefault:
        case LPPricing::kAuto:
          static_cast<void>(spx_->setIntParam(soplex::SoPlex::PRICER,
                                              soplex::SoPlex::PRICER_AUTO));
          break;
        case LPPricing::kFull:
          static_cast<void>(spx_->setIntParam(soplex::SoPlex::PRICER,
                                              soplex::SoPlex::PRICER_STEEP));
          break;
        case LPPricing::kPartial:
          static_cast<void>(spx_->setIntParam(soplex::SoPlex::PRICER,
                                              soplex::SoPlex::PRICER_PARMULT));
          break;
        case LPPricing::kSteep:
          static_cast<void>(spx_->setIntParam(soplex::SoPlex::PRICER,
                                              soplex::SoPlex::PRICER_STEEP));
          break;
        case LPPricing::kSteepQStart:
          static_cast<void>(spx_->setIntParam(
              soplex::SoPlex::PRICER, soplex::SoPlex::PRICER_QUICKSTEEP));
          break;
        case LPPricing::kDevex:
          static_cast<void>(spx_->setIntParam(soplex::SoPlex::PRICER,
                                              soplex::SoPlex::PRICER_DEVEX));
          break;
        default:
          return {absl::StatusCode::kInternal, "Error"};
      }
      break;
    case LPParameter::kScaling:
      CHECK(param_val >= 0 && param_val <= 2);
      if (param_val == 0) {
        static_cast<void>(spx_->setIntParam(soplex::SoPlex::SCALER,
                                            soplex::SoPlex::SCALER_OFF));
      } else if (param_val == 1) {
        static_cast<void>(spx_->setIntParam(soplex::SoPlex::SCALER,
                                            soplex::SoPlex::SCALER_BIEQUI));
      } else {
        static_cast<void>(spx_->setIntParam(soplex::SoPlex::SCALER,
                                            soplex::SoPlex::SCALER_LEASTSQ));
      }
      break;
    case LPParameter::kTiming:
      CHECK(param_val >= 0 && param_val < 3);
      static_cast<void>(spx_->setIntParam(soplex::SoPlex::TIMER, param_val));
      break;
    case LPParameter::kRandomSeed:
      spx_->setRandomSeed(
          static_cast<unsigned long>(static_cast<long>(param_val)));
      break;
    case LPParameter::kRefactor:
      DCHECK_GE(param_val, 0);
      static_cast<void>(
          spx_->setIntParam(soplex::SoPlex::FACTOR_UPDATE_MAX, param_val));
      break;
    default:
      return {absl::StatusCode::kInvalidArgument, "Parameter Unknown"};
  }

  return absl::OkStatus();
}

absl::StatusOr<double> LPSoplexInterface::GetRealParameter(
    LPParameter type) const {
  double param_val;
  VLOG(2) << "calling GetRealParameter().";

  switch (type) {
    case LPParameter::kFeasibilityTolerance:
      param_val = feasibility_tolerance();
      break;
    case LPParameter::kDualFeasibilityTolerance:
      param_val = optimality_tolerance();
      break;
    case LPParameter::kObjectiveLimit:
      if (spx_->intParam(soplex::SoPlex::OBJSENSE) ==
          soplex::SoPlex::OBJSENSE_MINIMIZE) {
        param_val = spx_->realParam(soplex::SoPlex::OBJLIMIT_UPPER);
      } else {
        param_val = spx_->realParam(soplex::SoPlex::OBJLIMIT_LOWER);
      }
      break;
    case LPParameter::kLPTimeLimit:
      param_val = spx_->realParam(soplex::SoPlex::TIMELIMIT);
      break;
    case LPParameter::kMarkowitz:
      param_val = spx_->realParam(soplex::SoPlex::MIN_MARKOWITZ);
      break;
    default:
      return absl::Status(absl::StatusCode::kInvalidArgument,
                          "Parameter Unknown");
  }

  return param_val;
}

absl::Status LPSoplexInterface::SetRealParameter(LPParameter type,
                                                 double param_val) {
  VLOG(2) << "calling SetRealParameter().";

  switch (type) {
    case LPParameter::kFeasibilityTolerance:
      DCHECK_GT(param_val, 0.0);
      // This needs to be a CHECK to trigger the function call.
      CHECK(spx_->setRealParam(soplex::SoPlex::FEASTOL, param_val))
          << "SoPlex: unsupported parameter value.\n";
      ;
      break;
    case LPParameter::kDualFeasibilityTolerance:
      DCHECK_GT(param_val, 0.0);
      // This needs to be a CHECK to trigger the function call.
      CHECK(spx_->setRealParam(soplex::SoPlex::OPTTOL, param_val))
          << "SoPlex: unsupported parameter value.\n";
      ;
      break;
    case LPParameter::kObjectiveLimit:
      // no restrictions on param_val
      if (spx_->intParam(soplex::SoPlex::OBJSENSE) ==
          soplex::SoPlex::OBJSENSE_MINIMIZE) {
        static_cast<void>(
            spx_->setRealParam(soplex::SoPlex::OBJLIMIT_UPPER, param_val));
      } else {
        static_cast<void>(
            spx_->setRealParam(soplex::SoPlex::OBJLIMIT_LOWER, param_val));
      }
      break;
    case LPParameter::kLPTimeLimit:
      DCHECK_GT(param_val, 0.0);
      // SoPleX requires 0 < param_val < DEFAULT_INFINITY (= 1e100).
      static_cast<void>(
          spx_->setRealParam(soplex::SoPlex::TIMELIMIT, param_val));
      break;
    case LPParameter::kMarkowitz:
      if (param_val < 1e-4) {
        param_val = 1e-4;
      } else if (param_val > 0.9999) {
        param_val = 0.9999;
      }
      static_cast<void>(
          spx_->setRealParam(soplex::SoPlex::MIN_MARKOWITZ, param_val));
      break;
    default:
      return {absl::StatusCode::kInvalidArgument, "Parameter Unknown"};
  }

  return absl::OkStatus();
}

// ==========================================================================
// Numerical methods.
// ==========================================================================

double LPSoplexInterface::Infinity() const {
  VLOG(2) << "calling Infinity().";

  return spx_->realParam(soplex::SoPlex::INFTY);
}

bool LPSoplexInterface::IsInfinity(double value) const {
  VLOG(2) << "calling IsInfinity().";

  return (value >= spx_->realParam(soplex::SoPlex::INFTY));
}

// ==========================================================================
// File interface methods.
// ==========================================================================

absl::Status LPSoplexInterface::ReadLPFromFile(const std::string& file_path) {
  VLOG(2) << "calling ReadLPFromFile().";

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

absl::StatusOr<std::string> LPSoplexInterface::WriteLPToFile(
    const std::string& file_path) const {
  VLOG(2) << "calling WriteLPToFile().";

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
