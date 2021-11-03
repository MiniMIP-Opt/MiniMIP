#include "src/lp_interface/lpi_glop.h"

#include <limits>

using operations_research::MPModelProto;
using operations_research::TimeLimit;
using operations_research::glop::BasisState;
using operations_research::glop::ColIndex;
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
using operations_research::glop::RowIndex;
using operations_research::glop::ScatteredColumn;
using operations_research::glop::ScatteredColumnIterator;
using operations_research::glop::ScatteredRow;
using operations_research::glop::ScatteredRowIterator;
using operations_research::glop::SparseColumn;
using operations_research::glop::SparseMatrix;
using operations_research::glop::VariableStatus;
using operations_research::glop::VariableStatusRow;

#define EPS 1e-6

// LP interface
namespace minimip {

LPGlopInterface::LPGlopInterface()
    : lp_modified_since_last_solve_(true),
      lp_time_limit_was_reached_(false),
      lp_info_(false),
      pricing_(LPPricing::kDefault),
      from_scratch_(false),
      num_threads_(0),
      timing_(0),
      niterations_(0LL),
      tmp_row_(new ScatteredRow()),
      tmp_column_(new ScatteredColumn()) {}

LPGlopInterface::~LPGlopInterface() {
  MiniMIPdebugMessage("LPGLopInterface Free\n");
}

// ==========================================================================
// LP model setters.
// ==========================================================================

// copies LP data with column matrix into LP solver
absl::Status LPGlopInterface::LoadSparseColumnLP(
    LPObjectiveSense obj_sense,  // objective sense
    const std::vector<double>&
        objective_values,  // objective function values of columns
    const std::vector<double>& lower_bounds,      // lower bounds of columns
    const std::vector<double>& upper_bounds,      // upper bounds of columns
    std::vector<std::string>& col_names,          // column names
    const std::vector<double>& left_hand_sides,   // left hand sides of rows
    const std::vector<double>& right_hand_sides,  // right hand sides of rows
    std::vector<std::string>& row_names,          // row names
    const std::vector<SparseVector>& cols         // sparse columns
) {
  linear_program_.Clear();
  MINIMIP_CALL(AddRows(std::vector<SparseVector>(), left_hand_sides,
                       right_hand_sides, row_names));
  MINIMIP_CALL(AddColumns(cols, objective_values, lower_bounds, upper_bounds,
                          col_names));
  MINIMIP_CALL(SetObjectiveSense(obj_sense));

  return absl::OkStatus();
}

// To be deleted when finalizing the unittests:
// absl::Status LPGlopInterface::LoadColumnLP(
//    LPObjectiveSense obj_sense,  // objective sense
//    int num_cols,                // number of columns
//    const std::vector<double>&
//        objective_values,  // objective function values of columns
//    const std::vector<double>& lower_bounds,      // lower bounds of columns
//    const std::vector<double>& upper_bounds,      // upper bounds of columns
//    std::vector<std::string>& col_names,          // column names
//    int num_rows,                                 // number of rows
//    const std::vector<double>& left_hand_sides,   // left hand sides of rows
//    const std::vector<double>& right_hand_sides,  // right hand sides of rows
//    std::vector<std::string>& row_names,          // row names
//    int num_non_zeros,  // number of non-zero elements in the constraint
//                        // matrix
//    const std::vector<int>& begin_cols,  // start index of each column in
//                                         // row_indices- and vals-array
//    const std::vector<int>&
//        row_indices,                 // row indices of constraint matrix
//        entries
//    const std::vector<double>& vals  // values of constraint matrix entries
//) {
//  linear_program_.Clear();
//  MINIMIP_CALL(AddRows(num_rows, left_hand_sides, right_hand_sides, row_names,
//                       0, std::vector<int>(), std::vector<int>(),
//                       std::vector<double>()));
//  //  MINIMIP_CALL(AddColumns(num_cols, objective_values, lower_bounds,
//  //                          upper_bounds, col_names, num_non_zeros,
//  //                          begin_cols, row_indices, vals));
//  MINIMIP_CALL(SetObjectiveSense(obj_sense));
//
//  return absl::OkStatus();
//}

// adds column to the LP
absl::Status LPGlopInterface::AddColumn(
    const SparseVector& col,  // column to be added
    double lower_bound,       // lower bound of new column
    double upper_bound,       // upper bound of new column
    double objective_value,   // objective function value of new columns
    std::string col_name      // column name
) {
  MiniMIPdebugMessage("add column.\n");

  // @todo add names
  if (!col.indices.empty()) {
#ifndef NDEBUG
    // perform check that no new rows are added
    RowIndex num_rows = linear_program_.num_constraints();
    for (size_t j = 0; j < col.indices.size(); ++j) {
      assert(0 <= col.indices[j] && col.indices[j] < num_rows.value());
      assert(col.values[j] != 0.0);
    }
#endif

    const ColIndex column = linear_program_.CreateNewVariable();
    linear_program_.SetVariableBounds(column, lower_bound, upper_bound);
    linear_program_.SetObjectiveCoefficient(column, objective_value);

    for (size_t j = 0; j < col.indices.size(); ++j) {
      linear_program_.SetCoefficient(RowIndex(col.indices[j]), column,
                                     col.values[j]);
    }
  } else {
    const ColIndex column = linear_program_.CreateNewVariable();
    linear_program_.SetVariableBounds(column, lower_bound, upper_bound);
    linear_program_.SetObjectiveCoefficient(column, objective_value);
  }
  lp_modified_since_last_solve_ = true;

  return absl::OkStatus();
}

// adds columns to the LP
absl::Status LPGlopInterface::AddColumns(
    const std::vector<SparseVector>& cols,    // column to be added
    const std::vector<double>& lower_bounds,  // lower bounds of new columns
    const std::vector<double>& upper_bounds,  // upper bounds of new columns
    const std::vector<double>&
        objective_values,  // objective function values of new columns
    std::vector<std::string>& col_names  // column names
) {
  for (size_t j = 0; j < cols.size(); j++) {
    MINIMIP_CALL(AddColumn(cols[j], lower_bounds[j], upper_bounds[j],
                           objective_values[j], col_names[j]));
  }
  return absl::OkStatus();
}

// To be deleted when finalizing the unittests:
//  absl::Status LPGlopInterface::AddColumns(
//      int num_cols,  // number of columns to be added
//      const std::vector<double>&
//          objective_values,  // objective function values of new columns
//      const std::vector<double>& lower_bounds,  // lower bounds of new columns
//      const std::vector<double>& upper_bounds,  // upper bounds of new columns
//      std::vector<std::string>& col_names,      // column names
//      int num_non_zeros,  // number of non-zero elements to be added to the
//                          // constraint matrix
//      const std::vector<int>&
//          begin_cols,  // start index of each column in indices- and
//          vals-array
//      const std::vector<int>&
//          indices,  // row indices of constraint matrix entries
//      const std::vector<double>& vals  // values of constraint matrix entries
//  ) {
//    MiniMIPdebugMessage("adding %d columns with %d nonzeros.\n", num_cols,
//                        num_non_zeros);
//
//    // @todo add names
//    if (num_non_zeros > 0) {
//      assert(num_cols > 0);
//
//#ifndef NDEBUG
//      // perform check that no new rows are added
//      RowIndex num_rows = linear_program_.num_constraints();
//      for (int j = 0; j < num_non_zeros; ++j) {
//        assert(0 <= indices[j] && indices[j] < num_rows.value());
//        assert(vals[j] != 0.0);
//      }
//#endif
//
//      int nz = 0;
//      for (int i = 0; i < num_cols; ++i) {
//        const ColIndex col = linear_program_.CreateNewVariable();
//        linear_program_.SetVariableBounds(col, lower_bounds[i],
//                                          upper_bounds[i]);
//        linear_program_.SetObjectiveCoefficient(col, objective_values[i]);
//        const int end = (num_non_zeros == 0 || i == num_cols - 1)
//                            ? num_non_zeros
//                            : begin_cols[i + 1];
//        while (nz < end) {
//          linear_program_.SetCoefficient(RowIndex(indices[nz]), col,
//          vals[nz]);
//          ++nz;
//        }
//      }
//      assert(nz == num_non_zeros);
//    } else {
//      for (int i = 0; i < num_cols; ++i) {
//        const ColIndex col = linear_program_.CreateNewVariable();
//        linear_program_.SetVariableBounds(col, lower_bounds[i],
//                                          upper_bounds[i]);
//        linear_program_.SetObjectiveCoefficient(col, objective_values[i]);
//      }
//    }
//
//    lp_modified_since_last_solve_ = true;
//
//    return absl::OkStatus();
//  }

// deletes all columns in the given range from LP
absl::Status LPGlopInterface::DeleteColumns(
    int first_col,  // first column to be deleted
    int last_col    // last column to be deleted
) {
  assert(0 <= first_col && first_col <= last_col &&
         last_col < linear_program_.num_variables());

  MiniMIPdebugMessage("deleting columns %d to %d.\n", first_col, last_col);

  const ColIndex num_cols = linear_program_.num_variables();
  DenseBooleanRow columns_to_delete(num_cols, false);
  for (int i = first_col; i <= last_col; ++i)
    columns_to_delete[ColIndex(i)] = true;

  linear_program_.DeleteColumns(columns_to_delete);
  lp_modified_since_last_solve_ = true;

  return absl::OkStatus();
}

// add row to the LP
absl::Status LPGlopInterface::AddRow(
    SparseVector row,        // row to be added
    double left_hand_side,   // left hand side of new row
    double right_hand_side,  // right hand side of new row
    std::string row_name     // row name
) {
  MiniMIPdebugMessage("adding row with %zu nonzeros.\n", row.indices.size());

  // @todo add names
  if (!row.indices.empty()) {
    assert(!row.values.empty());

#ifndef NDEBUG
    // perform check that no new columns are added - this is likely to be a
    // mistake
    const ColIndex num_cols = linear_program_.num_variables();
    for (size_t j = 0; j < row.indices.size(); ++j) {
      assert(row.values[j] != 0.0);
      assert(0 <= row.indices[j] && row.indices[j] < num_cols.value());
    }
#endif

    const RowIndex lprow = linear_program_.CreateNewConstraint();
    linear_program_.SetConstraintBounds(lprow, left_hand_side, right_hand_side);
    for (size_t j = 0; j < row.indices.size(); j++) {
      linear_program_.SetCoefficient(lprow, ColIndex(row.indices[j]),
                                     row.values[j]);
    }
  } else {
    const RowIndex lprow = linear_program_.CreateNewConstraint();
    linear_program_.SetConstraintBounds(lprow, left_hand_side, right_hand_side);
  }

  lp_modified_since_last_solve_ = true;

  return absl::OkStatus();
}

// add rows to the LP
absl::Status LPGlopInterface::AddRows(
    const std::vector<SparseVector>& rows,       // number of rows to be added
    const std::vector<double>& left_hand_sides,  // left hand sides of new rows
    const std::vector<double>&
        right_hand_sides,                // right hand sides of new rows
    std::vector<std::string>& row_names  // row names
) {
  MiniMIPdebugMessage("adding %zu rows.\n", rows.size());

  assert(left_hand_sides.size() == right_hand_sides.size());

  if (!rows.empty()) {
    for (size_t j = 0; j < rows.size(); j++) {
      assert(rows[j].indices.size() == rows[j].values.size());
      MINIMIP_CALL(AddRow(rows[j], left_hand_sides[j], right_hand_sides[j],
                          row_names[j]));
    }
  } else {
    for (size_t j = 0; j < left_hand_sides.size(); j++) {
      MINIMIP_CALL(AddRow(SparseVector(), left_hand_sides[j],
                          right_hand_sides[j], row_names[j]));
    }
  }
  return absl::OkStatus();
}

// To be deleted when finalizing the unittests:
// absl::Status LPGlopInterface::AddRows(
//     int num_rows,                                // number of rows to be
//     added const std::vector<double>& left_hand_sides,  // left hand sides
//     of new rows const std::vector<double>&
//         right_hand_sides,                 // right hand sides of new rows
//     std::vector<std::string>& row_names,  // row names
//     int num_non_zeros,  // number of non-zero elements to be added to the
//                         // constraint matrix
//     const std::vector<int>&
//         begin_rows,  // start index of each row in indices- and vals-array
//     const std::vector<int>&
//         indices,  // column indices of constraint matrix entries
//     const std::vector<double>& vals  // values of constraint matrix entries
//) {
//     MiniMIPdebugMessage("adding %d rows with %d nonzeros.\n", num_rows,
//                         num_non_zeros);
//
//     // @todo add names
//     if (num_non_zeros > 0) {
//       assert(num_rows > 0);
//
//  #ifndef NDEBUG
//       // perform check that no new columns are added - this is likely to be
//       a
//       // mistake
//       const ColIndex num_cols = linear_program_.num_variables();
//       for (int j = 0; j < num_non_zeros; ++j) {
//         assert(vals[j] != 0.0);
//         assert(0 <= indices[j] && indices[j] < num_cols.value());
//       }
//  #endif
//
//       int nz = 0;
//       for (int i = 0; i < num_rows; ++i) {
//         const RowIndex row = linear_program_.CreateNewConstraint();
//         linear_program_.SetConstraintBounds(row, left_hand_sides[i],
//                                             right_hand_sides[i]);
//         const int end = (num_non_zeros == 0 || i == num_rows - 1)
//                             ? num_non_zeros
//                             : begin_rows[i + 1];
//         while (nz < end) {
//           linear_program_.SetCoefficient(row, ColIndex(indices[nz]),
//           vals[nz]);
//           ++nz;
//         }
//       }
//       assert(nz == num_non_zeros);
//     } else {
//       for (int i = 0; i < num_rows; ++i) {
//         const RowIndex row = linear_program_.CreateNewConstraint();
//         linear_program_.SetConstraintBounds(row, left_hand_sides[i],
//                                             right_hand_sides[i]);
//       }
//     }
//
//     lp_modified_since_last_solve_ = true;
//
//     return absl::OkStatus();
// }

// delete rows from LP and update the current basis
void LPGlopInterface::DeleteRowsAndUpdateCurrentBasis(
    const DenseBooleanColumn&
        rows_to_delete  // array to mark rows that should be deleted
) {
  const RowIndex num_rows = linear_program_.num_constraints();
  const ColIndex num_cols = linear_program_.num_variables();

  // try to repair basis status if problem size has not changed before
  BasisState state = solver_.GetState();
  if (state.statuses.size() == num_cols.value() + num_rows.value()) {
    // Shift the status of the non-deleted rows. Note that if the deleted rows
    // where part of the basis (i.e., constraint not tight), then we should be
    // left with a correct basis afterward. This should be the most common use
    // case in MiniMIP.
    ColIndex new_size = num_cols;
    for (RowIndex row(0); row < num_rows; ++row) {
      if (rows_to_delete[row]) continue;
      state.statuses[new_size++] =
          state.statuses[num_cols + RowToColIndex(row)];
    }
    state.statuses.resize(new_size);
    solver_.LoadStateForNextSolve(state);
  }

  linear_program_.DeleteRows(rows_to_delete);
  lp_modified_since_last_solve_ = true;
}

// deletes all rows in the given range from LP
absl::Status LPGlopInterface::DeleteRows(
    int first_row,  // first row to be deleted
    int last_row    // last row to be deleted
) {
  assert(0 <= first_row && first_row <= last_row &&
         last_row < linear_program_.num_constraints());

  const RowIndex num_rows = linear_program_.num_constraints();
  DenseBooleanColumn rows_to_delete(num_rows, false);
  for (int i = first_row; i <= last_row; ++i)
    rows_to_delete[RowIndex(i)] = true;

  MiniMIPdebugMessage("deleting rows %d to %d.\n", first_row, last_row);
  DeleteRowsAndUpdateCurrentBasis(rows_to_delete);

  return absl::OkStatus();
}

// deletes rows from LP; the new position of a row must not be greater that
// its old position
absl::Status LPGlopInterface::DeleteRowSet(
    std::vector<bool>& deletion_status  // deletion status of rows
) {
  const RowIndex num_rows = linear_program_.num_constraints();
  DenseBooleanColumn rows_to_delete(num_rows, false);
  int new_index        = 0;
  int num_deleted_rows = 0;
  for (RowIndex row(0); row < num_rows; ++row) {
    int i = row.value();
    if (deletion_status[i] == 1) {
      rows_to_delete[row] = true;
      deletion_status[i]  = -1;
      ++num_deleted_rows;
    } else
      deletion_status[i] = new_index++;
  }

  MiniMIPdebugMessage("DeleteRowSet: deleting %d rows.\n", num_deleted_rows);
  DeleteRowsAndUpdateCurrentBasis(rows_to_delete);

  return absl::OkStatus();
}

// clears the whole LP
absl::Status LPGlopInterface::Clear() {
  MiniMIPdebugMessage("Clear\n");

  linear_program_.Clear();
  lp_modified_since_last_solve_ = true;

  return absl::OkStatus();
}

// clears current LPi state (like basis information) of the solver
absl::Status LPGlopInterface::ClearState() {
  solver_.ClearStateForNextSolve();

  return absl::OkStatus();
}

// changes lower and upper bounds of columns
absl::Status LPGlopInterface::SetColumnBounds(int col, double lower_bound,
                                              double upper_bound) {
  MiniMIPdebugMessage("set column bounds.\n");

  if (IsInfinity(lower_bound)) {
    MiniMIPerrorMessage(
        "LP Error: fixing lower bound for variable %d to infinity.\n", col);
    return absl::Status(absl::StatusCode::kInternal, "LP Error");
  }
  if (IsInfinity(-upper_bound)) {
    MiniMIPerrorMessage(
        "LP Error: fixing upper bound for variable %d to -infinity.\n", col);
    return absl::Status(absl::StatusCode::kInternal, "LP Error");
  }

  linear_program_.SetVariableBounds(ColIndex(col), lower_bound, upper_bound);

  lp_modified_since_last_solve_ = true;

  return absl::OkStatus();
}

// changes left- and right-hand sides of rows
absl::Status LPGlopInterface::SetRowSides(int row, double left_hand_side,
                                          double right_hand_side) {
  MiniMIPdebugMessage("set row sides\n");

  linear_program_.SetConstraintBounds(RowIndex(row), left_hand_side,
                                      right_hand_side);

  lp_modified_since_last_solve_ = true;

  return absl::OkStatus();
}

// changes the objective sense
absl::Status LPGlopInterface::SetObjectiveSense(
    LPObjectiveSense obj_sense  // new objective sense
) {
  switch (obj_sense) {
    case LPObjectiveSense::kMaximization:
      MiniMIPdebugMessage("changing objective sense to MAXIMIZE\n");
      linear_program_.SetMaximizationProblem(true);
      break;
    case LPObjectiveSense::kMinimization:
      MiniMIPdebugMessage("changing objective sense to MINIMIZE\n");
      linear_program_.SetMaximizationProblem(false);
      break;
  }
  lp_modified_since_last_solve_ = true;

  return absl::OkStatus();
}

// changes objective value of column in the LP
absl::Status LPGlopInterface::SetObjectiveCoefficient(
    int col, double objective_coefficient) {
  linear_program_.SetObjectiveCoefficient(ColIndex(col), objective_coefficient);
  lp_modified_since_last_solve_ = true;
  return absl::OkStatus();
}

// ==========================================================================
// LP model getters.
// ==========================================================================

// gets the number of rows in the LP
int LPGlopInterface::GetNumberOfRows() const {
  MiniMIPdebugMessage("getting number of rows.\n");

  return linear_program_.num_constraints().value();
}

// gets the number of columns in the LP
int LPGlopInterface::GetNumberOfColumns() const {
  MiniMIPdebugMessage("getting number of columns.\n");

  return linear_program_.num_variables().value();
}

// gets objective sense of the LP
LPObjectiveSense LPGlopInterface::GetObjectiveSense() const {
  MiniMIPdebugMessage("getting objective sense.\n");

  return linear_program_.IsMaximizationProblem()
             ? LPObjectiveSense::kMaximization
             : LPObjectiveSense::kMinimization;
}

// gets the number of nonzero elements in the LP constraint matrix
int LPGlopInterface::GetNumberOfNonZeros() const {
  MiniMIPdebugMessage("getting number of non-zeros.\n");

  return static_cast<int>(linear_program_.num_entries().value());
}

// gets columns from LP problem object
//
// Either both, lb and ub, have to be NULL, or both have to be non-NULL,
// either num_non_zeros, begin_cols, indices, and val have to be NULL, or all
// of them have to be non-NULL.
SparseVector LPGlopInterface::GetSparseColumnCoefficients(int col) const {
  assert(0 <= col && col < linear_program_.num_variables());

  SparseVector sparse_column;
  const SparseColumn& column = linear_program_.GetSparseColumn(ColIndex(col));

  for (const SparseColumn::Entry& entry : column) {
    const RowIndex row = entry.row();
    sparse_column.indices.push_back(row.value());
    sparse_column.values.push_back(entry.coefficient());
  }

  return sparse_column;
}

// gets rows from LP problem object
//
// Either both, left_hand_side and right_hand_side, have to be NULL, or both
// have to be non-NULL, either num_non_zeros, begin_rows, indices, and val
// have to be NULL, or all of them have to be non-NULL.
SparseVector LPGlopInterface::GetSparseRowCoefficients(int row) const {
  assert(0 <= row && row < linear_program_.num_constraints());

  const SparseMatrix& matrixtrans = linear_program_.GetTransposeSparseMatrix();
  const SparseColumn& column =
      matrixtrans.column(ColIndex(RowIndex(row).value()));
  SparseVector sparse_row;

  for (const SparseColumn::Entry& entry : column) {
    const RowIndex rowidx = entry.row();
    sparse_row.indices.push_back(rowidx.value());
    sparse_row.values.push_back(entry.coefficient());
  }

  return sparse_row;
}

// gets objective coefficient of column from LP problem object
double LPGlopInterface::GetObjectiveCoefficient(int col) const {
  MiniMIPdebugMessage("getting objective value of column %d\n", col);

  return linear_program_.objective_coefficients()[ColIndex(col)];
}

// gets current lower bound of column from LP problem object
double LPGlopInterface::GetLowerBound(int col) const {
  MiniMIPdebugMessage("getting lower bound of column %d\n", col);

  return linear_program_.variable_lower_bounds()[ColIndex(col)];
}

// gets current upper bound of column from LP problem object
double LPGlopInterface::GetUpperBound(int col) const {
  MiniMIPdebugMessage("getting upper bound of column %d\n", col);

  return linear_program_.variable_upper_bounds()[ColIndex(col)];
}

// gets current left hand side of row from LP problem object
double LPGlopInterface::GetLeftHandSide(int row) const {
  MiniMIPdebugMessage("getting left hand side of row %d\n", row);

  return linear_program_.constraint_lower_bounds()[RowIndex(row)];
}

// gets current right hand side of row from LP problem object
double LPGlopInterface::GetRightHandSide(int row) const {
  MiniMIPdebugMessage("getting right hand side of row %d\n", row);

  return linear_program_.constraint_upper_bounds()[RowIndex(row)];
}

// gets the matrix coefficient of column and row from LP problem object
double LPGlopInterface::GetMatrixCoefficient(
    int col,  // column number of coefficient
    int row   // row number of coefficient
) const {
  // quite slow method: possibly needs linear time if matrix is not sorted
  const SparseMatrix& matrix = linear_program_.GetSparseMatrix();

  return matrix.LookUpValue(RowIndex(row), ColIndex(col));
}

// ============================================================================
// Solving methods.
// ============================================================================

// update scaled linear program
void LPGlopInterface::updateScaledLP() {
  if (!lp_modified_since_last_solve_) return;

  scaled_lp_.PopulateFromLinearProgram(linear_program_);
  scaled_lp_.AddSlackVariablesWhereNecessary(false);

  // @todo: Avoid doing a copy if there is no scaling.
  // @todo: Avoid rescaling if not much changed.
  if (parameters_.use_scaling())
    scaler_.Scale(&scaled_lp_);
  else
    scaler_.Clear();
}

// check primal feasibility
bool LPGlopInterface::checkUnscaledPrimalFeasibility() const {
#if UNSCALEDFEAS_CHECK == 1
  // get unscaled solution
  const ColIndex num_cols = linear_program_.num_variables();
  DenseRow unscaledsol(num_cols);
  for (ColIndex col = ColIndex(0); col < num_cols; ++col)
    unscaledsol[col] =
        scaler_.UnscaleVariableValue(col, solver_.GetVariableValue(col));

  // if the solution is not feasible w.r.t. absolute tolerances, try to fix it
  // in the unscaled problem
  const double feastol = parameters_.primal_feasibility_tolerance();
  return linear_program_.SolutionIsLPFeasible(unscaledsol, feastol);

#elif UNSCALEDFEAS_CHECK == 2
  const double feastol = parameters_.primal_feasibility_tolerance();

  // check bounds of unscaled solution
  const ColIndex num_cols = linear_program_.num_variables();
  for (ColIndex col = ColIndex(0); col < num_cols; ++col) {
    const Fractional val =
        scaler_.UnscaleVariableValue(col, solver_.GetVariableValue(col));
    const Fractional lower_bound = linear_program_.variable_lower_bounds()[col];
    if (val < lower_bound - feastol) return false;
    const Fractional upper_bound = linear_program_.variable_upper_bounds()[col];
    if (val > upper_bound + feastol) return false;
  }

  // check activities of unscaled solution
  const RowIndex num_rows = linear_program_.num_constraints();
  for (RowIndex row(0); row < num_rows; ++row) {
    const Fractional val = scaler_.UnscaleConstraintActivity(
        row, solver_.GetConstraintActivity(row));
    const Fractional left_hand_side =
        linear_program_.constraint_lower_bounds()[row];
    if (val < left_hand_side - feastol) return false;
    const Fractional right_hand_side =
        linear_program_.constraint_upper_bounds()[row];
    if (val > right_hand_side + feastol) return false;
  }
#endif

  return true;
}

// common function between the two LPI Solve() functions
absl::Status LPGlopInterface::SolveInternal(
    bool recursive,                         // Is this a recursive call?
    std::unique_ptr<TimeLimit>& time_limit  // time limit
) {
  updateScaledLP();

  solver_.SetParameters(parameters_);
  lp_time_limit_was_reached_ = false;

  // possibly ignore warm start information for next solve
  if (from_scratch_) solver_.ClearStateForNextSolve();

  if (!solver_.Solve(scaled_lp_, time_limit.get()).ok()) {
    return absl::Status(absl::StatusCode::kInternal, "LP Error");
  }
  lp_time_limit_was_reached_ = time_limit->LimitReached();
  if (recursive)
    niterations_ += (std::int64_t)solver_.GetNumberOfIterations();
  else
    niterations_ = (std::int64_t)solver_.GetNumberOfIterations();

  MiniMIPdebugMessage(
      "status=%s  obj=%f  iterations=%ld.\n",
      GetProblemStatusString(solver_.GetProblemStatus()).c_str(),
      solver_.GetObjectiveValue(), solver_.GetNumberOfIterations());

  const ProblemStatus status = solver_.GetProblemStatus();
  if ((status == ProblemStatus::PRIMAL_FEASIBLE ||
       status == ProblemStatus::OPTIMAL) &&
      parameters_.use_scaling()) {
    if (!checkUnscaledPrimalFeasibility()) {
      MiniMIPdebugMessage(
          "Solution not feasible w.r.t. absolute tolerance %g -> "
          "reoptimize.\n",
          parameters_.primal_feasibility_tolerance());

      // Re-solve without scaling to try to fix the infeasibility.
      parameters_.set_use_scaling(false);
      lp_modified_since_last_solve_ = true;
      // inherit time limit, so used time is not reset;
      // do not change iteration limit for resolve
      MINIMIP_CALL(SolveInternal(true, time_limit));
      parameters_.set_use_scaling(true);
    }
  }

  lp_modified_since_last_solve_ = false;

  return absl::OkStatus();
}

// ==========================================================================
// Solving methods.
// ==========================================================================

// calls primal simplex to solve the LP
absl::Status LPGlopInterface::SolveLPWithPrimalSimplex() {
  MiniMIPdebugMessage("SolvePrimal: %d rows, %d cols.\n",
                      linear_program_.num_constraints().value(),
                      linear_program_.num_variables().value());
  std::unique_ptr<TimeLimit> time_limit =
      TimeLimit::FromParameters(parameters_);
  niterations_ = 0;

  parameters_.set_use_dual_simplex(false);
  return SolveInternal(false, time_limit);
}

// calls dual simplex to solve the LP
absl::Status LPGlopInterface::SolveLPWithDualSimplex() {
  MiniMIPdebugMessage("SolveDual: %d rows, %d cols.\n",
                      linear_program_.num_constraints().value(),
                      linear_program_.num_variables().value());
  std::unique_ptr<TimeLimit> time_limit =
      TimeLimit::FromParameters(parameters_);
  niterations_ = 0;

  parameters_.set_use_dual_simplex(true);
  return SolveInternal(false, time_limit);
}

// start strong branching
absl::Status LPGlopInterface::StartStrongBranching() {
  updateScaledLP();

  // @todo Save state and do all the branching from there.
  return absl::OkStatus();
}

// end strong branching
absl::Status LPGlopInterface::EndStrongBranching() {
  // @todo Restore the saved state.
  return absl::OkStatus();
}
// determine whether the dual bound is valid
bool LPGlopInterface::IsDualBoundValid(
    ProblemStatus status  // status to be checked
) const {
  return status == ProblemStatus::OPTIMAL ||
         status == ProblemStatus::DUAL_FEASIBLE ||
         status == ProblemStatus::DUAL_UNBOUNDED;
}

// performs strong branching iterations
absl::Status LPGlopInterface::strongbranch(
    int col_index,        // column to apply strong branching on
    double primal_sol,    // fractional current primal solution value of column
    int iteration_limit,  // iteration limit for strong branchings
    double& dual_bound_down_branch,  // stores dual bound after branching
                                     // column down
    double&
        dual_bound_up_branch,  // stores dual bound after branching column up
    bool& down_valid,  // stores whether the returned down value is a valid
                       // dual bound; otherwise, it can only be used as an
                       // estimate value
    bool& up_valid,    // stores whether the returned up value is a valid dual
                       // bound; otherwise, it can only be used as an estimate
                       // value
    int& iterations    // stores total number of strong branching iterations, or
                       // -1;
) {
  MiniMIPdebugMessage(
      "calling strongbranching on variable %d (%d iterations)\n", col_index,
      iteration_limit);

  // We work on the scaled problem.
  const ColIndex col(col_index);
  const Fractional lower_bound = scaled_lp_.variable_lower_bounds()[col];
  const Fractional upper_bound = scaled_lp_.variable_upper_bounds()[col];
  const double value = primal_sol * scaler_.VariableScalingFactor(col);

  // Configure solver.

  // @todo Use the iteration limit once glop supports incrementality.
  int num_iterations = 0;
  parameters_.set_use_dual_simplex(true);

  solver_.SetParameters(parameters_);
  const Fractional eps = parameters_.primal_feasibility_tolerance();

  std::unique_ptr<TimeLimit> time_limit =
      TimeLimit::FromParameters(parameters_);

  // Down branch.
  const Fractional new_upper_bound = EPSCEIL(value - 1.0, eps);
  if (new_upper_bound >= lower_bound - 0.5) {
    scaled_lp_.SetVariableBounds(col, lower_bound, new_upper_bound);

    if (solver_.Solve(scaled_lp_, time_limit.get()).ok()) {
      num_iterations += static_cast<int>(solver_.GetNumberOfIterations());
      dual_bound_down_branch = solver_.GetObjectiveValue();
      down_valid             = IsDualBoundValid(solver_.GetProblemStatus());

      MiniMIPdebugMessage(
          "dual_bound_down_branch: iteration_limit=%d col=%d [%f,%f] obj=%f "
          "status=%d iterations=%ld.\n",
          iteration_limit, col_index, lower_bound, EPSCEIL(value - 1.0, eps),
          solver_.GetObjectiveValue(),
          static_cast<int>(solver_.GetProblemStatus()),
          solver_.GetNumberOfIterations());
    } else {
      MiniMIPerrorMessage("error during solve");
      dual_bound_down_branch = 0.0;
      down_valid             = false;
    }
  } else {
    if (linear_program_.IsMaximizationProblem())
      dual_bound_down_branch = parameters_.objective_lower_limit();
    else
      dual_bound_down_branch = parameters_.objective_upper_limit();
    down_valid = true;
  }

  // Up branch.
  const Fractional new_lower_bound = EPSFLOOR(value + 1.0, eps);
  if (new_lower_bound <= upper_bound + 0.5) {
    scaled_lp_.SetVariableBounds(col, new_lower_bound, upper_bound);

    if (solver_.Solve(scaled_lp_, time_limit.get()).ok()) {
      num_iterations += static_cast<int>(solver_.GetNumberOfIterations());
      dual_bound_up_branch = solver_.GetObjectiveValue();
      up_valid             = IsDualBoundValid(solver_.GetProblemStatus());

      MiniMIPdebugMessage(
          "dual_bound_up_branch: iteration_limit=%d col=%d [%f,%f] obj=%f "
          "status=%d iterations=%ld.\n",
          iteration_limit, col_index, EPSFLOOR(value + 1.0, eps), upper_bound,
          solver_.GetObjectiveValue(),
          static_cast<int>(solver_.GetProblemStatus()),
          solver_.GetNumberOfIterations());
    } else {
      MiniMIPerrorMessage("error during solve");
      dual_bound_up_branch = 0.0;
      up_valid             = false;
    }
  } else {
    if (linear_program_.IsMaximizationProblem())
      dual_bound_up_branch = parameters_.objective_lower_limit();
    else
      dual_bound_up_branch = parameters_.objective_upper_limit();
    up_valid = true;
  }

  //  Restore bound.
  scaled_lp_.SetVariableBounds(col, lower_bound, upper_bound);
  if (iterations > 0) iterations = num_iterations;

  return absl::OkStatus();
}

// performs strong branching iterations on one @b fractional candidate
absl::StatusOr<LPGlopInterface::StrongBranchResult>
LPGlopInterface::StrongBranchValue(
    int col,             // column to apply strong branching on
    double primal_sol,   // current primal solution value of column
    int iteration_limit  // iteration limit for strong branchings
) {
  MiniMIPdebugMessage(
      "calling strong branching on variable %d (%d iterations)\n", col,
      iteration_limit);

  LPInterface::StrongBranchResult strong_branch_result;

  //@TODO: strongbranch() always returns absl::OkStatus() currently
  absl::Status absl_status_code = strongbranch(
      col, primal_sol, iteration_limit,
      strong_branch_result.dual_bound_down_branch,
      strong_branch_result.dual_bound_up_branch,
      strong_branch_result.down_valid, strong_branch_result.up_valid,
      strong_branch_result.iterations);

  if (absl_status_code != absl::OkStatus())
    return absl::Status(absl::StatusCode::kInternal, "This will never happen");
  else
    return strong_branch_result;
}

// ==========================================================================
// Solution information getters.
// ==========================================================================

// returns whether a solve method was called after the last modification of
// the LP
bool LPGlopInterface::IsSolved() const {
  // @todo Track this to avoid uneeded resolving.
  return (!lp_modified_since_last_solve_);
}

// returns true if LP is proven to have a primal unbounded ray (but not
// necessary a primal feasible point); this does not necessarily mean that the
// solver knows and can return the primal ray
bool LPGlopInterface::ExistsPrimalRay() const {
  return solver_.GetProblemStatus() == ProblemStatus::PRIMAL_UNBOUNDED;
}

// returns true if LP is proven to have a primal unbounded ray (but not
// necessary a primal feasible point), and the solver knows and can return the
// primal ray
bool LPGlopInterface::HasPrimalRay() const {
  return solver_.GetProblemStatus() == ProblemStatus::PRIMAL_UNBOUNDED;
}

// returns true if LP is proven to be primal unbounded
bool LPGlopInterface::IsPrimalUnbounded() const {
  return solver_.GetProblemStatus() == ProblemStatus::PRIMAL_UNBOUNDED;
}

// returns true if LP is proven to be primal infeasible
bool LPGlopInterface::IsPrimalInfeasible() const {
  const ProblemStatus status = solver_.GetProblemStatus();

  return status == ProblemStatus::DUAL_UNBOUNDED ||
         status == ProblemStatus::PRIMAL_INFEASIBLE;
}

// returns true if LP is proven to be primal feasible
bool LPGlopInterface::IsPrimalFeasible() const {
  const ProblemStatus status = solver_.GetProblemStatus();

  return status == ProblemStatus::PRIMAL_FEASIBLE ||
         status == ProblemStatus::OPTIMAL;
}

// returns true if LP is proven to have a dual unbounded ray (but not
// necessary a dual feasible point); this does not necessarily mean that the
// solver knows and can return the dual ray
bool LPGlopInterface::ExistsDualRay() const {
  const ProblemStatus status = solver_.GetProblemStatus();

  return status == ProblemStatus::DUAL_UNBOUNDED;
}

// returns true if LP is proven to have a dual unbounded ray (but not
// necessary a dual feasible point), and the solver knows and can return the
// dual ray
bool LPGlopInterface::HasDualRay() const {
  const ProblemStatus status = solver_.GetProblemStatus();

  return status == ProblemStatus::DUAL_UNBOUNDED;
}

// returns true if LP is proven to be dual unbounded
bool LPGlopInterface::IsDualUnbounded() const {
  const ProblemStatus status = solver_.GetProblemStatus();
  return status == ProblemStatus::DUAL_UNBOUNDED;
}

// returns true if LP is proven to be dual infeasible
bool LPGlopInterface::IsDualInfeasible() const {
  const ProblemStatus status = solver_.GetProblemStatus();
  return status == ProblemStatus::PRIMAL_UNBOUNDED ||
         status == ProblemStatus::DUAL_INFEASIBLE;
}

// returns true if LP is proven to be dual feasible
bool LPGlopInterface::IsDualFeasible() const {
  const ProblemStatus status = solver_.GetProblemStatus();

  return status == ProblemStatus::DUAL_FEASIBLE ||
         status == ProblemStatus::OPTIMAL;
}

// returns true if LP was solved to optimality
bool LPGlopInterface::IsOptimal() const {
  return solver_.GetProblemStatus() == ProblemStatus::OPTIMAL;
}

// returns true if current LP solution is stable
//
// This function should return true if the solution is reliable, i.e.,
// feasible and optimal (or proven infeasible/unbounded) with respect to the
// original problem. The optimality status might be with respect to a scaled
// version of the problem, but the solution might not be feasible to the
// unscaled original problem; in this case, IsStable() should return false.
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
    MiniMIPdebugMessage("OPTIMAL not reached and no limit: unstable.\n");
    return false;
  }

  if (status == ProblemStatus::ABNORMAL ||
      status == ProblemStatus::INVALID_PROBLEM ||
      status == ProblemStatus::IMPRECISE)
    return false;
  return true;
}  // @TODO: Case that neither if happens?

// returns true if the objective limit was reached
bool LPGlopInterface::ObjectiveLimitIsExceeded() const {
  return solver_.objective_limit_reached();
}

// returns true if the iteration limit was reached
bool LPGlopInterface::IterationLimitIsExceeded() const {
  assert(niterations_ >= static_cast<int>(solver_.GetNumberOfIterations()));

  int maxiter = static_cast<int>(parameters_.max_number_of_iterations());
  return maxiter >= 0 && niterations_ >= maxiter;
}

// returns true if the time limit was reached
bool LPGlopInterface::TimeLimitIsExceeded() const {
  return lp_time_limit_was_reached_;
}

// gets objective value of solution
double LPGlopInterface::GetObjectiveValue() {
  return solver_.GetObjectiveValue();
}

// Before calling this function, the caller must ensure that the LP has been
// solved to optimality, i.e., that IsOptimal() returns true.

// gets primal solution vector
absl::StatusOr<std::vector<double>> LPGlopInterface::GetPrimalSolution() const {
  MiniMIPdebugMessage("GetPrimalSolution\n");
  std::vector<double> primal_sol;

  const ColIndex num_cols = linear_program_.num_variables();
  for (ColIndex col(0); col < num_cols; ++col) {
    primal_sol.push_back(
        scaler_.UnscaleVariableValue(col, solver_.GetVariableValue(col)));
  }

  return primal_sol;
}

// Before calling this function, the caller must ensure that the LP has been
// solved to optimality, i.e., that IsOptimal() returns true.

// gets dual solution vector
absl::StatusOr<std::vector<double>> LPGlopInterface::GetDualSolution() const {
  MiniMIPdebugMessage("GetDualSolution\n");
  std::vector<double> dual_sol;

  const RowIndex num_rows = linear_program_.num_constraints();
  for (RowIndex row(0); row < num_rows; ++row) {
    dual_sol.push_back(
        scaler_.UnscaleDualValue(row, solver_.GetDualValue(row)));
  }

  return dual_sol;
}

// Before calling this function, the caller must ensure that the LP has been
// solved to optimality, i.e., that IsOptimal() returns true.

// gets row activity vector
absl::StatusOr<std::vector<double>> LPGlopInterface::GetRowActivity() const {
  MiniMIPdebugMessage("GetRowActivity\n");
  std::vector<double> activity;

  const RowIndex num_rows = linear_program_.num_constraints();
  for (RowIndex row(0); row < num_rows; ++row) {
    activity.push_back(scaler_.UnscaleConstraintActivity(
        row, solver_.GetConstraintActivity(row)));
  }

  return activity;
}

// Before calling this function, the caller must ensure that the LP has been
// solved to optimality, i.e., that IsOptimal() returns true.

// gets reduced cost vector
absl::StatusOr<std::vector<double>> LPGlopInterface::GetReducedCost() const {
  MiniMIPdebugMessage("GetReducedCost\n");
  std::vector<double> reduced_cost;

  const ColIndex num_cols = linear_program_.num_variables();
  for (ColIndex col(0); col < num_cols; ++col) {
    reduced_cost.push_back(
        scaler_.UnscaleReducedCost(col, solver_.GetReducedCost(col)));
  }

  return reduced_cost;
}

// gets primal ray for unbounded LPs
absl::StatusOr<std::vector<double>> LPGlopInterface::GetPrimalRay() const {
  MiniMIPdebugMessage("GetPrimalRay\n");
  std::vector<double> primal_ray;

  const ColIndex num_cols           = linear_program_.num_variables();
  const DenseRow& primal_ray_solver = solver_.GetPrimalRay();
  for (ColIndex col(0); col < num_cols; ++col)
    primal_ray.push_back(
        scaler_.UnscaleVariableValue(col, primal_ray_solver[col]));

  return primal_ray;
}

// gets dual Farkas proof for infeasibility
absl::StatusOr<std::vector<double>> LPGlopInterface::GetDualFarkasMultiplier()
    const {
  MiniMIPdebugMessage("GetDualFarkasMultiplier\n");
  std::vector<double> dual_farkas_multiplier;

  const RowIndex num_rows     = linear_program_.num_constraints();
  const DenseColumn& dual_ray = solver_.GetDualRay();
  for (RowIndex row(0); row < num_rows; ++row)
    dual_farkas_multiplier.push_back(
        -scaler_.UnscaleDualValue(row, dual_ray[row]));  // reverse sign

  return dual_farkas_multiplier;
}

// gets the number of LP iterations of the last solve call
int LPGlopInterface::GetIterations() const {
  return static_cast<int>(niterations_);
}

// LP Basis Methods

// @name LP Basis Methods
// @{

// convert Glop variable basis status to MiniMIP status
LPBasisStatus LPGlopInterface::ConvertGlopVariableStatus(
    VariableStatus status,   // variable status
    Fractional reduced_cost  // reduced cost of variable
) const {
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
      MiniMIPerrorMessage("invalid Glop basis status.\n");
      std::abort();
  }
}

// convert Glop constraint basis status to MiniMIP status
LPBasisStatus LPGlopInterface::ConvertGlopConstraintStatus(
    ConstraintStatus status,  // constraint status
    Fractional dual_value     // dual variable value
) const {
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
      MiniMIPerrorMessage("invalid Glop basis status.\n");
      std::abort();
  }
}

// Convert MiniMIP variable status to Glop status
VariableStatus LPGlopInterface::ConvertMiniMIPVariableStatus(
    LPBasisStatus status  // MiniMIP variable status
) const {
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
      MiniMIPerrorMessage("invalid MiniMIP basis status.\n");
      std::abort();
  }
}

// Convert a MiniMIP constraint status to its corresponding Glop slack
// VariableStatus.
//
// Note that we swap the upper/lower bounds.
VariableStatus LPGlopInterface::ConvertMiniMIPConstraintStatusToSlackStatus(
    LPBasisStatus status  // MiniMIP constraint status
) const {
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
      MiniMIPerrorMessage("invalid MiniMIP basis status.\n");
      std::abort();
  }
}

// ==========================================================================
// Getters and setters of the basis.
// ==========================================================================

// gets current basis status for columns and rows
absl::StatusOr<std::vector<LPBasisStatus>>
LPGlopInterface::GetColumnBasisStatus() const {
  MiniMIPdebugMessage("GetColumnBasisStatus\n");
  std::vector<LPBasisStatus> column_basis_status;

  assert(solver_.GetProblemStatus() == ProblemStatus::OPTIMAL);
  const ColIndex num_cols = linear_program_.num_variables();
  for (ColIndex col(0); col < num_cols; ++col) {
    column_basis_status.push_back((LPBasisStatus)ConvertGlopVariableStatus(
        solver_.GetVariableStatus(col), solver_.GetReducedCost(col)));
  }
  return column_basis_status;
}
// gets current basis status for columns and rows
absl::StatusOr<std::vector<LPBasisStatus>> LPGlopInterface::GetRowBasisStatus()
    const {
  MiniMIPdebugMessage("GetRowBasisStatus\n");
  std::vector<LPBasisStatus> row_basis_status;

  const RowIndex num_rows = linear_program_.num_constraints();
  for (RowIndex row(0); row < num_rows; ++row) {
    row_basis_status.push_back((LPBasisStatus)ConvertGlopConstraintStatus(
        solver_.GetConstraintStatus(row), solver_.GetDualValue(row)));
  }

  return row_basis_status;
}

// sets current basis status for columns and rows
absl::Status LPGlopInterface::SetBasisStatus(
    const std::vector<LPBasisStatus>&
        column_basis_status,  // array with column basis status
    const std::vector<LPBasisStatus>&
        row_basis_status  // array with row basis status
) {
  const ColIndex num_cols = linear_program_.num_variables();
  const RowIndex num_rows = linear_program_.num_constraints();

  MiniMIPdebugMessage("SetBasisStatus\n");

  BasisState state;
  state.statuses.reserve(ColIndex(num_cols.value() + num_rows.value()));

  for (ColIndex col(0); col < num_cols; ++col)
    state.statuses[col] =
        ConvertMiniMIPVariableStatus(column_basis_status[col.value()]);

  for (RowIndex row(0); row < num_rows; ++row)
    state.statuses[num_cols + RowToColIndex(row)] =
        ConvertMiniMIPConstraintStatusToSlackStatus(
            column_basis_status[row.value()]);

  solver_.LoadStateForNextSolve(state);

  return absl::OkStatus();
}

// returns the indices of the basic columns and rows; basic column n gives
// value n, basic row m gives value -1-m
std::vector<int> LPGlopInterface::GetBasisIndices() const {
  std::vector<int> basis_indices;

  MiniMIPdebugMessage("GetBasisIndices\n");

  // the order is important!
  const ColIndex num_cols = linear_program_.num_variables();
  const RowIndex num_rows = linear_program_.num_constraints();
  for (RowIndex row(0); row < num_rows; ++row) {
    const ColIndex col = solver_.GetBasis(row);
    if (col < num_cols)
      basis_indices[row.value()] = col.value();
    else {
      assert(col < num_cols.value() + num_rows.value());
      basis_indices[row.value()] = -1 - (col - num_cols).value();
    }
  }

  return basis_indices;
}

// ==========================================================================
// Getters of vectors in the inverted basis matrix.
// ==========================================================================

// get row of inverse basis matrix B^-1
//
// NOTE: The LP interface defines slack variables to have coefficient +1. This
// means that if, internally, the LP solver
//       uses a -1 coefficient, then rows associated with slacks variables
//       whose coefficient is -1, should be negated; see also the explanation
//       in lpi.h.
absl::StatusOr<SparseVector> LPGlopInterface::GetSparseRowOfBInverted(
    int row_number) const {
  SparseVector sparse_row;

  solver_.GetBasisFactorization().LeftSolveForUnitRow(ColIndex(row_number),
                                                      tmp_row_);
  scaler_.UnscaleUnitRowLeftSolve(solver_.GetBasis(RowIndex(row_number)),
                                  tmp_row_);

  const ColIndex size = tmp_row_->values.size();
  assert(size.value() == linear_program_.num_constraints());

  // Vectors in Glop might be stored in dense or sparse format dep
  //
  // ending on the values. If non_zeros are given, we
  // can directly loop over the non_zeros, otherwise we have to collect the
  // nonzeros.
  if (!tmp_row_->non_zeros.empty()) {
    ScatteredRowIterator end = tmp_row_->end();
    for (ScatteredRowIterator iter = tmp_row_->begin(); iter != end; ++iter) {
      int idx = (*iter).column().value();
      assert(0 <= idx && idx < linear_program_.num_constraints());
      sparse_row.values.push_back((*iter).coefficient());
      sparse_row.indices.push_back(idx);
    }
  }
  return sparse_row;
}

// get column of inverse basis matrix B^-1
//
// NOTE: The LP interface defines slack variables to have coefficient +1. This
// means that if, internally, the LP solver
//       uses a -1 coefficient, then rows associated with slacks variables
//       whose coefficient is -1, should be negated; see also the explanation
//       in lpi.h.
//
// column number of B^-1; this is NOT the number of the
// column in the LP; you have to call
// minimip::LPInterface.GetBasisIndices() to get the array
// which links the B^-1 column numbers to the row and
// column numbers of the LP! c must be between 0 and
// num_rows-1, since the basis has the size num_rows *
// num_rows
absl::StatusOr<SparseVector> LPGlopInterface::GetSparseColumnOfBInverted(
    int col_number) const {
  SparseVector sparse_column;
  // we need to loop through the rows to extract the values for column
  // col_number
  const ColIndex col(col_number);
  const RowIndex num_rows = linear_program_.num_constraints();

  const Fractional eps = parameters_.primal_feasibility_tolerance();

  for (int row = 0; row < num_rows; ++row) {
    solver_.GetBasisFactorization().LeftSolveForUnitRow(ColIndex(row),
                                                        tmp_row_);
    scaler_.UnscaleUnitRowLeftSolve(solver_.GetBasis(RowIndex(row)), tmp_row_);

    double value = (*tmp_row_)[col];
    if (fabs(value) >= eps) {
      sparse_column.values.push_back(value);
      sparse_column.indices.push_back(row);
    }
  }
  return sparse_column;
}

// get row of inverse basis matrix times constraint matrix B^-1 * A
//
// NOTE: The LP interface defines slack variables to have coefficient +1. This
// means that if, internally, the LP solver
//       uses a -1 coefficient, then rows associated with slacks variables
//       whose coefficient is -1, should be negated; see also the explanation
//       in lpi.h.
absl::StatusOr<SparseVector> LPGlopInterface::GetSparseRowOfBInvertedTimesA(
    int row_number) const {
  SparseVector sparse_row;
  // get row of basis inverse, loop through columns and muliply with matrix
  solver_.GetBasisFactorization().LeftSolveForUnitRow(ColIndex(row_number),
                                                      tmp_row_);
  scaler_.UnscaleUnitRowLeftSolve(solver_.GetBasis(RowIndex(row_number)),
                                  tmp_row_);

  const ColIndex num_cols = linear_program_.num_variables();

  const Fractional eps = parameters_.primal_feasibility_tolerance();

  for (ColIndex col(0); col < num_cols; ++col) {
    double value = operations_research::glop::ScalarProduct(
        tmp_row_->values, linear_program_.GetSparseColumn(col));
    if (fabs(value) >= eps) {
      sparse_row.values.push_back(value);
      sparse_row.indices.push_back(col.value());
    }
  }
  return sparse_row;
}

// get column of inverse basis matrix times constraint matrix B^-1 * A
//
// NOTE: The LP interface defines slack variables to have coefficient +1. This
// means that if, internally, the LP solver
//       uses a -1 coefficient, then rows associated with slacks variables
//       whose coefficient is -1, should be negated; see also the explanation
//       in lpi.h.
absl::StatusOr<SparseVector> LPGlopInterface::GetSparseColumnOfBInvertedTimesA(
    int col_number) const {
  SparseVector sparse_column;
  solver_.GetBasisFactorization().RightSolveForProblemColumn(
      ColIndex(col_number), tmp_column_);
  this->scaler_.UnscaleColumnRightSolve(solver_.GetBasisVector(),
                                        ColIndex(col_number), tmp_column_);

  const RowIndex num_rows = tmp_column_->values.size();

  // Vectors in Glop might be stored in dense or sparse format depending on
  // the values. If non_zeros are given, we can directly loop over the
  // non_zeros, otherwise we have to collect the nonzeros.
  ScatteredColumnIterator end = tmp_column_->end();
  for (ScatteredColumnIterator iter = tmp_column_->begin(); iter != end;
       ++iter) {
    int idx = (*iter).row().value();
    assert(0 <= idx && idx < num_rows);
    sparse_column.values.push_back((*iter).coefficient());
    sparse_column.indices.push_back(idx);
  }
  return sparse_column;
}

// ==========================================================================
// Getters and setters of the parameters.
// ==========================================================================

// gets integer parameter of LP
absl::StatusOr<int> LPGlopInterface::GetIntegerParameter(
    LPParameter type  // parameter number
) const {
  int param_val;
  switch (type) {
    case LPParameter::kFromScratch:
      param_val = (int)from_scratch_;
      MiniMIPdebugMessage(
          "GetIntegerParameter: LPParameter::kFromScratch = %d.\n", param_val);
      break;
    case LPParameter::kLPInfo:
      param_val = (int)lp_info_;
      MiniMIPdebugMessage("GetIntegerParameter: LPParameter::kLPInfo = %d.\n",
                          param_val);
      break;
    case LPParameter::kLPIterationLimit:
      param_val = (int)parameters_.max_number_of_iterations();
      MiniMIPdebugMessage(
          "GetIntegerParameter: LPParameter::kLPIterationLimit = %d.\n",
          param_val);
      break;
    case LPParameter::kPresolving:
      param_val = parameters_.use_preprocessing();
      MiniMIPdebugMessage(
          "GetIntegerParameter: LPParameter::kPresolving = %d.\n", param_val);
      break;
    case LPParameter::kPricing:
      param_val = (int)pricing_;
      MiniMIPdebugMessage("GetIntegerParameter: LPParameter::kPricing = %d.\n",
                          param_val);
      break;
#ifndef NOSCALING
    case LPParameter::kScaling:
      param_val = parameters_.use_scaling();
      MiniMIPdebugMessage("GetIntegerParameter: LPParameter::kScaling = %d.\n",
                          param_val);
      break;
#endif
    case LPParameter::kThreads:
      param_val = num_threads_;
      MiniMIPdebugMessage("GetIntegerParameter: LPParameter::kThreads = %d.\n",
                          param_val);
      break;
    case LPParameter::kTiming:
      param_val = timing_;
      MiniMIPdebugMessage("GetIntegerParameter: LPParameter::kTiming = %d.\n",
                          param_val);
      break;
    case LPParameter::kRandomSeed:
      param_val = (int)parameters_.random_seed();
      MiniMIPdebugMessage(
          "GetIntegerParameter: LPParameter::kRandomSeed = %d.\n", param_val);
      break;
    default:
      return absl::Status(absl::StatusCode::kInvalidArgument,
                          "Parameter Unknown");
  }

  return param_val;
}

// sets integer parameter of LP
absl::Status LPGlopInterface::SetIntegerParameter(
    LPParameter type,  // parameter number
    int param_val      // parameter value
) {
  switch (type) {
    case LPParameter::kFromScratch:
      MiniMIPdebugMessage(
          "SetIntegerParameter: LPParameter::kFromScratch -> %d.\n", param_val);
      from_scratch_ = static_cast<bool>(param_val);
      break;
    case LPParameter::kLPInfo:
      MiniMIPdebugMessage("SetIntegerParameter: LPParameter::kLPInfo -> %d.\n",
                          param_val);
      if (param_val == 0) {
        static_cast<void>(google::SetVLOGLevel("*", google::GLOG_INFO));
        lp_info_ = false;
      } else {
        static_cast<void>(google::SetVLOGLevel("*", google::GLOG_ERROR));
        lp_info_ = true;
      }
      break;
    case LPParameter::kLPIterationLimit:
      MiniMIPdebugMessage(
          "SetIntegerParameter: LPParameter::kLPIterationLimit -> %d.\n",
          param_val);
      parameters_.set_max_number_of_iterations(param_val);
      break;
    case LPParameter::kPresolving:
      MiniMIPdebugMessage(
          "SetIntegerParameter: LPParameter::kPresolving -> %d.\n", param_val);
      parameters_.set_use_preprocessing(param_val);
      break;
    case LPParameter::kPricing:
      MiniMIPdebugMessage("SetIntegerParameter: LPParameter::kPricing -> %d.\n",
                          param_val);
      pricing_ = (LPPricing)param_val;
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
          return absl::Status(absl::StatusCode::kInvalidArgument,
                              "Parameter Unknown");
      }
      break;
#ifndef NOSCALING
    case LPParameter::kScaling:
      MiniMIPdebugMessage("SetIntegerParameter: LPParameter::kScaling -> %d.\n",
                          param_val);
      parameters_.set_use_scaling(param_val);
      break;
#endif
    case LPParameter::kThreads:
      MiniMIPdebugMessage("SetIntegerParameter: LPParameter::kThreads -> %d.\n",
                          param_val);
      num_threads_ = param_val;
      if (param_val == 0)
        parameters_.set_num_omp_threads(1);
      else
        parameters_.set_num_omp_threads(param_val);
      break;
    case LPParameter::kTiming:
      MiniMIPdebugMessage("SetIntegerParameter: LPParameter::kTiming -> %d.\n",
                          param_val);
      assert(param_val <= 2);
      timing_ = param_val;
      if (param_val == 1)
        absl::SetFlag(&FLAGS_time_limit_use_usertime, true);
      else
        absl::SetFlag(&FLAGS_time_limit_use_usertime, false);
      break;
    case LPParameter::kRandomSeed:
      MiniMIPdebugMessage(
          "SetIntegerParameter: LPParameter::kRandomSeed -> %d.\n", param_val);
      parameters_.set_random_seed(param_val);
      break;
    default:
      return absl::Status(absl::StatusCode::kInvalidArgument,
                          "Parameter Unknown");
  }

  return absl::OkStatus();
}

// gets floating point parameter of LP
absl::StatusOr<double> LPGlopInterface::GetRealParameter(
    LPParameter type  // parameter number
) const {
  double param_val;

  switch (type) {
    case LPParameter::kFeasibilityTolerance:
      param_val = parameters_.primal_feasibility_tolerance();
      MiniMIPdebugMessage(
          "GetRealParameter: LPParameter::kFeasibilityTolerance = %g.\n",
          param_val);
      break;
    case LPParameter::kDualFeasibilityTolerance:
      param_val = parameters_.dual_feasibility_tolerance();
      MiniMIPdebugMessage(
          "GetRealParameter: LPParameter::kDualFeasibilityTolerance = %g.\n",
          param_val);
      break;
    case LPParameter::kObjectiveLimit:
      if (linear_program_.IsMaximizationProblem())
        param_val = parameters_.objective_lower_limit();
      else
        param_val = parameters_.objective_upper_limit();
      MiniMIPdebugMessage(
          "GetRealParameter: LPParameter::kObjectiveLimit = %f.\n", param_val);
      break;
    case LPParameter::kLPTimeLimit:
      if (absl::GetFlag(FLAGS_time_limit_use_usertime))
        param_val = parameters_.max_time_in_seconds();
      else
        param_val = parameters_.max_deterministic_time();
      MiniMIPdebugMessage("GetRealParameter: LPParameter::kLPTimeLimit = %f.\n",
                          param_val);
      break;

    default:
      return absl::Status(absl::StatusCode::kInvalidArgument,
                          "Parameter Unknown");
  }
  return param_val;
}

// sets floating point parameter of LP
absl::Status LPGlopInterface::SetRealParameter(
    LPParameter type,  // parameter number
    double param_val   // parameter value
) {
  switch (type) {
    case LPParameter::kFeasibilityTolerance:
      MiniMIPdebugMessage(
          "SetRealParameter: LPParameter::kFeasibilityTolerance -> %g.\n",
          param_val);
      parameters_.set_primal_feasibility_tolerance(param_val);
      break;
    case LPParameter::kDualFeasibilityTolerance:
      MiniMIPdebugMessage(
          "SetRealParameter: LPParameter::kDualFeasibilityTolerance -> %g.\n",
          param_val);
      parameters_.set_dual_feasibility_tolerance(param_val);
      break;
    case LPParameter::kObjectiveLimit:
      MiniMIPdebugMessage(
          "SetRealParameter: LPParameter::kObjectiveLimit -> %f.\n", param_val);
      if (linear_program_.IsMaximizationProblem())
        parameters_.set_objective_lower_limit(param_val);
      else
        parameters_.set_objective_upper_limit(param_val);
      break;
    case LPParameter::kLPTimeLimit:
      MiniMIPdebugMessage(
          "SetRealParameter: LPParameter::kLPTimeLimit -> %f.\n", param_val);
      if (absl::GetFlag(FLAGS_time_limit_use_usertime))
        parameters_.set_max_time_in_seconds(param_val);
      else
        parameters_.set_max_deterministic_time(param_val);
      break;
    default:
      return absl::Status(absl::StatusCode::kInvalidArgument,
                          "Parameter Unknown");
  }

  return absl::OkStatus();
}

// ==========================================================================
// Numerical methods.
// ==========================================================================

// returns value treated as infinity in the LP solver
double LPGlopInterface::Infinity() const {
  return std::numeric_limits<double>::infinity();
}

// checks if given value is treated as infinity in the LP solver
bool LPGlopInterface::IsInfinity(
    double value  // value to be checked for infinity
) const {
  return value == std::numeric_limits<double>::infinity();
}

// ==========================================================================
// File interface methods.
// ==========================================================================

// reads LP from a file
absl::Status LPGlopInterface::ReadLP(const char* file_name  // file name
) {
  assert(file_name != nullptr);

  const char* filespec(file_name);
  MPModelProto proto;
  if (!ReadFileToProto(filespec, &proto)) {
    MiniMIPerrorMessage("Could not read <%s>\n", file_name);
    return absl::Status(absl::StatusCode::kInternal, "Read Errror");
  }
  linear_program_.Clear();
  MPModelProtoToLinearProgram(proto, &linear_program_);

  return absl::OkStatus();
}

// writes LP to a file
absl::Status LPGlopInterface::WriteLP(const char* file_name  // file name
) const {
  assert(file_name != nullptr);

  MPModelProto proto;
  LinearProgramToMPModelProto(linear_program_, &proto);
  const char* filespec(file_name);
  if (!WriteProtoToFile(filespec, proto,
                        operations_research::ProtoWriteFormat::kProtoText,
                        true)) {
    MiniMIPerrorMessage("Could not write <%s>\n", file_name);
    return absl::Status(absl::StatusCode::kInternal, "Write Errror");
  }

  return absl::OkStatus();
}

}  // namespace minimip
