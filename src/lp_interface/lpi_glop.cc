#include "src/lp_interface/lpi_glop.h"

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
// constructor
LPGlopInterface::LPGlopInterface() : lp_modified_since_last_solve_(true),
                                     lp_time_limit_was_reached_(false),
                                     lp_info_(false),
                                     pricing_(LPPricing::kDefault),
                                     from_scratch_(false),
                                     num_threads_(0),
                                     timing_(0),
                                     niterations_(0LL),
                                     tmp_row_(new ScatteredRow()),
                                     tmp_column_(new ScatteredColumn()) {}

// Destructor default
LPGlopInterface::~LPGlopInterface() {
  MiniMIPdebugMessage("LPGLopInterface Free\n");
}

// Modification Methods

// @name Modification Methods
// @{

// copies LP data with column matrix into LP solver
RetCode LPGlopInterface::LoadColumnLP(
  LPObjectiveSense obj_sense,           // objective sense
  int num_cols,                       // number of columns
  const std::vector<double>& objective_values, // objective function values of columns
  const std::vector<double>& lower_bounds,     // lower bounds of columns
  const std::vector<double>& upper_bounds,     // upper bounds of columns
  std::vector<std::string>& col_names,               // column names
  int num_rows,                       // number of rows
  const std::vector<double>& left_hand_sides,  // left hand sides of rows
  const std::vector<double>& right_hand_sides, // right hand sides of rows
  std::vector<std::string>& row_names,               // row names
  int num_non_zeros,                  // number of non-zero elements in the constraint matrix
  const std::vector<int>& begin_cols,       // start index of each column in row_indices- and vals-array
  const std::vector<int>& row_indices,      // row indices of constraint matrix entries
  const std::vector<double>& vals              // values of constraint matrix entries
) {

  linear_program_.Clear();
  AddRows(num_rows, left_hand_sides, right_hand_sides, row_names, 0, std::vector<int>(), std::vector<int>(), std::vector<double>());
  AddColumns(num_cols, objective_values, lower_bounds, upper_bounds, col_names, num_non_zeros, begin_cols, row_indices, vals);
  ChangeObjectiveSense(obj_sense);

  return RetCode::kOkay;
}

// adds columns to the LP
RetCode LPGlopInterface::AddColumns(
  int num_cols,                       // number of columns to be added
  const std::vector<double>& objective_values, // objective function values of new columns
  const std::vector<double>& lower_bounds,     // lower bounds of new columns
  const std::vector<double>& upper_bounds,     // upper bounds of new columns
  std::vector<std::string>& col_names,               // column names
  int num_non_zeros,                  // number of non-zero elements to be added to the constraint matrix
  const std::vector<int>& begin_cols,       // start index of each column in indices- and vals-array
  const std::vector<int>& indices,          // row indices of constraint matrix entries
  const std::vector<double>& vals              // values of constraint matrix entries
) {

  MiniMIPdebugMessage("adding %d columns with %d nonzeros.\n", num_cols, num_non_zeros);

  // @todo add names
  if (num_non_zeros > 0) {
    assert(num_cols > 0);

#ifndef NDEBUG
    // perform check that no new rows are added
    RowIndex num_rows = linear_program_.num_constraints();
    for (int j = 0; j < num_non_zeros; ++j) {
      assert(0 <= indices[j] && indices[j] < num_rows.value());
      assert(vals[j] != 0.0);
    }
#endif

    int nz = 0;
    for (int i = 0; i < num_cols; ++i) {
      const ColIndex col = linear_program_.CreateNewVariable();
      linear_program_.SetVariableBounds(col, lower_bounds[i], upper_bounds[i]);
      linear_program_.SetObjectiveCoefficient(col, objective_values[i]);
      const int end = (num_non_zeros == 0 || i == num_cols - 1) ? num_non_zeros : begin_cols[i + 1];
      while (nz < end) {
        linear_program_.SetCoefficient(RowIndex(indices[nz]), col, vals[nz]);
        ++nz;
      }
    }
    assert(nz == num_non_zeros);
  } else {
    for (int i = 0; i < num_cols; ++i) {
      const ColIndex col = linear_program_.CreateNewVariable();
      linear_program_.SetVariableBounds(col, lower_bounds[i], upper_bounds[i]);
      linear_program_.SetObjectiveCoefficient(col, objective_values[i]);
    }
  }

  lp_modified_since_last_solve_ = true;

  return RetCode::kOkay;
}

// deletes all columns in the given range from LP
RetCode LPGlopInterface::DeleteColumns(
  int first_col, // first column to be deleted
  int last_col   // last column to be deleted
) {
  assert(0 <= first_col && first_col <= last_col && last_col < linear_program_.num_variables());

  MiniMIPdebugMessage("deleting columns %d to %d.\n", first_col, last_col);

  const ColIndex num_cols = linear_program_.num_variables();
  DenseBooleanRow columns_to_delete(num_cols, false);
  for (int i = first_col; i <= last_col; ++i)
    columns_to_delete[ColIndex(i)] = true;

  linear_program_.DeleteColumns(columns_to_delete);
  lp_modified_since_last_solve_ = true;

  return RetCode::kOkay;
}

// deletes columns from MiniMIP_LP; the new position of a column must not be greater that its old position
RetCode LPGlopInterface::DeleteColumnSet(
  std::vector<bool>& deletion_status // deletion status of columns
) {

  const ColIndex num_cols = linear_program_.num_variables();
  DenseBooleanRow columns_to_delete(num_cols, false);
  int new_index = 0;
  int num_deleted_columns = 0;
  for (ColIndex col(0); col < num_cols; ++col) {
    int i = col.value();
    if (deletion_status[i] == 1) {
      columns_to_delete[col] = true;
      deletion_status[i] = -1;
      ++num_deleted_columns;
    } else
      deletion_status[i] = new_index++;
  }
  MiniMIPdebugMessage("DeleteColumnset: deleting %d columns.\n", num_deleted_columns);
  linear_program_.DeleteColumns(columns_to_delete);
  lp_modified_since_last_solve_ = true;

  return RetCode::kOkay;
}

// adds rows to the LP
RetCode LPGlopInterface::AddRows(
  int num_rows,                       // number of rows to be added
  const std::vector<double>& left_hand_sides,  // left hand sides of new rows
  const std::vector<double>& right_hand_sides, // right hand sides of new rows
  std::vector<std::string>& row_names,               // row names
  int num_non_zeros,                  // number of non-zero elements to be added to the constraint matrix
  const std::vector<int>& begin_rows,       // start index of each row in indices- and vals-array
  const std::vector<int>& indices,          // column indices of constraint matrix entries
  const std::vector<double>& vals              // values of constraint matrix entries
) {

  MiniMIPdebugMessage("adding %d rows with %d nonzeros.\n", num_rows, num_non_zeros);

  // @todo add names
  if (num_non_zeros > 0) {
    assert(num_rows > 0);

#ifndef NDEBUG
    // perform check that no new columns are added - this is likely to be a mistake
    const ColIndex num_cols = linear_program_.num_variables();
    for (int j = 0; j < num_non_zeros; ++j) {
      assert(vals[j] != 0.0);
      assert(0 <= indices[j] && indices[j] < num_cols.value());
    }
#endif

    int nz = 0;
    for (int i = 0; i < num_rows; ++i) {
      const RowIndex row = linear_program_.CreateNewConstraint();
      linear_program_.SetConstraintBounds(row, left_hand_sides[i], right_hand_sides[i]);
      const int end = (num_non_zeros == 0 || i == num_rows - 1) ? num_non_zeros : begin_rows[i + 1];
      while (nz < end) {
        linear_program_.SetCoefficient(row, ColIndex(indices[nz]), vals[nz]);
        ++nz;
      }
    }
    assert(nz == num_non_zeros);
  } else {
    for (int i = 0; i < num_rows; ++i) {
      const RowIndex row = linear_program_.CreateNewConstraint();
      linear_program_.SetConstraintBounds(row, left_hand_sides[i], right_hand_sides[i]);
    }
  }

  lp_modified_since_last_solve_ = true;

  return RetCode::kOkay;
}

// delete rows from LP and update the current basis
void LPGlopInterface::DeleteRowsAndUpdateCurrentBasis(
  const DenseBooleanColumn& rows_to_delete // array to mark rows that should be deleted
) {
  const RowIndex num_rows = linear_program_.num_constraints();
  const ColIndex num_cols = linear_program_.num_variables();

  // try to repair basis status if problem size has not changed before
  BasisState state = solver_.GetState();
  if (state.statuses.size() == num_cols.value() + num_rows.value()) {
    // Shift the status of the non-deleted rows. Note that if the deleted rows where part of the basis (i.e., constraint
    // not tight), then we should be left with a correct basis afterward. This should be the most common use case in MiniMIP.
    ColIndex new_size = num_cols;
    for (RowIndex row(0); row < num_rows; ++row) {
      if (rows_to_delete[row])
        continue;
      state.statuses[new_size++] = state.statuses[num_cols + RowToColIndex(row)];
    }
    state.statuses.resize(new_size);
    solver_.LoadStateForNextSolve(state);
  }

  linear_program_.DeleteRows(rows_to_delete);
  lp_modified_since_last_solve_ = true;
}

// deletes all rows in the given range from LP
RetCode LPGlopInterface::DeleteRows(
  int first_row, // first row to be deleted
  int last_row   // last row to be deleted
) {
  assert(0 <= first_row && first_row <= last_row && last_row < linear_program_.num_constraints());

  const RowIndex num_rows = linear_program_.num_constraints();
  DenseBooleanColumn rows_to_delete(num_rows, false);
  for (int i = first_row; i <= last_row; ++i)
    rows_to_delete[RowIndex(i)] = true;

  MiniMIPdebugMessage("deleting rows %d to %d.\n", first_row, last_row);
  DeleteRowsAndUpdateCurrentBasis(rows_to_delete);

  return RetCode::kOkay;
}

// deletes rows from LP; the new position of a row must not be greater that its old position
RetCode LPGlopInterface::DeleteRowSet(
  std::vector<bool>& deletion_status // deletion status of rows
) {
  const RowIndex num_rows = linear_program_.num_constraints();
  DenseBooleanColumn rows_to_delete(num_rows, false);
  int new_index = 0;
  int num_deleted_rows = 0;
  for (RowIndex row(0); row < num_rows; ++row) {
    int i = row.value();
    if (deletion_status[i] == 1) {
      rows_to_delete[row] = true;
      deletion_status[i] = -1;
      ++num_deleted_rows;
    } else
      deletion_status[i] = new_index++;
  }

  MiniMIPdebugMessage("DeleteRowSet: deleting %d rows.\n", num_deleted_rows);
  DeleteRowsAndUpdateCurrentBasis(rows_to_delete);

  return RetCode::kOkay;
}

// clears the whole LP
RetCode LPGlopInterface::Clear() {

  MiniMIPdebugMessage("Clear\n");

  linear_program_.Clear();
  lp_modified_since_last_solve_ = true;

  return RetCode::kOkay;
}

// clears current LPi state (like basis information) of the solver
RetCode LPGlopInterface::ClearState() {

  solver_.ClearStateForNextSolve();

  return RetCode::kOkay;
}

// changes lower and upper bounds of columns
RetCode LPGlopInterface::ChangeBounds(
  int num_cols,                   // number of columns to change bounds for
  const std::vector<int>& indices,      // column indices
  const std::vector<double>& lower_bounds, // values for the new lower bounds
  const std::vector<double>& upper_bounds  // values for the new upper bounds
) {

  MiniMIPdebugMessage("changing %d bounds.\n", num_cols);

  for (int i = 0; i < num_cols; ++i) {
    MiniMIPdebugMessage("  col %d: [%g,%g]\n", indices[i], lower_bounds[i], upper_bounds[i]);

    if (IsInfinity(lower_bounds[i])) {
      MiniMIPerrorMessage("LP Error: fixing lower bound for variable %d to infinity.\n", indices[i]);
      return RetCode::kLPError;
    }
    if (IsInfinity(-upper_bounds[i])) {
      MiniMIPerrorMessage("LP Error: fixing upper bound for variable %d to -infinity.\n", indices[i]);
      return RetCode::kLPError;
    }

    linear_program_.SetVariableBounds(ColIndex(indices[i]), lower_bounds[i], upper_bounds[i]);
  }
  lp_modified_since_last_solve_ = true;

  return RetCode::kOkay;
}

// changes left and right hand sides of rows
RetCode LPGlopInterface::ChangeSides(
  int num_rows,                      // number of rows to change sides for
  const std::vector<int>& indices,         // row indices
  const std::vector<double>& left_hand_sides, // new values for left hand sides
  const std::vector<double>& right_hand_sides // new values for right hand sides
) {

  MiniMIPdebugMessage("changing %d sides\n", num_rows);

  for (int i = 0; i < num_rows; ++i)
    linear_program_.SetConstraintBounds(RowIndex(indices[i]), left_hand_sides[i], right_hand_sides[i]);

  lp_modified_since_last_solve_ = true;

  return RetCode::kOkay;
}

// changes the objective sense
RetCode LPGlopInterface::ChangeObjectiveSense(
  LPObjectiveSense obj_sense // new objective sense
) {

  switch (obj_sense) {
    case LPObjectiveSense::kMaximize:
      MiniMIPdebugMessage("changing objective sense to MAXIMIZE\n");
      linear_program_.SetMaximizationProblem(true);
      break;
    case LPObjectiveSense::kMinimize:
      MiniMIPdebugMessage("changing objective sense to MINIMIZE\n");
      linear_program_.SetMaximizationProblem(false);
      break;
  }
  lp_modified_since_last_solve_ = true;

  return RetCode::kOkay;
}

// changes objective values of columns in the LP
RetCode LPGlopInterface::ChangeObjective(
  int num_cols,                  // number of columns to change objective value for
  const std::vector<int>& indices,     // column indices to change objective value for
  const std::vector<double>& new_obj_vals // new objective values for columns
) {

  MiniMIPdebugMessage("changing %d objective values\n", num_cols);

  for (int i = 0; i < num_cols; ++i)
    linear_program_.SetObjectiveCoefficient(ColIndex(indices[i]), new_obj_vals[i]);

  lp_modified_since_last_solve_ = true;

  return RetCode::kOkay;
}


// Data Accessing Methods

// @name Data Accessing Methods
// @{

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

  return linear_program_.IsMaximizationProblem() ? LPObjectiveSense::kMaximize : LPObjectiveSense::kMinimize;
}

// gets the number of nonzero elements in the LP constraint matrix
int LPGlopInterface::GetNumberOfNonZeros() const {

  MiniMIPdebugMessage("getting number of non-zeros.\n");

  return static_cast<int>(linear_program_.num_entries().value());
}

// gets columns from LP problem object
//
// Either both, lb and ub, have to be NULL, or both have to be non-NULL,
// either num_non_zeros, begin_cols, indices, and val have to be NULL, or all of them have to be non-NULL.
RetCode LPGlopInterface::GetColumns(
  int first_col,            // first column to get from LP
  int last_col,             // last column to get from LP
  std::vector<double>& lower_bounds, // array to store the lower bound vector
  std::vector<double>& upper_bounds, // array to store the upper bound vector
  int& num_non_zeros,       // store the number of non-zero elements
  std::vector<int>& begin_cols,   // array to store start index of each column in indices- and vals-array
  std::vector<int>& indices,      // array to store row indices of constraint matrix entries
  std::vector<double>& vals          // array to store values of constraint matrix entries
) const {
  assert(0 <= first_col && first_col <= last_col && last_col < linear_program_.num_variables());

  const DenseRow& tmp_lower_bound = linear_program_.variable_lower_bounds();
  const DenseRow& tmp_upper_bound = linear_program_.variable_upper_bounds();

  if (num_non_zeros >= 0) {
    num_non_zeros = 0;
    int index = 0;
    for (ColIndex col(first_col); col <= ColIndex(last_col); ++col, ++index) {
      lower_bounds[index] = tmp_lower_bound[col];
      upper_bounds[index] = tmp_upper_bound[col];

      begin_cols[index] = num_non_zeros;
      const SparseColumn& column = linear_program_.GetSparseColumn(col);
      for (const SparseColumn::Entry& entry : column) {
        const RowIndex row = entry.row();
        indices[num_non_zeros] = row.value();
        vals[num_non_zeros] = entry.coefficient();
        ++(num_non_zeros);
      }
    }
  } else {
    int index = 0;
    for (ColIndex col(first_col); col <= ColIndex(last_col); ++col, ++index) {
      lower_bounds[index] = tmp_lower_bound[col];
      upper_bounds[index] = tmp_upper_bound[col];
    }
  }

  return RetCode::kOkay;
}

// gets rows from LP problem object
//
// Either both, left_hand_side and right_hand_side, have to be NULL, or both have to be non-NULL,
// either num_non_zeros, begin_rows, indices, and val have to be NULL, or all of them have to be non-NULL.
RetCode LPGlopInterface::GetRows(
  int first_row,                // first row to get from LP
  int last_row,                 // last row to get from LP
  std::vector<double>& left_hand_sides,  // array to store left hand side vector
  std::vector<double>& right_hand_sides, // array to store right hand side vector
  int& num_non_zeros,           // store the number of non-zero elements
  std::vector<int>& begin_rows,       // array to store start index of each row in indices- and vals-array
  std::vector<int>& indices,          // array to store column indices of constraint matrix entries
  std::vector<double>& vals              // array to store values of constraint matrix entries
) const {
  assert(0 <= first_row && first_row <= last_row && last_row < linear_program_.num_constraints());

  const DenseColumn& tmplhs = linear_program_.constraint_lower_bounds();
  const DenseColumn& tmprhs = linear_program_.constraint_upper_bounds();

  const SparseMatrix& matrixtrans = linear_program_.GetTransposeSparseMatrix();

  num_non_zeros = 0;
  int index = 0;
  for (RowIndex row(first_row); row <= RowIndex(last_row); ++row, ++index) {
    left_hand_sides[index] = tmplhs[row];
    right_hand_sides[index] = tmprhs[row];

    begin_rows[index] = num_non_zeros;
    const SparseColumn& column = matrixtrans.column(ColIndex(row.value()));
    int count = 0;
    for (const SparseColumn::Entry& entry : column) {
      count++;
    }
    for (const SparseColumn::Entry& entry : column) {
      const RowIndex rowidx = entry.row();
      indices[num_non_zeros] = rowidx.value();
      vals[num_non_zeros] = entry.coefficient();
      ++num_non_zeros;
    }
  }

  return RetCode::kOkay;
}

// gets objective coefficients from LP problem object
RetCode LPGlopInterface::GetObjective(
  int first_col,         // first column to get objective coefficient for
  int last_col,          // last column to get objective coefficient for
  std::vector<double>& obj_coeffs // array to store objective coefficients
) const {
  assert(first_col <= last_col);

  MiniMIPdebugMessage("getting objective values %d to %d\n", first_col, last_col);

  int index = 0;
  for (ColIndex col(first_col); col <= ColIndex(last_col); ++col) {
    obj_coeffs[index] = linear_program_.objective_coefficients()[col];
    ++index;
  }

  return RetCode::kOkay;
}

// gets current bounds from LP problem object
RetCode LPGlopInterface::GetBounds(
  int first_col,            // first column to get bounds for
  int last_col,             // last column to get bounds for
  std::vector<double>& lower_bounds, // array to store lower bound values
  std::vector<double>& upper_bounds  // array to store upper bound values
) const {
  assert(first_col <= last_col);

  MiniMIPdebugMessage("getting bounds %d to %d\n", first_col, last_col);

  int index = 0;
  for (ColIndex col(first_col); col <= ColIndex(last_col); ++col) {
    lower_bounds[index] = linear_program_.variable_lower_bounds()[col];

    upper_bounds[index] = linear_program_.variable_upper_bounds()[col];

    ++index;
  }

  return RetCode::kOkay;
}

// gets current row sides from LP problem object
RetCode LPGlopInterface::GetSides(
  int first_row,               // first row to get sides for
  int last_row,                // last row to get sides for
  std::vector<double>& left_hand_sides, // array to store left hand side values
  std::vector<double>& right_hand_sides // array to store right hand side values
) const {
  assert(first_row <= last_row);

  MiniMIPdebugMessage("getting row sides %d to %d\n", first_row, last_row);

  int index = 0;
  for (RowIndex row(first_row); row <= RowIndex(last_row); ++row) {
    left_hand_sides[index] = linear_program_.constraint_lower_bounds()[row];

    right_hand_sides[index] = linear_program_.constraint_upper_bounds()[row];

    ++index;
  }

  return RetCode::kOkay;
}

// gets a single coefficient
RetCode LPGlopInterface::GetCoefficient(
  int row,         // row number of coefficient
  int col_index, // column number of coefficient
  double& val       // array to store the value of the coefficient
) const {

  // quite slow method: possibly needs linear time if matrix is not sorted
  const SparseMatrix& matrix = linear_program_.GetSparseMatrix();
  val = matrix.LookUpValue(RowIndex(row), ColIndex(col_index));

  return RetCode::kOkay;
}

// @}

// Solving Methods

// @name Solving Methods
// @{

// update scaled linear program
void LPGlopInterface::updateScaledLP() {
  if (!lp_modified_since_last_solve_)
    return;

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
    unscaledsol[col] = scaler_.UnscaleVariableValue(col, solver_.GetVariableValue(col));

  // if the solution is not feasible w.r.t. absolute tolerances, try to fix it in the unscaled problem
  const double feastol = parameters_.primal_feasibility_tolerance();
  return linear_program_.SolutionIsLPFeasible(unscaledsol, feastol);

#elif UNSCALEDFEAS_CHECK == 2
  const double feastol = parameters_.primal_feasibility_tolerance();

  // check bounds of unscaled solution
  const ColIndex num_cols = linear_program_.num_variables();
  for (ColIndex col = ColIndex(0); col < num_cols; ++col) {
    const Fractional val = scaler_.UnscaleVariableValue(col, solver_.GetVariableValue(col));
    const Fractional lower_bound = linear_program_.variable_lower_bounds()[col];
    if (val < lower_bound - feastol)
      return false;
    const Fractional upper_bound = linear_program_.variable_upper_bounds()[col];
    if (val > upper_bound + feastol)
      return false;
  }

  // check activities of unscaled solution
  const RowIndex num_rows = linear_program_.num_constraints();
  for (RowIndex row(0); row < num_rows; ++row) {
    const Fractional val = scaler_.UnscaleConstraintActivity(row, solver_.GetConstraintActivity(row));
    const Fractional left_hand_side = linear_program_.constraint_lower_bounds()[row];
    if (val < left_hand_side - feastol)
      return false;
    const Fractional right_hand_side = linear_program_.constraint_upper_bounds()[row];
    if (val > right_hand_side + feastol)
      return false;
  }
#endif

  return true;
}

// common function between the two LPI Solve() functions
RetCode LPGlopInterface::SolveInternal(
  bool recursive,                        // Is this a recursive call?
  std::unique_ptr<TimeLimit>& time_limit // time limit
) {

  updateScaledLP();

  solver_.SetParameters(parameters_);
  lp_time_limit_was_reached_ = false;

  // possibly ignore warm start information for next solve
  if (from_scratch_)
    solver_.ClearStateForNextSolve();

  if (!solver_.Solve(scaled_lp_, time_limit.get()).ok()) {
    return RetCode::kLPError;
  }
  lp_time_limit_was_reached_ = time_limit->LimitReached();
  if (recursive)
    niterations_ += (long long) solver_.GetNumberOfIterations();
  else
    niterations_ = (long long) solver_.GetNumberOfIterations();

  MiniMIPdebugMessage("status=%s  obj=%f  iterations=%ld.\n", GetProblemStatusString(solver_.GetProblemStatus()).c_str(),
                      solver_.GetObjectiveValue(), solver_.GetNumberOfIterations());

  const ProblemStatus status = solver_.GetProblemStatus();
  if ((status == ProblemStatus::PRIMAL_FEASIBLE || status == ProblemStatus::OPTIMAL) && parameters_.use_scaling()) {
    if (!checkUnscaledPrimalFeasibility()) {
      MiniMIPdebugMessage("Solution not feasible w.r.t. absolute tolerance %g -> reoptimize.\n", parameters_.primal_feasibility_tolerance());

      // Re-solve without scaling to try to fix the infeasibility.
      parameters_.set_use_scaling(false);
      lp_modified_since_last_solve_ = true;
      SolveInternal(true, time_limit); // inherit time limit, so used time is not reset; do not change iteration limit for resolve
      parameters_.set_use_scaling(true);
    }
  }

  lp_modified_since_last_solve_ = false;

  return RetCode::kOkay;
}

// calls primal simplex to solve the LP
RetCode LPGlopInterface::SolvePrimal() {

  MiniMIPdebugMessage("SolvePrimal: %d rows, %d cols.\n", linear_program_.num_constraints().value(), linear_program_.num_variables().value());
  std::unique_ptr<TimeLimit> time_limit = TimeLimit::FromParameters(parameters_);
  niterations_ = 0;

  parameters_.set_use_dual_simplex(false);
  return SolveInternal(false, time_limit);
}

// calls dual simplex to solve the LP
RetCode LPGlopInterface::SolveDual() {

  MiniMIPdebugMessage("SolveDual: %d rows, %d cols.\n", linear_program_.num_constraints().value(), linear_program_.num_variables().value());
  std::unique_ptr<TimeLimit> time_limit = TimeLimit::FromParameters(parameters_);
  niterations_ = 0;

  parameters_.set_use_dual_simplex(true);
  return SolveInternal(false, time_limit);
}

// start strong branching
RetCode LPGlopInterface::StartStrongbranch() { 

  updateScaledLP();

  // @todo Save state and do all the branching from there.
  return RetCode::kOkay;
}

// end strong branching
RetCode LPGlopInterface::EndStrongbranch() { 

  // @todo Restore the saved state.
  return RetCode::kOkay;
}
// determine whether the dual bound is valid
bool LPGlopInterface::IsDualBoundValid(
  ProblemStatus status // status to be checked
) const {
  return status == ProblemStatus::OPTIMAL || status == ProblemStatus::DUAL_FEASIBLE || status == ProblemStatus::DUAL_UNBOUNDED;
}

// performs strong branching iterations
RetCode LPGlopInterface::strongbranch(
  int col_index,               // column to apply strong branching on
  double primal_sol,              // fractional current primal solution value of column
  int iteration_limit,           // iteration limit for strong branchings
  double& dual_bound_down_branch, // stores dual bound after branching column down
  double& dual_bound_up_branch,   // stores dual bound after branching column up
  bool& down_valid,                // stores whether the returned down value is a valid dual bound;
                                   // otherwise, it can only be used as an estimate value
  bool& up_valid,                  // stores whether the returned up value is a valid dual bound;
                                   // otherwise, it can only be used as an estimate value
  int& iterations                // stores total number of strong branching iterations, or -1;
) {

  MiniMIPdebugMessage("calling strongbranching on variable %d (%d iterations)\n", col_index, iteration_limit);

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

  std::unique_ptr<TimeLimit> time_limit = TimeLimit::FromParameters(parameters_);

  // Down branch.
  const Fractional new_upper_bound = EPSCEIL(value - 1.0, eps);
  if (new_upper_bound >= lower_bound - 0.5) {
    scaled_lp_.SetVariableBounds(col, lower_bound, new_upper_bound);

    if (solver_.Solve(scaled_lp_, time_limit.get()).ok()) {
      num_iterations += static_cast<int>(solver_.GetNumberOfIterations());
      dual_bound_down_branch = solver_.GetObjectiveValue();
      down_valid = IsDualBoundValid(solver_.GetProblemStatus());

      MiniMIPdebugMessage("dual_bound_down_branch: iteration_limit=%d col=%d [%f,%f] obj=%f status=%d iterations=%ld.\n", iteration_limit, col_index, lower_bound, EPSCEIL(value - 1.0, eps),
                          solver_.GetObjectiveValue(), static_cast<int>(solver_.GetProblemStatus()), solver_.GetNumberOfIterations());
    } else {
      MiniMIPerrorMessage("error during solve");
      dual_bound_down_branch = 0.0;
      down_valid = false;
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
      up_valid = IsDualBoundValid(solver_.GetProblemStatus());

      MiniMIPdebugMessage("dual_bound_up_branch: iteration_limit=%d col=%d [%f,%f] obj=%f status=%d iterations=%ld.\n", iteration_limit, col_index, EPSFLOOR(value + 1.0, eps), upper_bound,
                          solver_.GetObjectiveValue(), static_cast<int>(solver_.GetProblemStatus()), solver_.GetNumberOfIterations());
    } else {
      MiniMIPerrorMessage("error during solve");
      dual_bound_up_branch = 0.0;
      up_valid = false;
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
  if (iterations > 0)
    iterations = num_iterations;

  return RetCode::kOkay;
}

// performs strong branching iterations on one @b fractional candidate
RetCode LPGlopInterface::StrongbranchFractionalValue(
  int col,                     // column to apply strong branching on
  double primal_sol,              // fractional current primal solution value of column
  int iteration_limit,           // iteration limit for strong branchings
  double& dual_bound_down_branch, // stores dual bound after branching column down
  double& dual_bound_up_branch,   // stores dual bound after branching column up
  bool& down_valid,                // whether the returned down value is a valid dual bound; otherwise, it can only be used as an estimate value
  bool& up_valid,                  // whether the returned up value is a valid dual bound; otherwise, it can only be used as an estimate value
  int& iterations                // stores total number of strong branching iterations
) {

  MiniMIPdebugMessage("calling strong branching on fractional variable %d (%d iterations)\n", col, iteration_limit);

  MINIMIP_CALL(strongbranch(col, primal_sol, iteration_limit, dual_bound_down_branch, dual_bound_up_branch, down_valid, up_valid, iterations));

  return RetCode::kOkay;
}

// performs strong branching iterations on given @b fractional candidates
RetCode LPGlopInterface::StrongbranchFractionalValues(
  std::vector<int>& cols,                       // columns to apply strong branching on
  int num_cols,                         // number of columns
  std::vector<double>& primal_sols,              // fractional current primal solution values of columns
  int iteration_limit,                  // iteration limit for strong branchings
  std::vector<double>& dual_bound_down_branches, // stores dual bounds after branching columns down
  std::vector<double>& dual_bound_up_branches,   // stores dual bounds after branching columns up
  std::vector<bool>& down_valids,                 // stores whether the returned down values are valid dual bounds;
                                          // otherwise, they can only be used as an estimate values
  std::vector<bool>& up_valids,                   // stores whether the returned up values are a valid dual bounds;
                                          // otherwise, they can only be used as an estimate values
  int& iterations                       // stores total number of strong branching iterations
) {
  return RetCode::kNotImplemented;
}

// performs strong branching iterations on one candidate with @b integral value
RetCode LPGlopInterface::StrongbranchIntegerValue(
  int col,                     // column to apply strong branching on
  double primal_sol,              // current integral primal solution value of column
  int iteration_limit,           // iteration limit for strong branchings
  double& dual_bound_down_branch, // stores dual bound after branching column down
  double& dual_bound_up_branch,   // stores dual bound after branching column up
  bool& down_valid,                // stores whether the returned down value is a valid dual bound;
                                   // otherwise, it can only be used as an estimate value
  bool& up_valid,                  // stores whether the returned up value is a valid dual bound;
                                   // otherwise, it can only be used as an estimate value
  int& iterations                // stores total number of strong branching iterations
) {
  MiniMIPdebugMessage("calling strong branching on integer variable %d (%d iterations)\n", col, iteration_limit);

  MINIMIP_CALL(strongbranch(col, primal_sol, iteration_limit, dual_bound_down_branch, dual_bound_up_branch, down_valid, up_valid, iterations));

  return RetCode::kOkay;
}

// performs strong branching iterations on given candidates with @b integral values
RetCode LPGlopInterface::StrongbranchIntegerValues(
  std::vector<int>& cols,                       // columns to apply strong branching on
  int num_cols,                         // number of columns
  std::vector<double>& primal_sols,              // current integral primal solution values of columns
  int iteration_limit,                  // iteration limit for strong branchings
  std::vector<double>& dual_bound_down_branches, // stores dual bounds after branching columns down
  std::vector<double>& dual_bound_up_branches,   // stores dual bounds after branching columns up
  std::vector<bool>& down_valids,                 // stores whether the returned down values are valid dual bounds;
                                          // otherwise, they can only be used as an estimate values
  std::vector<bool>& up_valids,                   // stores whether the returned up values are a valid dual bounds;
                                          // otherwise, they can only be used as an estimate values
  int& iterations                       // stores total number of strong branching iterations
) {
  return RetCode::kNotImplemented;
}

// @}

// Solution Information Methods

// @name Solution Information Methods
// @{

// returns whether a solve method was called after the last modification of the LP
bool LPGlopInterface::IsSolved() const {

  // @todo Track this to avoid uneeded resolving.
  return (!lp_modified_since_last_solve_);
}

// returns true if LP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
// this does not necessarily mean that the solver knows and can return the primal ray
bool LPGlopInterface::ExistsPrimalRay() const {

  return solver_.GetProblemStatus() == ProblemStatus::PRIMAL_UNBOUNDED;
}

// returns true if LP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
// and the solver knows and can return the primal ray
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

  return status == ProblemStatus::DUAL_UNBOUNDED || status == ProblemStatus::PRIMAL_INFEASIBLE;
}

// returns true if LP is proven to be primal feasible
bool LPGlopInterface::IsPrimalFeasible() const {
  const ProblemStatus status = solver_.GetProblemStatus();

  return status == ProblemStatus::PRIMAL_FEASIBLE || status == ProblemStatus::OPTIMAL;
}

// returns true if LP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
// this does not necessarily mean that the solver knows and can return the dual ray
bool LPGlopInterface::ExistsDualRay() const {
  const ProblemStatus status = solver_.GetProblemStatus();

  return status == ProblemStatus::DUAL_UNBOUNDED;
}

// returns true if LP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
// and the solver knows and can return the dual ray
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
  return status == ProblemStatus::PRIMAL_UNBOUNDED || status == ProblemStatus::DUAL_INFEASIBLE;
}

// returns true if LP is proven to be dual feasible
bool LPGlopInterface::IsDualFeasible() const {
  const ProblemStatus status = solver_.GetProblemStatus();

  return status == ProblemStatus::DUAL_FEASIBLE || status == ProblemStatus::OPTIMAL;
}

// returns true if LP was solved to optimality
bool LPGlopInterface::IsOptimal() const {

  return solver_.GetProblemStatus() == ProblemStatus::OPTIMAL;
}

// returns true if current LP solution is stable
//
// This function should return true if the solution is reliable, i.e., feasible and optimal (or proven
// infeasible/unbounded) with respect to the original problem. The optimality status might be with respect to a scaled
// version of the problem, but the solution might not be feasible to the unscaled original problem; in this case,
// IsStable() should return false.
bool LPGlopInterface::IsStable() const {
  // For correctness, we need to report "unstable" if Glop was not able to prove optimality because of numerical
  // issues. Currently, Glop still reports primal/dual feasible if at the end, one status is within the tolerance but not
  // the other.
  const ProblemStatus status = solver_.GetProblemStatus();
  if ((status == ProblemStatus::PRIMAL_FEASIBLE || status == ProblemStatus::DUAL_FEASIBLE) && !ObjectiveLimitIsExceeded() && !IterationLimitIsExceeded() && !TimeLimitIsExceeded()) {
    MiniMIPdebugMessage("OPTIMAL not reached and no limit: unstable.\n");
    return false;
  }

  if (status == ProblemStatus::ABNORMAL || status == ProblemStatus::INVALID_PROBLEM || status == ProblemStatus::IMPRECISE)
    return false;
  return true;
} // @TODO: Case that neither if happens?

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
RetCode LPGlopInterface::GetObjectiveValue(
  double& obj_val // stores the objective value
) {

  obj_val = solver_.GetObjectiveValue();

  return RetCode::kOkay;
}

// gets primal and dual solution vectors for feasible LPs
//
// Before calling this function, the caller must ensure that the LP has been solved to optimality, i.e., that
// IsOptimal() returns true.
RetCode LPGlopInterface::GetSolution(
  double& obj_val,          // stores the objective value
  std::vector<double>& primal_sol,  // primal solution vector
  std::vector<double>& dual_sol,    // dual solution vector
  std::vector<double>& activity,    // row activity vector
  std::vector<double>& reduced_cost // reduced cost vector
) const {

  MiniMIPdebugMessage("GetSolution\n");
  obj_val = solver_.GetObjectiveValue();

  const ColIndex num_cols = linear_program_.num_variables();
  for (ColIndex col(0); col < num_cols; ++col) {
    int i = col.value();

    primal_sol[i] = scaler_.UnscaleVariableValue(col, solver_.GetVariableValue(col));

    reduced_cost[i] = scaler_.UnscaleReducedCost(col, solver_.GetReducedCost(col));
  }

  const RowIndex num_rows = linear_program_.num_constraints();
  for (RowIndex row(0); row < num_rows; ++row) {
    int j = row.value();

    dual_sol[j] = scaler_.UnscaleDualValue(row, solver_.GetDualValue(row));

    activity[j] = scaler_.UnscaleConstraintActivity(row, solver_.GetConstraintActivity(row));
  }

  return RetCode::kOkay;
}

// gets primal ray for unbounded LPs
RetCode LPGlopInterface::GetPrimalRay(
  std::vector<double>& primal_ray // primal ray
) const {

  MiniMIPdebugMessage("GetPrimalRay\n");

  const ColIndex num_cols = linear_program_.num_variables();
  const DenseRow& primal_ray_solver = solver_.GetPrimalRay();
  for (ColIndex col(0); col < num_cols; ++col)
    primal_ray[col.value()] = scaler_.UnscaleVariableValue(col, primal_ray_solver[col]);

  return RetCode::kOkay;
}

// gets dual Farkas proof for infeasibility
RetCode LPGlopInterface::GetDualFarkasMultiplier(
  std::vector<double>& dual_farkas_multiplier // dual Farkas row multipliers
) const {

  MiniMIPdebugMessage("GetDualFarkasMultiplier\n");

  const RowIndex num_rows = linear_program_.num_constraints();
  const DenseColumn& dual_ray = solver_.GetDualRay();
  for (RowIndex row(0); row < num_rows; ++row)
    dual_farkas_multiplier[row.value()] = -scaler_.UnscaleDualValue(row, dual_ray[row]); // reverse sign

  return RetCode::kOkay;
}

// gets the number of LP iterations of the last solve call
RetCode LPGlopInterface::GetIterations(
  int& iterations // number of iterations of the last solve call
) const {

  iterations = static_cast<int>(niterations_);

  return RetCode::kOkay;
}

// LP Basis Methods

// @name LP Basis Methods
// @{

// convert Glop variable basis status to MiniMIP status
LPBasisStatus LPGlopInterface::ConvertGlopVariableStatus(
  VariableStatus status,  // variable status
  Fractional reduced_cost // reduced cost of variable
) const {
  switch (status) {
    case VariableStatus::BASIC:
      return LPBasisStatus::kBasic;
    case VariableStatus::AT_UPPER_BOUND:
      return LPBasisStatus::kUpper;
    case VariableStatus::AT_LOWER_BOUND:
      return LPBasisStatus::kLower;
    case VariableStatus::FREE:
      return LPBasisStatus::kFree;
    case VariableStatus::FIXED_VALUE:
      return reduced_cost > 0.0 ? LPBasisStatus::kLower : LPBasisStatus::kUpper;
    default:
      MiniMIPerrorMessage("invalid Glop basis status.\n");
      std::abort();
  }
}

// convert Glop constraint basis status to MiniMIP status
LPBasisStatus LPGlopInterface::ConvertGlopConstraintStatus(
  ConstraintStatus status, // constraint status
  Fractional dual_value    // dual variable value
) const {
  switch (status) {
    case ConstraintStatus::BASIC:
      return LPBasisStatus::kBasic;
    case ConstraintStatus::AT_UPPER_BOUND:
      return LPBasisStatus::kUpper;
    case ConstraintStatus::AT_LOWER_BOUND:
      return LPBasisStatus::kLower;
    case ConstraintStatus::FREE:
      return LPBasisStatus::kFree;
    case ConstraintStatus::FIXED_VALUE:
      return dual_value > 0.0 ? LPBasisStatus::kLower : LPBasisStatus::kUpper;
    default:
      MiniMIPerrorMessage("invalid Glop basis status.\n");
      std::abort();
  }
}

// Convert MiniMIP variable status to Glop status
VariableStatus LPGlopInterface::ConvertMiniMIPVariableStatus(
  LPBasisStatus status // MiniMIP variable status
) const {
  switch (status) {
    case LPBasisStatus::kBasic:
      return VariableStatus::BASIC;
    case LPBasisStatus::kUpper:
      return VariableStatus::AT_UPPER_BOUND;
    case LPBasisStatus::kLower:
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

// Convert a MiniMIP constraint status to its corresponding Glop slack VariableStatus.
//
// Note that we swap the upper/lower bounds.
VariableStatus LPGlopInterface::ConvertMiniMIPConstraintStatusToSlackStatus(
  LPBasisStatus status // MiniMIP constraint status
) const {
  switch (status) {
    case LPBasisStatus::kBasic:
      return VariableStatus::BASIC;
    case LPBasisStatus::kUpper:
      return VariableStatus::AT_LOWER_BOUND;
    case LPBasisStatus::kLower:
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

// gets current basis status for columns and rows
RetCode LPGlopInterface::GetBase(
  std::vector<LPBasisStatus>& column_basis_status, // array to store column basis status, or NULL
  std::vector<LPBasisStatus>& row_basis_status     // array to store row basis status, or NULL
) const {
  MiniMIPdebugMessage("GetBase\n");

  assert(solver_.GetProblemStatus() == ProblemStatus::OPTIMAL);
  const ColIndex num_cols = linear_program_.num_variables();
  for (ColIndex col(0); col < num_cols; ++col) {
    int i = col.value();
    column_basis_status[i] = (LPBasisStatus) ConvertGlopVariableStatus(solver_.GetVariableStatus(col), solver_.GetReducedCost(col));
  }

  const RowIndex num_rows = linear_program_.num_constraints();
  for (RowIndex row(0); row < num_rows; ++row) {
    int i = row.value();
    row_basis_status[i] = (LPBasisStatus) ConvertGlopConstraintStatus(solver_.GetConstraintStatus(row), solver_.GetDualValue(row));
  }

  return RetCode::kOkay;
}

// sets current basis status for columns and rows
RetCode LPGlopInterface::SetBase(
  const std::vector<LPBasisStatus>& column_basis_status, // array with column basis status
  const std::vector<LPBasisStatus>& row_basis_status     // array with row basis status
) {

  const ColIndex num_cols = linear_program_.num_variables();
  const RowIndex num_rows = linear_program_.num_constraints();

  MiniMIPdebugMessage("SetBase\n");

  BasisState state;
  state.statuses.reserve(ColIndex(num_cols.value() + num_rows.value()));

  for (ColIndex col(0); col < num_cols; ++col)
    state.statuses[col] = ConvertMiniMIPVariableStatus(column_basis_status[col.value()]);

  for (RowIndex row(0); row < num_rows; ++row)
    state.statuses[num_cols + RowToColIndex(row)] = ConvertMiniMIPConstraintStatusToSlackStatus(column_basis_status[row.value()]);

  solver_.LoadStateForNextSolve(state);

  return RetCode::kOkay;
}

// returns the indices of the basic columns and rows; basic column n gives value n, basic row m gives value -1-m
RetCode LPGlopInterface::GetBasisIndices(
  std::vector<int>& basis_indices // array to store basis indices ready to keep number of rows entries
) const {
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

  return RetCode::kOkay;
}

// get row of inverse basis matrix B^-1
//
// @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
//       uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
//       see also the explanation in lpi.h.
RetCode LPGlopInterface::GetBInvertedRow(
  int row_number,         // row number
  std::vector<double>& row_coeffs, // array to store the coefficients of the row
  std::vector<int>& indices,    // array to store the non-zero indices
  int& num_indices          // the number of non-zero indices (-1: if we do not store sparsity information)
) const {

  solver_.GetBasisFactorization().LeftSolveForUnitRow(ColIndex(row_number), tmp_row_);
  scaler_.UnscaleUnitRowLeftSolve(solver_.GetBasis(RowIndex(row_number)), tmp_row_);

  const ColIndex size = tmp_row_->values.size();
  assert(size.value() == linear_program_.num_constraints());

  // if we want a sparse vector
  if (num_indices > 0 && !indices.empty()) {
    num_indices = 0;
    // Vectors in Glop might be stored in dense or sparse format dep
    //
    // ending on the values. If non_zeros are given, we
    // can directly loop over the non_zeros, otherwise we have to collect the nonzeros.
    if (!tmp_row_->non_zeros.empty()) {
      ScatteredRowIterator end = tmp_row_->end();
      for (ScatteredRowIterator iter = tmp_row_->begin(); iter != end; ++iter) {
        int idx = (*iter).column().value();
        assert(0 <= idx && idx < linear_program_.num_constraints());
        row_coeffs[idx] = (*iter).coefficient();
        indices[(num_indices)++] = idx;
      }
    } else {
      // use dense access to tmp_row_
      const Fractional eps = parameters_.primal_feasibility_tolerance();
      for (ColIndex col(0); col < size; ++col) {
        double val = (*tmp_row_)[col];
        if (fabs(val) >= eps) {
          row_coeffs[col.value()] = val;
          indices[(num_indices)++] = col.value();
        }
      }
    }
    return RetCode::kOkay;
  }

  // dense version
  for (ColIndex col(0); col < size; ++col)
    row_coeffs[col.value()] = (*tmp_row_)[col];

  if (num_indices >= 0)
    num_indices = -1;

  return RetCode::kOkay;
}

// get column of inverse basis matrix B^-1
//
// @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
//       uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
//       see also the explanation in lpi.h.
RetCode LPGlopInterface::GetBInvertedColumn(
  int col_number,         // column number of B^-1; this is NOT the number of the column in the LP;
                             // you have to call MiniMIP::LPInterface.GetBasisIndices() to get the array which links the
                             // B^-1 column numbers to the row and column numbers of the LP!
                             // c must be between 0 and num_rows-1, since the basis has the size
                             // num_rows * num_rows
  std::vector<double>& col_coeffs, // array to store the coefficients of the column
  std::vector<int>& indices,    // array to store the non-zero indices
  int& num_indices          // the number of non-zero indices (-1: if we do not store sparsity information)
) const {

  // we need to loop through the rows to extract the values for column col_number
  const ColIndex col(col_number);
  const RowIndex num_rows = linear_program_.num_constraints();

  // if we want a sparse vector
  if (num_indices > 0 && !indices.empty()) {
    const Fractional eps = parameters_.primal_feasibility_tolerance();

    num_indices = 0;
    for (int row = 0; row < num_rows; ++row) {
      solver_.GetBasisFactorization().LeftSolveForUnitRow(ColIndex(row), tmp_row_);
      scaler_.UnscaleUnitRowLeftSolve(solver_.GetBasis(RowIndex(row)), tmp_row_);

      double val = (*tmp_row_)[col];
      if (fabs(val) >= eps) {
        col_coeffs[row] = val;
        indices[(num_indices)++] = row;
      }
    }
    return RetCode::kOkay;
  }

  // dense version
  for (int row = 0; row < num_rows; ++row) {
    solver_.GetBasisFactorization().LeftSolveForUnitRow(ColIndex(row), tmp_row_);
    scaler_.UnscaleUnitRowLeftSolve(solver_.GetBasis(RowIndex(row)), tmp_row_);
    col_coeffs[row] = (*tmp_row_)[col];
  }

  if (num_indices >= 0)
    num_indices = -1;

  return RetCode::kOkay;
}

// get row of inverse basis matrix times constraint matrix B^-1 * A
//
// @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
//       uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
//       see also the explanation in lpi.h.
RetCode LPGlopInterface::GetBInvertedARow(
  int row_number,                   // row number
  const std::vector<double>& b_inverted_row, // row in (A_B)^-1 from prior call to MiniMIP::LPInterface.GetBInvRow()
  std::vector<double>& row_coeffs,           // array to store coefficients of the row
  std::vector<int>& indices,              // array to store the non-zero indices
  int& num_indices                    // thee number of non-zero indices (-1: if we do not store sparsity information)
) const {

  // get row of basis inverse, loop through columns and muliply with matrix
  solver_.GetBasisFactorization().LeftSolveForUnitRow(ColIndex(row_number), tmp_row_);
  scaler_.UnscaleUnitRowLeftSolve(solver_.GetBasis(RowIndex(row_number)), tmp_row_);

  const ColIndex num_cols = linear_program_.num_variables();

  // if we want a sparse vector
  if (num_indices > 0 && !indices.empty()) {
    const Fractional eps = parameters_.primal_feasibility_tolerance();

    num_indices = 0;
    for (ColIndex col(0); col < num_cols; ++col) {
      double val = operations_research::glop::ScalarProduct(tmp_row_->values, linear_program_.GetSparseColumn(col));
      if (fabs(val) >= eps) {
        row_coeffs[col.value()] = val;
        indices[num_indices++] = col.value();
      }
    }
    return RetCode::kOkay;
  }

  // dense version
  for (ColIndex col(0); col < num_cols; ++col) {
    double check = operations_research::glop::ScalarProduct(tmp_row_->values, linear_program_.GetSparseColumn(col));
    if (EPS <= fabs(check)) {
      row_coeffs[col.value()] = operations_research::glop::ScalarProduct(tmp_row_->values, linear_program_.GetSparseColumn(col));
    } else {
      row_coeffs[col.value()] = 0;
    }
  }
  num_indices = -1;

  return RetCode::kOkay;
}

// get column of inverse basis matrix times constraint matrix B^-1 * A
//
// @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
//       uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
//       see also the explanation in lpi.h.
RetCode LPGlopInterface::GetBInvertedAColumn(
  int col_number,         // column number
  std::vector<double>& col_coeffs, // array to store coefficients of the column
  std::vector<int>& indices,    // array to store the non-zero indices
  int& num_indices          // the number of non-zero indices (-1: if we do not store sparsity information)
) const {

  solver_.GetBasisFactorization().RightSolveForProblemColumn(ColIndex(col_number), tmp_column_);
  scaler_.UnscaleColumnRightSolve(solver_.GetBasisVector(), ColIndex(col_number), tmp_column_);

  const RowIndex num_rows = tmp_column_->values.size();

  // if we want a sparse vector
  if (num_indices > 0 && !indices.empty()) {
    num_indices = 0;
    // Vectors in Glop might be stored in dense or sparse format depending on the values. If non_zeros are given, we
    // can directly loop over the non_zeros, otherwise we have to collect the nonzeros.
    if (!tmp_column_->non_zeros.empty()) {
      ScatteredColumnIterator end = tmp_column_->end();
      for (ScatteredColumnIterator iter = tmp_column_->begin(); iter != end; ++iter) {
        int idx = (*iter).row().value();
        assert(0 <= idx && idx < num_rows);
        col_coeffs[idx] = (*iter).coefficient();
        indices[(num_indices)++] = idx;
      }
    } else {
      // use dense access to tmp_column_
      const Fractional eps = parameters_.primal_feasibility_tolerance();
      for (RowIndex row(0); row < num_rows; ++row) {
        double val = (*tmp_column_)[row];
        if (fabs(val) > eps) {
          col_coeffs[row.value()] = val;
          indices[(num_indices)++] = row.value();
        }
      }
    }
    return RetCode::kOkay;
  }

  // dense version
  for (RowIndex row(0); row < num_rows; ++row)
    col_coeffs[row.value()] = (*tmp_column_)[row];

  if (num_indices >= 0)
    num_indices = -1;

  return RetCode::kOkay;
}

// @}

// Parameter Methods

// @name Parameter Methods
// @{

// gets integer parameter of LP
RetCode LPGlopInterface::GetIntegerParameter(
  LPParameter type, // parameter number
  int& param_val  // buffer to store the parameter value
) const {

  switch (type) {
    case LPParameter::kFromScratch:
      param_val = (int) from_scratch_;
      MiniMIPdebugMessage("GetIntegerParameter: LPParameter::kFromScratch = %d.\n", param_val);
      break;
    case LPParameter::kLPInfo:
      param_val = (int) lp_info_;
      MiniMIPdebugMessage("GetIntegerParameter: LPParameter::kLPInfo = %d.\n", param_val);
      break;
    case LPParameter::kLPIterationLimit:
      param_val = (int) parameters_.max_number_of_iterations();
      MiniMIPdebugMessage("GetIntegerParameter: LPParameter::kLPIterationLimit = %d.\n", param_val);
      break;
    case LPParameter::kPresolving:
      param_val = parameters_.use_preprocessing();
      MiniMIPdebugMessage("GetIntegerParameter: LPParameter::kPresolving = %d.\n", param_val);
      break;
    case LPParameter::kPricing:
      param_val = (int) pricing_;
      MiniMIPdebugMessage("GetIntegerParameter: LPParameter::kPricing = %d.\n", param_val);
      break;
#ifndef NOSCALING
    case LPParameter::kScaling:
      param_val = parameters_.use_scaling();
      MiniMIPdebugMessage("GetIntegerParameter: LPParameter::kScaling = %d.\n", param_val);
      break;
#endif
    case LPParameter::kThreads:
      param_val = num_threads_;
      MiniMIPdebugMessage("GetIntegerParameter: LPParameter::kThreads = %d.\n", param_val);
      break;
    case LPParameter::kTiming:
      param_val = timing_;
      MiniMIPdebugMessage("GetIntegerParameter: LPParameter::kTiming = %d.\n", param_val);
      break;
    case LPParameter::kRandomSeed:
      param_val = (int) parameters_.random_seed();
      MiniMIPdebugMessage("GetIntegerParameter: LPParameter::kRandomSeed = %d.\n", param_val);
      break;
    default:
      return RetCode::kParameterUnknown;
  }

  return RetCode::kOkay;
}

// sets integer parameter of LP
RetCode LPGlopInterface::SetIntegerParameter(
  LPParameter type, // parameter number
  int param_val   // parameter value
) {

  switch (type) {
    case LPParameter::kFromScratch:
      MiniMIPdebugMessage("SetIntegerParameter: LPParameter::kFromScratch -> %d.\n", param_val);
      from_scratch_ = static_cast<bool>(param_val);
      break;
    case LPParameter::kLPInfo:
      MiniMIPdebugMessage("SetIntegerParameter: LPParameter::kLPInfo -> %d.\n", param_val);
      if (param_val == 0) {
        static_cast<void>(google::SetVLOGLevel("*", google::GLOG_INFO));
        lp_info_ = false;
      } else {
        static_cast<void>(google::SetVLOGLevel("*", google::GLOG_ERROR));
        lp_info_ = true;
      }
      break;
    case LPParameter::kLPIterationLimit:
      MiniMIPdebugMessage("SetIntegerParameter: LPParameter::kLPIterationLimit -> %d.\n", param_val);
      parameters_.set_max_number_of_iterations(param_val);
      break;
    case LPParameter::kPresolving:
      MiniMIPdebugMessage("SetIntegerParameter: LPParameter::kPresolving -> %d.\n", param_val);
      parameters_.set_use_preprocessing(param_val);
      break;
    case LPParameter::kPricing:
      MiniMIPdebugMessage("SetIntegerParameter: LPParameter::kPricing -> %d.\n", param_val);
      pricing_ = (LPPricing) param_val;
      switch (pricing_) {
        case LPPricing::kDefault:
        case LPPricing::kAuto:
        case LPPricing::kPartial:
        case LPPricing::kSteep:
        case LPPricing::kSteepQStart:
          parameters_.set_feasibility_rule(operations_research::glop::GlopParameters_PricingRule_STEEPEST_EDGE);
          break;
        case LPPricing::kFull:
          // Dantzig does not really fit, but use it anyway
          parameters_.set_feasibility_rule(operations_research::glop::GlopParameters_PricingRule_DANTZIG);
          break;
        case LPPricing::kDevex:
          parameters_.set_feasibility_rule(operations_research::glop::GlopParameters_PricingRule_DEVEX);
          break;
        default:
          return RetCode::kParameterUnknown;
      }
      break;
#ifndef NOSCALING
    case LPParameter::kScaling:
      MiniMIPdebugMessage("SetIntegerParameter: LPParameter::kScaling -> %d.\n", param_val);
      parameters_.set_use_scaling(param_val);
      break;
#endif
    case LPParameter::kThreads:
      MiniMIPdebugMessage("SetIntegerParameter: LPParameter::kThreads -> %d.\n", param_val);
      num_threads_ = param_val;
      if (param_val == 0)
        parameters_.set_num_omp_threads(1);
      else
        parameters_.set_num_omp_threads(param_val);
      break;
    case LPParameter::kTiming:
      MiniMIPdebugMessage("SetIntegerParameter: LPParameter::kTiming -> %d.\n", param_val);
      assert(param_val <= 2);
      timing_ = param_val;
      if (param_val == 1)
        absl::SetFlag(&FLAGS_time_limit_use_usertime, true);
      else
        absl::SetFlag(&FLAGS_time_limit_use_usertime, false);
      break;
    case LPParameter::kRandomSeed:
      MiniMIPdebugMessage("SetIntegerParameter: LPParameter::kRandomSeed -> %d.\n", param_val);
      parameters_.set_random_seed(param_val);
      break;
    default:
      return RetCode::kParameterUnknown;
  }

  return RetCode::kOkay;
}

// gets floating point parameter of LP
RetCode LPGlopInterface::GetRealParameter(
  LPParameter type,  // parameter number
  double& param_val // buffer to store the parameter value
) const {

  switch (type) {
    case LPParameter::kFeasibilityTolerance:
      param_val = parameters_.primal_feasibility_tolerance();
      MiniMIPdebugMessage("GetRealParameter: LPParameter::kFeasibilityTolerance = %g.\n", param_val);
      break;
    case LPParameter::kDualFeasibilityTolerance:
      param_val = parameters_.dual_feasibility_tolerance();
      MiniMIPdebugMessage("GetRealParameter: LPParameter::kDualFeasibilityTolerance = %g.\n", param_val);
      break;
    case LPParameter::kObjectiveLimit:
      if (linear_program_.IsMaximizationProblem())
        param_val = parameters_.objective_lower_limit();
      else
        param_val = parameters_.objective_upper_limit();
      MiniMIPdebugMessage("GetRealParameter: LPParameter::kObjectiveLimit = %f.\n", param_val);
      break;
    case LPParameter::kLPTimeLimit:
      if (absl::GetFlag(FLAGS_time_limit_use_usertime))
        param_val = parameters_.max_time_in_seconds();
      else
        param_val = parameters_.max_deterministic_time();
      MiniMIPdebugMessage("GetRealParameter: LPParameter::kLPTimeLimit = %f.\n", param_val);
      break;

    default:
      return RetCode::kParameterUnknown;
  }

  return RetCode::kOkay;
}

// sets floating point parameter of LP
RetCode LPGlopInterface::SetRealParameter(
  LPParameter type, // parameter number
  double param_val // parameter value
) {

  switch (type) {
    case LPParameter::kFeasibilityTolerance:
      MiniMIPdebugMessage("SetRealParameter: LPParameter::kFeasibilityTolerance -> %g.\n", param_val);
      parameters_.set_primal_feasibility_tolerance(param_val);
      break;
    case LPParameter::kDualFeasibilityTolerance:
      MiniMIPdebugMessage("SetRealParameter: LPParameter::kDualFeasibilityTolerance -> %g.\n", param_val);
      parameters_.set_dual_feasibility_tolerance(param_val);
      break;
    case LPParameter::kObjectiveLimit:
      MiniMIPdebugMessage("SetRealParameter: LPParameter::kObjectiveLimit -> %f.\n", param_val);
      if (linear_program_.IsMaximizationProblem())
        parameters_.set_objective_lower_limit(param_val);
      else
        parameters_.set_objective_upper_limit(param_val);
      break;
    case LPParameter::kLPTimeLimit:
      MiniMIPdebugMessage("SetRealParameter: LPParameter::kLPTimeLimit -> %f.\n", param_val);
      if (absl::GetFlag(FLAGS_time_limit_use_usertime))
        parameters_.set_max_time_in_seconds(param_val);
      else
        parameters_.set_max_deterministic_time(param_val);
      break;
    default:
      return RetCode::kParameterUnknown;
  }

  return RetCode::kOkay;
}

// @}

// Numerical Methods

// @name Numerical Methods
// @{

// returns value treated as infinity in the LP solver
double LPGlopInterface::Infinity() const {
  return std::numeric_limits<double>::infinity();
}

// checks if given value is treated as infinity in the LP solver
bool LPGlopInterface::IsInfinity(
  double val // value to be checked for infinity
) const {
  return val == std::numeric_limits<double>::infinity();
}

// @}

// File Interface Methods

// @name File Interface Methods
// @{

// reads LP from a file
RetCode LPGlopInterface::ReadLP(
  const char* file_name // file name
) {
  assert(file_name != nullptr);

  const char* filespec(file_name);
  MPModelProto proto;
  if (!ReadFileToProto(filespec, &proto)) {
    MiniMIPerrorMessage("Could not read <%s>\n", file_name);
    return RetCode::kReadError;
  }
  linear_program_.Clear();
  MPModelProtoToLinearProgram(proto, &linear_program_);

  return RetCode::kOkay;
}

// writes LP to a file
RetCode LPGlopInterface::WriteLP(
  const char* file_name // file name
) const {
  assert(file_name != nullptr);

  MPModelProto proto;
  LinearProgramToMPModelProto(linear_program_, &proto);
  const char* filespec(file_name);
  if (!WriteProtoToFile(filespec, proto, operations_research::ProtoWriteFormat::kProtoText, true)) {
    MiniMIPerrorMessage("Could not write <%s>\n", file_name);
    return RetCode::kReadError;
  }

  return RetCode::kOkay;
}

} // namespace minimip