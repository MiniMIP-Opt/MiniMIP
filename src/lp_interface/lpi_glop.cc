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

/** LP interface */
namespace minimip {
/* constructor */
LPGlopInterface::LPGlopInterface() : lp_modified_since_last_solve_(true),
                                     lp_time_limit_was_reached_(false),
                                     lp_info_(false),
                                     pricing_(LPPricing::DEFAULT),
                                     from_scratch_(false),
                                     num_threads_(0),
                                     timing_(0),
                                     niterations_(0LL),
                                     tmp_row_(new ScatteredRow()),
                                     tmp_column_(new ScatteredColumn()) {}

/* Destructor default */
LPGlopInterface::~LPGlopInterface() {
  MiniMIPdebugMessage("LPGLopInterface Free\n");
}

/*
* Modification Methods
*/

/**@name Modification Methods */
/**@{ */

/** copies LP data with column matrix into LP solver */
RetCode LPGlopInterface::LoadColumnLP(
  LPObjectiveSense obj_sense,           /**< objective sense */
  LPNum num_cols,                       /**< number of columns */
  const LPValueArray& objective_values, /**< objective function values of columns */
  const LPValueArray& lower_bounds,     /**< lower bounds of columns */
  const LPValueArray& upper_bounds,     /**< upper bounds of columns */
  StringArray& col_names,               /**< column names */
  LPNum num_rows,                       /**< number of rows */
  const LPValueArray& left_hand_sides,  /**< left hand sides of rows */
  const LPValueArray& right_hand_sides, /**< right hand sides of rows */
  StringArray& row_names,               /**< row names */
  LPNum num_non_zeros,                  /**< number of non-zero elements in the constraint matrix */
  const LPIndexArray& begin_cols,       /**< start index of each column in row_indices- and vals-array */
  const LPIndexArray& row_indices,      /**< row indices of constraint matrix entries */
  const LPValueArray& vals              /**< values of constraint matrix entries */
) {

  linear_program_.Clear();
  AddRows(num_rows, left_hand_sides, right_hand_sides, row_names, 0, LPIndexArray(), LPIndexArray(), LPValueArray());
  AddColumns(num_cols, objective_values, lower_bounds, upper_bounds, col_names, num_non_zeros, begin_cols, row_indices, vals);
  ChangeObjectiveSense(obj_sense);

  return RetCode::OKAY;
}

/** adds columns to the LP */
RetCode LPGlopInterface::AddColumns(
  LPNum num_cols,                       /**< number of columns to be added */
  const LPValueArray& objective_values, /**< objective function values of new columns */
  const LPValueArray& lower_bounds,     /**< lower bounds of new columns */
  const LPValueArray& upper_bounds,     /**< upper bounds of new columns */
  StringArray& col_names,               /**< column names */
  LPNum num_non_zeros,                  /**< number of non-zero elements to be added to the constraint matrix */
  const LPIndexArray& begin_cols,       /**< start index of each column in indices- and vals-array */
  const LPIndexArray& indices,          /**< row indices of constraint matrix entries */
  const LPValueArray& vals              /**< values of constraint matrix entries */
) {

  MiniMIPdebugMessage("adding %d columns with %d nonzeros.\n", num_cols, num_non_zeros);

  /* @todo add names */
  if (num_non_zeros > 0) {
    assert(num_cols > 0);

#ifndef NDEBUG
    /* perform check that no new rows are added */
    RowIndex num_rows = linear_program_.num_constraints();
    for (LPNum j = 0; j < num_non_zeros; ++j) {
      assert(0 <= indices[j] && static_cast<int>(indices[j]) < num_rows.value());
      assert(vals[j] != 0.0);
    }
#endif

    LPIndex nz = 0;
    for (LPNum i = 0; i < num_cols; ++i) {
      const ColIndex col = linear_program_.CreateNewVariable();
      linear_program_.SetVariableBounds(col, lower_bounds[i], upper_bounds[i]);
      linear_program_.SetObjectiveCoefficient(col, objective_values[i]);
      const LPIndex end = (num_non_zeros == 0 || i == num_cols - 1) ? num_non_zeros : begin_cols[i + 1];
      while (nz < end) {
        linear_program_.SetCoefficient(RowIndex(static_cast<int>(indices[nz])), col, vals[nz]);
        ++nz;
      }
    }
    assert(nz == num_non_zeros);
  } else {
    for (LPNum i = 0; i < num_cols; ++i) {
      const ColIndex col = linear_program_.CreateNewVariable();
      linear_program_.SetVariableBounds(col, lower_bounds[i], upper_bounds[i]);
      linear_program_.SetObjectiveCoefficient(col, objective_values[i]);
    }
  }

  lp_modified_since_last_solve_ = true;

  return RetCode::OKAY;
}

/** deletes all columns in the given range from LP */
RetCode LPGlopInterface::DeleteColumns(
  LPNum first_col, /**< first column to be deleted */
  LPNum last_col   /**< last column to be deleted */
) {
  assert(0 <= first_col && first_col <= last_col && last_col < linear_program_.num_variables());

  MiniMIPdebugMessage("deleting columns %d to %d.\n", first_col, last_col);

  const ColIndex num_cols = linear_program_.num_variables();
  DenseBooleanRow columns_to_delete(num_cols, false);
  for (LPNum i = first_col; i <= last_col; ++i)
    columns_to_delete[ColIndex(static_cast<int>(i))] = true;

  linear_program_.DeleteColumns(columns_to_delete);
  lp_modified_since_last_solve_ = true;

  return RetCode::OKAY;
}

/** deletes columns from MiniMIP_LP; the new position of a column must not be greater that its old position */
RetCode LPGlopInterface::DeleteColumnSet(
  BoolArray& deletion_status /**< deletion status of columns */
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

  return RetCode::OKAY;
}

/** adds rows to the LP */
RetCode LPGlopInterface::AddRows(
  LPNum num_rows,                       /**< number of rows to be added */
  const LPValueArray& left_hand_sides,  /**< left hand sides of new rows */
  const LPValueArray& right_hand_sides, /**< right hand sides of new rows */
  StringArray& row_names,               /**< row names */
  LPNum num_non_zeros,                  /**< number of non-zero elements to be added to the constraint matrix */
  const LPIndexArray& begin_rows,       /**< start index of each row in indices- and vals-array */
  const LPIndexArray& indices,          /**< column indices of constraint matrix entries */
  const LPValueArray& vals              /**< values of constraint matrix entries */
) {

  MiniMIPdebugMessage("adding %d rows with %d nonzeros.\n", num_rows, num_non_zeros);

  /* @todo add names */
  if (num_non_zeros > 0) {
    assert(num_rows > 0);

#ifndef NDEBUG
    /* perform check that no new columns are added - this is likely to be a mistake */
    const ColIndex num_cols = linear_program_.num_variables();
    for (LPNum j = 0; j < num_non_zeros; ++j) {
      assert(vals[j] != 0.0);
      assert(0 <= indices[j] && static_cast<int>(indices[j]) < num_cols.value());
    }
#endif

    LPNum nz = 0;
    for (LPNum i = 0; i < num_rows; ++i) {
      const RowIndex row = linear_program_.CreateNewConstraint();
      linear_program_.SetConstraintBounds(row, left_hand_sides[i], right_hand_sides[i]);
      const LPNum end = (num_non_zeros == 0 || i == num_rows - 1) ? num_non_zeros : begin_rows[i + 1];
      while (nz < end) {
        linear_program_.SetCoefficient(row, ColIndex(static_cast<int>(indices[nz])), vals[nz]);
        ++nz;
      }
    }
    assert(nz == num_non_zeros);
  } else {
    for (LPNum i = 0; i < num_rows; ++i) {
      const RowIndex row = linear_program_.CreateNewConstraint();
      linear_program_.SetConstraintBounds(row, left_hand_sides[i], right_hand_sides[i]);
    }
  }

  lp_modified_since_last_solve_ = true;

  return RetCode::OKAY;
}

/** delete rows from LP and update the current basis */
void LPGlopInterface::DeleteRowsAndUpdateCurrentBasis(
  const DenseBooleanColumn& rows_to_delete /**< array to mark rows that should be deleted */
) {
  const RowIndex num_rows = linear_program_.num_constraints();
  const ColIndex num_cols = linear_program_.num_variables();

  /* try to repair basis status if problem size has not changed before */
  BasisState state = solver_.GetState();
  if (state.statuses.size() == num_cols.value() + num_rows.value()) {
    /* Shift the status of the non-deleted rows. Note that if the deleted rows where part of the basis (i.e., constraint
     * not tight), then we should be left with a correct basis afterward. This should be the most common use case in MiniMIP. */
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

/** deletes all rows in the given range from LP */
RetCode LPGlopInterface::DeleteRows(
  LPNum first_row, /**< first row to be deleted */
  LPNum last_row   /**< last row to be deleted */
) {
  assert(0 <= first_row && first_row <= last_row && last_row < linear_program_.num_constraints());

  const RowIndex num_rows = linear_program_.num_constraints();
  DenseBooleanColumn rows_to_delete(num_rows, false);
  for (LPNum i = first_row; i <= last_row; ++i)
    rows_to_delete[RowIndex(static_cast<int>(i))] = true;

  MiniMIPdebugMessage("deleting rows %d to %d.\n", first_row, last_row);
  DeleteRowsAndUpdateCurrentBasis(rows_to_delete);

  return RetCode::OKAY;
}

/** deletes rows from LP; the new position of a row must not be greater that its old position */
RetCode LPGlopInterface::DeleteRowSet(
  BoolArray& deletion_status /**< deletion status of rows */
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

  return RetCode::OKAY;
}

/** clears the whole LP */
RetCode LPGlopInterface::Clear() {

  MiniMIPdebugMessage("Clear\n");

  linear_program_.Clear();
  lp_modified_since_last_solve_ = true;

  return RetCode::OKAY;
}

/** clears current LPi state (like basis information) of the solver */
RetCode LPGlopInterface::ClearState() {

  solver_.ClearStateForNextSolve();

  return RetCode::OKAY;
}

/** changes lower and upper bounds of columns */
RetCode LPGlopInterface::ChangeBounds(
  LPNum num_cols,                   /**< number of columns to change bounds for */
  const LPIndexArray& indices,      /**< column indices */
  const LPValueArray& lower_bounds, /**< values for the new lower bounds */
  const LPValueArray& upper_bounds  /**< values for the new upper bounds */
) {

  MiniMIPdebugMessage("changing %d bounds.\n", num_cols);

  for (LPNum i = 0; i < num_cols; ++i) {
    MiniMIPdebugMessage("  col %d: [%g,%g]\n", indices[i], lower_bounds[i], upper_bounds[i]);

    if (IsInfinity(lower_bounds[i])) {
      MiniMIPerrorMessage("LP Error: fixing lower bound for variable %d to infinity.\n", indices[i]);
      return RetCode::LP_ERROR;
    }
    if (IsInfinity(-upper_bounds[i])) {
      MiniMIPerrorMessage("LP Error: fixing upper bound for variable %d to -infinity.\n", indices[i]);
      return RetCode::LP_ERROR;
    }

    linear_program_.SetVariableBounds(ColIndex(static_cast<int>(indices[i])), lower_bounds[i], upper_bounds[i]);
  }
  lp_modified_since_last_solve_ = true;

  return RetCode::OKAY;
}

/** changes left and right hand sides of rows */
RetCode LPGlopInterface::ChangeSides(
  LPNum num_rows,                      /**< number of rows to change sides for */
  const LPIndexArray& indices,         /**< row indices */
  const LPValueArray& left_hand_sides, /**< new values for left hand sides */
  const LPValueArray& right_hand_sides /**< new values for right hand sides */
) {

  MiniMIPdebugMessage("changing %d sides\n", num_rows);

  for (LPNum i = 0; i < num_rows; ++i)
    linear_program_.SetConstraintBounds(RowIndex(static_cast<int>(indices[i])), left_hand_sides[i], right_hand_sides[i]);

  lp_modified_since_last_solve_ = true;

  return RetCode::OKAY;
}

/** changes the objective sense */
RetCode LPGlopInterface::ChangeObjectiveSense(
  LPObjectiveSense obj_sense /**< new objective sense */
) {

  switch (obj_sense) {
    case LPObjectiveSense::OBJ_SENSE_MAXIMIZE:
      MiniMIPdebugMessage("changing objective sense to MAXIMIZE\n");
      linear_program_.SetMaximizationProblem(true);
      break;
    case LPObjectiveSense::OBJ_SENSE_MINIMIZE:
      MiniMIPdebugMessage("changing objective sense to MINIMIZE\n");
      linear_program_.SetMaximizationProblem(false);
      break;
  }
  lp_modified_since_last_solve_ = true;

  return RetCode::OKAY;
}

/** changes objective values of columns in the LP */
RetCode LPGlopInterface::ChangeObjective(
  LPNum num_cols,                  /**< number of columns to change objective value for */
  const LPIndexArray& indices,     /**< column indices to change objective value for */
  const LPValueArray& new_obj_vals /**< new objective values for columns */
) {

  MiniMIPdebugMessage("changing %d objective values\n", num_cols);

  for (LPNum i = 0; i < num_cols; ++i)
    linear_program_.SetObjectiveCoefficient(ColIndex(static_cast<int>(indices[i])), new_obj_vals[i]);

  lp_modified_since_last_solve_ = true;

  return RetCode::OKAY;
}


/*
 * Data Accessing Methods
 */

/**@name Data Accessing Methods */
/**@{ */

/** gets the number of rows in the LP */
LPNum LPGlopInterface::GetNumberOfRows() {

  MiniMIPdebugMessage("getting number of rows.\n");

  return linear_program_.num_constraints().value();
}

/** gets the number of columns in the LP */
LPNum LPGlopInterface::GetNumberOfColumns() {

  MiniMIPdebugMessage("getting number of columns.\n");

  return linear_program_.num_variables().value();
}

/** gets objective sense of the LP */
LPObjectiveSense LPGlopInterface::GetObjectiveSense() {

  MiniMIPdebugMessage("getting objective sense.\n");

  return linear_program_.IsMaximizationProblem() ? LPObjectiveSense::OBJ_SENSE_MAXIMIZE : LPObjectiveSense::OBJ_SENSE_MINIMIZE;
}

/** gets the number of nonzero elements in the LP constraint matrix */
LPNum LPGlopInterface::GetNumberOfNonZeros() {

  MiniMIPdebugMessage("getting number of non-zeros.\n");

  return static_cast<int>(linear_program_.num_entries().value());
}

/** gets columns from LP problem object
 *
 *  Either both, lb and ub, have to be NULL, or both have to be non-NULL,
 *  either num_non_zeros, begin_cols, indices, and val have to be NULL, or all of them have to be non-NULL.
 */
RetCode LPGlopInterface::GetColumns(
  LPNum first_col,            /**< first column to get from LP */
  LPNum last_col,             /**< last column to get from LP */
  LPValueArray& lower_bounds, /**< array to store the lower bound vector */
  LPValueArray& upper_bounds, /**< array to store the upper bound vector */
  LPNum& num_non_zeros,       /**< store the number of non-zero elements */
  LPIndexArray& begin_cols,   /**< array to store start index of each column in indices- and vals-array */
  LPIndexArray& indices,      /**< array to store row indices of constraint matrix entries */
  LPValueArray& vals          /**< array to store values of constraint matrix entries */
) {
  assert(0 <= first_col && first_col <= last_col && last_col < linear_program_.num_variables());

  const DenseRow& tmp_lower_bound = linear_program_.variable_lower_bounds();
  const DenseRow& tmp_upper_bound = linear_program_.variable_upper_bounds();

  if (num_non_zeros >= 0) {
    num_non_zeros = 0;
    int index = 0;
    for (ColIndex col(static_cast<int>(first_col)); col <= ColIndex(static_cast<int>(last_col)); ++col, ++index) {
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
    for (ColIndex col(static_cast<int>(first_col)); col <= ColIndex(static_cast<int>(last_col)); ++col, ++index) {
      lower_bounds[index] = tmp_lower_bound[col];
      upper_bounds[index] = tmp_upper_bound[col];
    }
  }

  return RetCode::OKAY;
}

/** gets rows from LP problem object
 *
 *  Either both, left_hand_side and right_hand_side, have to be NULL, or both have to be non-NULL,
 *  either num_non_zeros, begin_rows, indices, and val have to be NULL, or all of them have to be non-NULL.
 */
RetCode LPGlopInterface::GetRows(
  LPNum first_row,                /**< first row to get from LP */
  LPNum last_row,                 /**< last row to get from LP */
  LPValueArray& left_hand_sides,  /**< array to store left hand side vector */
  LPValueArray& right_hand_sides, /**< array to store right hand side vector */
  LPNum& num_non_zeros,           /**< store the number of non-zero elements */
  LPIndexArray& begin_rows,       /**< array to store start index of each row in indices- and vals-array */
  LPIndexArray& indices,          /**< array to store column indices of constraint matrix entries */
  LPValueArray& vals              /**< array to store values of constraint matrix entries */
) {
  assert(0 <= first_row && first_row <= last_row && last_row < linear_program_.num_constraints());

  const DenseColumn& tmplhs = linear_program_.constraint_lower_bounds();
  const DenseColumn& tmprhs = linear_program_.constraint_upper_bounds();

  const SparseMatrix& matrixtrans = linear_program_.GetTransposeSparseMatrix();

  num_non_zeros = 0;
  int index = 0;
  for (RowIndex row(static_cast<int>(first_row)); row <= RowIndex(static_cast<int>(last_row)); ++row, ++index) {
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

  return RetCode::OKAY;
}

/** gets objective coefficients from LP problem object */
RetCode LPGlopInterface::GetObjective(
  LPNum first_col,         /**< first column to get objective coefficient for */
  LPNum last_col,          /**< last column to get objective coefficient for */
  LPValueArray& obj_coeffs /**< array to store objective coefficients */
) {
  assert(first_col <= last_col);

  MiniMIPdebugMessage("getting objective values %d to %d\n", first_col, last_col);

  int index = 0;
  for (ColIndex col(static_cast<int>(first_col)); col <= ColIndex(static_cast<int>(last_col)); ++col) {
    obj_coeffs[index] = linear_program_.objective_coefficients()[col];
    ++index;
  }

  return RetCode::OKAY;
}

/** gets current bounds from LP problem object */
RetCode LPGlopInterface::GetBounds(
  LPNum first_col,            /**< first column to get bounds for */
  LPNum last_col,             /**< last column to get bounds for */
  LPValueArray& lower_bounds, /**< array to store lower bound values */
  LPValueArray& upper_bounds  /**< array to store upper bound values */
) {
  assert(first_col <= last_col);

  MiniMIPdebugMessage("getting bounds %d to %d\n", first_col, last_col);

  int index = 0;
  for (ColIndex col(static_cast<int>(first_col)); col <= ColIndex(static_cast<int>(last_col)); ++col) {
    lower_bounds[index] = linear_program_.variable_lower_bounds()[col];

    upper_bounds[index] = linear_program_.variable_upper_bounds()[col];

    ++index;
  }

  return RetCode::OKAY;
}

/** gets current row sides from LP problem object */
RetCode LPGlopInterface::GetSides(
  LPNum first_row,               /**< first row to get sides for */
  LPNum last_row,                /**< last row to get sides for */
  LPValueArray& left_hand_sides, /**< array to store left hand side values */
  LPValueArray& right_hand_sides /**< array to store right hand side values */
) {
  assert(first_row <= last_row);

  MiniMIPdebugMessage("getting row sides %d to %d\n", first_row, last_row);

  LPIndex index = 0;
  for (RowIndex row(static_cast<int>(first_row)); row <= RowIndex(static_cast<int>(last_row)); ++row) {
    left_hand_sides[index] = linear_program_.constraint_lower_bounds()[row];

    right_hand_sides[index] = linear_program_.constraint_upper_bounds()[row];

    ++index;
  }

  return RetCode::OKAY;
}

/** gets a single coefficient */
RetCode LPGlopInterface::GetCoefficient(
  LPNum row,         /**< row number of coefficient */
  LPIndex col_index, /**< column number of coefficient */
  LPValue& val       /**< array to store the value of the coefficient */
) {

  /* quite slow method: possibly needs linear time if matrix is not sorted */
  const SparseMatrix& matrix = linear_program_.GetSparseMatrix();
  val = matrix.LookUpValue(RowIndex(static_cast<int>(row)), ColIndex(static_cast<int>(col_index)));

  return RetCode::OKAY;
}

/**@} */

/*
* Solving Methods
*/

/**@name Solving Methods */
/**@{ */

/** update scaled linear program */
void LPGlopInterface::updateScaledLP() {
  if (!lp_modified_since_last_solve_)
    return;

  scaled_lp_.PopulateFromLinearProgram(linear_program_);
  scaled_lp_.AddSlackVariablesWhereNecessary(false);

  /* @todo: Avoid doing a copy if there is no scaling. */
  /* @todo: Avoid rescaling if not much changed. */
  if (parameters_.use_scaling())
    scaler_.Scale(&scaled_lp_);
  else
    scaler_.Clear();
}

/** check primal feasibility */
bool LPGlopInterface::checkUnscaledPrimalFeasibility() {

#if UNSCALEDFEAS_CHECK == 1
  /* get unscaled solution */
  const ColIndex num_cols = linear_program_.num_variables();
  DenseRow unscaledsol(num_cols);
  for (ColIndex col = ColIndex(0); col < num_cols; ++col)
    unscaledsol[col] = scaler_.UnscaleVariableValue(col, solver_.GetVariableValue(col));

  /* if the solution is not feasible w.r.t. absolute tolerances, try to fix it in the unscaled problem */
  const LPValue feastol = parameters_.primal_feasibility_tolerance();
  return linear_program_.SolutionIsLPFeasible(unscaledsol, feastol);

#elif UNSCALEDFEAS_CHECK == 2
  const LPValue feastol = parameters_.primal_feasibility_tolerance();

  /* check bounds of unscaled solution */
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

  /* check activities of unscaled solution */
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

/** common function between the two LPI Solve() functions */
RetCode LPGlopInterface::SolveInternal(
  bool recursive,                        /**< Is this a recursive call? */
  std::unique_ptr<TimeLimit>& time_limit /**< time limit */
) {

  updateScaledLP();

  solver_.SetParameters(parameters_);
  lp_time_limit_was_reached_ = false;

  /* possibly ignore warm start information for next solve */
  if (from_scratch_)
    solver_.ClearStateForNextSolve();

  if (!solver_.Solve(scaled_lp_, time_limit.get()).ok()) {
    return RetCode::LP_ERROR;
  }
  lp_time_limit_was_reached_ = time_limit->LimitReached();
  if (recursive)
    niterations_ += (LPLongInt) solver_.GetNumberOfIterations();
  else
    niterations_ = (LPLongInt) solver_.GetNumberOfIterations();

  MiniMIPdebugMessage("status=%s  obj=%f  iterations=%ld.\n", GetProblemStatusString(solver_.GetProblemStatus()).c_str(),
                      solver_.GetObjectiveValue(), solver_.GetNumberOfIterations());

  const ProblemStatus status = solver_.GetProblemStatus();
  if ((status == ProblemStatus::PRIMAL_FEASIBLE || status == ProblemStatus::OPTIMAL) && parameters_.use_scaling()) {
    if (!checkUnscaledPrimalFeasibility()) {
      MiniMIPdebugMessage("Solution not feasible w.r.t. absolute tolerance %g -> reoptimize.\n", parameters_.primal_feasibility_tolerance());

      /* Re-solve without scaling to try to fix the infeasibility. */
      parameters_.set_use_scaling(false);
      lp_modified_since_last_solve_ = true;
      SolveInternal(true, time_limit); /* inherit time limit, so used time is not reset; do not change iteration limit for resolve */
      parameters_.set_use_scaling(true);
    }
  }

  lp_modified_since_last_solve_ = false;

  return RetCode::OKAY;
}

/** calls primal simplex to solve the LP */
RetCode LPGlopInterface::SolvePrimal() {

  MiniMIPdebugMessage("SolvePrimal: %d rows, %d cols.\n", linear_program_.num_constraints().value(), linear_program_.num_variables().value());
  std::unique_ptr<TimeLimit> time_limit = TimeLimit::FromParameters(parameters_);
  niterations_ = 0;

  parameters_.set_use_dual_simplex(false);
  return SolveInternal(false, time_limit);
}

/** calls dual simplex to solve the LP */
RetCode LPGlopInterface::SolveDual() {

  MiniMIPdebugMessage("SolveDual: %d rows, %d cols.\n", linear_program_.num_constraints().value(), linear_program_.num_variables().value());
  std::unique_ptr<TimeLimit> time_limit = TimeLimit::FromParameters(parameters_);
  niterations_ = 0;

  parameters_.set_use_dual_simplex(true);
  return SolveInternal(false, time_limit);
}

/** start strong branching */
RetCode LPGlopInterface::StartStrongbranch() { /*lint --e{715}*/

  updateScaledLP();

  /* @todo Save state and do all the branching from there. */
  return RetCode::OKAY;
}

/** end strong branching */
RetCode LPGlopInterface::EndStrongbranch() { /*lint --e{715}*/

  /* @todo Restore the saved state. */
  return RetCode::OKAY;
}
/** determine whether the dual bound is valid */
bool LPGlopInterface::IsDualBoundValid(
  ProblemStatus status /**< status to be checked */
) {
  return status == ProblemStatus::OPTIMAL || status == ProblemStatus::DUAL_FEASIBLE || status == ProblemStatus::DUAL_UNBOUNDED;
}

/** performs strong branching iterations */
RetCode LPGlopInterface::strongbranch(
  LPIndex col_index,               /**< column to apply strong branching on */
  LPValue primal_sol,              /**< fractional current primal solution value of column */
  LPNum iteration_limit,           /**< iteration limit for strong branchings */
  LPValue& dual_bound_down_branch, /**< stores dual bound after branching column down */
  LPValue& dual_bound_up_branch,   /**< stores dual bound after branching column up */
  bool& down_valid,                /**< stores whether the returned down value is a valid dual bound;
                                    *   otherwise, it can only be used as an estimate value */
  bool& up_valid,                  /**< stores whether the returned up value is a valid dual bound;
                                    *   otherwise, it can only be used as an estimate value */
  LPNum& iterations                /**< stores total number of strong branching iterations, or -1; */
) {

  MiniMIPdebugMessage("calling strongbranching on variable %d (%d iterations)\n", col_index, iteration_limit);

  /* We work on the scaled problem. */
  const ColIndex col(static_cast<int>(col_index));
  const Fractional lower_bound = scaled_lp_.variable_lower_bounds()[col];
  const Fractional upper_bound = scaled_lp_.variable_upper_bounds()[col];
  const LPValue value = primal_sol * scaler_.VariableScalingFactor(col);

  /* Configure solver. */

  /* @todo Use the iteration limit once glop supports incrementality. */
  int num_iterations = 0;
  parameters_.set_use_dual_simplex(true);

  solver_.SetParameters(parameters_);
  const Fractional eps = parameters_.primal_feasibility_tolerance();

  std::unique_ptr<TimeLimit> time_limit = TimeLimit::FromParameters(parameters_);

  /* Down branch. */
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

  /* Up branch. */
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

  /*  Restore bound. */
  scaled_lp_.SetVariableBounds(col, lower_bound, upper_bound);
  if (iterations > 0)
    iterations = num_iterations;

  return RetCode::OKAY;
}

/** performs strong branching iterations on one @b fractional candidate */
RetCode LPGlopInterface::StrongbranchFractionalValue(
  LPIndex col,                     /**< column to apply strong branching on */
  LPValue primal_sol,              /**< fractional current primal solution value of column */
  LPNum iteration_limit,           /**< iteration limit for strong branchings */
  LPValue& dual_bound_down_branch, /**< stores dual bound after branching column down */
  LPValue& dual_bound_up_branch,   /**< stores dual bound after branching column up */
  bool& down_valid,                /**< whether the returned down value is a valid dual bound; otherwise, it can only be used as an estimate value */
  bool& up_valid,                  /**< whether the returned up value is a valid dual bound; otherwise, it can only be used as an estimate value */
  LPNum& iterations                /**< stores total number of strong branching iterations */
) {

  MiniMIPdebugMessage("calling strong branching on fractional variable %d (%d iterations)\n", col, iteration_limit);

  MINIMIP_CALL(strongbranch(col, primal_sol, iteration_limit, dual_bound_down_branch, dual_bound_up_branch, down_valid, up_valid, iterations));

  return RetCode::OKAY;
}

/** performs strong branching iterations on given @b fractional candidates */
RetCode LPGlopInterface::StrongbranchFractionalValues(
  LPNumArray& cols,                       /**< columns to apply strong branching on */
  LPNum num_cols,                         /**< number of columns */
  LPValueArray& primal_sols,              /**< fractional current primal solution values of columns */
  LPNum iteration_limit,                  /**< iteration limit for strong branchings */
  LPValueArray& dual_bound_down_branches, /**< stores dual bounds after branching columns down */
  LPValueArray& dual_bound_up_branches,   /**< stores dual bounds after branching columns up */
  BoolArray& down_valids,                 /**< stores whether the returned down values are valid dual bounds;
                                           *   otherwise, they can only be used as an estimate values */
  BoolArray& up_valids,                   /**< stores whether the returned up values are a valid dual bounds;
                                           *   otherwise, they can only be used as an estimate values */
  LPNum& iterations                       /**< stores total number of strong branching iterations */
) {
  return RetCode::NOT_IMPLEMENTED;
}

/** performs strong branching iterations on one candidate with @b integral value */
RetCode LPGlopInterface::StrongbranchIntegerValue(
  LPIndex col,                     /**< column to apply strong branching on */
  LPValue primal_sol,              /**< current integral primal solution value of column */
  LPNum iteration_limit,           /**< iteration limit for strong branchings */
  LPValue& dual_bound_down_branch, /**< stores dual bound after branching column down */
  LPValue& dual_bound_up_branch,   /**< stores dual bound after branching column up */
  bool& down_valid,                /**< stores whether the returned down value is a valid dual bound;
                                    *   otherwise, it can only be used as an estimate value */
  bool& up_valid,                  /**< stores whether the returned up value is a valid dual bound;
                                    *   otherwise, it can only be used as an estimate value */
  LPNum& iterations                /**< stores total number of strong branching iterations */
) {
  MiniMIPdebugMessage("calling strong branching on integer variable %d (%d iterations)\n", col, iteration_limit);

  MINIMIP_CALL(strongbranch(col, primal_sol, iteration_limit, dual_bound_down_branch, dual_bound_up_branch, down_valid, up_valid, iterations));

  return RetCode::OKAY;
}

/** performs strong branching iterations on given candidates with @b integral values */
RetCode LPGlopInterface::StrongbranchIntegerValues(
  LPNumArray& cols,                       /**< columns to apply strong branching on */
  LPNum num_cols,                         /**< number of columns */
  LPValueArray& primal_sols,              /**< current integral primal solution values of columns */
  LPNum iteration_limit,                  /**< iteration limit for strong branchings */
  LPValueArray& dual_bound_down_branches, /**< stores dual bounds after branching columns down */
  LPValueArray& dual_bound_up_branches,   /**< stores dual bounds after branching columns up */
  BoolArray& down_valids,                 /**< stores whether the returned down values are valid dual bounds;
                                           *   otherwise, they can only be used as an estimate values */
  BoolArray& up_valids,                   /**< stores whether the returned up values are a valid dual bounds;
                                           *   otherwise, they can only be used as an estimate values */
  LPNum& iterations                       /**< stores total number of strong branching iterations */
) {
  return RetCode::NOT_IMPLEMENTED;
}

/**@} */

/*
* Solution Information Methods
*/

/**@name Solution Information Methods */
/**@{ */

/** returns whether a solve method was called after the last modification of the LP */
bool LPGlopInterface::WasSolved() {

  /* @todo Track this to avoid uneeded resolving. */
  return (!lp_modified_since_last_solve_);
}

/** returns true if LP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
 *  this does not necessarily mean that the solver knows and can return the primal ray
 */
bool LPGlopInterface::ExistsPrimalRay() {

  return solver_.GetProblemStatus() == ProblemStatus::PRIMAL_UNBOUNDED;
}

/** returns true if LP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
 *  and the solver knows and can return the primal ray
 */
bool LPGlopInterface::HasPrimalRay() {

  return solver_.GetProblemStatus() == ProblemStatus::PRIMAL_UNBOUNDED;
}

/** returns true if LP is proven to be primal unbounded */
bool LPGlopInterface::IsPrimalUnbounded() {

  return solver_.GetProblemStatus() == ProblemStatus::PRIMAL_UNBOUNDED;
}

/** returns true if LP is proven to be primal infeasible */
bool LPGlopInterface::IsPrimalInfeasible() {

  const ProblemStatus status = solver_.GetProblemStatus();

  return status == ProblemStatus::DUAL_UNBOUNDED || status == ProblemStatus::PRIMAL_INFEASIBLE;
}

/** returns true if LP is proven to be primal feasible */
bool LPGlopInterface::IsPrimalFeasible() {
  const ProblemStatus status = solver_.GetProblemStatus();

  return status == ProblemStatus::PRIMAL_FEASIBLE || status == ProblemStatus::OPTIMAL;
}

/** returns true if LP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
 *  this does not necessarily mean that the solver knows and can return the dual ray
 */
bool LPGlopInterface::ExistsDualRay() {
  const ProblemStatus status = solver_.GetProblemStatus();

  return status == ProblemStatus::DUAL_UNBOUNDED;
}

/** returns true if LP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
 *  and the solver knows and can return the dual ray
 */
bool LPGlopInterface::HasDualRay() {

  const ProblemStatus status = solver_.GetProblemStatus();

  return status == ProblemStatus::DUAL_UNBOUNDED;
}

/** returns true if LP is proven to be dual unbounded */
bool LPGlopInterface::IsDualUnbounded() {

  const ProblemStatus status = solver_.GetProblemStatus();
  return status == ProblemStatus::DUAL_UNBOUNDED;
}

/** returns true if LP is proven to be dual infeasible */
bool LPGlopInterface::IsDualInfeasible() {

  const ProblemStatus status = solver_.GetProblemStatus();
  return status == ProblemStatus::PRIMAL_UNBOUNDED || status == ProblemStatus::DUAL_INFEASIBLE;
}

/** returns true if LP is proven to be dual feasible */
bool LPGlopInterface::IsDualFeasible() {
  const ProblemStatus status = solver_.GetProblemStatus();

  return status == ProblemStatus::DUAL_FEASIBLE || status == ProblemStatus::OPTIMAL;
}

/** returns true if LP was solved to optimality */
bool LPGlopInterface::IsOptimal() {

  return solver_.GetProblemStatus() == ProblemStatus::OPTIMAL;
}

/** returns true if current LP solution is stable
 *
 *  This function should return true if the solution is reliable, i.e., feasible and optimal (or proven
 *  infeasible/unbounded) with respect to the original problem. The optimality status might be with respect to a scaled
 *  version of the problem, but the solution might not be feasible to the unscaled original problem; in this case,
 *  IsStable() should return false.
 */
bool LPGlopInterface::IsStable() {
  /* For correctness, we need to report "unstable" if Glop was not able to prove optimality because of numerical
   * issues. Currently, Glop still reports primal/dual feasible if at the end, one status is within the tolerance but not
   * the other. */
  const ProblemStatus status = solver_.GetProblemStatus();
  if ((status == ProblemStatus::PRIMAL_FEASIBLE || status == ProblemStatus::DUAL_FEASIBLE) && !IsObjectiveLimitExceeded() && !IsIterationLimitExceeded() && !IsTimeLimitExceeded()) {
    MiniMIPdebugMessage("OPTIMAL not reached and no limit: unstable.\n");
    return false;
  }

  if (status == ProblemStatus::ABNORMAL || status == ProblemStatus::INVALID_PROBLEM || status == ProblemStatus::IMPRECISE)
    return false;
  return true;
} /* @TODO: Case that neither if happens? */

/** returns true if the objective limit was reached */
bool LPGlopInterface::IsObjectiveLimitExceeded() {

  return solver_.objective_limit_reached();
}

/** returns true if the iteration limit was reached */
bool LPGlopInterface::IsIterationLimitExceeded() {
  assert(niterations_ >= static_cast<int>(solver_.GetNumberOfIterations()));

  int maxiter = static_cast<int>(parameters_.max_number_of_iterations());
  return maxiter >= 0 && niterations_ >= maxiter;
}

/** returns true if the time limit was reached */
bool LPGlopInterface::IsTimeLimitExceeded() {

  return lp_time_limit_was_reached_;
}

/** gets objective value of solution */
RetCode LPGlopInterface::GetObjectiveValue(
  LPValue& obj_val /**< stores the objective value */
) {

  obj_val = solver_.GetObjectiveValue();

  return RetCode::OKAY;
}

/** gets primal and dual solution vectors for feasible LPs
 *
 *  Before calling this function, the caller must ensure that the LP has been solved to optimality, i.e., that
 *  IsOptimal() returns true.
 */
RetCode LPGlopInterface::GetSolution(
  LPValue& obj_val,          /**< stores the objective value */
  LPValueArray& primal_sol,  /**< primal solution vector */
  LPValueArray& dual_sol,    /**< dual solution vector */
  LPValueArray& activity,    /**< row activity vector */
  LPValueArray& reduced_cost /**< reduced cost vector */
) {

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

  return RetCode::OKAY;
}

/** gets primal ray for unbounded LPs */
RetCode LPGlopInterface::GetPrimalRay(
  LPValueArray& primal_ray /**< primal ray */
) {

  MiniMIPdebugMessage("GetPrimalRay\n");

  const ColIndex num_cols = linear_program_.num_variables();
  const DenseRow& primal_ray_solver = solver_.GetPrimalRay();
  for (ColIndex col(0); col < num_cols; ++col)
    primal_ray[col.value()] = scaler_.UnscaleVariableValue(col, primal_ray_solver[col]);

  return RetCode::OKAY;
}

/** gets dual Farkas proof for infeasibility */
RetCode LPGlopInterface::GetDualFarkasMultiplier(
  LPValueArray& dual_farkas_multiplier /**< dual Farkas row multipliers */
) {

  MiniMIPdebugMessage("GetDualFarkasMultiplier\n");

  const RowIndex num_rows = linear_program_.num_constraints();
  const DenseColumn& dual_ray = solver_.GetDualRay();
  for (RowIndex row(0); row < num_rows; ++row)
    dual_farkas_multiplier[row.value()] = -scaler_.UnscaleDualValue(row, dual_ray[row]); /* reverse sign */

  return RetCode::OKAY;
}

/** gets the number of LP iterations of the last solve call */
RetCode LPGlopInterface::GetIterations(
  LPNum& iterations /**< number of iterations of the last solve call */
) {

  iterations = static_cast<int>(niterations_);

  return RetCode::OKAY;
}

/*
* LP Basis Methods
*/

/**@name LP Basis Methods */
/**@{ */

/** convert Glop variable basis status to MiniMIP status */
LPBaseStat LPGlopInterface::ConvertGlopVariableStatus(
  VariableStatus status,  /**< variable status */
  Fractional reduced_cost /**< reduced cost of variable */
) {
  switch (status) {
    case VariableStatus::BASIC:
      return LPBaseStat::BASESTAT_BASIC;
    case VariableStatus::AT_UPPER_BOUND:
      return LPBaseStat::BASESTAT_UPPER;
    case VariableStatus::AT_LOWER_BOUND:
      return LPBaseStat::BASESTAT_LOWER;
    case VariableStatus::FREE:
      return LPBaseStat::BASESTAT_ZERO;
    case VariableStatus::FIXED_VALUE:
      return reduced_cost > 0.0 ? LPBaseStat::BASESTAT_LOWER : LPBaseStat::BASESTAT_UPPER;
    default:
      MiniMIPerrorMessage("invalid Glop basis status.\n");
      std::abort();
  }
}

/** convert Glop constraint basis status to MiniMIP status */
LPBaseStat LPGlopInterface::ConvertGlopConstraintStatus(
  ConstraintStatus status, /**< constraint status */
  Fractional dual_value    /**< dual variable value */
) {
  switch (status) {
    case ConstraintStatus::BASIC:
      return LPBaseStat::BASESTAT_BASIC;
    case ConstraintStatus::AT_UPPER_BOUND:
      return LPBaseStat::BASESTAT_UPPER;
    case ConstraintStatus::AT_LOWER_BOUND:
      return LPBaseStat::BASESTAT_LOWER;
    case ConstraintStatus::FREE:
      return LPBaseStat::BASESTAT_ZERO;
    case ConstraintStatus::FIXED_VALUE:
      return dual_value > 0.0 ? LPBaseStat::BASESTAT_LOWER : LPBaseStat::BASESTAT_UPPER;
    default:
      MiniMIPerrorMessage("invalid Glop basis status.\n");
      std::abort();
  }
}

/** Convert MiniMIP variable status to Glop status */
VariableStatus LPGlopInterface::ConvertMiniMIPVariableStatus(
  LPBaseStat status /**< MiniMIP variable status */
) {
  switch (status) {
    case LPBaseStat::BASESTAT_BASIC:
      return VariableStatus::BASIC;
    case LPBaseStat::BASESTAT_UPPER:
      return VariableStatus::AT_UPPER_BOUND;
    case LPBaseStat::BASESTAT_LOWER:
      return VariableStatus::AT_LOWER_BOUND;
    case LPBaseStat::BASESTAT_ZERO:
      return VariableStatus::FREE;
    default:
      MiniMIPerrorMessage("invalid MiniMIP basis status.\n");
      std::abort();
  }
}

/** Convert a MiniMIP constraint status to its corresponding Glop slack VariableStatus.
 *
 *  Note that we swap the upper/lower bounds.
 */
VariableStatus LPGlopInterface::ConvertMiniMIPConstraintStatusToSlackStatus(
  LPBaseStat status /**< MiniMIP constraint status */
) {
  switch (status) {
    case LPBaseStat::BASESTAT_BASIC:
      return VariableStatus::BASIC;
    case LPBaseStat::BASESTAT_UPPER:
      return VariableStatus::AT_LOWER_BOUND;
    case LPBaseStat::BASESTAT_LOWER:
      return VariableStatus::AT_UPPER_BOUND;
    case LPBaseStat::BASESTAT_ZERO:
      return VariableStatus::FREE;
    default:
      MiniMIPerrorMessage("invalid MiniMIP basis status.\n");
      std::abort();
  }
}

/** gets current basis status for columns and rows */
RetCode LPGlopInterface::GetBase(
  LPBaseStatArray& column_basis_status, /**< array to store column basis status, or NULL */
  LPBaseStatArray& row_basis_status     /**< array to store row basis status, or NULL */
) {
  MiniMIPdebugMessage("GetBase\n");

  assert(solver_.GetProblemStatus() == ProblemStatus::OPTIMAL);
  const ColIndex num_cols = linear_program_.num_variables();
  for (ColIndex col(0); col < num_cols; ++col) {
    int i = col.value();
    column_basis_status[i] = (LPBaseStat) ConvertGlopVariableStatus(solver_.GetVariableStatus(col), solver_.GetReducedCost(col));
  }

  const RowIndex num_rows = linear_program_.num_constraints();
  for (RowIndex row(0); row < num_rows; ++row) {
    int i = row.value();
    row_basis_status[i] = (LPBaseStat) ConvertGlopConstraintStatus(solver_.GetConstraintStatus(row), solver_.GetDualValue(row));
  }

  return RetCode::OKAY;
}

/** sets current basis status for columns and rows */
RetCode LPGlopInterface::SetBase(
  const LPBaseStatArray& column_basis_status, /**< array with column basis status */
  const LPBaseStatArray& row_basis_status     /**< array with row basis status */
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

  return RetCode::OKAY;
}

/** returns the indices of the basic columns and rows; basic column n gives value n, basic row m gives value -1-m */
RetCode LPGlopInterface::GetBasisIndices(
  IntArray& basis_indices /**< array to store basis indices ready to keep number of rows entries */
) {
  MiniMIPdebugMessage("GetBasisIndices\n");

  /* the order is important! */
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

  return RetCode::OKAY;
}

/** get row of inverse basis matrix B^-1
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 */
RetCode LPGlopInterface::GetBInvertedRow(
  LPNum row_number,         /**< row number */
  LPValueArray& row_coeffs, /**< array to store the coefficients of the row */
  LPIndexArray& indices,    /**< array to store the non-zero indices */
  int& num_indices          /**< the number of non-zero indices (-1: if we do not store sparsity information) */
) {

  solver_.GetBasisFactorization().LeftSolveForUnitRow(ColIndex(static_cast<int>(row_number)), tmp_row_);
  scaler_.UnscaleUnitRowLeftSolve(solver_.GetBasis(RowIndex(static_cast<int>(row_number))), tmp_row_);

  const ColIndex size = tmp_row_->values.size();
  assert(size.value() == linear_program_.num_constraints());

  /* if we want a sparse vector */
  if (num_indices > 0 && !indices.empty()) {
    num_indices = 0;
    /* Vectors in Glop might be stored in dense or sparse format dep
     *
     * ending on the values. If non_zeros are given, we
     * can directly loop over the non_zeros, otherwise we have to collect the nonzeros. */
    if (!tmp_row_->non_zeros.empty()) {
      ScatteredRowIterator end = tmp_row_->end();
      for (ScatteredRowIterator iter = tmp_row_->begin(); iter != end; ++iter) {
        int idx = (*iter).column().value();
        assert(0 <= idx && idx < linear_program_.num_constraints());
        row_coeffs[idx] = (*iter).coefficient();
        indices[(num_indices)++] = idx;
      }
    } else {
      /* use dense access to tmp_row_ */
      const Fractional eps = parameters_.primal_feasibility_tolerance();
      for (ColIndex col(0); col < size; ++col) {
        LPValue val = (*tmp_row_)[col];
        if (fabs(val) >= eps) {
          row_coeffs[col.value()] = val;
          indices[(num_indices)++] = col.value();
        }
      }
    }
    return RetCode::OKAY;
  }

  /* dense version */
  for (ColIndex col(0); col < size; ++col)
    row_coeffs[col.value()] = (*tmp_row_)[col];

  if (num_indices >= 0)
    num_indices = -1;

  return RetCode::OKAY;
}

/** get column of inverse basis matrix B^-1
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 */
RetCode LPGlopInterface::GetBInvertedColumn(
  LPNum col_number,         /**< column number of B^-1; this is NOT the number of the column in the LP;
                             *   you have to call MiniMIP::LPInterface.GetBasisIndices() to get the array which links the
                             *   B^-1 column numbers to the row and column numbers of the LP!
                             *   c must be between 0 and num_rows-1, since the basis has the size
                             *   num_rows * num_rows */
  LPValueArray& col_coeffs, /**< array to store the coefficients of the column */
  LPIndexArray& indices,    /**< array to store the non-zero indices */
  int& num_indices          /**< the number of non-zero indices (-1: if we do not store sparsity information) */
) {

  /* we need to loop through the rows to extract the values for column col_number */
  const ColIndex col(static_cast<int>(col_number));
  const RowIndex num_rows = linear_program_.num_constraints();

  /* if we want a sparse vector */
  if (num_indices > 0 && !indices.empty()) {
    const Fractional eps = parameters_.primal_feasibility_tolerance();

    num_indices = 0;
    for (int row = 0; row < num_rows; ++row) {
      solver_.GetBasisFactorization().LeftSolveForUnitRow(ColIndex(row), tmp_row_);
      scaler_.UnscaleUnitRowLeftSolve(solver_.GetBasis(RowIndex(row)), tmp_row_);

      LPValue val = (*tmp_row_)[col];
      if (fabs(val) >= eps) {
        col_coeffs[row] = val;
        indices[(num_indices)++] = row;
      }
    }
    return RetCode::OKAY;
  }

  /* dense version */
  for (int row = 0; row < num_rows; ++row) {
    solver_.GetBasisFactorization().LeftSolveForUnitRow(ColIndex(row), tmp_row_);
    scaler_.UnscaleUnitRowLeftSolve(solver_.GetBasis(RowIndex(row)), tmp_row_);
    col_coeffs[row] = (*tmp_row_)[col];
  }

  if (num_indices >= 0)
    num_indices = -1;

  return RetCode::OKAY;
}

/** get row of inverse basis matrix times constraint matrix B^-1 * A
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 */
RetCode LPGlopInterface::GetBInvertedARow(
  LPNum row_number,                   /**< row number */
  const LPValueArray& b_inverted_row, /**< row in (A_B)^-1 from prior call to MiniMIP::LPInterface.GetBInvRow() */
  LPValueArray& row_coeffs,           /**< array to store coefficients of the row */
  LPIndexArray& indices,              /**< array to store the non-zero indices */
  int& num_indices                    /**< thee number of non-zero indices (-1: if we do not store sparsity information) */
) {

  /* get row of basis inverse, loop through columns and muliply with matrix */
  solver_.GetBasisFactorization().LeftSolveForUnitRow(ColIndex(static_cast<int>(row_number)), tmp_row_);
  scaler_.UnscaleUnitRowLeftSolve(solver_.GetBasis(RowIndex(static_cast<int>(row_number))), tmp_row_);

  const ColIndex num_cols = linear_program_.num_variables();

  /* if we want a sparse vector */
  if (num_indices > 0 && !indices.empty()) {
    const Fractional eps = parameters_.primal_feasibility_tolerance();

    num_indices = 0;
    for (ColIndex col(0); col < num_cols; ++col) {
      LPValue val = operations_research::glop::ScalarProduct(tmp_row_->values, linear_program_.GetSparseColumn(col));
      if (fabs(val) >= eps) {
        row_coeffs[col.value()] = val;
        indices[num_indices++] = col.value();
      }
    }
    return RetCode::OKAY;
  }

  /* dense version */
  for (ColIndex col(0); col < num_cols; ++col) {
    LPValue check = operations_research::glop::ScalarProduct(tmp_row_->values, linear_program_.GetSparseColumn(col));
    if (EPS <= fabs(check)) {
      row_coeffs[col.value()] = operations_research::glop::ScalarProduct(tmp_row_->values, linear_program_.GetSparseColumn(col));
    } else {
      row_coeffs[col.value()] = 0;
    }
  }
  num_indices = -1;

  return RetCode::OKAY;
}

/** get column of inverse basis matrix times constraint matrix B^-1 * A
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 */
RetCode LPGlopInterface::GetBInvertedAColumn(
  LPNum col_number,         /**< column number */
  LPValueArray& col_coeffs, /**< array to store coefficients of the column */
  LPIndexArray& indices,    /**< array to store the non-zero indices */
  int& num_indices          /**< the number of non-zero indices (-1: if we do not store sparsity information) */
) {

  solver_.GetBasisFactorization().RightSolveForProblemColumn(ColIndex(static_cast<int>(col_number)), tmp_column_);
  scaler_.UnscaleColumnRightSolve(solver_.GetBasisVector(), ColIndex(static_cast<int>(col_number)), tmp_column_);

  const RowIndex num_rows = tmp_column_->values.size();

  /* if we want a sparse vector */
  if (num_indices > 0 && !indices.empty()) {
    num_indices = 0;
    /* Vectors in Glop might be stored in dense or sparse format depending on the values. If non_zeros are given, we
     * can directly loop over the non_zeros, otherwise we have to collect the nonzeros. */
    if (!tmp_column_->non_zeros.empty()) {
      ScatteredColumnIterator end = tmp_column_->end();
      for (ScatteredColumnIterator iter = tmp_column_->begin(); iter != end; ++iter) {
        int idx = (*iter).row().value();
        assert(0 <= idx && idx < num_rows);
        col_coeffs[idx] = (*iter).coefficient();
        indices[(num_indices)++] = idx;
      }
    } else {
      /* use dense access to tmp_column_ */
      const Fractional eps = parameters_.primal_feasibility_tolerance();
      for (RowIndex row(0); row < num_rows; ++row) {
        LPValue val = (*tmp_column_)[row];
        if (fabs(val) > eps) {
          col_coeffs[row.value()] = val;
          indices[(num_indices)++] = row.value();
        }
      }
    }
    return RetCode::OKAY;
  }

  /* dense version */
  for (RowIndex row(0); row < num_rows; ++row)
    col_coeffs[row.value()] = (*tmp_column_)[row];

  if (num_indices >= 0)
    num_indices = -1;

  return RetCode::OKAY;
}

/**@} */

/*
* Parameter Methods
*/

/**@name Parameter Methods */
/**@{ */

/** gets integer parameter of LP */
RetCode LPGlopInterface::GetIntegerParameter(
  LPParameter type, /**< parameter number */
  LPNum& param_val  /**< buffer to store the parameter value */
) {

  switch (type) {
    case LPParameter::FROM_SCRATCH:
      param_val = (LPNum) from_scratch_;
      MiniMIPdebugMessage("GetIntegerParameter: LPParameter::FROM_SCRATCH = %d.\n", param_val);
      break;
    case LPParameter::LP_INFO:
      param_val = (LPNum) lp_info_;
      MiniMIPdebugMessage("GetIntegerParameter: LPParameter::LP_INFO = %d.\n", param_val);
      break;
    case LPParameter::LP_ITERATION_LIMIT:
      param_val = (LPNum) parameters_.max_number_of_iterations();
      MiniMIPdebugMessage("GetIntegerParameter: LPParameter::LP_ITERATION_LIMIT = %d.\n", param_val);
      break;
    case LPParameter::PRESOLVING:
      param_val = parameters_.use_preprocessing();
      MiniMIPdebugMessage("GetIntegerParameter: LPParameter::PRESOLVING = %d.\n", param_val);
      break;
    case LPParameter::PRICING:
      param_val = (LPNum) pricing_;
      MiniMIPdebugMessage("GetIntegerParameter: LPParameter::PRICING = %d.\n", param_val);
      break;
#ifndef NOSCALING
    case LPParameter::SCALING:
      param_val = parameters_.use_scaling();
      MiniMIPdebugMessage("GetIntegerParameter: LPParameter::SCALING = %d.\n", param_val);
      break;
#endif
    case LPParameter::THREADS:
      param_val = num_threads_;
      MiniMIPdebugMessage("GetIntegerParameter: LPParameter::THREADS = %d.\n", param_val);
      break;
    case LPParameter::TIMING:
      param_val = timing_;
      MiniMIPdebugMessage("GetIntegerParameter: LPParameter::TIMING = %d.\n", param_val);
      break;
    case LPParameter::RANDOMSEED:
      param_val = (LPNum) parameters_.random_seed();
      MiniMIPdebugMessage("GetIntegerParameter: LPParameter::RANDOMSEED = %d.\n", param_val);
      break;
    default:
      return RetCode::PARAMETER_UNKNOWN;
  }

  return RetCode::OKAY;
}

/** sets integer parameter of LP */
RetCode LPGlopInterface::SetIntegerParameter(
  LPParameter type, /**< parameter number */
  LPNum param_val   /**< parameter value */
) {

  switch (type) {
    case LPParameter::FROM_SCRATCH:
      MiniMIPdebugMessage("SetIntegerParameter: LPParameter::FROM_SCRATCH -> %d.\n", param_val);
      from_scratch_ = static_cast<bool>(param_val);
      break;
    case LPParameter::LP_INFO:
      MiniMIPdebugMessage("SetIntegerParameter: LPParameter::LP_INFO -> %d.\n", param_val);
      if (param_val == 0) {
        static_cast<void>(google::SetVLOGLevel("*", google::GLOG_INFO));
        lp_info_ = false;
      } else {
        static_cast<void>(google::SetVLOGLevel("*", google::GLOG_ERROR));
        lp_info_ = true;
      }
      break;
    case LPParameter::LP_ITERATION_LIMIT:
      MiniMIPdebugMessage("SetIntegerParameter: LPParameter::LP_ITERATION_LIMIT -> %d.\n", param_val);
      parameters_.set_max_number_of_iterations(param_val);
      break;
    case LPParameter::PRESOLVING:
      MiniMIPdebugMessage("SetIntegerParameter: LPParameter::PRESOLVING -> %d.\n", param_val);
      parameters_.set_use_preprocessing(param_val);
      break;
    case LPParameter::PRICING:
      MiniMIPdebugMessage("SetIntegerParameter: LPParameter::PRICING -> %d.\n", param_val);
      pricing_ = (LPPricing) param_val;
      switch (pricing_) {
        case LPPricing::DEFAULT:
        case LPPricing::AUTO:
        case LPPricing::PARTIAL:
        case LPPricing::STEEP:
        case LPPricing::STEEPQSTART:
          parameters_.set_feasibility_rule(operations_research::glop::GlopParameters_PricingRule_STEEPEST_EDGE);
          break;
        case LPPricing::FULL:
          /* Dantzig does not really fit, but use it anyway */
          parameters_.set_feasibility_rule(operations_research::glop::GlopParameters_PricingRule_DANTZIG);
          break;
        case LPPricing::DEVEX:
          parameters_.set_feasibility_rule(operations_research::glop::GlopParameters_PricingRule_DEVEX);
          break;
        default:
          return RetCode::PARAMETER_UNKNOWN;
      }
      break;
#ifndef NOSCALING
    case LPParameter::SCALING:
      MiniMIPdebugMessage("SetIntegerParameter: LPParameter::SCALING -> %d.\n", param_val);
      parameters_.set_use_scaling(param_val);
      break;
#endif
    case LPParameter::THREADS:
      MiniMIPdebugMessage("SetIntegerParameter: LPParameter::THREADS -> %d.\n", param_val);
      num_threads_ = param_val;
      if (param_val == 0)
        parameters_.set_num_omp_threads(1);
      else
        parameters_.set_num_omp_threads(static_cast<int>(param_val));
      break;
    case LPParameter::TIMING:
      MiniMIPdebugMessage("SetIntegerParameter: LPParameter::TIMING -> %d.\n", param_val);
      assert(param_val <= 2);
      timing_ = param_val;
      if (param_val == 1)
        absl::SetFlag(&FLAGS_time_limit_use_usertime, true);
      else
        absl::SetFlag(&FLAGS_time_limit_use_usertime, false);
      break;
    case LPParameter::RANDOMSEED:
      MiniMIPdebugMessage("SetIntegerParameter: LPParameter::RANDOMSEED -> %d.\n", param_val);
      parameters_.set_random_seed(static_cast<int>(param_val));
      break;
    default:
      return RetCode::PARAMETER_UNKNOWN;
  }

  return RetCode::OKAY;
}

/** gets floating point parameter of LP */
RetCode LPGlopInterface::GetRealParameter(
  LPParameter type,  /**< parameter number */
  LPValue& param_val /**< buffer to store the parameter value */
) {

  switch (type) {
    case LPParameter::FEASIBLITY_TOLERANCE:
      param_val = parameters_.primal_feasibility_tolerance();
      MiniMIPdebugMessage("GetRealParameter: LPParameter::FEASIBLITY_TOLERANCE = %g.\n", param_val);
      break;
    case LPParameter::DUAL_FEASIBILITY_TOLERANCE:
      param_val = parameters_.dual_feasibility_tolerance();
      MiniMIPdebugMessage("GetRealParameter: LPParameter::DUAL_FEASIBILITY_TOLERANCE = %g.\n", param_val);
      break;
    case LPParameter::OBJECTIVE_LIMIT:
      if (linear_program_.IsMaximizationProblem())
        param_val = parameters_.objective_lower_limit();
      else
        param_val = parameters_.objective_upper_limit();
      MiniMIPdebugMessage("GetRealParameter: LPParameter::OBJECTIVE_LIMIT = %f.\n", param_val);
      break;
    case LPParameter::LP_TIME_LIMIT:
      if (absl::GetFlag(FLAGS_time_limit_use_usertime))
        param_val = parameters_.max_time_in_seconds();
      else
        param_val = parameters_.max_deterministic_time();
      MiniMIPdebugMessage("GetRealParameter: LPParameter::LP_TIME_LIMIT = %f.\n", param_val);
      break;

    default:
      return RetCode::PARAMETER_UNKNOWN;
  }

  return RetCode::OKAY;
}

/** sets floating point parameter of LP */
RetCode LPGlopInterface::SetRealParameter(
  LPParameter type, /**< parameter number */
  LPValue param_val /**< parameter value */
) {

  switch (type) {
    case LPParameter::FEASIBLITY_TOLERANCE:
      MiniMIPdebugMessage("SetRealParameter: LPParameter::FEASIBLITY_TOLERANCE -> %g.\n", param_val);
      parameters_.set_primal_feasibility_tolerance(param_val);
      break;
    case LPParameter::DUAL_FEASIBILITY_TOLERANCE:
      MiniMIPdebugMessage("SetRealParameter: LPParameter::DUAL_FEASIBILITY_TOLERANCE -> %g.\n", param_val);
      parameters_.set_dual_feasibility_tolerance(param_val);
      break;
    case LPParameter::OBJECTIVE_LIMIT:
      MiniMIPdebugMessage("SetRealParameter: LPParameter::OBJECTIVE_LIMIT -> %f.\n", param_val);
      if (linear_program_.IsMaximizationProblem())
        parameters_.set_objective_lower_limit(param_val);
      else
        parameters_.set_objective_upper_limit(param_val);
      break;
    case LPParameter::LP_TIME_LIMIT:
      MiniMIPdebugMessage("SetRealParameter: LPParameter::LP_TIME_LIMIT -> %f.\n", param_val);
      if (absl::GetFlag(FLAGS_time_limit_use_usertime))
        parameters_.set_max_time_in_seconds(param_val);
      else
        parameters_.set_max_deterministic_time(param_val);
      break;
    default:
      return RetCode::PARAMETER_UNKNOWN;
  }

  return RetCode::OKAY;
}

/**@} */

/*
* Numerical Methods
*/

/**@name Numerical Methods */
/**@{ */

/** returns value treated as infinity in the LP solver */
LPValue LPGlopInterface::Infinity() {
  return std::numeric_limits<LPValue>::infinity();
}

/** checks if given value is treated as infinity in the LP solver */
bool LPGlopInterface::IsInfinity(
  LPValue val /**< value to be checked for infinity */
) {
  return val == std::numeric_limits<LPValue>::infinity();
}

/**@} */

/*
* File Interface Methods
*/

/**@name File Interface Methods */
/**@{ */

/** reads LP from a file */
RetCode LPGlopInterface::ReadLP(
  const char* file_name /**< file name */
) {
  assert(file_name != nullptr);

  const char* filespec(file_name);
  MPModelProto proto;
  if (!ReadFileToProto(filespec, &proto)) {
    MiniMIPerrorMessage("Could not read <%s>\n", file_name);
    return RetCode::READ_ERROR;
  }
  linear_program_.Clear();
  MPModelProtoToLinearProgram(proto, &linear_program_);

  return RetCode::OKAY;
}

/** writes LP to a file */
RetCode LPGlopInterface::WriteLP(
  const char* file_name /**< file name */
) {
  assert(file_name != nullptr);

  MPModelProto proto;
  LinearProgramToMPModelProto(linear_program_, &proto);
  const char* filespec(file_name);
  if (!WriteProtoToFile(filespec, proto, operations_research::ProtoWriteFormat::kProtoText, true)) {
    MiniMIPerrorMessage("Could not write <%s>\n", file_name);
    return RetCode::READ_ERROR;
  }

  return RetCode::OKAY;
}

} /* namespace minimip*/
