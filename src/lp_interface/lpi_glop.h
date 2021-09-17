#ifndef MINIMIP_SRC_LP_INTERFACE_LPI_GLOP_H_
#define MINIMIP_SRC_LP_INTERFACE_LPI_GLOP_H_

#include "src/lp_interface/lp_types.h"
#include "src/lp_interface/lpi.h"
#include "src/minimip/minimip_def.h"

/* turn off some warnings from includes */
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wignored-qualifiers"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#pragma GCC diagnostic ignored "-Wctor-dtor-privacy"
#pragma GCC diagnostic ignored "-Woverflow"

#include "ortools/lp_data/lp_data_utils.h"
#include "ortools/lp_data/lp_print_utils.h"
#include "ortools/lp_data/proto_utils.h"
#include "ortools/util/file_util.h"
#include "ortools/util/stats.h"

/* turn warnings on again */
#pragma GCC diagnostic warning "-Wsign-compare"
#pragma GCC diagnostic warning "-Wpedantic"
#pragma GCC diagnostic warning "-Wignored-qualifiers"
#pragma GCC diagnostic warning "-Wshadow"
#pragma GCC diagnostic warning "-Wnon-virtual-dtor"
#pragma GCC diagnostic warning "-Wctor-dtor-privacy"
#pragma GCC diagnostic warning "-Woverflow"

#include <cassert>
#include <string>
#include <vector>

using operations_research::TimeLimit;
using operations_research::glop::ConstraintStatus;
using operations_research::glop::DenseBooleanColumn;
using operations_research::glop::Fractional;
using operations_research::glop::GlopParameters;
using operations_research::glop::LinearProgram;
using operations_research::glop::LpScalingHelper;
using operations_research::glop::ProblemStatus;
using operations_research::glop::RevisedSimplex;
using operations_research::glop::ScatteredColumn;
using operations_research::glop::ScatteredRow;
using operations_research::glop::VariableStatus;

namespace minimip {

/** LPInterface base class */
class LPGlopInterface : public LPInterface {
 private:
  LinearProgram linear_program_; /**< the linear program */
  LinearProgram scaled_lp_;      /**< scaled linear program */
  RevisedSimplex solver_;        /**< direct reference to the revised simplex, not passing through lp_solver */
  GlopParameters parameters_;    /**< parameters */
  LpScalingHelper scaler_;       /**< scaler auxiliary class */

  /* the following is used by WasSolved() */
  bool lp_modified_since_last_solve_;
  bool lp_time_limit_was_reached_;

  /* store the values of some parameters in order to be able to return them */
  bool lp_info_;      /**< whether additional output is turned on */
  LPPricing pricing_; /**< MiniMIP pricing setting  */
  bool from_scratch_; /**< store whether basis is ignored for next solving call */
  LPNum num_threads_; /**< number of threads used to solve the LP (0 = automatic) */
  LPNum timing_;      /**< type of timer (1 - cpu, 2 - wallclock, 0 - off) */

  /* other data */
  LPLongInt niterations_; /**< number of iterations used */

  /* Temporary vectors allocated here for speed. This gain is non-negligible
   * because in many situations, only a few entries of these vectors are
   * inspected (hypersparsity) and allocating them is in O(num_rows) or
   * O(num_cols) instead of O(num_non_zeros) to read/clear them. */
  ScatteredRow* tmp_row_;       /**< temporary vector */
  ScatteredColumn* tmp_column_; /**< temporary vector */

  /** delete rows from LP and update the current basis */
  void DeleteRowsAndUpdateCurrentBasis(
    const DenseBooleanColumn& rows_to_delete /**< array to mark rows that should be deleted */
  );

  /** update scaled linear program */
  void updateScaledLP();

  /** check primal feasibility */
  bool checkUnscaledPrimalFeasibility();

  /** common function between the two LPI Solve() functions */
  RetCode SolveInternal(
    bool recursive,                        /**< Is this a recursive call? */
    std::unique_ptr<TimeLimit>& time_limit /**< time limit */
  );

  /** determine whether the dual bound is valid */
  bool IsDualBoundValid(
    ProblemStatus status /**< status to be checked */
  );

  /** performs strong branching iterations */
  RetCode strongbranch(
    LPIndex col_index,               /**< column to apply strong branching on */
    LPValue primal_sol,              /**< fractional current primal solution value of column */
    LPNum iteration_limit,           /**< iteration limit for strong branchings */
    LPValue& dual_bound_down_branch, /**< stores dual bound after branching column down */
    LPValue& dual_bound_up_branch,   /**< stores dual bound after branching column up */
    bool& down_valid,                /**< stores whether the returned down value is a valid dual bound;
                                      *   otherwise, it can only be used as an estimate value */
    bool& up_valid,                  /**< stores whether the returned up value is a valid dual bound;
                                      *   otherwise, it can only be used as an estimate value */
    LPNum& iterations                /**< stores total number of strong branching iterations, or -1; may be NULL */
  );

  /*
   * LP Basis Methods
   */

  /**@name LP Basis Methods */
  /**@{ */

  /** convert Glop variable basis status to MiniMIP status */
  LPBaseStat ConvertGlopVariableStatus(
    VariableStatus status,  /**< variable status */
    Fractional reduced_cost /**< reduced cost of variable */
  );

  /** convert Glop constraint basis status to MiniMIP status */
  LPBaseStat ConvertGlopConstraintStatus(
    ConstraintStatus status, /**< constraint status */
    Fractional dual_value    /**< dual variable value */
  );

  /** Convert MiniMIP variable status to Glop status */
  VariableStatus ConvertMiniMIPVariableStatus(
    LPBaseStat status /**< MiniMIP variable status */
  );

  /** Convert a MiniMIP constraint status to its corresponding Glop slack VariableStatus.
  *
  *  Note that we swap the upper/lower bounds.
  */
  VariableStatus ConvertMiniMIPConstraintStatusToSlackStatus(
    LPBaseStat status /**< MiniMIP constraint status */
  );

 public:
  /* constructor */
  LPGlopInterface();

  /* copy constructor */
  // implicitly declared - default copy constructor given by compiler

  /* Destructor */
  ~LPGlopInterface() override;

  /* Methods */
  /**@name Modification Methods */
  /**@{ */

  /** copies LP data with column matrix into LP solver */
  RetCode LoadColumnLP(
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
    ) override;

  /** adds columns to the LP
    *
    *  @note The indices array is not checked for duplicates, problems may appear if indices are added more than once.
    */
  RetCode AddColumns(
    LPNum num_cols,                       /**< number of columns to be added */
    const LPValueArray& objective_values, /**< objective function values of new columns */
    const LPValueArray& lower_bounds,     /**< lower bounds of new columns */
    const LPValueArray& upper_bounds,     /**< upper bounds of new columns */
    StringArray& col_names,               /**< column names */
    LPNum num_non_zeros,                  /**< number of non-zero elements to be added to the constraint matrix */
    const LPIndexArray& begin_cols,       /**< start index of each column in indices- and vals-array */
    const LPIndexArray& indices,          /**< row indices of constraint matrix entries */
    const LPValueArray& vals              /**< values of constraint matrix entries */
    ) override;

  /** deletes all columns in the given range from LP */
  RetCode DeleteColumns(
    LPNum first_col, /**< first column to be deleted */
    LPNum last_col   /**< last column to be deleted */
    ) override;

  /** deletes columns from LP; the new position of a column must not be greater than its old position */
  RetCode DeleteColumnSet(
    BoolArray& deletion_status /**< deletion status of columns */
    ) override;

  /** adds rows to the LP
    *
    *  @note The indices array is not checked for duplicates, problems may appear if indices are added more than once.
    */
  RetCode AddRows(
    LPNum num_rows,                       /**< number of rows to be added */
    const LPValueArray& left_hand_sides,  /**< left hand sides of new rows */
    const LPValueArray& right_hand_sides, /**< right hand sides of new rows */
    StringArray& row_names,               /**< row names */
    LPNum num_non_zeros,                  /**< number of non-zero elements to be added to the constraint matrix */
    const LPIndexArray& begin_rows,       /**< start index of each row in indices- and vals-array */
    const LPIndexArray& indices,          /**< column indices of constraint matrix entries */
    const LPValueArray& vals              /**< values of constraint matrix entries */
    ) override;

  /** deletes all rows in the given range from LP */
  RetCode DeleteRows(
    LPNum first_row, /**< first row to be deleted */
    LPNum last_row   /**< last row to be deleted */
    ) override;

  /** deletes rows from LP; the new position of a row must not be greater that its old position */
  RetCode DeleteRowSet(
    BoolArray& deletion_status /**< deletion status of rows */
    ) override;

  /** clears the whole LP */
  RetCode Clear() override;

  /** clears current LPInterface state (like basis information) of the solver */
  RetCode ClearState() override;

  /** changes lower and upper bounds of columns */
  RetCode ChangeBounds(
    LPNum num_cols,                   /**< number of columns to change bounds for */
    const LPIndexArray& indices,      /**< column indices */
    const LPValueArray& lower_bounds, /**< values for the new lower bounds */
    const LPValueArray& upper_bounds  /**< values for the new upper bounds */
    ) override;

  /** changes left and right hand sides of rows */
  RetCode ChangeSides(
    LPNum num_rows,                      /**< number of rows to change sides for */
    const LPIndexArray& indices,         /**< row indices */
    const LPValueArray& left_hand_sides, /**< new values for left hand sides */
    const LPValueArray& right_hand_sides /**< new values for right hand sides */
    ) override;

  /** changes the objective sense */
  RetCode ChangeObjectiveSense(
    LPObjectiveSense obj_sense /**< new objective sense */
    ) override;

  /** changes objective values of columns in the LP */
  RetCode ChangeObjective(
    LPNum num_cols,                  /**< number of columns to change objective value for */
    const LPIndexArray& indices,     /**< column indices to change objective value for */
    const LPValueArray& new_obj_vals /**< new objective values for columns */
    ) override;

  /**@} */

  /**@name Data Accessing Methods */
  /**@{ */

  /** gets the number of rows in the LP */
   LPNum GetNumberOfRows() override;

  /** gets the number of columns in the LP */
   LPNum GetNumberOfColumns() override;

  /** gets the number of non-zero elements in the LP constraint matrix */
   LPNum GetNumberOfNonZeros() override;

  /** gets the objective sense of the LP */
   LPObjectiveSense GetObjectiveSense() override;

  /** gets columns from LP problem object; the arrays have to be large enough to store all values; */
  RetCode GetColumns(
    LPNum first_col,            /**< first column to get from LP */
    LPNum last_col,             /**< last column to get from LP */
    LPValueArray& lower_bounds, /**< array to store the lower bound vector */
    LPValueArray& upper_bounds, /**< array to store the upper bound vector */
    LPNum& num_non_zeros,       /**< store the number of non-zero elements */
    LPIndexArray& begin_cols,   /**< array to store start index of each column in indices- and vals-array */
    LPIndexArray& indices,      /**< array to store row indices of constraint matrix entries */
    LPValueArray& vals          /**< array to store values of constraint matrix entries */
    ) override;

  /** gets rows from LP problem object; the arrays have to be large enough to store all values. */
  RetCode GetRows(
    LPNum first_row,                /**< first row to get from LP */
    LPNum last_row,                 /**< last row to get from LP */
    LPValueArray& left_hand_sides,  /**< array to store left hand side vector */
    LPValueArray& right_hand_sides, /**< array to store right hand side vector */
    LPNum& num_non_zeros,           /**< store the number of non-zero elements */
    LPIndexArray& begin_rows,       /**< array to store start index of each row in indices- and vals-array */
    LPIndexArray& indices,          /**< array to store column indices of constraint matrix entries */
    LPValueArray& vals              /**< array to store values of constraint matrix entries */
    ) override;

  /** gets objective coefficients from LP problem object */
  RetCode GetObjective(
    LPNum first_col,         /**< first column to get objective coefficient for */
    LPNum last_col,          /**< last column to get objective coefficient for */
    LPValueArray& obj_coeffs /**< array to store objective coefficients */
    ) override;

  /** gets current bounds from LP problem object */
  RetCode GetBounds(
    LPNum first_col,            /**< first column to get bounds for */
    LPNum last_col,             /**< last column to get bounds for */
    LPValueArray& lower_bounds, /**< array to store lower bound values */
    LPValueArray& upper_bounds  /**< array to store upper bound values */
    ) override;

  /** gets current row sides from LP problem object */
  RetCode GetSides(
    LPNum first_row,               /**< first row to get sides for */
    LPNum last_row,                /**< last row to get sides for */
    LPValueArray& left_hand_sides, /**< array to store left hand side values */
    LPValueArray& right_hand_sides /**< array to store right hand side values */
    ) override;

  /** gets a single coefficient */
  RetCode GetCoefficient(
    LPNum row,   /**< row number of coefficient */
    LPNum col,   /**< column number of coefficient */
    LPValue& val /**< array to store the value of the coefficient */
    ) override;

  /**@} */

  /*
  * Solving Methods
  */

  /**@name Solving Methods */
  /**@{ */

  /** calls primal simplex to solve the LP */
  RetCode SolvePrimal() override;

  /** calls dual simplex to solve the LP */
  RetCode SolveDual() override;

  /** start strong branching - call before any strong branching */
  RetCode StartStrongbranch() override;

  /** end strong branching - call after any strong branching */
  RetCode EndStrongbranch() override;

  /** performs strong branching iterations on one @b fractional candidate */
  RetCode StrongbranchFractionalValue(
    LPNum col,                 /**< column to apply strong branching on */
    LPValue primal_sol,              /**< fractional current primal solution value of column */
    LPNum iteration_limit,           /**< iteration limit for strong branchings */
    LPValue& dual_bound_down_branch, /**< stores dual bound after branching column down */
    LPValue& dual_bound_up_branch,   /**< stores dual bound after branching column up */
    bool& down_valid,                /**< whether the returned down value is a valid dual bound; otherwise, it can only be used as an estimate value */
    bool& up_valid,                  /**< whether the returned up value is a valid dual bound; otherwise, it can only be used as an estimate value */
    LPNum& iterations                /**< stores total number of strong branching iterations */
    ) override;

  /** performs strong branching iterations on given @b fractional candidates */
  RetCode StrongbranchFractionalValues(
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
    ) override;

  /** performs strong branching iterations on one candidate with @b integral value */
  RetCode StrongbranchIntegerValue(
    LPNum col,                       /**< column to apply strong branching on */
    LPValue primal_sol,              /**< current integral primal solution value of column */
    LPNum iteration_limit,           /**< iteration limit for strong branchings */
    LPValue& dual_bound_down_branch, /**< stores dual bound after branching column down */
    LPValue& dual_bound_up_branch,   /**< stores dual bound after branching column up */
    bool& down_valid,                /**< stores whether the returned down value is a valid dual bound;
                                      *   otherwise, it can only be used as an estimate value */
    bool& up_valid,                  /**< stores whether the returned up value is a valid dual bound;
                                      *   otherwise, it can only be used as an estimate value */
    LPNum& iterations                /**< stores total number of strong branching iterations */
    ) override;

  /** performs strong branching iterations on given candidates with @b integral values */
  RetCode StrongbranchIntegerValues(
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
    ) override;
  /**@} */

  /**@name Solution Information Methods */
  /**@{ */

  /** returns whether a solve method was called after the last modification of the LP */
  bool WasSolved() override;

  /** returns true if LP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
    *  this does not necessarily mean that the solver knows and can return the primal ray
    */
  bool ExistsPrimalRay() override;

  /** returns true if LP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
    *  and the solver knows and can return the primal ray
    */
  bool HasPrimalRay() override;

  /** returns true if LP is proven to be primal unbounded */
  bool IsPrimalUnbounded() override;

  /** returns TRUE if LP is proven to be primal infeasible */
  bool IsPrimalInfeasible() override;

  /** returns TRUE if LP is proven to be primal feasible */
  bool IsPrimalFeasible() override;

  /** returns TRUE if LP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
   *  this does not necessarily mean that the solver knows and can return the dual ray
   */
  bool ExistsDualRay() override;

  /** returns TRUE if LP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
   *  and the solver knows and can return the dual ray
   */
  bool HasDualRay() override;

  /** returns TRUE if LP is proven to be dual unbounded */
  bool IsDualUnbounded() override;

  /** returns TRUE if LP is proven to be dual infeasible */
  bool IsDualInfeasible() override;

  /** returns TRUE if LP is proven to be dual feasible */
  bool IsDualFeasible() override;

  /** returns TRUE if LP was solved to optimality */
  bool IsOptimal() override;

  /** returns TRUE if current LP solution is stable
   *
   *  This function should return true if the solution is reliable, i.e., feasible and optimal (or proven
   *  infeasible/unbounded) with respect to the original problem. The optimality status might be with respect to a scaled
   *  version of the problem, but the solution might not be feasible to the unscaled original problem; in this case,
   *  MiniMIP::LPInterface.IsStable() should return false.
   */
  bool IsStable() override;

  /** returns TRUE if the objective limit was reached */
  bool IsObjectiveLimitExceeded() override;

  /** returns TRUE if the iteration limit was reached */
  bool IsIterationLimitExceeded() override;

  /** returns TRUE if the time limit was reached */
  bool IsTimeLimitExceeded() override;

  /** gets objective value of solution */
  RetCode GetObjectiveValue(
    LPValue& obj_val /**< the objective value */
    ) override;

  /** gets primal and dual solution vectors for feasible LPs
   *
   *  Before calling this function, the caller must ensure that the LP has been solved to optimality, i.e., that
   *  MiniMIP::LPInterface.IsOptimal() returns true.
   */
  RetCode GetSolution(
    LPValue& obj_val,          /**< stores the objective value */
    LPValueArray& primal_sol,  /**< primal solution vector */
    LPValueArray& dual_sol,    /**< dual solution vector */
    LPValueArray& activity,    /**< row activity vector */
    LPValueArray& reduced_cost /**< reduced cost vector */
    ) override;

  /** gets primal ray for unbounded LPs */
  RetCode GetPrimalRay(
    LPValueArray& primal_ray /**< primal ray */
    ) override;

  /** gets dual Farkas proof for infeasibility */
  RetCode GetDualFarkasMultiplier(
    LPValueArray& dual_farkas_multiplier /**< dual Farkas row multipliers */
    ) override;

  /** gets the number of LP iterations of the last solve call */
  RetCode GetIterations(
    LPNum& iterations /**< number of iterations of the last solve call */
    ) override;

  /**@} */

  /**@name LP Basis Methods */
  /**@{ */

  /** gets current basis status for columns and rows */
  RetCode GetBase(
    LPBaseStatArray& column_basis_status, /**< array to store column basis status */
    LPBaseStatArray& row_basis_status     /**< array to store row basis status */
    ) override;

  /** sets current basis status for columns and rows */
  RetCode SetBase(
    const LPBaseStatArray& column_basis_status, /**< array with column basis status */
    const LPBaseStatArray& row_basis_status     /**< array with row basis status */
    ) override;

  /** returns the indices of the basic columns and rows; basic column n gives value n, basic row m gives value -1-m */
  RetCode GetBasisIndices(
    IntArray& basis_indices /**< array to store basis indices ready to keep number of rows entries */
    ) override;

  /** get row of inverse basis matrix B^-1
   *
   *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
   *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
   *        see also the explanation in lpi.h.
   */
  RetCode GetBInvertedRow(
    LPNum row_number,         /**< row number */
    LPValueArray& row_coeffs, /**< array to store the coefficients of the row */
    LPIndexArray& indices,    /**< array to store the non-zero indices */
    int& num_indices          /**< the number of non-zero indices (-1: if we do not store sparsity information) */
    ) override;

  /** get column of inverse basis matrix B^-1
   *
   *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
   *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated
   */
  RetCode GetBInvertedColumn(
    LPNum col_number,         /**< column number of B^-1; this is NOT the number of the column in the LP;
                               *   you have to call MiniMIP::LPInterface.GetBasisIndices() to get the array which links the
                               *   B^-1 column numbers to the row and column numbers of the LP!
                               *   c must be between 0 and num_rows-1, since the basis has the size
                               *   num_rows * num_rows */
    LPValueArray& col_coeffs, /**< array to store the coefficients of the column */
    LPIndexArray& indices,    /**< array to store the non-zero indices */
    int& num_indices          /**< the number of non-zero indices (-1: if we do not store sparsity information) */
    ) override;

  /** get row of inverse basis matrix times constraint matrix B^-1 * A
   *
   *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
   *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
   *        see also the explanation in lpi.h.
   */
  RetCode GetBInvertedARow(
    LPNum row_number,                   /**< row number */
    const LPValueArray& b_inverted_row, /**< row in (A_B)^-1 from prior call to MiniMIP::LPInterface.GetBInvRow() */
    LPValueArray& row_coeffs,           /**< array to store coefficients of the row */
    LPIndexArray& indices,              /**< array to store the non-zero indices */
    int& num_indices                    /**< thee number of non-zero indices (-1: if we do not store sparsity information) */
    ) override;

  /** get column of inverse basis matrix times constraint matrix B^-1 * A
   *
   *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
   *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
   *        see also the explanation in lpi.h.
   */
  RetCode GetBInvertedAColumn(
    LPNum col_number,         /**< column number */
    LPValueArray& col_coeffs, /**< array to store coefficients of the column */
    LPIndexArray& indices,    /**< array to store the non-zero indices */
    int& num_indices          /**< the number of non-zero indices (-1: if we do not store sparsity information) */
    ) override;

  /**@} */

  /**@name Parameter Methods */
  /**@{ */

  /** gets integer parameter of LP */
  RetCode GetIntegerParameter(
    LPParameter type, /**< parameter number */
    LPNum& param_val  /**< returns the parameter value */
    ) override;

  /** sets integer parameter of LP */
  RetCode SetIntegerParameter(
    LPParameter type, /**< parameter number */
    LPNum param_val   /**< parameter value */
    ) override;

  /** gets floating point parameter of LP */
  RetCode GetRealParameter(
    LPParameter type,  /**< parameter number */
    LPValue& param_val /**< returns the parameter value */
    ) override;

  /** sets floating point parameter of LP */
  RetCode SetRealParameter(
    LPParameter type, /**< parameter number */
    LPValue param_val /**< parameter value */
    ) override;

  /**@} */

  /**@name Numerical Methods */
  /**@{ */

  /** returns value treated as infinity in the LP solver */
  LPValue Infinity() override;

  /** checks if given value is treated as infinity in the LP solver */
  bool IsInfinity(
    LPValue val /**< value to be checked for infinity */
    ) override;

  /**@} */

  /**@name File Interface Methods */
  /**@{ */

  /** reads LP from a file */
  RetCode ReadLP(
    const char* file_name /**< file name */
    ) override;

  /** writes LP to a file */
  RetCode WriteLP(
    const char* file_name /**< file name */
    ) override;

  /**@} */

};// class LPGlopInterface
} /* namespace minimip*/

#endif//MINIMIP_SRC_LP_INTERFACE_LPI_GLOP_H_
