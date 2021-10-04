#ifndef MINIMIP_SRC_LP_INTERFACE_LPI_H
#define MINIMIP_SRC_LP_INTERFACE_LPI_H

#include "src/lp_interface/lp_types.h"
#include "src/messagehandler/message_handler.h"

#include <string>
#include <vector>

namespace minimip {

// LPInterface base class
//
// Just testing how the change appears.
class LPInterface : private messagehandler {

 public:
  // This is yet another change in a separate commit.
  // Methods
  virtual ~LPInterface() = default;

  // @name Modification Methods
  // @{
  // copies LP data with column matrix into LP solver
  virtual RetCode LoadColumnLP(
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
    ) = 0;

  // adds columns to the LP
  //
  //  @note The indices array is not checked for duplicates, problems may appear if indices are added more than once.
  virtual RetCode AddColumns(
    int num_cols,                       // number of columns to be added
    const std::vector<double>& objective_values, // objective function values of new columns
    const std::vector<double>& lower_bounds,     // lower bounds of new columns
    const std::vector<double>& upper_bounds,     // upper bounds of new columns
    std::vector<std::string>& col_names,               // column names
    int num_non_zeros,                  // number of non-zero elements to be added to the constraint matrix
    const std::vector<int>& begin_cols,       // start index of each column in indices- and vals-array
    const std::vector<int>& indices,          // row indices of constraint matrix entries
    const std::vector<double>& vals              // values of constraint matrix entries
    ) = 0;

  // deletes all columns in the given range from LP
  virtual RetCode DeleteColumns(
    int first_col, // first column to be deleted
    int last_col   // last column to be deleted
    ) = 0;

  // deletes columns from LP; the new position of a column must not be greater than its old position
  virtual RetCode DeleteColumnSet(
    std::vector<bool>& deletion_status // deletion status of columns
    ) = 0;

  // adds rows to the LP
  //
  //  @note The indices array is not checked for duplicates, problems may appear if indices are added more than once.
  virtual RetCode AddRows(
    int num_rows,                       // number of rows to be added
    const std::vector<double>& left_hand_sides,  // left hand sides of new rows
    const std::vector<double>& right_hand_sides, // right hand sides of new rows
    std::vector<std::string>& row_names,               // row names
    int num_non_zeros,                  // number of non-zero elements to be added to the constraint matrix
    const std::vector<int>& begin_rows,       // start index of each row in indices- and vals-array
    const std::vector<int>& indices,          // column indices of constraint matrix entries
    const std::vector<double>& vals              // values of constraint matrix entries
    ) = 0;

  // deletes all rows in the given range from LP
  virtual RetCode DeleteRows(
    int first_row, // first row to be deleted
    int last_row   // last row to be deleted
    ) = 0;

  // deletes rows from LP; the new position of a row must not be greater that its old position
  virtual RetCode DeleteRowSet(
    std::vector<bool>& deletion_status // deletion status of rows
    ) = 0;

  // clears the whole LP
  virtual RetCode Clear() = 0;

  // clears current LPInterface state (like basis information) of the solver
  virtual RetCode ClearState() = 0;

  // changes lower and upper bounds of columns
  virtual RetCode ChangeBounds(
    int num_cols,                   // number of columns to change bounds for
    const std::vector<int>& indices,      // column indices
    const std::vector<double>& lower_bounds, // values for the new lower bounds
    const std::vector<double>& upper_bounds  // values for the new upper bounds
    ) = 0;

  // changes left and right hand sides of rows
  virtual RetCode ChangeSides(
    int num_rows,                      // number of rows to change sides for
    const std::vector<int>& indices,         // row indices
    const std::vector<double>& left_hand_sides, // new values for left hand sides
    const std::vector<double>& right_hand_sides // new values for right hand sides
    ) = 0;

  // changes the objective sense
  virtual RetCode ChangeObjectiveSense(
    LPObjectiveSense obj_sense // new objective sense
    ) = 0;

  // changes objective values of columns in the LP
  virtual RetCode ChangeObjective(
    int num_cols,                  // number of columns to change objective value for
    const std::vector<int>& indices,     // column indices to change objective value for
    const std::vector<double>& new_obj_vals // new objective values for columns
    ) = 0;

  // @}

  // @name Data Accessing Methods
  // @{

  // gets the number of rows in the LP
  virtual int GetNumberOfRows() = 0;

  // gets the number of columns in the LP
  virtual int GetNumberOfColumns() = 0;

  // gets the number of non-zero elements in the LP constraint matrix
  virtual int GetNumberOfNonZeros() = 0;

  // gets the objective sense of the LP
  virtual LPObjectiveSense GetObjectiveSense() = 0;

  // gets columns from LP problem object
  virtual RetCode GetColumns(
    int first_col,            // first column to get from LP
    int last_col,             // last column to get from LP
    std::vector<double>& lower_bounds, // array to store the lower bound vector
    std::vector<double>& upper_bounds, // array to store the upper bound vector
    int& num_non_zeros,       // store the number of non-zero elements
    std::vector<int>& begin_cols,   // array to store start index of each column in indices- and vals-array
    std::vector<int>& indices,      // array to store row indices of constraint matrix entries
    std::vector<double>& vals          // array to store values of constraint matrix entries
    ) = 0;

  // gets rows from LP problem object
  virtual RetCode GetRows(
    int first_row,                // first row to get from LP
    int last_row,                 // last row to get from LP
    std::vector<double>& left_hand_sides,  // array to store left hand side vector
    std::vector<double>& right_hand_sides, // array to store right hand side vector
    int& num_non_zeros,           // store the number of non-zero elements
    std::vector<int>& begin_rows,       // array to store start index of each row in indices- and vals-array
    std::vector<int>& indices,          // array to store column indices of constraint matrix entries
    std::vector<double>& vals              // array to store values of constraint matrix entries
    ) = 0;

  // gets objective coefficients from LP problem object
  virtual RetCode GetObjective(
    int first_col,         // first column to get objective coefficient for
    int last_col,          // last column to get objective coefficient for
    std::vector<double>& obj_coeffs // array to store objective coefficients
    ) = 0;

  // gets current bounds from LP problem object
  virtual RetCode GetBounds(
    int first_col,            // first column to get bounds for
    int last_col,             // last column to get bounds for
    std::vector<double>& lower_bounds, // array to store lower bound values
    std::vector<double>& upper_bounds  // array to store upper bound values
    ) = 0;

  // gets current row sides from LP problem object
  virtual RetCode GetSides(
    int first_row,               // first row to get sides for
    int last_row,                // last row to get sides for
    std::vector<double>& left_hand_sides, // array to store left hand side values
    std::vector<double>& right_hand_sides // array to store right hand side values
    ) = 0;

  // gets a single coefficient
  virtual RetCode GetCoefficient(
    int row,   // row number of coefficient
    int col,   // column number of coefficient
    double& val // array to store the value of the coefficient
    ) = 0;

  // @}

  // @name Solving Methods
  // @{

  // calls primal simplex to solve the LP
  virtual RetCode SolvePrimal() = 0;

  // calls dual simplex to solve the LP
  virtual RetCode SolveDual() = 0;

  // start strong branching - call before any strong branching
  virtual RetCode StartStrongbranch() = 0;

  // end strong branching - call after any strong branching
  virtual RetCode EndStrongbranch() = 0;

  // performs strong branching iterations on one @b fractional candidate
  virtual RetCode StrongbranchFractionalValue(
    int col,                       // column to apply strong branching on
    double primal_sol,              // fractional current primal solution value of column
    int iteration_limit,           // iteration limit for strong branchings
    double& dual_bound_down_branch, // stores dual bound after branching column down
    double& dual_bound_up_branch,   // stores dual bound after branching column up
    bool& down_valid,                // whether the returned down value is a valid dual bound; otherwise, it can only be used as an estimate value
    bool& up_valid,                  // whether the returned up value is a valid dual bound; otherwise, it can only be used as an estimate value
    int& iterations                // stores total number of strong branching iterations
    ) = 0;

  // performs strong branching iterations on given @b fractional candidates
  virtual RetCode StrongbranchFractionalValues(
    std::vector<int>& cols,                       // columns to apply strong branching on
    int num_cols,                         // number of columns
    std::vector<double>& primal_sols,              // fractional current primal solution values of columns
    int iteration_limit,                  // iteration limit for strong branchings
    std::vector<double>& dual_bound_down_branches, // stores dual bounds after branching columns down
    std::vector<double>& dual_bound_up_branches,   // stores dual bounds after branching columns up
    std::vector<bool>& down_valids,                 // stores whether the returned down values are valid dual bounds;
                                            //  otherwise, they can only be used as an estimate values
    std::vector<bool>& up_valids,                   // stores whether the returned up values are a valid dual bounds;
                                            //  otherwise, they can only be used as an estimate values
    int& iterations                       // stores total number of strong branching iterations
    ) = 0;

  // performs strong branching iterations on one candidate with @b integral value
  virtual RetCode StrongbranchIntegerValue(
    int col,                       // column to apply strong branching on
    double primal_sol,              // current integral primal solution value of column
    int iteration_limit,           // iteration limit for strong branchings
    double& dual_bound_down_branch, // stores dual bound after branching column down
    double& dual_bound_up_branch,   // stores dual bound after branching column up
    bool& down_valid,                // stores whether the returned down value is a valid dual bound;
                                     //     *   otherwise, it can only be used as an estimate value
    bool& up_valid,                  // stores whether the returned up value is a valid dual bound;
                                     //     *   otherwise, it can only be used as an estimate value
    int& iterations                // stores total number of strong branching iterations
    ) = 0;

  // performs strong branching iterations on given candidates with @b integral values
  virtual RetCode StrongbranchIntegerValues(
    std::vector<int>& cols,                       // columns to apply strong branching on
    int num_cols,                         // number of columns
    std::vector<double>& primal_sols,              // current integral primal solution values of columns
    int iteration_limit,                  // iteration limit for strong branchings
    std::vector<double>& dual_bound_down_branches, // stores dual bounds after branching columns down
    std::vector<double>& dual_bound_up_branches,   // stores dual bounds after branching columns up
    std::vector<bool>& down_valids,                 // stores whether the returned down values are valid dual bounds;
                                            //  otherwise, they can only be used as an estimate values
    std::vector<bool>& up_valids,                   // stores whether the returned up values are a valid dual bounds;
                                            //  otherwise, they can only be used as an estimate values
    int& iterations                       // stores total number of strong branching iterations
    ) = 0;
  // @}

  // @name Solution Information Methods
  // @{

  // returns whether a solve method was called after the last modification of the LP
  virtual bool WasSolved() = 0;

  // returns true if LP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
  //  this does not necessarily mean that the solver knows and can return the primal ray
  virtual bool ExistsPrimalRay() = 0;

  // returns true if LP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
  //  and the solver knows and can return the primal ray
  virtual bool HasPrimalRay() = 0;

  // returns true if LP is proven to be primal unbounded
  virtual bool IsPrimalUnbounded() = 0;

  // returns true if LP is proven to be primal infeasible
  virtual bool IsPrimalInfeasible() = 0;

  // returns true if LP is proven to be primal feasible
  virtual bool IsPrimalFeasible() = 0;

  // returns true if LP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
  // this does not necessarily mean that the solver knows and can return the dual ray
  virtual bool ExistsDualRay() = 0;

  // returns true if LP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
  // and the solver knows and can return the dual ray
  virtual bool HasDualRay() = 0;

  // returns true if LP is proven to be dual unbounded
  virtual bool IsDualUnbounded() = 0;

  // returns true if LP is proven to be dual infeasible
  virtual bool IsDualInfeasible() = 0;

  // returns true if LP is proven to be dual feasible
  virtual bool IsDualFeasible() = 0;

  // returns true if LP was solved to optimality
  virtual bool IsOptimal() = 0;

  // returns true if current LP solution is stable
  //
  // This function should return true if the solution is reliable, i.e., feasible and optimal (or proven
  // infeasible/unbounded) with respect to the original problem. The optimality status might be with respect to a scaled
  // version of the problem, but the solution might not be feasible to the unscaled original problem; in this case,
  // MiniMIP::LPInterface.IsStable() should return false.
  virtual bool IsStable() = 0;

  // returns true if the objective limit was reached
  virtual bool IsObjectiveLimitExceeded() = 0;

  // returns true if the iteration limit was reached
  virtual bool IsIterationLimitExceeded() = 0;

  // returns true if the time limit was reached
  virtual bool IsTimeLimitExceeded() = 0;

  // gets objective value of solution
  virtual RetCode GetObjectiveValue(
    double& obj_val // the objective value
    ) = 0;

  // gets primal and dual solution vectors for feasible LPs
  //
  // Before calling this function, the caller must ensure that the LP has been solved to optimality, i.e., that
  // MiniMIP::LPInterface.IsOptimal() returns true.
  virtual RetCode GetSolution(
    double& obj_val,          // stores the objective value
    std::vector<double>& primal_sol,  // primal solution vector
    std::vector<double>& dual_sol,    // dual solution vector
    std::vector<double>& activity,    // row activity vector
    std::vector<double>& reduced_cost // reduced cost vector
    ) = 0;

  // gets primal ray for unbounded LPs
  virtual RetCode GetPrimalRay(
    std::vector<double>& primal_ray // primal ray
    ) = 0;

  // gets dual Farkas proof for infeasibility
  virtual RetCode GetDualFarkasMultiplier(
    std::vector<double>& dual_farkas_multiplier // dual Farkas row multipliers
    ) = 0;

  // gets the number of LP iterations of the last solve call
  virtual RetCode GetIterations(
    int& iterations // number of iterations of the last solve call
    ) = 0;

  // @}

  // @name LP Basis Methods
  // @{

  // gets current basis status for columns and rows
  virtual RetCode GetBase(
    std::vector<LPBasisStatus>& column_basis_status, // array to store column basis status
    std::vector<LPBasisStatus>& row_basis_status     // array to store row basis status
    ) = 0;

  // sets current basis status for columns and rows
  virtual RetCode SetBase(
    const std::vector<LPBasisStatus>& column_basis_status, // array with column basis status
    const std::vector<LPBasisStatus>& row_basis_status     // array with row basis status
    ) = 0;

  // returns the indices of the basic columns and rows; basic column n gives value n, basic row m gives value -1-m
  virtual RetCode GetBasisIndices(
    std::vector<int>& basis_indices // array to store basis indices ready to keep number of rows entries
    ) = 0;

  // get row of inverse basis matrix B^-1
  //
  // @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
  //       uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
  //       see also the explanation in lpi.h.
  virtual RetCode GetBInvertedRow(
    int row_number,         // row number
    std::vector<double>& row_coeffs, // array to store the coefficients of the row
    std::vector<int>& indices,    // array to store the non-zero indices
    int& num_indices          // the number of non-zero indices (-1: if we do not store sparsity information)
    ) = 0;

  // get column of inverse basis matrix B^-1
  //
  // @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
  //       uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated
  virtual RetCode GetBInvertedColumn(
    int col_number,         // column number of B^-1; this is NOT the number of the column in the LP;
                               // you have to call MiniMIP::LPInterface.GetBasisIndices() to get the array which links the
                               // B^-1 column numbers to the row and column numbers of the LP!
                               // c must be between 0 and num_rows-1, since the basis has the size
                               // num_rows * num_rows
    std::vector<double>& col_coeffs, // array to store the coefficients of the column
    std::vector<int>& indices,    // array to store the non-zero indices
    int& num_indices          // the number of non-zero indices (-1: if we do not store sparsity information)
    ) = 0;

  // get row of inverse basis matrix times constraint matrix B^-1 * A
  //
  // @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
  //       uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
  //       see also the explanation in lpi.h.
  virtual RetCode GetBInvertedARow(
    int row_number,                   // row number
    const std::vector<double>& b_inverted_row, // row in (A_B)^-1 from prior call to MiniMIP::LPInterface.GetBInvRow()
    std::vector<double>& row_coeffs,           // array to store coefficients of the row
    std::vector<int>& indices,              // array to store the non-zero indices
    int& num_indices                    // thee number of non-zero indices (-1: if we do not store sparsity information)
    ) = 0;

  // get column of inverse basis matrix times constraint matrix B^-1 * A
  //
  // @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
  //       uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
  //       see also the explanation in lpi.h.
  virtual RetCode GetBInvertedAColumn(
    int col_number,         // column number
    std::vector<double>& col_coeffs, // array to store coefficients of the column
    std::vector<int>& indices,    // array to store the non-zero indices
    int& num_indices          // the number of non-zero indices (-1: if we do not store sparsity information)
    ) = 0;

  // @}

  // @name Parameter Methods
  // @{

  // gets integer parameter of LP
  virtual RetCode GetIntegerParameter(
    LPParameter type, // parameter number
    int& param_val  // returns the parameter value
    ) = 0;

  // sets integer parameter of LP
  virtual RetCode SetIntegerParameter(
    LPParameter type, // parameter number
    int param_val   // parameter value
    ) = 0;

  // gets floating point parameter of LP
  virtual RetCode GetRealParameter(
    LPParameter type,  // parameter number
    double& param_val // returns the parameter value
    ) = 0;

  // sets floating point parameter of LP
  virtual RetCode SetRealParameter(
    LPParameter type, // parameter number
    double param_val // parameter value
    ) = 0;

  // @}

  // @name Numerical Methods
  // @{

  // returns value treated as infinity in the LP solver
  virtual double Infinity() = 0;

  // checks if given value is treated as infinity in the LP solver
  virtual bool IsInfinity(
    double val // value to be checked for infinity
    ) = 0;

  // @}

  // @name File Interface Methods
  // @{

  // reads LP from a file
  virtual RetCode ReadLP(
    const char* file_name // file name
    ) = 0;

  // writes LP to a file
  virtual RetCode WriteLP(
    const char* file_name // file name
    ) = 0;

  // @}
};

} // namespace minimip
#endif // MINIMIP_SRC_LP_INTERFACE_LPI_H