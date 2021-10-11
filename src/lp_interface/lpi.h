#ifndef MINIMIP_SRC_LP_INTERFACE_LPI_H
#define MINIMIP_SRC_LP_INTERFACE_LPI_H

#include "src/lp_interface/lp_types.h"
#include "src/messagehandler/message_handler.h"

#include "absl/status/status.h"
// #include "absl/status/statusor.h"

#include <string>
#include <vector>

namespace minimip {

// LPInterface abstract class. This is used by MiniMIP to communicate with an
// underlying LP solver.
class LPInterface : private messagehandler {

 public:
  virtual ~LPInterface() = default;

  // ==========================================================================
  // LP model setters.
  // ==========================================================================

  // copies LP data with column matrix into LP solver
  virtual absl::Status LoadColumnLP(
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
  // NOTE: The indices array is not checked for duplicates, problems may appear if indices are added more than once.
  virtual absl::Status AddColumns(
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
  virtual absl::Status DeleteColumns(
    int first_col, // first column to be deleted
    int last_col   // last column to be deleted
    ) = 0;

  // adds rows to the LP
  //
  // NOTE: The indices array is not checked for duplicates, problems may appear if indices are added more than once.
  virtual absl::Status AddRows(
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
  virtual absl::Status DeleteRows(
    int first_row, // first row to be deleted
    int last_row   // last row to be deleted
    ) = 0;

  // deletes rows from LP; the new position of a row must not be greater that its old position
  virtual absl::Status DeleteRowSet(
    std::vector<bool>& deletion_status // deletion status of rows
    ) = 0;

  // clears the whole LP
  virtual absl::Status Clear() = 0;

  // clears current LPInterface state (like basis information) of the solver
  virtual absl::Status ClearState() = 0;

  // changes lower and upper bounds of columns
  virtual absl::Status SetColumnBounds(
    int num_cols,                   // number of columns to change bounds for
    const std::vector<int>& indices,      // column indices
    const std::vector<double>& lower_bounds, // values for the new lower bounds
    const std::vector<double>& upper_bounds  // values for the new upper bounds
    ) = 0;

  // changes left and right hand sides of rows
  virtual absl::Status SetRowSides(
    int num_rows,                      // number of rows to change sides for
    const std::vector<int>& indices,         // row indices
    const std::vector<double>& left_hand_sides, // new values for left hand sides
    const std::vector<double>& right_hand_sides // new values for right hand sides
    ) = 0;

  // changes the objective sense
  virtual absl::Status SetObjectiveSense(
    LPObjectiveSense obj_sense // new objective sense
    ) = 0;

  // changes objective values of columns in the LP
  virtual absl::Status SetObjectiveCoefficients(
    int num_cols,                  // number of columns to change objective value for
    const std::vector<int>& indices,     // column indices to change objective value for
    const std::vector<double>& objective_coefficients // new objective values for columns
    ) = 0;

  // ==========================================================================
  // LP model getters.
  // ==========================================================================

  // gets the number of rows in the LP
  virtual int GetNumberOfRows() const = 0;

  // gets the number of columns in the LP
  virtual int GetNumberOfColumns() const = 0;

  // gets the number of non-zero elements in the LP constraint matrix
  virtual int GetNumberOfNonZeros() const = 0;

  // gets the objective sense of the LP
  virtual LPObjectiveSense GetObjectiveSense() const = 0;

  // gets columns from LP problem object
  virtual absl::Status GetColumns(
    int first_col,            // first column to get from LP
    int last_col,             // last column to get from LP
    std::vector<double>& lower_bounds, // array to store the lower bound vector
    std::vector<double>& upper_bounds, // array to store the upper bound vector
    int& num_non_zeros,       // store the number of non-zero elements
    std::vector<int>& begin_cols,   // array to store start index of each column in indices- and vals-array
    std::vector<int>& indices,      // array to store row indices of constraint matrix entries
    std::vector<double>& vals          // array to store values of constraint matrix entries
    ) const = 0;

  // gets rows from LP problem object
  virtual absl::Status GetRows(
    int first_row,                // first row to get from LP
    int last_row,                 // last row to get from LP
    std::vector<double>& left_hand_sides,  // array to store left hand side vector
    std::vector<double>& right_hand_sides, // array to store right hand side vector
    int& num_non_zeros,           // store the number of non-zero elements
    std::vector<int>& begin_rows,       // array to store start index of each row in indices- and vals-array
    std::vector<int>& indices,          // array to store column indices of constraint matrix entries
    std::vector<double>& vals              // array to store values of constraint matrix entries
    ) const = 0;

  // gets objective coefficients from LP problem object
  virtual absl::Status GetObjective(
    int first_col,         // first column to get objective coefficient for
    int last_col,          // last column to get objective coefficient for
    std::vector<double>& obj_coeffs // array to store objective coefficients
    ) const = 0;

  // gets current bounds from LP problem object
  virtual absl::Status GetBounds(
    int first_col,            // first column to get bounds for
    int last_col,             // last column to get bounds for
    std::vector<double>& lower_bounds, // array to store lower bound values
    std::vector<double>& upper_bounds  // array to store upper bound values
    ) const = 0;

  // gets current row sides from LP problem object
  virtual absl::Status GetSides(
    int first_row,               // first row to get sides for
    int last_row,                // last row to get sides for
    std::vector<double>& left_hand_sides, // array to store left hand side values
    std::vector<double>& right_hand_sides // array to store right hand side values
    ) const = 0;

  // gets a single coefficient
  virtual double GetMatrixCoefficient(
    int col,   // column number of coefficient
    int row   // row number of coefficient
    ) const = 0;

  // ==========================================================================
  // Solving methods.
  // ==========================================================================

  // calls primal simplex to solve the LP
  virtual absl::Status SolveLpWithPrimalSimplex() = 0;

  // calls dual simplex to solve the LP
  virtual absl::Status SolveLpWithDualSimplex() = 0;

  // start strong branching - call before any strong branching
  virtual absl::Status StartStrongbranch() = 0;

  // end strong branching - call after any strong branching
  virtual absl::Status EndStrongbranch() = 0;

  // performs strong branching iterations on one @b fractional candidate
  virtual absl::Status StrongbranchFractionalValue(
    int col,                       // column to apply strong branching on
    double primal_sol,              // fractional current primal solution value of column
    int iteration_limit,           // iteration limit for strong branchings
    double& dual_bound_down_branch, // stores dual bound after branching column down
    double& dual_bound_up_branch,   // stores dual bound after branching column up
    bool& down_valid,                // whether the returned down value is a valid dual bound; otherwise, it can only be used as an estimate value
    bool& up_valid,                  // whether the returned up value is a valid dual bound; otherwise, it can only be used as an estimate value
    int& iterations                // stores total number of strong branching iterations
    ) = 0;

  // performs strong branching iterations on one candidate with @b integral value
  virtual absl::Status StrongbranchIntegerValue(
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

  // ==========================================================================
  // Solution information getters.
  // ==========================================================================

  // returns whether a solve method was called after the last modification of the LP
  virtual bool IsSolved() const = 0;

  // returns true if current LP solution is stable
  //
  // This function should return true if the solution is reliable, i.e., feasible and optimal (or proven
  // infeasible/unbounded) with respect to the original problem. The optimality status might be with respect to a scaled
  // version of the problem, but the solution might not be feasible to the unscaled original problem; in this case,
  // minimip::LPInterface.IsStable() should return false.
  virtual bool IsStable() const = 0;

  // returns true if LP was solved to optimality
  virtual bool IsOptimal() const = 0;

  // returns true if LP is proven to be primal feasible
  virtual bool IsPrimalFeasible() const = 0;

  // returns true if LP is proven to be primal infeasible
  virtual bool IsPrimalInfeasible() const = 0;

  // returns true if LP is proven to be primal unbounded
  virtual bool IsPrimalUnbounded() const = 0;

  // returns true if LP is proven to be dual feasible
  virtual bool IsDualFeasible() const = 0;

  // returns true if LP is proven to be dual infeasible
  virtual bool IsDualInfeasible() const = 0;

  // returns true if LP is proven to be dual unbounded
  virtual bool IsDualUnbounded() const = 0;

  // returns true if LP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
  //  this does not necessarily mean that the solver knows and can return the primal ray
  virtual bool ExistsPrimalRay() const = 0;

  // returns true if LP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
  //  and the solver knows and can return the primal ray
  virtual bool HasPrimalRay() const = 0;
  
  // returns true if LP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
  // this does not necessarily mean that the solver knows and can return the dual ray
  virtual bool ExistsDualRay() const = 0;

  // returns true if LP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
  // and the solver knows and can return the dual ray
  virtual bool HasDualRay() const = 0;

  // returns true if the objective limit was reached
  virtual bool ObjectiveLimitIsExceeded() const = 0;

  // returns true if the iteration limit was reached
  virtual bool IterationLimitIsExceeded() const = 0;

  // returns true if the time limit was reached
  virtual bool TimeLimitIsExceeded() const = 0;

  // gets objective value of solution
  virtual absl::Status GetObjectiveValue(
    double& obj_val // the objective value
    ) = 0;

  // gets primal and dual solution vectors for feasible LPs
  //
  // Before calling this function, the caller must ensure that the LP has been solved to optimality, i.e., that
  // minimip::LPInterface.IsOptimal() returns true.
  virtual absl::Status GetSolution(
    double& obj_val,          // stores the objective value
    std::vector<double>& primal_sol,  // primal solution vector
    std::vector<double>& dual_sol,    // dual solution vector
    std::vector<double>& activity,    // row activity vector
    std::vector<double>& reduced_cost // reduced cost vector
    ) const = 0;

  // gets primal ray for unbounded LPs
  virtual absl::Status GetPrimalRay(
    std::vector<double>& primal_ray // primal ray
    ) const = 0;

  // gets dual Farkas proof for infeasibility
  virtual absl::Status GetDualFarkasMultiplier(
    std::vector<double>& dual_farkas_multiplier // dual Farkas row multipliers
    ) const = 0;

  // gets the number of LP iterations of the last solve call
  virtual int GetIterations() const = 0;

  // ==========================================================================
  // Getters and setters of the basis.
  // ==========================================================================

  // gets current basis status for columns and rows
  virtual absl::Status GetBase(
    std::vector<LPBasisStatus>& column_basis_status, // array to store column basis status
    std::vector<LPBasisStatus>& row_basis_status     // array to store row basis status
    ) const = 0;

  // sets current basis status for columns and rows
  virtual absl::Status SetBase(
    const std::vector<LPBasisStatus>& column_basis_status, // array with column basis status
    const std::vector<LPBasisStatus>& row_basis_status     // array with row basis status
    ) = 0;

  // returns the indices of the basic columns and rows; basic column n gives value n, basic row m gives value -1-m
  virtual absl::Status GetBasisIndices(
    std::vector<int>& basis_indices // array to store basis indices ready to keep number of rows entries
    ) const = 0;


  // ==========================================================================
  // Getters of vectors in the inverted basis matrix.
  // ==========================================================================

  // get row of inverse basis matrix B^-1
  //
  // NOTE: The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
  //       uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
  //       see also the explanation in lpi.h.
  virtual absl::Status GetBInvertedRow(
    int row_number,         // row number
    std::vector<double>& row_coeffs, // array to store the coefficients of the row
    std::vector<int>& indices,    // array to store the non-zero indices
    int& num_indices          // the number of non-zero indices (-1: if we do not store sparsity information)
    ) const = 0;

  // get column of inverse basis matrix B^-1
  //
  // NOTE: The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
  //       uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated
  virtual absl::Status GetBInvertedColumn(
    int col_number,         // column number of B^-1; this is NOT the number of the column in the LP;
                               // you have to call minimip::LPInterface.GetBasisIndices() to get the array which links the
                               // B^-1 column numbers to the row and column numbers of the LP!
                               // c must be between 0 and num_rows-1, since the basis has the size
                               // num_rows * num_rows
    std::vector<double>& col_coeffs, // array to store the coefficients of the column
    std::vector<int>& indices,    // array to store the non-zero indices
    int& num_indices          // the number of non-zero indices (-1: if we do not store sparsity information)
    ) const = 0;

  // get row of inverse basis matrix times constraint matrix B^-1 * A
  //
  // NOTE: The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
  //       uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
  //       see also the explanation in lpi.h.
  virtual absl::Status GetBInvertedARow(
    int row_number,                   // row number
    const std::vector<double>& b_inverted_row, // row in (A_B)^-1 from prior call to minimip::LPInterface.GetBInvRow()
    std::vector<double>& row_coeffs,           // array to store coefficients of the row
    std::vector<int>& indices,              // array to store the non-zero indices
    int& num_indices                    // thee number of non-zero indices (-1: if we do not store sparsity information)
    ) const = 0;

  // get column of inverse basis matrix times constraint matrix B^-1 * A
  //
  // NOTE: The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
  //       uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
  //       see also the explanation in lpi.h.
  virtual absl::Status GetBInvertedAColumn(
    int col_number,         // column number
    std::vector<double>& col_coeffs, // array to store coefficients of the column
    std::vector<int>& indices,    // array to store the non-zero indices
    int& num_indices          // the number of non-zero indices (-1: if we do not store sparsity information)
    ) const = 0;

  // ==========================================================================
  // Getters and setters of the parameters.
  // ==========================================================================


  // gets integer parameter of LP
  virtual absl::Status GetIntegerParameter(
    LPParameter type, // parameter number
    int& param_val  // returns the parameter value
    ) const = 0;

  // sets integer parameter of LP
  virtual absl::Status SetIntegerParameter(
    LPParameter type, // parameter number
    int param_val   // parameter value
    ) = 0;

  // gets floating point parameter of LP
  virtual absl::Status GetRealParameter(
    LPParameter type,  // parameter number
    double& param_val // returns the parameter value
    ) const = 0;

  // sets floating point parameter of LP
  virtual absl::Status SetRealParameter(
    LPParameter type, // parameter number
    double param_val // parameter value
    ) = 0;

  // ==========================================================================
  // Numerical methods.
  // ==========================================================================

  // returns value treated as infinity in the LP solver
  virtual double Infinity() const = 0;

  // checks if given value is treated as infinity in the LP solver
  virtual bool IsInfinity(
    double value // value to be checked for infinity
    ) const = 0;

  // ==========================================================================
  // File interface methods.
  // ==========================================================================

  // reads LP from a file
  virtual absl::Status ReadLP(
    const char* file_name // file name
    ) = 0;

  // writes LP to a file
  virtual absl::Status WriteLP(
    const char* file_name // file name
    ) const = 0;

};

} // namespace minimip
#endif // MINIMIP_SRC_LP_INTERFACE_LPI_H