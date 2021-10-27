#ifndef SRC_LP_INTERFACE_LPI_GLOP_H_
#define SRC_LP_INTERFACE_LPI_GLOP_H_

#include "src/lp_interface/lp_types.h"
#include "src/lp_interface/lpi.h"
#include "src/minimip/minimip_def.h"

// turn off some warnings from includes
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

// turn warnings on again
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
using operations_research::glop::LPScalingHelper;
using operations_research::glop::ProblemStatus;
using operations_research::glop::RevisedSimplex;
using operations_research::glop::ScatteredColumn;
using operations_research::glop::ScatteredRow;
using operations_research::glop::VariableStatus;

namespace minimip {

// LPInterface base class
class LPGlopInterface : public LPInterface {
 private:
  LinearProgram linear_program_;  // the linear program
  LinearProgram scaled_lp_;       // scaled linear program
  RevisedSimplex solver_;      // direct reference to the revised simplex, not
                               // passing through lp_solver
  GlopParameters parameters_;  // parameters
  LPScalingHelper scaler_;     // scaler auxiliary class

  // the following is used by IsSolved()
  bool lp_modified_since_last_solve_;
  bool lp_time_limit_was_reached_;

  // store the values of some parameters in order to be able to return them
  bool lp_info_;       // whether additional output is turned on
  LPPricing pricing_;  // MiniMIP pricing setting
  bool from_scratch_;  // store whether basis is ignored for next solving call
  int num_threads_;    // number of threads used to solve the LP (0 = automatic)
  int timing_;         // type of timer (1 - cpu, 2 - wallclock, 0 - off)

  // other data
  long long niterations_;  // number of iterations used

  // Temporary vectors allocated here for speed. This gain is non-negligible
  // because in many situations, only a few entries of these vectors are
  // inspected (hypersparsity) and allocating them is in O(num_rows) or
  // O(num_cols) instead of O(num_non_zeros) to read/clear them.
  ScatteredRow* tmp_row_;        // temporary vector
  ScatteredColumn* tmp_column_;  // temporary vector

  // delete rows from LP and update the current basis
  void DeleteRowsAndUpdateCurrentBasis(
      const DenseBooleanColumn&
          rows_to_delete  // array to mark rows that should be deleted
  );

  // update scaled linear program
  void updateScaledLP();

  // check primal feasibility
  bool checkUnscaledPrimalFeasibility() const;

  // common function between the two LPI Solve() functions
  absl::Status SolveInternal(
      bool recursive,                         // Is this a recursive call?
      std::unique_ptr<TimeLimit>& time_limit  // time limit
  );

  // determine whether the dual bound is valid
  bool IsDualBoundValid(ProblemStatus status  // status to be checked
  ) const;

  // performs strong branching iterations
  absl::Status strongbranch(
      int col,      // column to apply strong branching on
      double primal_sol,  // fractional current primal solution value of column
      int iteration_limit,             // iteration limit for strong branchings
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
      int& iterations  // stores total number of strong branching iterations, or
                       // -1; may be NULL
  );

  // LP Basis Methods

  // @name LP Basis Methods
  // @{

  // convert Glop variable basis status to MiniMIP status
  LPBasisStatus ConvertGlopVariableStatus(
      VariableStatus status,   // variable status
      Fractional reduced_cost  // reduced cost of variable
  ) const;

  // convert Glop constraint basis status to MiniMIP status
  LPBasisStatus ConvertGlopConstraintStatus(
      ConstraintStatus status,  // constraint status
      Fractional dual_value     // dual variable value
  ) const;

  // Convert MiniMIP variable status to Glop status
  VariableStatus ConvertMiniMIPVariableStatus(
      LPBasisStatus status  // MiniMIP variable status
  ) const;

  // Convert a MiniMIP constraint status to its corresponding Glop slack
  // VariableStatus.
  //
  // Note that we swap the upper/lower bounds.
  VariableStatus ConvertMiniMIPConstraintStatusToSlackStatus(
      LPBasisStatus status  // MiniMIP constraint status
  ) const;

 public:
  LPGlopInterface();

  ~LPGlopInterface() override;

  // ==========================================================================
  // LP model setters.
  // ==========================================================================

  // copies LP data with column matrix into LP solver
  absl::Status LoadColumnLP(
      LPObjectiveSense obj_sense,  // objective sense
      int num_cols,                // number of columns
      const std::vector<double>&
          objective_values,  // objective function values of columns
      const std::vector<double>& lower_bounds,      // lower bounds of columns
      const std::vector<double>& upper_bounds,      // upper bounds of columns
      std::vector<std::string>& col_names,          // column names
      int num_rows,                                 // number of rows
      const std::vector<double>& left_hand_sides,   // left hand sides of rows
      const std::vector<double>& right_hand_sides,  // right hand sides of rows
      std::vector<std::string>& row_names,          // row names
      int num_non_zeros,  // number of non-zero elements in the constraint
                          // matrix
      const std::vector<int>& begin_cols,  // start index of each column in
                                           // row_indices- and vals-array
      const std::vector<int>&
          row_indices,  // row indices of constraint matrix entries
      const std::vector<double>& vals  // values of constraint matrix entries
      ) override;

  // adds columns to the LP
  //
  // NOTE: The indices array is not checked for duplicates, problems may appear
  // if indices are added more than once.
  absl::Status AddColumns(
      int num_cols,  // number of columns to be added
      const std::vector<double>&
          objective_values,  // objective function values of new columns
      const std::vector<double>& lower_bounds,  // lower bounds of new columns
      const std::vector<double>& upper_bounds,  // upper bounds of new columns
      std::vector<std::string>& col_names,      // column names
      int num_non_zeros,  // number of non-zero elements to be added to the
                          // constraint matrix
      const std::vector<int>&
          begin_cols,  // start index of each column in indices- and vals-array
      const std::vector<int>&
          indices,  // row indices of constraint matrix entries
      const std::vector<double>& vals  // values of constraint matrix entries
      ) override;

  // deletes all columns in the given range from LP
  absl::Status DeleteColumns(int first_col,  // first column to be deleted
                             int last_col    // last column to be deleted
                             ) override;

  // adds rows to the LP
  //
  // NOTE: The indices array is not checked for duplicates, problems may appear
  // if indices are added more than once.
  absl::Status AddRows(
      int num_rows,  // number of rows to be added
      const std::vector<double>&
          left_hand_sides,  // left hand sides of new rows
      const std::vector<double>&
          right_hand_sides,                 // right hand sides of new rows
      std::vector<std::string>& row_names,  // row names
      int num_non_zeros,  // number of non-zero elements to be added to the
                          // constraint matrix
      const std::vector<int>&
          begin_rows,  // start index of each row in indices- and vals-array
      const std::vector<int>&
          indices,  // column indices of constraint matrix entries
      const std::vector<double>& vals  // values of constraint matrix entries
      ) override;

  // deletes all rows in the given range from LP
  absl::Status DeleteRows(int first_row,  // first row to be deleted
                          int last_row    // last row to be deleted
                          ) override;

  // deletes rows from LP; the new position of a row must not be greater that
  // its old position
  absl::Status DeleteRowSet(
      std::vector<bool>& deletion_status  // deletion status of rows
      ) override;

  // clears the whole LP
  absl::Status Clear() override;

  // clears current LPInterface state (like basis information) of the solver
  absl::Status ClearState() override;

  // changes lower and upper bounds of columns
  absl::Status SetColumnBounds(
      int num_cols,  // number of columns to change bounds for
      const std::vector<int>& indices,  // column indices
      const std::vector<double>&
          lower_bounds,  // values for the new lower bounds
      const std::vector<double>&
          upper_bounds  // values for the new upper bounds
      ) override;

  // changes left and right hand sides of rows
  absl::Status SetRowSides(
      int num_rows,                     // number of rows to change sides for
      const std::vector<int>& indices,  // row indices
      const std::vector<double>&
          left_hand_sides,  // new values for left hand sides
      const std::vector<double>&
          right_hand_sides  // new values for right hand sides
      ) override;

  // changes the objective sense
  absl::Status SetObjectiveSense(
      LPObjectiveSense obj_sense  // new objective sense
      ) override;

  // changes objective values of columns in the LP
  absl::Status SetObjectiveCoefficients(
      int num_cols,  // number of columns to change objective value for
      const std::vector<int>&
          indices,  // column indices to change objective value for
      const std::vector<double>&
          objective_coefficients  // new objective values for columns
      ) override;

  // ==========================================================================
  // LP model getters.
  // ==========================================================================

  // gets the number of rows in the LP
  int GetNumberOfRows() const override;

  // gets the number of columns in the LP
  int GetNumberOfColumns() const override;

  // gets the number of non-zero elements in the LP constraint matrix
  int GetNumberOfNonZeros() const override;

  // gets the objective sense of the LP
  LPObjectiveSense GetObjectiveSense() const override;

  // gets the sparse coefficients of the column from LP problem object
  SparseVector GetSparseColumnCoefficients(int col) const override;

  // gets the sparse coefficients of the row from LP problem object
  SparseVector GetSparseRowCoefficients(int row) const override;

  // gets objective coefficient of column from LP problem object
  double GetObjectiveCoefficient(int col) const override;

  // gets current lower bound of column from LP problem object
  double GetLowerBound(int col) const override;

  // gets current upper bound of column from LP problem object
  double GetUpperBound(int col) const override;

  // gets current left hand side of row from LP problem object
  double GetLeftHandSide(int row) const override;

  // gets current right hand side of row from LP problem object
  double GetRightHandSide(int row) const override;

  // gets the matrix coefficient of column and row from LP problem object
  double GetMatrixCoefficient(int col,  // column number of coefficient
                                      int row   // row number of coefficient
  ) const override;

  // ==========================================================================
  // Solving methods.
  // ==========================================================================

  // calls primal simplex to solve the LP
  absl::Status SolveLPWithPrimalSimplex() override;

  // calls dual simplex to solve the LP
  absl::Status SolveLPWithDualSimplex() override;

  // start strong branching - call before any strong branching
  absl::Status StartStrongBranching() override;

  // end strong branching - call after any strong branching
  absl::Status EndStrongBranching() override;

  // performs strong branching iterations on one branching candidate
  absl::StatusOr<StrongBranchResult> StrongBranchValue(
      int col,              // column to apply strong branching on
      double primal_sol,    // current primal solution value of column
      int iteration_limit  // iteration limit for strong branchings
      ) override;

  // ==========================================================================
  // Solution information getters.
  // ==========================================================================

  // returns whether a solve method was called after the last modification of
  // the LP
  bool IsSolved() const override;

  // returns true if current LP solution is stable
  //
  // This function should return true if the solution is reliable, i.e.,
  // feasible and optimal (or proven infeasible/unbounded) with respect to the
  // original problem. The optimality status might be with respect to a scaled
  // version of the problem, but the solution might not be feasible to the
  // unscaled original problem; in this case, minimip::LPInterface.IsStable()
  // should return false.
  bool IsStable() const override;

  // returns true if LP was solved to optimality
  bool IsOptimal() const override;

  // returns true if LP is proven to be primal feasible
  bool IsPrimalFeasible() const override;

  // returns true if LP is proven to be primal infeasible
  bool IsPrimalInfeasible() const override;

  // returns true if LP is proven to be primal unbounded
  bool IsPrimalUnbounded() const override;

  // returns true if LP is proven to be dual feasible
  bool IsDualFeasible() const override;

  // returns true if LP is proven to be dual infeasible
  bool IsDualInfeasible() const override;

  // returns true if LP is proven to be dual unbounded
  bool IsDualUnbounded() const override;

  // returns true if LP is proven to have a primal unbounded ray (but not
  // necessary a primal feasible point);
  //  this does not necessarily mean that the solver knows and can return the
  //  primal ray
  bool ExistsPrimalRay() const override;

  // returns true if LP is proven to have a primal unbounded ray (but not
  // necessary a primal feasible point),
  //  and the solver knows and can return the primal ray
  bool HasPrimalRay() const override;

  // returns true if LP is proven to have a dual unbounded ray (but not
  // necessary a dual feasible point); this does not necessarily mean that the
  // solver knows and can return the dual ray
  bool ExistsDualRay() const override;

  // returns true if LP is proven to have a dual unbounded ray (but not
  // necessary a dual feasible point), and the solver knows and can return the
  // dual ray
  bool HasDualRay() const override;

  // returns true if the objective limit was reached
  bool ObjectiveLimitIsExceeded() const override;

  // returns true if the iteration limit was reached
  bool IterationLimitIsExceeded() const override;

  // returns true if the time limit was reached
  bool TimeLimitIsExceeded() const override;

  // gets objective value of solution
  absl::Status GetObjectiveValue(double& obj_val  // the objective value
                                 ) override;

  // gets primal and dual solution vectors for feasible LPs
  //
  // Before calling this function, the caller must ensure that the LP has been
  // solved to optimality, i.e., that minimip::LPInterface.IsOptimal() returns
  // true.
  absl::Status GetSolution(
      double& obj_val,                   // stores the objective value
      std::vector<double>& primal_sol,   // primal solution vector
      std::vector<double>& dual_sol,     // dual solution vector
      std::vector<double>& activity,     // row activity vector
      std::vector<double>& reduced_cost  // reduced cost vector
  ) const override;

  // gets primal ray for unbounded LPs
  absl::Status GetPrimalRay(std::vector<double>& primal_ray  // primal ray
  ) const override;

  // gets dual Farkas proof for infeasibility
  absl::Status GetDualFarkasMultiplier(
      std::vector<double>&
          dual_farkas_multiplier  // dual Farkas row multipliers
  ) const override;

  // gets the number of LP iterations of the last solve call
  int GetIterations() const override;

  // ==========================================================================
  // Getters and setters of the basis.
  // ==========================================================================

  // gets current basis status for columns and rows
  absl::Status GetBase(
      std::vector<LPBasisStatus>&
          column_basis_status,  // array to store column basis status
      std::vector<LPBasisStatus>&
          row_basis_status  // array to store row basis status
  ) const override;

  // sets current basis status for columns and rows
  absl::Status SetBase(
      const std::vector<LPBasisStatus>&
          column_basis_status,  // array with column basis status
      const std::vector<LPBasisStatus>&
          row_basis_status  // array with row basis status
      ) override;

  // returns the indices of the basic columns and rows; basic column n gives
  // value n, basic row m gives value -1-m
  absl::Status GetBasisIndices(
      std::vector<int>& basis_indices  // array to store basis indices ready to
                                       // keep number of rows entries
  ) const override;

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
  absl::Status GetRowOfBInverted(
      int row_number,  // row number
      std::vector<double>&
          row_coeffs,             // array to store the coefficients of the row
      std::vector<int>& indices,  // array to store the non-zero indices
      int& num_indices  // the number of non-zero indices (-1: if we do not
                        // store sparsity information)
  ) const override;

  // get column of inverse basis matrix B^-1
  //
  // NOTE: The LP interface defines slack variables to have coefficient +1. This
  // means that if, internally, the LP solver
  //       uses a -1 coefficient, then rows associated with slacks variables
  //       whose coefficient is -1, should be negated
  absl::Status GetColumnOfBInverted(
      int col_number,  // column number of B^-1; this is NOT the number of the
                       // column in the LP; you have to call
                       // minimip::LPInterface.GetBasisIndices() to get the
                       // array which links the B^-1 column numbers to the row
                       // and column numbers of the LP! c must be between 0 and
                       // num_rows-1, since the basis has the size num_rows *
                       // num_rows
      std::vector<double>&
          col_coeffs,  // array to store the coefficients of the column
      std::vector<int>& indices,  // array to store the non-zero indices
      int& num_indices  // the number of non-zero indices (-1: if we do not
                        // store sparsity information)
  ) const override;

  // get row of inverse basis matrix times constraint matrix B^-1 * A
  //
  // NOTE: The LP interface defines slack variables to have coefficient +1. This
  // means that if, internally, the LP solver
  //       uses a -1 coefficient, then rows associated with slacks variables
  //       whose coefficient is -1, should be negated; see also the explanation
  //       in lpi.h.
  absl::Status GetRowOfBInvertedTimesA(
      int row_number,  // row number
      std::vector<double>&
          row_coeffs,             // array to store coefficients of the row
      std::vector<int>& indices,  // array to store the non-zero indices
      int& num_indices  // thee number of non-zero indices (-1: if we do not
                        // store sparsity information)
  ) const override;

  // get column of inverse basis matrix times constraint matrix B^-1 * A
  //
  // NOTE: The LP interface defines slack variables to have coefficient +1. This
  // means that if, internally, the LP solver
  //       uses a -1 coefficient, then rows associated with slacks variables
  //       whose coefficient is -1, should be negated; see also the explanation
  //       in lpi.h.
  absl::Status GetColumnOfBInvertedTimesA(
      int col_number,  // column number
      std::vector<double>&
          col_coeffs,             // array to store coefficients of the column
      std::vector<int>& indices,  // array to store the non-zero indices
      int& num_indices  // the number of non-zero indices (-1: if we do not
                        // store sparsity information)
  ) const override;

  // ==========================================================================
  // Getters and setters of the parameters.
  // ==========================================================================

  // gets integer parameter of LP
  absl::Status GetIntegerParameter(
      LPParameter type,  // parameter number
      int& param_val     // returns the parameter value
  ) const override;

  // sets integer parameter of LP
  absl::Status SetIntegerParameter(LPParameter type,  // parameter number
                                   int param_val      // parameter value
                                   ) override;

  // gets floating point parameter of LP
  absl::Status GetRealParameter(
      LPParameter type,  // parameter number
      double& param_val  // returns the parameter value
  ) const override;

  // sets floating point parameter of LP
  absl::Status SetRealParameter(LPParameter type,  // parameter number
                                double param_val   // parameter value
                                ) override;

  // ==========================================================================
  // Numerical methods.
  // ==========================================================================

  // returns value treated as infinity in the LP solver
  double Infinity() const override;

  // checks if given value is treated as infinity in the LP solver
  bool IsInfinity(double value  // value to be checked for infinity
  ) const override;

  // ==========================================================================
  // File interface methods.
  // ==========================================================================

  // reads LP from a file
  absl::Status ReadLP(const char* file_name  // file name
                      ) override;

  // writes LP to a file
  absl::Status WriteLP(const char* file_name  // file name
  ) const override;
};  // class LPGlopInterface

}  // namespace minimip

#endif  // SRC_LP_INTERFACE_LPI_GLOP_H_
