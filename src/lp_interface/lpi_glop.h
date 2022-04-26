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

#ifndef SRC_LP_INTERFACE_LPI_GLOP_H_
#define SRC_LP_INTERFACE_LPI_GLOP_H_

#include <cinttypes>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "src/lp_interface/lp_types.h"
#include "src/lp_interface/lpi.h"
#include "src/minimip/minimip_def.h"
#include "src/minimip/sparse_types.h"

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

// LPInterface base class
class LPGlopInterface : public LPInterface {
 private:
  LinearProgram linear_program_;  // the linear program
  LinearProgram scaled_lp_;       // scaled linear program
  RevisedSimplex solver_;      // direct reference to the revised simplex, not
                               // passing through lp_solver
  GlopParameters parameters_;  // parameters
  LpScalingHelper scaler_;     // scaler auxiliary class

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
  std::int64_t niterations_;  // number of iterations used

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
      int col,            // column to apply strong branching on
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
  absl::Status LoadSparseColumnLP(
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
      ) override;

  // add column to the LP
  absl::Status AddColumn(
      const AbstractSparseVector& col,  // column to be added
      double lower_bound,               // lower bound of new column
      double upper_bound,               // upper bound of new column
      double objective_value,      // objective function value of new column
      const std::string& col_name  // column name
      ) override;

  // add columns to the LP
  absl::Status AddColumns(
      const std::vector<SparseVector>& cols,    // columns to be added
      const std::vector<double>& lower_bounds,  // lower bounds of new columns
      const std::vector<double>& upper_bounds,  // upper bounds of new columns
      const std::vector<double>&
          objective_values,  // objective function values of new columns
      const std::vector<std::string>& col_names  // column names
      ) override;

  // deletes all columns in the given range from LP
  absl::Status DeleteColumns(int first_col,  // first column to be deleted
                             int last_col    // last column to be deleted
                             ) override;

  // add row to the LP
  absl::Status AddRow(const AbstractSparseVector& row,  // row to be added
                      double left_hand_side,       // left hand side of new row
                      double right_hand_side,      // right hand side of new row
                      const std::string& row_name  // row name
                      ) override;

  // adds rows to the LP
  absl::Status AddRows(
      const std::vector<SparseVector>& rows,  // number of rows to be added
      const std::vector<double>&
          left_hand_sides,  // left hand sides of new rows
      const std::vector<double>&
          right_hand_sides,                // right hand sides of new rows
      std::vector<std::string>& row_names  // row names
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

  // change lower bound and upper bound of column
  absl::Status SetColumnBounds(int col, double lower_bound,
                               double upper_bound) override;

  // change left- and right-hand side of row
  absl::Status SetRowSides(int row, double left_hand_side,
                           double right_hand_side) override;

  // changes the objective sense
  absl::Status SetObjectiveSense(
      LPObjectiveSense obj_sense  // new objective sense
      ) override;

  // changes objective value of column in the LP
  absl::Status SetObjectiveCoefficient(int col,
                                       double objective_coefficient) override;

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
      int col,             // column to apply strong branching on
      double primal_sol,   // current primal solution value of column
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
  double GetObjectiveValue() override;

  // gets primal and dual solution vectors for feasible LPs
  //
  // Before calling this function, the caller must ensure that the LP has been
  // solved to optimality, i.e., that minimip::LPInterface.IsOptimal() returns
  // true.

  // gets primal solution vector
  absl::StatusOr<std::vector<double>> GetPrimalSolution() const override;

  // gets row activity vector
  absl::StatusOr<std::vector<double>> GetRowActivity() const override;

  // gets dual solution vector
  absl::StatusOr<std::vector<double>> GetDualSolution() const override;

  // gets reduced cost vector
  absl::StatusOr<std::vector<double>> GetReducedCost() const override;

  // gets primal ray for unbounded LPs
  absl::StatusOr<std::vector<double>> GetPrimalRay() const override;

  // gets dual Farkas proof for infeasibility
  absl::StatusOr<std::vector<double>> GetDualFarkasMultiplier() const override;

  // gets the number of LP iterations of the last solve call
  int GetIterations() const override;

  // ==========================================================================
  // Getters and setters of the basis.
  // ==========================================================================

  // gets current basis status for columns and rows

  // gets current basis status for columns and rows
  absl::StatusOr<std::vector<LPBasisStatus>> GetColumnBasisStatus()
      const override;
  absl::StatusOr<std::vector<LPBasisStatus>> GetRowBasisStatus() const override;

  // sets current basis status for columns and rows
  absl::Status SetBasisStatus(
      const std::vector<LPBasisStatus>& column_basis_status,
      const std::vector<LPBasisStatus>& row_basis_status) override;

  // returns the indices of the basic columns and rows; basic column n gives
  // value n, basic row m gives value -1-m
  std::vector<int> GetBasisIndices() const override;

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
  absl::StatusOr<SparseVector> GetSparseRowOfBInverted(
      int row_number) const override;

  // get column of inverse basis matrix B^-1
  //
  // NOTE: The LP interface defines slack variables to have coefficient +1. This
  // means that if, internally, the LP solver
  //       uses a -1 coefficient, then rows associated with slacks variables
  //       whose coefficient is -1, should be negated
  //
  // column number of B^-1; this is NOT the number of the
  // column in the LP; you have to call
  // minimip::LPInterface.GetBasisIndices() to get the
  // array which links the B^-1 column numbers to the row
  // and column numbers of the LP! c must be between 0 and
  // num_rows-1, since the basis has the size num_rows *
  // num_rows
  absl::StatusOr<SparseVector> GetSparseColumnOfBInverted(
      int col_number) const override;

  // get row of inverse basis matrix times constraint matrix B^-1 * A
  //
  // NOTE: The LP interface defines slack variables to have coefficient +1. This
  // means that if, internally, the LP solver
  //       uses a -1 coefficient, then rows associated with slacks variables
  //       whose coefficient is -1, should be negated; see also the explanation
  //       in lpi.h.
  absl::StatusOr<SparseVector> GetSparseRowOfBInvertedTimesA(
      int row_number) const override;

  // get column of inverse basis matrix times constraint matrix B^-1 * A
  //
  // NOTE: The LP interface defines slack variables to have coefficient +1. This
  // means that if, internally, the LP solver
  //       uses a -1 coefficient, then rows associated with slacks variables
  //       whose coefficient is -1, should be negated; see also the explanation
  //       in lpi.h.
  absl::StatusOr<SparseVector> GetSparseColumnOfBInvertedTimesA(
      int col_number) const override;

  // ==========================================================================
  // Getters and setters of the parameters.
  // ==========================================================================

  // gets integer parameter of LP
  absl::StatusOr<int> GetIntegerParameter(LPParameter type  // parameter number
  ) const override;

  // sets integer parameter of LP
  absl::Status SetIntegerParameter(LPParameter type,  // parameter number
                                   int param_val      // parameter value
                                   ) override;

  // gets floating point parameter of LP
  absl::StatusOr<double> GetRealParameter(LPParameter type  // parameter number
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
