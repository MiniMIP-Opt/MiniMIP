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

// Linear Programming Solver interface (aka LP interface, or lpi) in MiniMip.
//
// Any external LP solver that implements this interface can - in principle -
// be used by MiniMip. This interface is strongly typed, i.e., we use
// `minimip::RowIndex` and `minimip::ColIndex` to index rows (constraints)
// and columns (variables), respectively. Wherever appropriate we use dense
// `absl::StrongVector`, or sparse `minimip::SparseCol` and
// `minimip::SparseRow`, all of which are strongly typed.
//
// If a function can fail, we use `absl::Status` and `absl::StatusOr` to
// indicate the error (note, MiniMip is exception free).

#ifndef SRC_LP_INTERFACE_LPI_H_
#define SRC_LP_INTERFACE_LPI_H_

#include <string>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "src/data_structures/mip_data.h"
#include "src/data_structures/strong_sparse_vector.h"
#include "src/lp_interface/lp_types.h"

namespace minimip {

// This represents an element of the basis which may be either a column or
// a row of the original problem.
class ColOrRowIndex {
 public:
  explicit ColOrRowIndex(ColIndex only_col)
      : col_(only_col), row_(kInvalidRow) {}
  explicit ColOrRowIndex(RowIndex only_row)
      : col_(kInvalidCol), row_(only_row) {}

  ColIndex col() const { return col_; }
  RowIndex row() const { return row_; }

  bool operator==(const ColOrRowIndex& o) const {
    return col_ == o.col_ && row_ == o.row_;
  }

 private:
  ColIndex col_;
  RowIndex row_;
};

// LPInterface is used by MiniMIP to communicate with an underlying LP solver.
class LPInterface {
 public:
  virtual ~LPInterface() = default;

  // ==========================================================================
  // LP model setters.
  // ==========================================================================

  // Sets the entire problem: variables, constraints, objective, bounds, sides,
  // and names.
  virtual absl::Status PopulateFromMipData(const MipData& mip_data) = 0;

  // Adds a new column (aka "variable") to the LP problem stored in the
  // underlying LP solver. The newly created variable gets the next available
  // ColIndex, i.e., in the snippet below the newly created variable can be
  // referenced by `next_free_index`.
  //
  //   const ColIndex next_free_index = lp_interface.GetNumberOfColumns();
  //   lp_interface.AddColumn(...);
  //
  // The `column` is simply appended to the coefficient matrix from the right.
  // Use `lower_bound = -Infinity()` and `upper_bound = Infinity()` for imposing
  // no lower and upper bounds on the added variable, respectively. Use
  // `objective_coefficient=0.0` if the variable is not present in the
  // objective. `name` can be empty.
  virtual absl::Status AddColumn(const SparseCol& column, double lower_bound,
                                 double upper_bound,
                                 double objective_coefficient,
                                 const std::string& name) = 0;

  // Like `AddColumn()`, but adds multiple columns at once. This is useful if
  // the underlying LP solver supports batched operations.
  virtual absl::Status AddColumns(
      const StrongSparseMatrix& matrix,
      const absl::StrongVector<ColIndex, double>& lower_bounds,
      const absl::StrongVector<ColIndex, double>& upper_bounds,
      const absl::StrongVector<ColIndex, double>& objective_coefficients,
      const absl::StrongVector<ColIndex, std::string>& names) = 0;

  // Deletes all columns from `first_col` to `last_col` inclusive on both ends.
  // The column indices larger than `last_col` are decremented by `last_col -
  // first_col + 1`. E.g., assume an LP with columns from 0 to 7, then
  // `DeleteColumns(2, 5)` will result in re-indexing column 6 as column 2, and
  // column 7 as column 3.
  virtual absl::Status DeleteColumns(ColIndex first_col, ColIndex last_col) = 0;

  // Adds a new row (aka constraint) to the LP problem stored in the underlying
  // LP solver. The newly created constraints get the next available RowIndex,
  // i.e., in the snippet below the newly created constraint can be referenced
  // by `next_free_index`.
  //
  // const RowIndex next_free_index = lp_interface.GetNumberOfRows();
  // lp_interface.AddRow(...);
  //
  // The `row` is simply appended to the coefficient matrix from the bottom. Use
  // `lower_bound = -Infinity()` and `upper_bound = Infinity()` for imposing no
  // limit on the left and right-hand sides, respectively, of the added
  // constraint. `name` can be empty.
  virtual absl::Status AddRow(const SparseRow& row, double left_hand_side,
                              double right_hand_side,
                              const std::string& name) = 0;

  // Like `AddRow()`, but adds multiple rows at once. This is useful if
  // the underlying LP solver supports batched operations.
  virtual absl::Status AddRows(
      const absl::StrongVector<RowIndex, SparseRow>& rows,
      const absl::StrongVector<RowIndex, double>& left_hand_sides,
      const absl::StrongVector<RowIndex, double>& right_hand_sides,
      const absl::StrongVector<RowIndex, std::string>& names) = 0;

  // Deletes all rows from `first_row` to `last_row` inclusive on both ends.
  // The rowindices larger than `last_row` are decremented by `last_row -
  // first_row + 1`. E.g., assume an LP with rows from 0 to 7, then
  // `DeleteRows(2, 5)` will result in re-indexing row 6 as row 2, and
  // row 7 as row 3.
  virtual absl::Status DeleteRows(RowIndex first_row, RowIndex last_row) = 0;

  // Deletes all rows `r` for which `rows_to_delete[r]` is true. Computes and
  // returns the "row_mapping":
  //   row_mapping[old_index] == kInvalidRow if a row at `old_index` was deleted
  //   row_mapping[old_index] == new_index, otherwise.
  // Note, it is guaranteed that always new_index <= old_index (rows are not
  // reshuffled on deletion, but the matrix gets "compacted").
  virtual absl::StatusOr<absl::StrongVector<RowIndex, RowIndex>> DeleteRowSet(
      const absl::StrongVector<RowIndex, bool>& rows_to_delete) = 0;

  // Clears the whole LP solver (including the loaded model).
  virtual absl::Status Clear() = 0;

  // Clears only the current state (e.g., basis information) of the LP solver,
  // but not the loaded model.
  virtual absl::Status ClearState() = 0;

  // Sets lower and upper bounds for `col`.
  virtual absl::Status SetColumnBounds(ColIndex col, double lower_bound,
                                       double upper_bound) = 0;

  // Sets left and right-hand sides for `row`, assuming lower-or-equal
  // relationship. I.e., left_hand_side <= constraint[row] <= right_hand_side.
  // Use '-Infinity()` and `Infinity()` to impose no left and right-hand side
  // limit, respectively.
  virtual absl::Status SetRowSides(RowIndex row, double left_hand_side,
                                   double right_hand_side) = 0;

  // Sets the objective sense.
  virtual absl::Status SetObjectiveSense(bool is_maximization) = 0;

  // Sets the objective coefficient for `col`.
  virtual absl::Status SetObjectiveCoefficient(
      ColIndex col, double objective_coefficient) = 0;

  // ==========================================================================
  // LP model getters.
  //
  // In all calls taking `col` and/or `row` as arguments, these indices must be
  // valid, i.e.,
  // `col` must be a valid column index, i.e., 0 <= col < GetNumberOfColumns()
  // `row` must be a valid row index, i.e., 0 <= row < GetNumberOfRows()
  // ==========================================================================

  // Returns the number of rows stored in the LP solver.
  virtual RowIndex GetNumberOfRows() const = 0;

  // Returns the number of columns stored in the LP solver.
  virtual ColIndex GetNumberOfColumns() const = 0;

  // Returns the total number of non-zeros in the constraint matrix stored in
  // the LP solver (i.e., across all rows and columns).
  virtual int64_t GetNumberOfNonZeros() const = 0;

  // Returns the objective sense.
  virtual bool IsMaximization() const = 0;

  // Returns the coefficients from the constraint matrix for `col`.
  virtual SparseCol GetSparseColumnCoefficients(ColIndex col) const = 0;

  // Returns the coefficients from the constraint matrix for `row`.
  virtual SparseRow GetSparseRowCoefficients(RowIndex row) const = 0;

  // Returns the objective coefficient for `col`. Returns 0.0 if the column
  // (variable) is not present in the objective.
  virtual double GetObjectiveCoefficient(ColIndex col) const = 0;

  // Returns the lower bound for `col` (aka variable lower bound).
  virtual double GetLowerBound(ColIndex col) const = 0;

  // Returns the upper bound for `col` (aka variable upper bound).
  virtual double GetUpperBound(ColIndex col) const = 0;

  // Returns the left-hand side for `row` (aka constraint lower bound).
  virtual double GetLeftHandSide(RowIndex row) const = 0;

  // Returns the right-hand side for `row` (aka constraint upper bound).
  virtual double GetRightHandSide(RowIndex row) const = 0;

  // Returns a single coefficient from the constraint matrix for `col` and
  // `row`. Depending on the underlying LP solver this may be "relatively
  // inefficient" and should not be used for batch accesses. Instead, iterate
  // over entries of SparseRow (resp. SparseCol) obtained with
  // `GetSparseRowCoefficients()` (resp. `GetSparseColumnCoefficients()`).
  virtual double GetMatrixCoefficient(ColIndex col, RowIndex row) const = 0;

  // ==========================================================================
  // Solving methods.
  // ==========================================================================

  // Calls primal simplex to solve the currently loaded LP.
  virtual absl::Status SolveLPWithPrimalSimplex() = 0;

  // Calls dual simplex to solve the currently loaded LP.
  virtual absl::Status SolveLPWithDualSimplex() = 0;

  // Informs the LP solver to enter "strong branching" mode. This may entail
  // setting loose precision requirements, stringent iteration / time limits,
  // etc. Call this before performing any strong branching, otherwise strong
  // branching will be inefficient.
  virtual absl::Status StartStrongBranching() = 0;

  // Informs the LP solver to quit "strong branching" mode. Call this once
  // strong branching is over, and the LP solver is used to solve "regular"
  // nodes' LP relaxations, otherwise LP relaxation of the nodes in the MIP
  // search tree may be inaccurate or remain not solved at all.
  virtual absl::Status EndStrongBranching() = 0;

  // Result of strong branching on a single variable.
  struct StrongBranchResult {
    // The objective value of the LP relaxation after branching down.
    double dual_bound_down_branch;

    // The objective value of the LP relaxation after branching up.
    double dual_bound_up_branch;

    // Whether `dual_bound_down_branch` is precise. If not, it can only be used
    // as an estimate (which may still be fine, depending on the use case).
    bool down_valid;

    // Whether `dual_bound_up_branch` is precise. If not, it can only be used
    // as an estimate (which may still be fine, depending on the use case).
    bool up_valid;

    // The total number of strong branching iterations. Must be non-negative.
    int64_t iterations;
  };

  // Performs strong branching iterations on a single strong branching
  // candidate.
  virtual absl::StatusOr<StrongBranchResult> SolveDownAndUpStrongBranch(
      ColIndex col, double primal_value, int iteration_limit) = 0;

  // ==========================================================================
  // Solution information getters.
  // ==========================================================================

  // Whether a solve method was called after the last modification of the LP.
  // This does not mean the solve was successful, only that it was called.
  virtual bool IsSolved() const = 0;

  // Whether current LP solution is stable.
  //
  // This function should return true if the solution is reliable, i.e.,
  // feasible and optimal (or proven infeasible/unbounded) with respect to the
  // original problem. The optimality status might be with respect to a scaled
  // version of the problem, but the solution might not be feasible to the
  // unscaled original problem; in this case, IsStable() should return false.
  //
  // TODO(lpawel): Explain more precisely what stability means.
  virtual bool IsStable() const = 0;

  // LP solve statuses. After solving a feasible bounded model to completion,
  // these should represent that actual model status. If the problem is primal
  // infeasible and solved with primal simplex, the solver isn't required to
  // prove dual unboundedness/infeasibility, and similarly for dual infeasible
  // models when using dual simplex.
  virtual bool IsOptimal() const = 0;
  virtual bool IsPrimalFeasible() const = 0;
  virtual bool IsPrimalInfeasible() const = 0;
  virtual bool IsPrimalUnbounded() const = 0;
  virtual bool IsDualFeasible() const = 0;
  virtual bool IsDualInfeasible() const = 0;
  virtual bool IsDualUnbounded() const = 0;

  // Returns true if the LP is proven to have a primal unbounded ray (but not
  // necessary a primal feasible point). This does not mean that the solver
  // knows and can return the primal ray.
  virtual bool ExistsPrimalRay() const = 0;

  // Returns true if LP is proven to have a primal unbounded ray (but not
  // necessary a primal feasible point). The LP solver knows and can return the
  // primal ray.
  virtual bool HasPrimalRay() const = 0;

  // Returns true if LP is proven to have a dual unbounded ray (but not
  // necessary a dual feasible point). This does not necessarily mean that the
  // solver knows and can return the dual ray.
  virtual bool ExistsDualRay() const = 0;

  // Returns true if LP is proven to have a dual unbounded ray (but not
  // necessary a dual feasible point), and the solver knows and can return the
  // dual ray.
  virtual bool HasDualRay() const = 0;

  // Returns true if the objective limit was reached.
  virtual bool ObjectiveLimitIsExceeded() const = 0;

  // Returns true if the time limit was reached.
  virtual bool TimeLimitIsExceeded() const = 0;

  // Returns true if the iteration limit was reached.
  virtual bool IterationLimitIsExceeded() const = 0;

  // Gets the number of LP iterations of the last solve call.
  virtual int64_t GetNumIterations() const = 0;

  // Returns the objective value corresponding to `GetPrimalValues()`.
  // `IsOptimal()` must be true when calling this.
  virtual double GetObjectiveValue() const = 0;

  // Returns the primal values for all columns (i.e., the solution).
  // `IsOptimal()` must be true when calling this.
  virtual absl::StatusOr<absl::StrongVector<ColIndex, double>> GetPrimalValues()
      const = 0;

  // Returns dual values for all rows. `IsOptimal() must be true when calling
  // this.
  virtual absl::StatusOr<absl::StrongVector<RowIndex, double>> GetDualValues()
      const = 0;

  // Returns reduced costs for all columns. `IsOptimal() must be true when
  // calling this.
  virtual absl::StatusOr<absl::StrongVector<ColIndex, double>> GetReducedCosts()
      const = 0;

  // Returns constraint activities corresponding to `GetPrimalValues()`.
  // `IsOptimal()` must be true when calling this.
  virtual absl::StatusOr<absl::StrongVector<RowIndex, double>>
  GetRowActivities() const = 0;

  // Returns primal and dual rays. `HasPrimalRay()` and `HasDualRay() must be
  // true, respectively, before calling these.
  virtual absl::StatusOr<absl::StrongVector<ColIndex, double>> GetPrimalRay()
      const = 0;
  virtual absl::StatusOr<absl::StrongVector<RowIndex, double>> GetDualRay()
      const = 0;

  // ==========================================================================
  // Getters and setters of the basis.
  // ==========================================================================

  // Returns the basis status for all variables.
  virtual absl::StatusOr<absl::StrongVector<ColIndex, LPBasisStatus>>
  GetBasisStatusForColumns() const = 0;

  // Returns the basis status for all constraints.
  virtual absl::StatusOr<absl::StrongVector<RowIndex, LPBasisStatus>>
  GetBasisStatusForRows() const = 0;

  virtual absl::Status SetBasisStatusForColumnsAndRows(
      const absl::StrongVector<ColIndex, LPBasisStatus>& col_statuses,
      const absl::StrongVector<RowIndex, LPBasisStatus>& row_statuses) = 0;

  // Returns the indices of the basic columns and rows. The size of the returned
  // vector is always `GetNumberOfRows()`.
  virtual std::vector<ColOrRowIndex> GetColumnsAndRowsInBasis() const = 0;

  // ==========================================================================
  // Getters of vectors in the inverted basis matrix.
  //
  // Note 1, all these functions take `col_in_basis` or `row_in_basis` as
  // arguments. These are *not* columns and rows of the coefficient matrix
  // (i.e., A). These are columns and rows of the inverted basis. To get a
  // column or row index in A corresponding to `*_in_basis` index use
  // `GetColsAndRowsInBasis()[*_in_basis.value()`.]
  //
  // Note 2, the LP interface assumes slack variables were added with +1.0
  // coefficient. If the underlying LP solver adds slacks with -1.0 coefficients
  // then rows associated with the slack variables must be negated before
  // returning!
  //
  // ==========================================================================

  // Gets a row of the inverse basis matrix, i.e., B^-1[row_in_basis][*]).
  virtual absl::StatusOr<SparseRow> GetSparseRowOfBInverted(
      RowIndex row_in_basis) const = 0;

  // Gets a column of the inverse basis matrix, i.e., B^-1[*][col_in_basis].
  virtual absl::StatusOr<SparseCol> GetSparseColumnOfBInverted(
      ColIndex col_in_basis) const = 0;

  // Gets a row of the inverse basis matrix multiplied by the constraint matrix,
  // i.e., (B^-1 * A)[row_in_basis][*].
  virtual absl::StatusOr<SparseRow> GetSparseRowOfBInvertedTimesA(
      RowIndex row_in_basis) const = 0;

  // Gets a column of the inverse basis matrix multiplied by constraint matrix,
  // i.e., (B^-1 * A)[*][col_in_basis].
  virtual absl::StatusOr<SparseCol> GetSparseColumnOfBInvertedTimesA(
      ColIndex col_in_basis) const = 0;

  // ==========================================================================
  // Getters and setters of the parameters.
  // ==========================================================================

  // Gets an integer parameter of LP.
  // TODO(lpawel): Setting boolean parameters is currently done via "integer"
  // parameters. Fix this.
  virtual absl::StatusOr<int> GetIntegerParameter(LPParameter type) const = 0;

  // Sets an integer parameter of LP.
  virtual absl::Status SetIntegerParameter(LPParameter type, int param_val) = 0;

  // Gets a floating point parameter of LP.
  virtual absl::StatusOr<double> GetRealParameter(LPParameter type) const = 0;

  // Sets a floating point parameter of LP.
  virtual absl::Status SetRealParameter(LPParameter type, double param_val) = 0;

  // ==========================================================================
  // Numerical methods.
  // ==========================================================================

  // The "canonical" value interpreted as positive infinity in the LP solver.
  virtual double Infinity() const = 0;

  // Checks if LP solver treats `value` as positive infinity.
  // Note 1: The LP solver may interpret more values than just `Infinity()` as
  //         infinity. Thus, to detect whether a value is treated as infinity
  //         use `IsInfinity(value)` and *not* `value == Infinity()`.
  // Note 2: To check whether the LP solver interprets a value as negative
  //         infinity use `IsInfinity(-value)`.
  virtual bool IsInfinity(double value) const = 0;

  // ==========================================================================
  // File interface methods.
  // ==========================================================================

  // Reads an LP from a file in the format supported by the underlying LP solver
  // (may differ between implementations). After this call, the LP model should
  // *essentially* the same as before calling `WriteLPToFile`. By "essentially
  // the same", we mean that the models should have the same solution space and
  // objective values, but e.g. the order of variables may be different, and
  // double-sided constraints may be split to two single-sided constraints. This
  // should be documented for each implementation.
  virtual absl::Status ReadLPFromFile(const std::string& file_path) = 0;

  // Writes an LP to a file in the format supported by the underlying LP solver
  // (may differ between implementations). Returns the path to the written file.
  // This may differ from the provided path in the file extension only.
  virtual absl::StatusOr<std::string> WriteLPToFile(
      const std::string& file_path) const = 0;
};

}  // namespace minimip
#endif  // SRC_LP_INTERFACE_LPI_H_
