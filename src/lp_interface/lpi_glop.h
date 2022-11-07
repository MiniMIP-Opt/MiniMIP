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

#include <cassert>
#include <cinttypes>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "ortools/base/logging.h"
#include "ortools/base/strong_vector.h"
#include "ortools/lp_data/lp_data_utils.h"
#include "ortools/lp_data/lp_print_utils.h"
#include "ortools/lp_data/proto_utils.h"
#include "ortools/util/file_util.h"
#include "ortools/util/stats.h"
#include "ortools/util/strong_integers.h"
#include "src/lp_interface/lp_types.h"
#include "src/lp_interface/lpi.h"

namespace minimip {

class LPGlopInterface : public LPInterface {
 public:
  LPGlopInterface();

  // ==========================================================================
  // LP model setters.
  // ==========================================================================

  // Sets the entire problem: variables, constraints, objective, bounds, sides,
  // and names.
  absl::Status PopulateFromMipData(const MipData& mip_data) final;

  absl::Status AddColumn(const SparseCol& col, double lower_bound,
                         double upper_bound, double objective_coefficient,
                         const std::string& name) final;

  absl::Status AddColumns(
      const StrongSparseMatrix& matrix,
      const absl::StrongVector<ColIndex, double>& lower_bounds,
      const absl::StrongVector<ColIndex, double>& upper_bounds,
      const absl::StrongVector<ColIndex, double>& objective_coefficients,
      const absl::StrongVector<ColIndex, std::string>& names) final;

  absl::Status DeleteColumns(ColIndex first_col, ColIndex last_col) final;

  absl::Status AddRow(const SparseRow& row, double left_hand_side,
                      double right_hand_side, const std::string& name) final;

  absl::Status AddRows(
      const absl::StrongVector<RowIndex, SparseRow>& rows,
      const absl::StrongVector<RowIndex, double>& left_hand_sides,
      const absl::StrongVector<RowIndex, double>& right_hand_sides,
      const absl::StrongVector<RowIndex, std::string>& names) final;

  absl::Status DeleteRows(RowIndex first_row, RowIndex last_row) final;

  absl::StatusOr<absl::StrongVector<RowIndex, RowIndex>> DeleteRowSet(
      const absl::StrongVector<RowIndex, bool>& rows_to_delete) final;

  absl::Status Clear() final;

  absl::Status ClearState() final;

  absl::Status SetColumnBounds(ColIndex col, double lower_bound,
                               double upper_bound) final;

  absl::Status SetRowSides(RowIndex row, double left_hand_side,
                           double right_hand_side) final;

  absl::Status SetObjectiveSense(bool is_maximization) final;

  absl::Status SetObjectiveCoefficient(ColIndex col,
                                       double objective_coefficient) final;

  // ==========================================================================
  // LP model getters.
  // ==========================================================================

  RowIndex GetNumberOfRows() const final;
  ColIndex GetNumberOfColumns() const final;
  int64_t GetNumberOfNonZeros() const final;

  bool IsMaximization() const final;
  SparseCol GetSparseColumnCoefficients(ColIndex col) const final;
  SparseRow GetSparseRowCoefficients(RowIndex row) const final;

  double GetObjectiveCoefficient(ColIndex col) const final;
  double GetLowerBound(ColIndex col) const final;
  double GetUpperBound(ColIndex col) const final;
  double GetLeftHandSide(RowIndex row) const final;
  double GetRightHandSide(RowIndex row) const final;
  double GetMatrixCoefficient(ColIndex col, RowIndex row) const final;

  // ==========================================================================
  // Solving methods.
  // ==========================================================================

  absl::Status SolveLPWithPrimalSimplex() final;
  absl::Status SolveLPWithDualSimplex() final;
  absl::Status StartStrongBranching() final;
  absl::Status EndStrongBranching() final;

  absl::StatusOr<StrongBranchResult> SolveDownAndUpStrongBranch(
      ColIndex col, double primal_value, int iteration_limit) final;

  // ==========================================================================
  // Solution information getters.
  // ==========================================================================

  bool IsSolved() const final;
  bool IsStable() const final;
  bool IsOptimal() const final;
  bool IsPrimalFeasible() const final;
  bool IsPrimalInfeasible() const final;
  bool IsPrimalUnbounded() const final;
  bool IsDualFeasible() const final;
  bool IsDualInfeasible() const final;
  bool IsDualUnbounded() const final;

  bool ExistsPrimalRay() const final;
  bool HasPrimalRay() const final;
  bool ExistsDualRay() const final;
  bool HasDualRay() const final;

  bool ObjectiveLimitIsExceeded() const final;
  bool TimeLimitIsExceeded() const final;
  bool IterationLimitIsExceeded() const final;
  int64_t GetNumIterations() const final;

  double GetObjectiveValue() final;

  absl::StatusOr<absl::StrongVector<ColIndex, double>> GetPrimalValues()
      const final;
  absl::StatusOr<absl::StrongVector<RowIndex, double>> GetDualValues()
      const final;

  absl::StatusOr<absl::StrongVector<ColIndex, double>> GetReducedCosts()
      const final;
  absl::StatusOr<absl::StrongVector<RowIndex, double>> GetRowActivities()
      const final;

  absl::StatusOr<absl::StrongVector<ColIndex, double>> GetPrimalRay()
      const final;
  absl::StatusOr<absl::StrongVector<RowIndex, double>> GetDualRay() const final;

  // ==========================================================================
  // Getters and setters of the basis.
  // ==========================================================================

  absl::StatusOr<absl::StrongVector<ColIndex, LPBasisStatus>>
  GetBasisStatusForColumns() const final;
  absl::StatusOr<absl::StrongVector<RowIndex, LPBasisStatus>>
  GetBasisStatusForRows() const final;

  absl::Status SetBasisStatusForColumnsAndRows(
      const absl::StrongVector<ColIndex, LPBasisStatus>& column_basis_statuses,
      const absl::StrongVector<RowIndex, LPBasisStatus>& row_basis_statuses)
      final;

  std::vector<ColOrRowIndex> GetColumnsAndRowsInBasis() const final;

  // ==========================================================================
  // Getters of vectors in the inverted basis matrix.
  // ==========================================================================
  absl::StatusOr<SparseRow> GetSparseRowOfBInverted(
      RowIndex row_in_basis) const final;

  absl::StatusOr<SparseCol> GetSparseColumnOfBInverted(
      ColIndex col_in_basis) const final;

  absl::StatusOr<SparseRow> GetSparseRowOfBInvertedTimesA(
      RowIndex row_in_basis) const final;

  absl::StatusOr<SparseCol> GetSparseColumnOfBInvertedTimesA(
      ColIndex col_in_basis) const final;

  // ==========================================================================
  // Getters and setters of the parameters.
  // ==========================================================================

  absl::StatusOr<int> GetIntegerParameter(LPParameter type) const final;

  absl::Status SetIntegerParameter(LPParameter type, int param_value) final;

  absl::StatusOr<double> GetRealParameter(LPParameter type) const final;

  absl::Status SetRealParameter(LPParameter type, double param_value) final;

  // ==========================================================================
  // Numerical methods.
  // ==========================================================================

  double Infinity() const final;
  bool IsInfinity(double value) const final;

  // ==========================================================================
  // File interface methods.
  // ==========================================================================

  absl::Status ReadLPFromFile(const std::string& file_path) final;

  absl::StatusOr<std::string> WriteLPToFile(
      const std::string& file_path) const final;

 private:
  // Helper function to delete rows from LP and update the current basis.
  void DeleteRowsAndUpdateCurrentBasis(
      const operations_research::glop::DenseBooleanColumn& rows_to_delete);

  // Glop's `Solve()` happens inside this function. With this we can call
  // `SolveInternal()` from within `SolveInternal()` (in case of numerical
  // issues and resolving from scratch or unscaled).
  absl::Status SolveInternal(bool recursive,
                             operations_research::TimeLimit* time_limit);

  // Contains the current lp to be solved in the original form.
  operations_research::glop::LinearProgram lp_;

  // Contains the scaled lp.
  operations_research::glop::LinearProgram scaled_lp_;

  // Revised primal and dual simplex solver from Glop. We use `RevisedSimplex`
  // directly (and not `glop::LpSolver`) to bypass the extra layers for speed.
  operations_research::glop::RevisedSimplex solver_;

  // Glop current parameters.
  operations_research::glop::GlopParameters parameters_;

  // Scaler to compute `scaled_lp_` from `lp_`.
  operations_research::glop::LpScalingHelper scaler_;

  // Indicator whether we need to recompute `scaled_lp_` from `lp_` on next
  // solve.
  bool lp_modified_since_last_solve_;

  // Indicator whether the last solve hit a time limit.
  bool lp_time_limit_was_reached_;

  // Number of simplex iterations executed since the last solve.
  int64_t num_iterations_of_last_solve_;

  // Store the values of some parameters in order to be able to return them
  bool lp_info_;       // whether additional output is turned on
  LPPricing pricing_;  // MiniMIP pricing setting
  bool from_scratch_;  // store whether basis is ignored for next solving call
  int num_threads_;    // number of threads used to solve the LP (0 = automatic)
  int timing_;         // type of timer (1 - cpu, 2 - wallclock, 0 - off)

  // Temporary vectors allocated here for speed. This gain is non-negligible
  // because in many situations, only a few entries of these vectors are
  // inspected (hypersparsity) and allocating them is in O(num_rows) or
  // O(num_cols) instead of O(num_non_zeros) to read/clear them.
  std::unique_ptr<operations_research::glop::ScatteredRow> tmp_row_;
  std::unique_ptr<operations_research::glop::ScatteredColumn> tmp_column_;

};  // class LPGlopInterface

}  // namespace minimip

#endif  // SRC_LP_INTERFACE_LPI_GLOP_H_
