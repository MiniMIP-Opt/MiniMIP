// Copyright 2024 the MiniMIP Project
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

#ifndef MINIMIP_LP_INTERFACE_LPI_SOPLEX_H_
#define MINIMIP_LP_INTERFACE_LPI_SOPLEX_H_

// include SoPlex solver
#include <algorithm>
#include <string>
#include <vector>

#include "minimip/lp_interface/lpi.h"
#include "minimip/parameters.pb.h"
#include "ortools/util/file_util.h"
#include "soplex.h"

namespace minimip {

class LpSoplexInterface : public LpInterface {
 public:
  LpSoplexInterface();

  ~LpSoplexInterface() override;

  LpParameters::SolverType GetSolverType() const override {
    return LpParameters::LP_SOPLEX;
  }
  // ==========================================================================
  // LP model setters.
  // ==========================================================================

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

  absl::Status SolveLpWithPrimalSimplex() final;
  absl::Status SolveLpWithDualSimplex() final;
  absl::Status StartStrongBranching() final;
  absl::Status EndStrongBranching() final;

  // Strongbranching is applied to the given column, with the corresponding
  // current primal solution value. The double references are used to store the
  // dual bound after branching up and down. Additionally, the validity of both
  // bounds is stored, if one bound is not valid it can be used as an estimate.
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

  double GetObjectiveValue() const final;

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

  absl::StatusOr<absl::StrongVector<ColIndex, LpBasisStatus>>
  GetBasisStatusForColumns() const final;
  absl::StatusOr<absl::StrongVector<RowIndex, LpBasisStatus>>
  GetBasisStatusForRows() const final;

  absl::Status SetBasisStatusForColumnsAndRows(
      const absl::StrongVector<ColIndex, LpBasisStatus>& column_basis_statuses,
      const absl::StrongVector<RowIndex, LpBasisStatus>& row_basis_statuses)
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

  LpParameters GetLpParameters() const final;
  absl::Status SetLpParameters(const LpParameters& lp_parameters) final;

  // ==========================================================================
  // Numerical methods.
  // ==========================================================================

  double Infinity() const final;
  bool IsInfinity(double value) const final;

  // ==========================================================================
  // File interface methods.
  // ==========================================================================

  absl::Status ReadLpFromFile(const std::string& file_path) final;

  // Note that Soplex splits double-sided constraints in the stored file. I.e.
  // l_j <= sum(a_ijx_i) <= u_j
  // becomes
  // l_j <= sum(a_ijx_i)
  //        sum(a_ijx_i) <= u_j
  absl::StatusOr<std::string> WriteLpToFile(
      const std::string& file_path) const final;

  // ==========================================================================
  // Private interface methods.
  // ==========================================================================

 private:
  double feasibility_tolerance() const {
    DCHECK(spx_ != nullptr);
    return spx_->realParam(soplex::SoPlex::FEASTOL);
  }

  double objective_limit() const {
    DCHECK(spx_ != nullptr);
    return (spx_->intParam(soplex::SoPlex::OBJSENSE) ==
            soplex::SoPlex::OBJSENSE_MINIMIZE)
               ? spx_->realParam(soplex::SoPlex::OBJLIMIT_UPPER)
               : spx_->realParam(soplex::SoPlex::OBJLIMIT_LOWER);
  }

  // Marks the current LP to be unsolved.
  void InvalidateSolution();

  bool PreStrongBranchingBasisFreed() const;

  void FreePreStrongBranchingBasis();

  void RestorePreStrongbranchingBasis();

  void SavePreStrongbranchingBasis();

  bool CheckConsistentBounds() const;

  bool CheckConsistentSides() const;

  soplex::DataArray<soplex::SPxSolver::VarStatus>& RowsBasisStatus();

  soplex::DataArray<soplex::SPxSolver::VarStatus>& ColumnsBasisStatus();

  soplex::SPxSolver::Status LpSolve(bool print_warning);

  // solves LP -- used for both, primal and dual simplex, because SoPlex doesn't
  // distinct the two cases
  absl::Status SoPlexSolve();

  // ==========================================================================
  // Member variables of the SoplexLpInterface.
  // ==========================================================================

  // SoPlex solver implementing primal and dual simplex algorithm.
  std::unique_ptr<soplex::SoPlex> spx_;

  // Whether to start solving from scratch (i.e., don't leverage
  // incrementalism).
  bool solve_from_scratch_;

  // Whether the LP that is currently loaded inside `spx_` has already been
  // solved.
  bool is_solved_;

  // The column basis status used for strong branching.
  soplex::DataArray<soplex::SPxSolver::VarStatus> col_basis_status_;

  // The row basis status used for strong branching.
  soplex::DataArray<soplex::SPxSolver::VarStatus> row_basis_status_;
};

}  // namespace minimip

#endif  // MINIMIP_LP_INTERFACE_LPI_SOPLEX_H_
