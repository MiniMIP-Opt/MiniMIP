// Copyright 2023 the MiniMIP Project
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

#include "src/cutting_interface/aggregating_separator.h"

#include <limits>

#include "cuts_selector.h"
#include "cuts_separator.h"
#include "src/solver.h"

namespace minimip {

namespace {

// Compute the fractionality of a real value.
double Fractionality(const Solver& solver, double v) {
  return solver.IsIntegerWithinTolerance(v) ? 0.0 : v - std::floor(v);
}

// Check if a slack variable must by integer valued.
bool HasIntegralSlackVariable(const Solver& solver, RowIndex row) {
  return solver.mip_data().is_integral_constraint()[row];
}

// Round the coefficient of a single non-negative integer variable using Mixed
// Integer Rounding.
double MIRRoundInteger(const Solver& solver, double coefficient, double f0) {
  const double fj = Fractionality(solver, coefficient);
  return fj <= f0
             ? solver.FloorWithTolerance(coefficient)
             : solver.FloorWithTolerance(coefficient) + (fj - f0) / (1 - f0);
}

// Round the coefficient of a single non-negative continuous variable using
// Mixed Integer Rounding.
double MIRRoundContinuous(double coefficient, double f0) {
  return coefficient >= 0.0 ? 0.0 : coefficient / (1 - f0);
}

// Round the coefficient of a single non-negative integer variable using strong
// Chvatal-Gomory rounding.
double CGRoundInteger(const Solver& solver, double coefficient, double k,
                      double f0) {
  const double fj = Fractionality(solver, coefficient);
  const double p = fj <= f0 ? 0 : std::ceil(k * (fj - f0) / (1 - f0));
  if (fj <= f0) DCHECK_EQ(p, 0.0);
  if (fj > f0) {
    DCHECK_LT(f0 + 1 / k * (p - 1) * (1 - f0), fj);
    DCHECK_LE(fj, f0 + 1 / k * p * (1 - f0));
  }
  return solver.FloorWithTolerance(coefficient) + p / (k + 1);
}

// Creates a row aggregation by summing rows with the given weights. See the
// top-level comment in the header file for details.
AggregatedRow AggregateByWeight(
    const Solver& solver, const SparseCol& weights,
    const absl::StrongVector<RowIndex, bool>& use_right_hand_side) {
  AggregatedRow aggregated_row;
  for (auto [row, weight] : weights.entries()) {
    const double side_value = use_right_hand_side[row]
                                  ? solver.lpi()->GetRightHandSide(row)
                                  : solver.lpi()->GetLeftHandSide(row);
    const double slack_sign = use_right_hand_side[row] ? 1 : -1;

    aggregated_row.variable_coefficients.AddMultipleOfVector(
        weight, solver.lpi()->GetSparseRowCoefficients(row));
    aggregated_row.slack_coefficients.AddEntry(row, slack_sign * weight);
    aggregated_row.slack_signs.AddEntry(row, slack_sign);
    aggregated_row.side_values.AddEntry(row, side_value);
    aggregated_row.right_hand_side += weight * side_value;
  }
  aggregated_row.slack_coefficients.CleanUpIfNeeded();
  aggregated_row.slack_signs.CleanUpIfNeeded();
  aggregated_row.side_values.CleanUpIfNeeded();
  return aggregated_row;
}

// Transforms the row to a space where all variables are non-negative. Returns
// false if the transformation failed. Note that after calling this function,
// the coefficients in variable_coefficients refers to auxiliary variables, not
// the original ones. Hence, it is important to call
// `TransformBackToOriginalVariables` before the resulting cut is used or slack
// variables are substituted.
[[nodiscard]] bool TransformToNonNegativeVariables(
    const Solver& solver, AggregatedRow& aggregated_row) {
  bool success = true;
  aggregated_row.variable_coefficients.Transform([&solver, &aggregated_row,
                                                  &success](
                                                     ColIndex column,
                                                     double coefficient) {
    const double lb = solver.lpi()->GetLowerBound(column);
    const double ub = solver.lpi()->GetUpperBound(column);

    // Variable is already positive, no transformation needed.
    if (lb >= 0.0) return coefficient;

    if (!solver.lpi()->IsInfinity(-lb)) {
      // We have x[j] >= lb <-> x[j] - lb >= 0 and can define the new
      // variable
      //
      // x[j]' = x[j] - lb <-> x[j] = x[j]' + lb
      //
      // where x[j]' >= 0. Inserted into our contraint, we get
      //
      // c[j] * x[j] = c[j] * (x[j]' + lb) <= r
      // <->
      // c[j] * x[j]' <= r - c[j] * lb
      aggregated_row.right_hand_side -= coefficient * lb;
      return coefficient;
    }

    if (!solver.lpi()->IsInfinity(ub)) {
      // We have x[j] <= ub <-> ub - x[j] >= 0 and can define the new
      // variable
      //
      // x[j]' = ub - x[j] <-> x[j] = ub - x[j]'
      //
      // where x[j]' >= 0. Inserted into our contraint, we get
      //
      // c[j] * x[j] = c[j] * (ub - x[j]') <= r
      // <->
      // -c[j] * x[j]' <= r - c[j] * ub
      aggregated_row.right_hand_side -= coefficient * ub;
      return -coefficient;
    }

    VLOG(2) << "Variable " << column
            << " is unbounded and cannot be included in rounded constraint.";
    success = false;
    return coefficient;
  });
  return success;
}

// The inverse of `TransformToNonNegativeVariables`.
[[nodiscard]] bool TransformBackToOriginalVariables(
    const Solver& solver, AggregatedRow& aggregated_row) {
  bool success = true;
  aggregated_row.variable_coefficients.Transform(
      [&solver, &aggregated_row, &success](ColIndex column,
                                           double coefficient) {
        const double lb = solver.lpi()->GetLowerBound(column);
        const double ub = solver.lpi()->GetUpperBound(column);

        // Variable is already positive, no transformation needed.
        if (lb >= 0.0) return coefficient;

        if (!solver.lpi()->IsInfinity(-lb)) {
          // We used x[j]' = x[j] - lb, so transforming back we get
          //
          // c[j]' * x[j]' = c[j]' * (x[j] - lb) <= r'
          // <->
          // c[j]' * x[j] <= r' + c[j]' * lb
          aggregated_row.right_hand_side += coefficient * lb;
          return coefficient;
        }

        if (!solver.lpi()->IsInfinity(ub)) {
          // We used x[j]' = ub - x[j], so transforming back we get
          //
          // c[j]' * x[j]' = c[j]' * (ub - x[j]) <= r'
          // <->
          // -c[j]' * x[j] <= r' - c[j]' * ub
          aggregated_row.right_hand_side -= coefficient * ub;
          return -coefficient;
        }

        LOG(DFATAL) << "Variable " << column
                    << " is unbounded, which should have been detected in the "
                       "original transform.";
        success = false;
        return coefficient;
      });
  return success;
}

// Substitutes all slack variables using their definition. See the top-level
// comment in the header file for details.
void SubstituteSlackVariables(const Solver& solver,
                              AggregatedRow& aggregated_row) {
  for (const auto [row, slack_coefficient] :
       aggregated_row.slack_coefficients.entries()) {
    const double slack_sign = aggregated_row.slack_signs[row];
    const double side_value = aggregated_row.side_values[row];
    aggregated_row.variable_coefficients -=
        slack_sign * slack_coefficient *
        solver.lpi()->GetSparseRowCoefficients(row);
    aggregated_row.right_hand_side -=
        slack_sign * slack_coefficient * side_value;
  }
  aggregated_row.slack_coefficients.Clear();
  aggregated_row.slack_signs.Clear();
  aggregated_row.side_values.Clear();
}

// Compute which side to use for each row in the LP based on the current basis
// information. `true` means that the right hand side should be chosen, while
// `false` means the left hand side.
absl::StatusOr<absl::StrongVector<RowIndex, bool>> ChooseActiveSidesByBasis(
    const Solver& solver) {
  ASSIGN_OR_RETURN((const absl::StrongVector<RowIndex, double> activities),
                   solver.lpi()->GetRowActivities());

  absl::StrongVector<RowIndex, bool> use_right_hand_side(activities.size());
  for (RowIndex row(0); row < activities.size(); ++row) {
    const double rhs = solver.lpi()->GetRightHandSide(row);
    const double lhs = solver.lpi()->GetLeftHandSide(row);

    // If the constraint is one-sided, we don't have a choice.
    if (solver.lpi()->IsInfinity(rhs)) {
      use_right_hand_side[row] = false;
      continue;
    }
    if (solver.lpi()->IsInfinity(-lhs)) {
      use_right_hand_side[row] = true;
      continue;
    }

    // For equality constraints, the side doesn't matter.
    if (lhs == rhs) {
      use_right_hand_side[row] = true;
      continue;
    }
    DCHECK_LT(lhs, rhs);

    // Use the side that is closest to being violated. Note that this is a
    // heuristic choice.
    use_right_hand_side[row] = (activities[row] - lhs) / (rhs - lhs) >= 0.5;
  }
  return use_right_hand_side;
}

// Check if a cut removes the current lp optimum.
bool RemovesLPOptimum(const SparseRow& row, double right_hand_side,
                      const SparseRow& lp_optimum) {
  return row.DotProduct(lp_optimum) > right_hand_side;
}

}  // namespace

std::optional<AggregatedRow> MIRRounder::RoundAggregatedRow(
    const Solver& solver, AggregatedRow aggregated_row) const {
  const double f0 = Fractionality(solver, aggregated_row.right_hand_side);
  aggregated_row.variable_coefficients.Transform(
      [f0, &solver](ColIndex col, double coefficient) {
        return solver.mip_data().integer_variables().contains(col)
                   ? MIRRoundInteger(solver, coefficient, f0)
                   : MIRRoundContinuous(coefficient, f0);
      });
  aggregated_row.slack_coefficients.Transform(
      [f0, &solver](RowIndex row, double coefficient) {
        return HasIntegralSlackVariable(solver, row)
                   ? MIRRoundInteger(solver, coefficient, f0)
                   : MIRRoundContinuous(coefficient, f0);
      });
  aggregated_row.right_hand_side =
      solver.FloorWithTolerance(aggregated_row.right_hand_side);
  VLOG(3) << "MIR rounding succeeded.";
  return aggregated_row;
}

std::optional<AggregatedRow> StrongCGRounder::RoundAggregatedRow(
    const Solver& solver, AggregatedRow aggregated_row) const {
  // CG rounding cannot be applied if there are continuous variables with
  // negative coefficients.
  if (!std::all_of(aggregated_row.slack_coefficients.entries().begin(),
                   aggregated_row.slack_coefficients.entries().end(),
                   [&solver](const SparseEntry<RowIndex>& entry) {
                     return HasIntegralSlackVariable(solver, entry.index) ||
                            entry.value >= 0.0;
                   })) {
    VLOG(3)
        << "SKIP: Rows containing negative weight continuous slack variables "
           "cannot be rounded using CG.";
    return std::nullopt;
  }
  if (!std::all_of(aggregated_row.variable_coefficients.entries().begin(),
                   aggregated_row.variable_coefficients.entries().end(),
                   [&solver](const SparseEntry<ColIndex>& entry) {
                     return !solver.mip_data().integer_variables().contains(
                                entry.index) ||
                            entry.value >= 0.0;
                   })) {
    VLOG(3) << "SKIP: Rows containing negative weight continuous variables "
               "cannot be rounded using CG.";
    return std::nullopt;
  }

  const double f0 = Fractionality(solver, aggregated_row.right_hand_side);
  const double k = std::ceil(1 / f0) - 1;
  DCHECK_LE(1 / (k + 1), f0);
  DCHECK_LT(f0, 1 / k);
  aggregated_row.slack_coefficients.Transform(
      [&solver, k, f0](RowIndex row, double val) {
        if (HasIntegralSlackVariable(solver, row)) {
          return CGRoundInteger(solver, val, k, f0);
        }
        DCHECK(val >= 0.0);
        return 0.0;
      });
  aggregated_row.variable_coefficients.Transform(
      [&solver, k, f0](ColIndex col, double val) {
        if (solver.mip_data().integer_variables().contains(col)) {
          return CGRoundInteger(solver, val, k, f0);
        }
        DCHECK(val >= 0.0);
        return 0.0;
      });
  aggregated_row.right_hand_side =
      solver.FloorWithTolerance(aggregated_row.right_hand_side);
  VLOG(3) << "CG rounding succeeded.";
  return aggregated_row;
}

absl::StatusOr<std::vector<CutData>>
TableauRoundingSeparator::GenerateCuttingPlanes(const Solver& solver) {
  const int max_num_cuts = params_.max_num_cuts();
  if (!solver.lpi()->IsSolved()) {
    return absl::FailedPreconditionError(
        "Generator called without solving LP.");
  }
  if (!solver.lpi()->IsOptimal()) {
    return absl::FailedPreconditionError("Generator called on non-optimal LP.");
  }

  // Extract data from the LP solver.
  const RowIndex num_rows = solver.lpi()->GetNumberOfRows();
  ASSIGN_OR_RETURN(
      (const absl::StrongVector<RowIndex, bool> use_right_hand_side),
      ChooseActiveSidesByBasis(solver));
  const std::vector<ColOrRowIndex> col_or_row_in_basis =
      solver.lpi()->GetColumnsAndRowsInBasis();
  ASSIGN_OR_RETURN(
      const SparseRow lp_optimum, ([&solver]() -> absl::StatusOr<SparseRow> {
        ASSIGN_OR_RETURN(
            (const absl::StrongVector<ColIndex, double> primal_values),
            solver.lpi()->GetPrimalValues());
        return SparseRow(primal_values);
      }()));

  // Generate the cuts.
  std::vector<CutData> cutting_planes;
  cutting_planes.reserve(max_num_cuts);
  for (RowIndex tableau_row(0);
       tableau_row < num_rows && cutting_planes.size() < max_num_cuts;
       ++tableau_row) {
    // Gomory cuts are generated from rows of the simplex tableau where the
    // corresponding basic variable is integer but has a fractional value in
    // the LP solution.
    const ColIndex basic_column =
        col_or_row_in_basis[tableau_row.value()].col();
    if (basic_column == kInvalidCol) {
      VLOG(3) << "SKIP row " << tableau_row
              << ": corresponds to a slack variable";
      continue;
    }
    if (!solver.mip_data().integer_variables().contains(basic_column)) {
      VLOG(3) << "SKIP row " << tableau_row
              << ": corresponds to a continuous variable";
      continue;
    }
    if (solver.IsIntegerWithinTolerance(lp_optimum[basic_column])) {
      VLOG(3) << "SKIP row " << tableau_row
              << ": corresponds to an integer variable with integer value";
      continue;
    }
    VLOG(3) << "Using row " << tableau_row << " corresponding to column "
            << basic_column << " with fractional value "
            << lp_optimum[basic_column];
    ASSIGN_OR_RETURN(const SparseRow basis_row,
                     solver.lpi()->GetSparseRowOfBInverted(tableau_row));
    SparseCol row_weights;
    for (auto [index, value] : basis_row.entries()) {
      row_weights.AddEntry(RowIndex(index.value()), value);
    }
    AggregatedRow aggregated_row =
        AggregateByWeight(solver, row_weights, use_right_hand_side);

    // See if the rounders can generate a cut from this aggregated row.
    if (!TransformToNonNegativeVariables(solver, aggregated_row)) {
      VLOG(3) << "SKIP: Failed to transform variables for row " << tableau_row;
      continue;
    }

    const SparseRow& objective = solver.mip_data().objective();
    int total_cuts = 0;

    for (const std::unique_ptr<Rounder>& rounder : rounders_) {
      std::optional<AggregatedRow> rounded_row =
          rounder->RoundAggregatedRow(solver, aggregated_row);
      if (!rounded_row.has_value()) continue;
      if (!TransformBackToOriginalVariables(solver, *rounded_row)) {
        return absl::InternalError(
            "Failed to transform back to original variables, even though this "
            "should always be possible.");
      }
      SubstituteSlackVariables(solver, *rounded_row);

      SparseRow& row = rounded_row->variable_coefficients;

      if (!RemovesLPOptimum(row, rounded_row->right_hand_side, lp_optimum)) {
        VLOG(3) << "Final cut doesn't remove LP optimum.";
        continue;
      }

      // TODO(cgraczy): add a unique identifier to the name for multiple rounds.
      std::string cutname = rounder->GetName() + std::to_string(total_cuts++);

      const double objective_parallelism =
          row.DotProduct(objective) /
          std::sqrt(row.DotProduct(row) * objective.DotProduct(objective));

      bool is_forced = false;
      if (row.entries().size() == 1) {
        is_forced = true;
      }

      int number_of_non_zeros = row.entries().size();
      int number_of_integer_variables = row.entries().size();
      double efficacy =
          (row.DotProduct(lp_optimum) - rounded_row->right_hand_side) /
          sqrt(row.DotProduct(row));

      CutData cut(std::move(rounded_row->variable_coefficients),
                  rounded_row->right_hand_side, number_of_non_zeros,
                  number_of_integer_variables, objective_parallelism, efficacy,
                  cutname, is_forced);
      cutting_planes.push_back(std::move(cut));
    }
  }
  return cutting_planes;
}

}  // namespace minimip
