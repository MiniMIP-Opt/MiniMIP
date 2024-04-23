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

#include "mip_data.h"

namespace minimip {
namespace {

bool IsIntegralCoefficient(double coefficient) {
  return coefficient == std::round(coefficient);
}

}  // namespace
// ==========================================================================
// Constructor functionality
// ==========================================================================

// TODO(cgraczyk): Enforce minimization assumption on the internal
// representation.
MipData::MipData(const MiniMipProblem& problem)
    : solution_hints_(problem.hints.size()),
      variable_names_(problem.variables.size()),
      lower_bounds_(problem.variables.size()),
      upper_bounds_(problem.variables.size()),
      variable_types_(problem.variables.size()),
      constraint_names_(problem.constraints.size()),
      left_hand_sides_(problem.constraints.size()),
      right_hand_sides_(problem.constraints.size()),
      is_integral_constraint_(problem.constraints.size()),
      constraint_matrix_(ColIndex(problem.variables.size()),
                         RowIndex(problem.constraints.size())) {
  DCHECK(FindErrorInMiniMipProblem(problem).empty());

  problem_name_ = problem.name;
  is_maximization_ = problem.is_maximization;
  objective_offset_ = problem.objective_offset;

  for (ColIndex col_idx(0); col_idx < problem.variables.size(); ++col_idx) {
    const MiniMipVariable& variable = problem.variables[col_idx.value()];

    lower_bounds_[col_idx] = variable.lower_bound;
    upper_bounds_[col_idx] = variable.upper_bound;

    if (variable.is_integer) {
      variable_types_[col_idx] = VariableType::kInteger;
      integer_variables_.insert(col_idx);

      // TODO(lpawel): Move the logic that makes bounds integral to the
      // preprocessor.
      CHECK_LE(lower_bounds_[col_idx], upper_bounds_[col_idx]);
    } else {
      variable_types_[col_idx] = VariableType::kFractional;
    }
    variable_names_[col_idx] = variable.name;
    if (variable.objective_coefficient != 0) {
      objective_.AddEntry(ColIndex(col_idx), variable.objective_coefficient);
    }
  }
  objective_.CleanUpIfNeeded();

  for (int i = 0; i < problem.hints.size(); ++i) {
    solution_hints_[i] = problem.hints.at(i);
  }

  for (RowIndex row_idx(0); row_idx < problem.constraints.size(); ++row_idx) {
    const MiniMipConstraint& constraint = problem.constraints[row_idx.value()];
    left_hand_sides_[row_idx] = constraint.left_hand_side;
    right_hand_sides_[row_idx] = constraint.right_hand_side;
    constraint_names_[row_idx] = constraint.name;

    SparseRow sparse_constraint;
    for (int col_idx = 0; col_idx < constraint.var_indices.size(); ++col_idx) {
      sparse_constraint.AddEntry(ColIndex{constraint.var_indices[col_idx]},
                                 constraint.coefficients[col_idx]);
    }
    is_integral_constraint_[row_idx] = std::all_of(
        sparse_constraint.entries().begin(), sparse_constraint.entries().end(),
        [this](const SparseEntry<ColIndex>& entry) {
          return IsIntegralCoefficient(entry.value) &&
                 integer_variables_.contains(entry.index);
        });
    constraint_matrix_.PopulateRow(RowIndex(row_idx), sparse_constraint);
  }

  integer_variables_.reserve(problem.variables.size());

  for (ColIndex col_idx(0); col_idx < variable_types_.size(); ++col_idx) {
    if (variable_types_[col_idx] == VariableType::kInteger) {
      integer_variables_.insert(col_idx);
    }
  }
}

bool MipData::SolutionIsIntegral(const absl::StrongVector<ColIndex, double>& solution_values, double tolerance) const {
  return std::all_of(
      integer_variables_.begin(), integer_variables_.end(),
      [&solution_values, tolerance](ColIndex col) {
        return std::abs(solution_values[col] - std::round(solution_values[col])) <=
            tolerance;
      });
}

}  // namespace minimip
