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

#include "mip_data.h"

namespace minimip {
// ==========================================================================
// Constructor functionality
// ==========================================================================

MipData::MipData()
    : problem_name_(),
      is_maximization_(),
      objective_offset_(),
      objective_(SparseRow()),
      solution_hints_({}),
      lower_bounds_({}),
      upper_bounds_({}),
      left_hand_sides_({}),
      right_hand_sides_({}),
      constraint_matrix_(minimip::ColIndex(0), minimip::RowIndex(0)),
      variable_types_({}),
      variable_names_({}),
      constraint_names_({}) {}

// TODO(cgraczyk): Enforce minimization assumption on the internal
// representation.
MipData::MipData(const MiniMipProblem& problem)
    : objective_(SparseRow()),
      solution_hints_(problem.hints.size()),
      lower_bounds_(problem.variables.size()),
      upper_bounds_(problem.variables.size()),
      left_hand_sides_(problem.constraints.size()),
      right_hand_sides_(problem.constraints.size()),
      binary_variables_({}),
      integer_variables_({}),
      constraint_matrix_(ColIndex(problem.variables.size()),
                         RowIndex(problem.constraints.size())),
      variable_types_(problem.variables.size()),
      variable_names_(problem.variables.size()),
      constraint_names_(problem.constraints.size()) {
  DCHECK(FindErrorInMiniMipProblem(problem).empty());

  problem_name_ = problem.name;
  is_maximization_ = problem.is_maximization;
  objective_offset_ = problem.objective_offset;

  for (int col_idx = 0; col_idx < problem.variables.size(); ++col_idx) {
    const MiniMipVariable& variable = problem.variables[col_idx];

    lower_bounds_[col_idx] = variable.lower_bound;
    upper_bounds_[col_idx] = variable.upper_bound;

    if (variable.is_integer) {
      variable_types_[col_idx] = VariableType::kInteger;
    } else
      variable_types_[col_idx] = VariableType::kFractional;

    variable_names_[col_idx] = variable.name;

    if (variable.objective_coefficient != 0) {
      objective_.AddEntry(ColIndex(col_idx), variable.objective_coefficient);
    }
  }

  for (int i = 0; i < problem.hints.size(); ++i) {
    solution_hints_[i] = problem.hints.at(i);
  }

  for (int row_idx = 0; row_idx < problem.constraints.size(); ++row_idx) {
    SparseRow sparse_constraint;
    const MiniMipConstraint& constraint = problem.constraints[row_idx];

    left_hand_sides_[row_idx] = constraint.left_hand_side;
    right_hand_sides_[row_idx] = constraint.right_hand_side;
    constraint_names_[row_idx] = constraint.name;

    for (int col_idx = 0; col_idx < constraint.var_indices.size(); ++col_idx) {
      sparse_constraint.AddEntry((ColIndex)constraint.var_indices[col_idx],
                                 constraint.coefficients[col_idx]);
    }
    constraint_matrix_.PopulateRow(RowIndex(row_idx), sparse_constraint);
  }

  binary_variables_.reserve(problem.variables.size());
  integer_variables_.reserve(problem.variables.size());

  for (int col_idx = 0; col_idx < variable_types_.size(); ++col_idx) {
    if (variable_types_[col_idx] == VariableType::kBinary) {
      //@TODO(cgraczyk): This is always empty because variable type kBinary is
      // never set.
      binary_variables_.push_back(col_idx);
    }
    if (variable_types_[col_idx] == VariableType::kInteger) {
      integer_variables_.push_back(col_idx);
    }
  }
}

}  // namespace minimip
