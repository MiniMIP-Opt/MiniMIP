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

#include "problem.h"

#include <fstream>

namespace minimip {

// Returns the first error encountered while validating MiniMipProblem.
// Returns empty string if `problem` is valid.
std::string FindErrorInMiniMipProblem(const MiniMipProblem& problem) {
  for (const MiniMipVariable& variable : problem.variables) {
    if (variable.name.empty()) {
      return "No variable name given.";
    }
    if (variable.upper_bound < variable.lower_bound) {
      return absl::StrFormat(
          "Variable %s's lower variable bound is less than its upper variable "
          "bound.",
          variable.name);
    }
  }

  for (const MiniMipConstraint& constraint : problem.constraints) {
    if (constraint.name.empty()) {
      return "No constraint name given.";
    }
    absl::flat_hash_set<int> var_indices_set;
    for (int var_index : constraint.var_indices) {
      if (var_index < 0) {
        return absl::StrFormat("Negative variable index in constraint %s.",
                               constraint.name);
      }
      if (var_index > problem.variables.size()) {
        return absl::StrFormat("Variable index %d does not exist.", var_index);
      }
      if (!var_indices_set.insert(var_index).second) {
        return absl::StrFormat(
            "Variable indices in constraint %s are not unique.",
            constraint.name);
      }
    }
    if (constraint.var_indices.size() != constraint.coefficients.size()) {
      return absl::StrFormat(
          "Constraint %s has unequal number of indices and coefficients.",
          constraint.name);
    }
    if (constraint.right_hand_side < constraint.left_hand_side) {
      return absl::StrFormat(
          "Constraint %s has lesser left- than right-hand side.",
          constraint.name);
    }
  }
  return "";
}

}  // namespace minimip
