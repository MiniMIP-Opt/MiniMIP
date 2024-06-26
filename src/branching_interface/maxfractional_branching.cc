// Copyright (2024) the MiniMIP Project
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

#include "maxfractional_branching.h"

namespace minimip {

const absl::StatusOr<BranchingVariable>
MaxFractionalBranching::NextBranchingVariable(
    const SolverContextInterface& context) const {
  // Return a random variable.
  auto integer_variables = context.mip_data().integer_variables();

  BranchingVariable branching_variable;
  if (!integer_variables.empty()) {
    double max_fractional_part = 0.0;
    absl::StrongVector<ColIndex, double> primal_values =
        context.lpi()->GetPrimalValues().value();

    for (ColIndex col : integer_variables) {
      double value = primal_values[col];
      if (!context.IsIntegerWithinTolerance(value)) {
        double fractional_part = abs(value - std::floor(value));
        if (fractional_part > max_fractional_part) {
          max_fractional_part = fractional_part;
          branching_variable.index = col;
        }
      }
    }
    branching_variable.branching_up_first =
        (params_.branching_direction() ==
         BranchingParameters_BranchingDirection_UP);

    if (params_.branching_direction() ==
        BranchingParameters_BranchingDirection_UNSPECIFIED) {
      branching_variable.branching_up_first = (max_fractional_part <= 0.5);
    }
    return branching_variable;
  }

  return absl::InvalidArgumentError("No integer variables to branch on.");
}

}  // namespace minimip