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

#include "random_branching.h"

#include <random>

namespace minimip {

const absl::StatusOr<BranchingVariable> RandomBranching::NextBranchingVariable(
    const SolverContextInterface& context) const {
  // Return a random variable.
  int integer_variables = context.mip_data().integer_variables().size();
  BranchingVariable branching_variable;

  if (integer_variables > 0) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis(0, integer_variables - 1);
    branching_variable.index = ColIndex(dis(gen));
    branching_variable.branching_up_first =
        params_.branching_direction() ==
        BranchingParameters_BranchingDirection_UP;
    return branching_variable;
  }

  return absl::InvalidArgumentError("No integer variables to branch on.");
}

}  // namespace minimip