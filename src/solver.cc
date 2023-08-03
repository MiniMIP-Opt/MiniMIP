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

#include "src/solver.h"
namespace minimip {

absl::StatusOr<MiniMipResult> Solver::Solve() {
  MipTree tree = this->mip_tree_;
  LPInterface* lp = this->mutable_lpi();

  lp->PopulateFromMipData(this->mip_data_);

  lp->SolveLPWithDualSimplex();
  if (!lp->IsOptimal()) {
    return absl::InternalError("LP is infeasible");
  }

#ifdef CUTTING
  // Add CuttingPlane loop here to add cuts to the LP before branching.
#endif

  double objective_value = lp->GetObjectiveValue();
  absl::StrongVector<ColIndex, double> primal_values =
      lp->GetPrimalValues().value();
  absl::flat_hash_set<ColIndex> integer_variables =
      this->mip_data_.integer_variables();

  tree.SetLpRelaxationDataInNode(kRootNode, objective_value);

  // Check if the solution is integral
  absl::StrongVector<ColIndex, double> fractional_parts(primal_values.size());
  SparseEntry<ColIndex> max_fractional_part = SparseEntry(kInvalidCol, 0);

  for (ColIndex col : integer_variables) {
    if (!this->IsIntegerWithinTolerance(primal_values[col])) {
      double fractional =
          primal_values[col] - this->FloorWithTolerance(primal_values[col]);
      if (fractional > max_fractional_part.value) {
        max_fractional_part.value = fractional;
        max_fractional_part.index = col;
      }
      // Get all branching candidates and Branch on the variable with the
      // largest fractional part
      fractional_parts[col] = fractional;
    }
  }

  if (fractional_parts.empty()) {
    std::vector<double> solution_values(primal_values.size());
    std::transform(primal_values.begin(), primal_values.end(),
                   solution_values.begin(), [](double val) { return val; });

    // The solution is integral, thus primal optimal
    MiniMipSolution solution(solution_values, objective_value);
  } else {
    // The solution is not integral, thus branch on the variable with the
    // largest fractional part
    ColIndex branching_variable = max_fractional_part.index;
    double fractional_part = max_fractional_part.value;

    // Branch on the variable with the largest fractional part
    // Create two new nodes

    // Add the nodes to the tree
    bool branch_down = (fractional_part < 0.5);

    // for each candidate create new nodes??

    tree.AddNodeByBranchingFromParent(kRootNode, branching_variable, true,
                                      fractional_part);
    // Update the LP relaxation data in the nodes
    // Solve the LP relaxation in the nodes
    // If the LP relaxation is infeasible, then prune the node
    // If the LP relaxation is feasible, then check if the solution is integral
    // If the solution is integral, then update the incumbent
    // If the solution is not integral, then branch on a variable
  }

  return absl::UnimplementedError("Solve isn't implemented yet.");
}

}  // namespace minimip
