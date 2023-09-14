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

#include <queue>
namespace minimip {

absl::StatusOr<MiniMipResult> Solver::Solve() {
  MipTree tree;
  LPInterface* lp = this->mutable_lpi();
  MiniMipResult result;

  CHECK_OK(lp->PopulateFromMipData(this->mip_data_));
  CHECK_OK(lp->SolveLPWithDualSimplex());
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

  std::deque<NodeIndex> node_deque;
  node_deque.push_back(kRootNode);

  while (!node_deque.empty()) {
    NodeIndex current_node = node_deque.front();
    node_deque.pop_front();

    // Apply the branching decisions and solve the LP relaxation in the current
    // node
    NodeData current_node_data = tree.node(current_node);
    if (current_node_data.parent != kInvalidNode) {
      // Apply branching decisions here based on currentNodeData
      // For example, update the bounds of currentNodeData.branch_variable
    }

#ifdef CUTTING
    // Add CuttingPlane loop here to add cuts to the LP before branching.
#endif
    CHECK_OK(lp->SolveLPWithDualSimplex());  // Solve the LP

    if (!lp->IsOptimal()) {
      // The node is infeasible, prune it
      continue;
    }

    // Update the objective value for the current node
    objective_value = lp->GetObjectiveValue();
    current_node_data.lp_objective_value = objective_value;

    // Check if the solution is integral
    bool is_integral = true;
    ColIndex branching_variable;
    double max_fractional_part = 0;

    for (ColIndex col : integer_variables) {
      if (!this->IsIntegerWithinTolerance(primal_values[col])) {
        is_integral = false;
        double fractional =
            primal_values[col] - this->FloorWithTolerance(primal_values[col]);
        if (fractional > max_fractional_part) {
          max_fractional_part = fractional;
          branching_variable = col;
        }
      }
    }

    if (is_integral) {
      // Update the incumbent solution if it's better
      if (objective_value < result.best_solution.objective_value) {
        result.best_solution.objective_value = objective_value;
        result.best_solution.variable_values =
            std::vector<double>(primal_values.begin(), primal_values.end());
      }
    } else {
      // Branch on the variable with the largest fractional part
      NodeIndex left_child = tree.AddNodeByBranchingFromParent(
          current_node, branching_variable, true, max_fractional_part);
      NodeIndex right_child = tree.AddNodeByBranchingFromParent(
          current_node, branching_variable, false, max_fractional_part);

      // Update node data
      NodeData left_node_data = tree.node(left_child);
      NodeData right_node_data = tree.node(right_child);

      left_node_data.branch_variable = branching_variable;
      left_node_data.branch_down = true;
      left_node_data.branch_primal_value_in_parent_lp = max_fractional_part;

      right_node_data.branch_variable = branching_variable;
      right_node_data.branch_down = false;
      right_node_data.branch_primal_value_in_parent_lp = max_fractional_part;

      // Add the new nodes to the deque for further processing
      node_deque.push_front(left_child);
      node_deque.push_back(right_child);
    }
  }

  if (result.solve_status != MiniMipSolveStatus::kOptimal) {
    return absl::InternalError("No feasible solution found");
  }

  return result;
}

}  // namespace minimip
