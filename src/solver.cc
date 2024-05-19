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

#include "src/solver.h"

#include <queue>
#include <vector>

namespace minimip {
namespace {
absl::Status SolveCurrentNode() {
  return absl::OkStatus();
}
}  // namespace

absl::Status Solver::Solve() {
  std::deque<NodeIndex> node_queue = {kRootNode};
  LpInterface* lp = mutable_lpi();

  bool backjump = true;
  int num_consecutive_infeasible_nodes = 0;

  while (true) {
    if (node_queue.empty()) {
      LOG(INFO) << "Explored all nodes. Terminating branch-and-bound.";
      break;
    }

    // ==========================================================================
    // Identify next node to explore.
    // ==========================================================================

    // Get the current node from the queue.
    // TODO(lpawel): Implement a dedicated node queue class.
    const NodeIndex current_node =
        backjump ? node_queue.back() : node_queue.front();

    if (backjump) {
      node_queue.pop_back();
    } else {
      node_queue.pop_front();
    }
    backjump = false;

    // ==========================================================================
    // Solve the lp relaxation check infeasibility
    // ==========================================================================

    // Set the LP relaxation bounds for the current node
    const DenseRow current_lower_bounds =
        mip_tree_.RetrieveLowerBounds(current_node, mip_data_.lower_bounds());
    const DenseRow current_upper_bounds =
        mip_tree_.RetrieveUpperBounds(current_node, mip_data_.upper_bounds());

    for (ColIndex col : mip_data_.integer_variables()) {
      CHECK_OK(lp->SetColumnBounds(col, current_lower_bounds[col],
                                     current_upper_bounds[col]));
    }

    // TODO(CG): First run the cutting loop to add cuts to the LP before setting the node data
    CHECK_OK(lp->SolveLpWithDualSimplex());

    if (lp->IsOptimal()) {
      CHECK_OK(cut_runner_->SeparateCurrentLPSolution(*this));
    }

    if (!lp->IsOptimal()) {
      if (current_node == kRootNode) {
        // If the root node's lpi_ solution is not optimal, return the problem
        // status TODO: set result.solve_status etc.
        return absl::InternalError("The problem is infeasible or unbounded");
      }
      mip_tree_.CloseNodeAndReclaimNodesUpToRootIfPossible(current_node);

      ++num_consecutive_infeasible_nodes;
      if (num_consecutive_infeasible_nodes > params_.backtrack_limit()) {
        num_consecutive_infeasible_nodes = 0;
        backjump = true;
      }
      continue;
    }

    // The node was not infeasible -- reset the counter of infeasible nodes.
    num_consecutive_infeasible_nodes = 0;

    // ==========================================================================
    // Set lp relaxation data in the current node and check for integrality.
    // ==========================================================================

    const double objective_value = lp->GetObjectiveValue();
    mip_tree_.SetLpRelaxationDataInNode(current_node, objective_value);

    absl::StrongVector<ColIndex, double> primal_values = lp->GetPrimalValues().value();
    const bool solution_is_incumbent = mip_data_.SolutionIsIntegral(
        primal_values, params_.integrality_tolerance());

    // The node is integral, there will be no branching from that node.
    if (solution_is_incumbent) {
      const MiniMipSolution solution = {
          std::vector<double>(primal_values.begin(), primal_values.end()),
          objective_value};

      CHECK_OK(result_.AddSolution(solution));

      // If the root node is MIP optimal, return this solution
      if (current_node == kRootNode) {
        break;
      }
      backjump = true;
      continue;
    }

    // ==========================================================================
    // Branching
    // ==========================================================================

    // Find the variable with the maximum fractional part to branch on
    ColIndex branching_variable = branching_interface_->NextBranchingVariable(*this).value();

    // Create binary child nodes for chosen branching variable
    const NodeIndex down_child = mip_tree_.AddNodeByBranchingFromParent(
        current_node, branching_variable, true,
        primal_values[branching_variable]);
    const NodeIndex up_child = mip_tree_.AddNodeByBranchingFromParent(
        current_node, branching_variable, false,
        primal_values[branching_variable]);

    // Add logic to process child nodes...
    node_queue.push_front(down_child);
    node_queue.push_back(up_child);
  }
  result_.solve_status = MiniMipSolveStatus::kOptimal;
  return absl::OkStatus();
}
}  // namespace minimip
