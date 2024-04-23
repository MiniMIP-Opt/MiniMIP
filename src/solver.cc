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
#include <random>
#include <vector>

namespace minimip {

absl::Status Solver::Solve() {

  std::deque<NodeIndex> node_queue = {kRootNode};
  DenseRow root_lower_bounds = mip_data_.lower_bounds();
  DenseRow root_upper_bounds = mip_data_.upper_bounds();

  while (!node_queue.empty()) {
    // Get the current node from the queue
    NodeIndex current_node = node_queue.front();
    node_queue.pop_front();
    NodeData current_node_data = mip_tree_.node(current_node);

    // Set the LP relaxation bounds for the current node
    DenseRow current_lower_bounds =
        mip_tree_.RetrieveLowerBounds(current_node, root_lower_bounds);
    DenseRow current_upper_bounds =
        mip_tree_.RetrieveUpperBounds(current_node, root_upper_bounds);

    for (ColIndex col : mip_data_.integer_variables()) {
      CHECK_OK(lpi_->SetColumnBounds(col, current_lower_bounds[col],
                                   current_upper_bounds[col]));
    }

// First run the cutting loop to add cuts to the LP before setting the node data
#ifdef CUTTING
    // Add CuttingPlane loop here to add cuts to the LP before branching.
#endif

    // Solve and check current LP solution
    CHECK_OK(lpi_->SolveLpWithDualSimplex());

    if (!lpi_->IsOptimal()) {
      if (current_node == kRootNode) {
        // If the root node's lpi_ solution is not optimal, return the problem
        // status TODO: set result.solve_status etc.
        return absl::InternalError("The problem is infeasible or unbounded");
      }
      mip_tree_.CloseNodeAndReclaimNodesUpToRootIfPossible(current_node);
      continue;
    }

    // Set the lpi_ relaxation data in the current node
    double objective_value =
        mip_data_.is_maximization()
            ? lpi_->GetObjectiveValue()
            : -lpi_->GetObjectiveValue();  // TODO: make general
    mip_tree_.SetLpRelaxationDataInNode(current_node, objective_value);

    // Check if solution is MIP feasible
    absl::StrongVector<ColIndex, double> primal_values =
        lpi_->GetPrimalValues().value();
    bool solution_is_incumbent =
        mip_data_.SolutionIsIntegral(primal_values,params_.integrality_tolerance());

    // if the root node is MIP feasible, we can skip the branch and bound
    // process and return the solution
    if (solution_is_incumbent) {
      MiniMipSolution solution = {
          std::vector<double>(primal_values.begin(), primal_values.end()),
          objective_value};

      if ((mip_data_.is_maximization() and
           result_.best_solution.objective_value < objective_value) or
          (!mip_data_.is_maximization() and
           result_.best_solution.objective_value > objective_value)) {
        result_.best_solution = solution;
      }
      result_.additional_solutions.push_back(solution);

      // If the root node is MIP optimal, return this solution
      if (current_node == kRootNode or node_queue.empty()) {
        result_.solve_status = MiniMipSolveStatus::kOptimal;
        return absl::OkStatus();
      }

    } else {
      // Find the variable with the maximum fractional part to branch on
      ColIndex branching_variable;
      double max_fractional_part = 0.0;

      for (ColIndex col : mip_data_.integer_variables()) {
        double value = primal_values[col];
        double fractional_part = value - std::floor(value);
        if (fractional_part > max_fractional_part) {
          max_fractional_part = fractional_part;
          branching_variable = col;
        }
      }

      // Assume AddNodeByBranchingFromParent is a function that handles the
      // logic of creating child nodes
      NodeIndex left_child = mip_tree_.AddNodeByBranchingFromParent(
          current_node, branching_variable, true,
          primal_values[branching_variable]);
      NodeIndex right_child = mip_tree_.AddNodeByBranchingFromParent(
          current_node, branching_variable, false,
          primal_values[branching_variable]);

      // Add logic to process child nodes...
      node_queue.push_front(left_child);
      node_queue.push_back(right_child);
    }
  }
  return absl::OkStatus();
}
}  // namespace minimip
