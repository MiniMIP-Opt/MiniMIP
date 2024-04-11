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

bool CheckMIPFeasibility(const absl::StrongVector<ColIndex, double> &primal_values, const absl::flat_hash_set<ColIndex> &integer_variables, double tolerance) {
  return std::all_of(integer_variables.begin(), integer_variables.end(),
                     [&primal_values, tolerance](ColIndex col) {
                       return std::abs(primal_values[col] - std::round(primal_values[col])) > tolerance;
                     });
}

absl::StatusOr<MiniMipResult> Solver::Solve() {
  // Initialize the mip tree and LP interface
  MipTree tree = this->mutable_mip_tree();
  LpInterface *lp = this->mutable_lpi();
  MiniMipResult result; // TODO: Initialize this with the best known solution
  std::deque<NodeIndex> node_queue = {kRootNode};

  result.best_solution.objective_value = std::numeric_limits<double>::infinity();
  result.best_solution.objective_value *= this->mip_data_.is_maximization() ? -1.0 : 1.0; // TODO: make this more general

  absl::flat_hash_set<ColIndex> integer_variables = this->mip_data_.integer_variables();
  DenseRow root_lower_bounds = this->mip_data_.lower_bounds();
  DenseRow root_upper_bounds = this->mip_data_.upper_bounds();


  while (!node_queue.empty()) {
    // Get the current node from the queue
    NodeIndex current_node = node_queue.front();
    node_queue.pop_front();
    NodeData current_node_data = tree.node(current_node);

    DenseRow current_lower_bounds = tree.RetrieveLowerBounds(current_node, root_lower_bounds);
    DenseRow current_upper_bounds = tree.RetrieveUpperBounds(current_node, root_upper_bounds);

    for (ColIndex col : integer_variables ) {
      CHECK_OK(lp->SetColumnBounds(col, current_lower_bounds[col], current_upper_bounds[col]));
    }

    // First run the cutting loop to add cuts to the LP before setting the node data
    #ifdef CUTTING
        // Add CuttingPlane loop here to add cuts to the LP before branching.
    #endif

    // Solve and check current LP solution
    CHECK_OK(lp->SolveLpWithDualSimplex());

    if (!lp->IsOptimal()){
      if (current_node == kRootNode) {
        // If the root node's LP solution is not optimal, return the problem status TODO: set result.solve_status etc.
        return absl::InternalError("The problem is infeasible or unbounded");
      }
      tree.CloseNodeAndReclaimNodesUpToRootIfPossible(current_node);
      continue;
    }

    // Check if solution is MIP feasible
    double objective_value = this->mip_data_.is_maximization() ? lp->GetObjectiveValue() : -lp->GetObjectiveValue(); // TODO: make general
    absl::StrongVector<ColIndex, double> primal_values = lp->GetPrimalValues().value();

    bool solution_is_incumbent = std::all_of(integer_variables.begin(), integer_variables.end(),
                                   [&primal_values, tolerance = params_.integrality_tolerance()](ColIndex col) {
                                     return std::abs(primal_values[col] - std::round(primal_values[col])) <= tolerance;
                                   });

    // if the root node is MIP feasible, we can skip the branch and bound process and return the solution
    if (solution_is_incumbent) {
      MiniMipSolution solution = {std::vector<double>(primal_values.begin(), primal_values.end()), objective_value};

      if( (this->mip_data_.is_maximization() and result.best_solution.objective_value < objective_value )
          or ( !this->mip_data_.is_maximization() and result.best_solution.objective_value > objective_value)  ) {
          result.best_solution = solution;
        }
        result.additional_solutions.push_back(solution);

        // If the root node is MIP optimal, return this solution
        if (current_node == kRootNode or node_queue.empty()) {
          result.solve_status = MiniMipSolveStatus::kOptimal;
          return result;
        }

    }
    else {

      // Setup branch and bound node data
      tree.SetLpRelaxationDataInNode(current_node, objective_value);

      // Find the variable with the maximum fractional part to branch on
      ColIndex branching_variable;
      double max_fractional_part = 0.0;

      for (ColIndex col : integer_variables) {
        double value = primal_values[col];
        double fractional_part = value - std::floor(value);
        if (fractional_part > max_fractional_part) {
          max_fractional_part = fractional_part;
          branching_variable = col;
        }
      }

      // Branch on the variable with the maximum fractional part
      //ASSIGN_OR_RETURN(LpInterface::StrongBranchResult strong_branch_result,
      //                 lp->SolveDownAndUpStrongBranch(branching_variable, primal_values[branching_variable], 0));

      //result.best_dual_bound = lp->GetObjectiveValue();

        // Assume AddNodeByBranchingFromParent is a function that handles the logic of creating child nodes
        NodeIndex left_child = tree.AddNodeByBranchingFromParent(current_node,
                                                                 branching_variable,
                                                                 true,
                                                                 primal_values[branching_variable]);
        NodeIndex right_child = tree.AddNodeByBranchingFromParent(current_node,
                                                                    branching_variable,
                                                                    false,
                                                                    primal_values[branching_variable]);

        // Add logic to process child nodes...
        node_queue.push_front(left_child);
        node_queue.push_back(right_child);

    }

  }
}
}  // namespace minimip
