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

#include "default_runner.h"

namespace minimip {

bool DefaultRunner::CutCondition(const SolverContextInterface& context) {
  VLOG(10) << "calling CutCondition().";
  // max cutrounds per node

  // max cuts
  /* if the LP is infeasible = cutoff , unbounded, exceeded the objective limit
   * or a global performance limit was reached, we don't need to separate cuts
   * (the global limits are only checked at the root node in order to not query
   * system time too often)
   */
  if (context.lpi()->IsPrimalInfeasible() ||
      context.lpi()->IsPrimalUnbounded() || context.lpi()->IsDualInfeasible() ||
      context.lpi()->IsDualUnbounded() ||
      context.lpi()->TimeLimitIsExceeded() ||
      context.lpi()->IterationLimitIsExceeded() ||
      context.lpi()->ObjectiveLimitIsExceeded()) {
    VLOG(3) << "CutRunner: LP is infeasible, unbounded, or reached a limit.";
    return false;
  }

  if (!context.lpi()->IsSolved()) {
    VLOG(3) << "CutRunner: No cuts were added to LP in the last round.";
    return false;
  }

  if (num_cutrounds >= params_.max_iterations()) {
    VLOG(3) << "CutRunner: Maximum number of cut rounds.";
    return false;
  }

  if (num_of_cuts_added_since_last_run >= params_.max_num_cuts_at_node()) {
    VLOG(3) << "CutRunner: Maximum number of cuts added in this node.";
    return false;
  }
  if (context.cut_registry().active_cuts().size() >= params_.max_num_cuts()) {
    VLOG(3) << "CutRunner: Maximum number of active cuts in the LP.";
    return false;
  }

  // Continue the cut separation.
  return true;
}

absl::Status DefaultRunner::SeparateCurrentLPSolution(
    SolverContextInterface& context) {
  VLOG(10) << "calling SeparateCurrentLPSolution().";
  LpInterface* mutable_lpi = context.mutable_lpi();

  bool should_separate = CutCondition(context);

  while (should_separate) {
    std::vector<int> new_cut_indices;

    for (const std::unique_ptr<CutGeneratorInterface>& generator :
         generators_) {
      ASSIGN_OR_RETURN(std::vector<CutData> cuts,
                       generator->GenerateCuttingPlanes(context));
      if (cuts.empty()) {
        continue;
      }

      ASSIGN_OR_RETURN(std::vector<CutData> filtered_cuts,
                       selector_->SelectCuttingPlanes(context, cuts));

      if (filtered_cuts.empty()) {
        continue;
      }

      for (CutData& cut : filtered_cuts) {
        new_cut_indices.push_back(
            context.mutable_cut_registry().AddCut(std::move(cut)));
      }
    }

    // If no new cuts have been added to the LP, we can stop the cut separation.
    if (new_cut_indices.empty()) {
      break;
    }

    for (int cut_index : new_cut_indices) {
      const CutData& cut = context.mutable_cut_registry().GetCut(cut_index);
      RETURN_IF_ERROR(mutable_lpi->AddRow(cut.row(), -mutable_lpi->Infinity(),
                                          cut.right_hand_side(), cut.name()));
    }

    RETURN_IF_ERROR(mutable_lpi->SolveLpWithDualSimplex());
    num_cutrounds++;
    should_separate = CutCondition(context);
  }
  return absl::OkStatus();
};

}  // namespace minimip
