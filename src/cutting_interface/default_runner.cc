// Copyright 2023 the MiniMIP Project
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

bool DefaultRunner::CutCondition(const Solver& solver) {
  // max cutrounds per node

  // max cuts
  /* if the LP is infeasible = cutoff , unbounded, exceeded the objective limit
   * or a global performance limit was reached, we don't need to separate cuts
   * (the global limits are only checked at the root node in order to not query
   * system time too often)
   */

  return false;
}

absl::Status DefaultRunner::SeparateCurrentLPSolution(
    const Solver& solver, LPInterface* mutable_lpi,
    CutRegistry& mutable_cut_registry) {
  int i = 0;

  while (i < 1) {  // CutCondition(solver)
    std::vector<int> new_cut_indices;

    for (const std::unique_ptr<CutGenerator>& generator : generators_) {
      absl::StatusOr<std::vector<CutData>> cuts =
          generator->GenerateCuttingPlanes(solver);
      if (cuts.status() != absl::OkStatus()) {
        return cuts.status();
      }

      absl::StatusOr<std::vector<CutData>> filtered_cuts =
          selector_->SelectCuttingPlanes(solver, cuts.value());

      if (filtered_cuts.status() != absl::OkStatus()) {
        return filtered_cuts.status();
      }

      for (CutData& cut : filtered_cuts.value()) {
        new_cut_indices.push_back(mutable_cut_registry.AddCut(std::move(cut)));
      }
    }

    for (int cut_index : new_cut_indices) {
      const CutData& cut = mutable_cut_registry.GetCut(cut_index);
      RETURN_IF_ERROR(mutable_lpi->AddRow(cut.row(), -mutable_lpi->Infinity(),
                                          cut.right_hand_side(), cut.name()));
    }
    RETURN_IF_ERROR(mutable_lpi->SolveLPWithDualSimplex());
    i++;
  }
  return absl::OkStatus();
};

}  // namespace minimip
