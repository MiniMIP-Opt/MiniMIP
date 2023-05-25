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

#include "cuts_runner.h"

namespace minimip {

//
absl::Status CutRunner::SeparateCurrentLPSolution(
    const Solver& solver, CutStorage& mutable_cut_storage) {
  for (const std::unique_ptr<Separator>& separator : separators_) {
    absl::StatusOr<std::vector<CutData>> cuts =
        separator->GenerateCuttingPlanes(solver);
    if (cuts.status() != absl::OkStatus()) {
      return cuts.status();
    }

    absl::StatusOr<std::vector<CutData>> filtered_cuts =
        selector_->SelectCuttingPlanes(solver, cuts.value());
    for (CutData& cut : filtered_cuts.value()) {
      mutable_cut_storage.AddCut(std::move(cut));
    }
  }
  return absl::OkStatus();
};

}  // namespace minimip
