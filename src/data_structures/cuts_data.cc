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

#include "cuts_data.h"

#include <utility>

namespace minimip {

// ==========================================================================
// Constructors
// ==========================================================================

CutRegistry::CutRegistry()
    : cuts_({}), active_cut_indices_({}), total_number_of_cuts_found_(0) {}

// Initialize CutRegistry from initial separation round.
CutRegistry::CutRegistry(std::vector<CutData> cuts,
                         std::vector<int> cut_indices)
    : cuts_(std::move(cuts)),
      active_cut_indices_(std::move(cut_indices)),
      total_number_of_cuts_found_(cuts_.size()) {
  DCHECK_LE(active_cut_indices_.size(), cuts_.size());
}

// Use this to initialize by deep copy from another matrix `m`. Under-the-hood
// we just use a private copy / move constructor and assignment operator.
void CutRegistry::PopulateFromCutRegistry(CutRegistry cut_registry) {
  *this = std::move(cut_registry);
}

CutData CutRegistry::CreateCut(const MipData& mip_data,
                               const SparseRow& lp_optimum, SparseRow row,
                               double right_hand_side, std::string name,
                               bool is_forced) const {
  // Calculate characteristics of the cut.
  const SparseRow& objective = mip_data.objective();
  const int number_of_non_zeros = row.entries().size();

  int number_of_integer_variables = 0;

  for (const auto& entry : row.entries()) {
    if (mip_data.integer_variables().contains(entry.index)) {
      number_of_integer_variables++;
    }
  }

  const double objective_parallelism =
      row.DotProduct(objective) /
      std::sqrt(row.DotProduct(row) * objective.DotProduct(objective));

  double efficacy = ComputeEfficacy(row, right_hand_side, lp_optimum);

  // Create CutData object
  CutData cut(std::move(row), right_hand_side, number_of_non_zeros,
              number_of_integer_variables, objective_parallelism, efficacy,
              std::move(name), is_forced);

  return cut;
}

// ==========================================================================
// Methods for managing the cuts in the cut registry.
// ==========================================================================

// Add a cut to registry.
// TODO(cgraczy): Add checks to enforce correct initialization on cuts.
int CutRegistry::AddCut(CutData&& cut_data) {
  cut_data.SetIndex(cuts_.size());
  cuts_.push_back(std::move(cut_data));
  return cuts_.size() - 1;
}

// Activate stored cut.
void CutRegistry::ActivateCut(int cut_index) {
  DCHECK_GE(cut_index, 0);
  DCHECK_LT(cut_index, cuts_.size());
  active_cut_indices_.push_back(cut_index);
}

// Remove a single cut from registry.
void CutRegistry::RemoveCut(const int& cut_index) {
  DCHECK_GE(cut_index, 0);
  DCHECK_LT(cut_index, cuts_.size());
  auto it = cuts_.begin();
  it = it + cut_index - 1;
  cuts_.erase(it);

  for (int i = cut_index; i < cuts_.size(); i++) {
    cuts_[i].SetIndex(cuts_[i].index() - 1);
  }
}

// Getter for individual cuts.
const CutData& CutRegistry::GetCut(const int& cut_index) const {
  DCHECK_GE(cut_index, 0);
  DCHECK_LT(cut_index, cuts_.size());
  return cuts_[cut_index];
}

void CutRegistry::UpdateCutEfficacy(const minimip::SparseRow& lp_optimum) {
  for (CutData& cut : cuts_) {
    cut.SetEfficacy(
        ComputeEfficacy(cut.row(), cut.right_hand_side(), lp_optimum));
  }
}

}  // namespace minimip
