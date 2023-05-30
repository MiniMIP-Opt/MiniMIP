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

CutStorage::CutStorage()
    : cuts_({}), active_cut_indices_({}), total_number_of_cuts_found_(0) {}

// Initialize CutStorage from initial separation round.
CutStorage::CutStorage(std::vector<CutData> cuts, std::vector<int> cut_indices)
    : cuts_(std::move(cuts)),
      active_cut_indices_(std::move(cut_indices)),
      total_number_of_cuts_found_(cuts_.size()) {
  DCHECK_LE(active_cut_indices_.size(), cuts_.size());
}

// Use this to initialize by deep copy from another matrix `m`. Under-the-hood
// we just use a private copy / move constructor and assignment operator.
void CutStorage::PopulateFromCutStorage(CutStorage cut_storage) {
  *this = std::move(cut_storage);
}

// ==========================================================================
// Methods for managing the cuts in the cut storage.
// ==========================================================================

// Add a cut to storage.
// TODO(cgraczy): Add checks to enforce correct initialization on cuts.
int CutStorage::AddCut(CutData&& cut_data) {
  DCHECK_EQ(cut_data.index(), 0);
  cuts_.push_back(cut_data);
  total_number_of_cuts_found_ += 1;
  return cut_data.index();
}

// Activate stored cut.
void CutStorage::ActivateCut(int cut_index) {
  DCHECK_GE(cut_index, 0);
  DCHECK_LT(cut_index, cuts_.size());
  active_cut_indices_.push_back(cut_index);
}

// Remove a single cut from storage.
void CutStorage::RemoveCut(const int& cut_index) {
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
const CutData& CutStorage::GetCut(const int& cut_index) const {
  DCHECK_GE(cut_index, 0);
  DCHECK_LT(cut_index, cuts_.size());
  return cuts_[cut_index];
}

}  // namespace minimip
