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

#ifndef SRC_DATA_STRUCTURES_CUTS_DATA_H_
#define SRC_DATA_STRUCTURES_CUTS_DATA_H_

#include "src/data_structures/mip_data.h"

#include <utility>

#include "src/data_structures/strong_sparse_matrix.h"

namespace minimip {

// ==========================================================================
// API Input Datastructures
// ==========================================================================

// This is the cutting plane object used in the cutting plane interface to
// add more information to the generated rows from the cut generators.
struct CuttingPlane {
  SparseRow row;
  bool is_active;
  bool added_at_root;
  bool forced;
  unsigned int from_sepa_round_n;
  unsigned int number_of_non_zeros;
  unsigned int cut_position;
  unsigned int integer_variable_support;
  double original_score;
  double current_score;
  double efficacy;
  std::string name;
};
using Cut = CuttingPlane;

// The cut storage is used to track the globally valid cuts generated while
// solving the Mixed Integer Problem. The storage functions as a register and
// provides access to cuts for selection and the lp.
class CutStorage {
  CutStorage() {}

  // Initialize CutStorage from initial separation round.
  CutStorage(std::vector<Cut> cuts,std::vector<unsigned int> cut_positions)
    : cuts_(std::move(cuts)),
      active_cuts_positions_(std::move(cut_positions)),
      current_number_of_cuts_(cuts.size()),
      total_number_of_cuts_found_(cuts.size()),
      current_number_of_active_cuts_(cut_positions.size())
      {
        DCHECK_LE(cut_positions.size(),cuts.size());
      }


  // CutStorage is not copyable to make sure a copy will not be
  // triggered by accident (copy constructor and assign operator are private).
  // CutStorage is (no-throw) moveable.
  CutStorage(CutStorage&&) noexcept = default;
  CutStorage& operator = (CutStorage&&) noexcept = default;

  // Use this to initialize by deep copy from another matrix `m`. Under-the-hood
  // we just use a private copy / move constructor and assignment operator.
  void PopulateFromCutStorage(CutStorage cut_storage){
    *this = std::move(cut_storage);
  }

  // Add a cut to storage.
  void AddCut(Cut& cut) {
    cuts_.push_back(cut);
    current_number_of_cuts_ += 1;
    total_number_of_cuts_found_ += 1;
    if (cut.forced) {
      ActivateCut(cut);
    }
  }

  // Add a vector of cuts to storage.
  void AddCuts(std::vector<Cut> cuts){
    cuts_.insert(cuts_.end(),cuts.begin(),cuts.end());
    current_number_of_cuts_ += cuts.size();
    total_number_of_cuts_found_+=cuts.size();
  }

  // Activate stored cut.
  void ActivateCut(const Cut& cut){
    active_cuts_positions_.push_back(cut.cut_position);
    current_number_of_active_cuts_ += 1;
  }

  // Activate cuts.
  void ActivateCuts(std::vector<unsigned int>& active_cuts){
    DCHECK_LE(active_cuts.size(),cuts_.size());
    current_number_of_active_cuts_ = active_cuts.size();
    active_cuts_positions_ = std::move(active_cuts);
  }

  // Remove a single cut from storage.
  void RemoveCut(const unsigned int& cut_position){
    auto it = cuts_.begin();
    it = it+cut_position-1;
    cuts_.erase(it);

    for(unsigned int i = cut_position; i < cuts_.size(); i++)
      cuts_[i].cut_position -= 1;

    current_number_of_cuts_ -= 1;
  }

  // Getter for individual cuts.
  const Cut& GetCut(const unsigned int& cut_position) const {
    DCHECK_LT(cut_position, current_number_of_cuts_);
    return cuts_[cut_position];
  }

  // Getter for some subset of cuts.
  std::vector<Cut> GetCuts(const std::vector<unsigned int> cut_indices) const {
    std::vector<Cut> CutSubset;
    CutSubset.reserve(cut_indices.size());
    for(unsigned int cut_indice : cut_indices)
      CutSubset.push_back(cuts_[cut_indice]);
    return CutSubset;
  }

  // Getter for all cuts.
  const std::vector<Cut>& cuts() const {
    return cuts_;
  }

  // Getter for all active cuts.
  std::vector<Cut> active_cuts() const {
    return GetCuts(active_cuts_positions_);
  }

  // Getter for active cut positions.
  const std::vector<unsigned int>& active_cut_positions() const {
    return active_cuts_positions_;
  }

  // Getter for number of currently active cuts.
  const unsigned int& current_number_of_active_cuts() const {
    return current_number_of_active_cuts_;
  }

  // Getter for number of current cuts in the storage.
  const unsigned int& current_number_of_cuts() const {
    return current_number_of_cuts_;
  }

  // Getter for total number of cuts in the storage.
  const unsigned int& total_number_of_cuts_found() const {
    return total_number_of_cuts_found_;
  }

 private:
  std::vector<Cut> cuts_;
  std::vector<unsigned int> active_cuts_positions_;
  unsigned int current_number_of_cuts_;
  unsigned int total_number_of_cuts_found_;
  unsigned int current_number_of_active_cuts_;
};

} // namespace
#endif //SRC_DATA_STRUCTURES_CUTS_DATA_H_
