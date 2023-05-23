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

#ifndef SRC_DATA_STRUCTURES_CUTS_DATA_H_
#define SRC_DATA_STRUCTURES_CUTS_DATA_H_

#include <utility>

#include "src/data_structures/mip_data.h"
#include "src/data_structures/strong_sparse_vector.h"
#include "src/lp_interface/lpi.h"

namespace minimip {

// ============================================================================
// The CutData struct stores all information of a cutting plane d <= a^Tx <= b.
// ============================================================================
struct CutData {
  // The coefficients of the cutting plane are stored as a SparseRow
  SparseRow row;
  double right_hand_side = std::numeric_limits<double>::infinity();
  double left_hand_side = -std::numeric_limits<double>::infinity();

  // If the cut will always be selected and activated. This is necessary for
  // cut selection (e.g., bound changes returned from separators).
  bool is_forced = false;

  // The node at which the cut is added to the storage. A value of zero means
  // the cut was added at the root, otherwise the node of the tree is given.
  int added_at_node = -1;

  // The origin of a cut is specified to the separation round in the node it
  // originated from. Useful for cut selection between rounds.
  int from_separation_round_n = -1;

  // Important characteristics of a cutting plane.
  int number_of_non_zeros = -1;
  int number_of_integer_variables = -1;

  // The relative parallelism to the objective function.
  double objective_parallelism = -std::numeric_limits<double>::infinity();

  // The cut name is set by the separator it originates from.
  std::string name;
};

class Cut {
 public:
  // Constructor that takes a CutData object
  explicit Cut(CutData& data) : data_(data) {
    DCHECK(!data_.row.IsClean());
    DCHECK_NE(data_.right_hand_side, -data_.left_hand_side);
    DCHECK_NE(data_.added_at_node, -1);
    DCHECK_NE(data_.from_separation_round_n, -1);
    DCHECK_NE(data_.number_of_non_zeros, -1);
    DCHECK_NE(data_.number_of_integer_variables, -1);
    DCHECK_NE(data_.objective_parallelism,
              -std::numeric_limits<double>::infinity());
    DCHECK(!data_.name.empty());
  }

  // Setter methods for the changing data
  void setEfficacy(double efficacy) { efficacy_ = efficacy; }
  void setScore(double score) { current_score_ = score; }
  void setIndex(int index) { cut_index_ = index; }
  void setActive(bool is_active) { is_active_ = is_active; }

  // Getter methods for the changing data
  double efficacy() const { return efficacy_; }
  double score() const { return current_score_; }
  int index() const { return cut_index_; }
  bool isActive() const { return is_active_; }

  // Getter method for the constant data
  const CutData& getConstData() const { return data_; }

  // Method to create a new Cut with updated CutData and moved fields.
  Cut recreateWithUpdatedData(CutData& newData) const {
    Cut new_cut(newData);
    new_cut.efficacy_ = efficacy_;
    new_cut.current_score_ = current_score_;
    new_cut.cut_index_ = cut_index_;
    new_cut.is_active_ = is_active_;
    return new_cut;
  }

 private:
  // The CutData struct stores all constant information of a cutting plane d <=
  // a^Tx <= b.
  CutData data_;

  // The index of the cut in the storage, for easy reference.
  // The index is set to the current number of cuts in storage once the cut is
  // added to the storage (i.e. cut_index += cuts().size()).
  int cut_index_ = 0;

  // If the cut is currently applied to the problem, is_active is true.
  bool is_active_ = false;

  // The efficacy of the cut is the orthogonal projection of the LP optimum onto
  // the cutting plane.
  double efficacy_ = 0.0;

  // The current score of the cut is the score the cut was given when first
  // selected.
  double current_score_ = 0.0;
};

// ============================================================================
// CutStorage contains the generated cutting planes and all relevant meta-data
// for cutting plane management(e.g., it tracks which cuts are currently
// active).
//
// NOTE: As of 2022/11/09, MiniMip uses only globally valid cuts
// (though, a cutting plane might have been generated in the inner nodes).
// ============================================================================

class CutStorage {
 public:
  // ==========================================================================
  // Constructors
  // ==========================================================================

  CutStorage();

  // Initialize CutStorage from initial separation round.
  CutStorage(std::vector<Cut> cuts, std::vector<int> cut_indices);

  // CutStorage is not copyable to make sure a copy will not be
  // triggered by accident (copy constructor and assign operator are private).
  // CutStorage is (no-throw) movable.
  CutStorage(CutStorage&&) noexcept = default;
  CutStorage& operator=(CutStorage&&) noexcept = default;

  // Use this to initialize by deep copy from another cut storage.
  // Under-the-hood we just use a private copy / move constructor and
  // assignment operator.
  void PopulateFromCutStorage(CutStorage cut_storage);

  // ==========================================================================
  // Methods for managing the cuts in the cut storage.
  // ==========================================================================

  // Add a cut to storage.
  int AddCut(Cut&& cut);

  // Activate stored cut.
  void ActivateCut(int cut_index);

  // Remove a single cut from storage.
  void RemoveCut(const int& cut_index);

  // Getter for individual cuts.
  const Cut& GetCut(const int& cut_index) const;

  // Getter for all cuts.
  const std::vector<Cut>& cuts() const { return cuts_; }

  // Getter for active cut indices.
  const std::vector<int>& active_cuts() const { return active_cut_indices_; }

  // Getter for total number of cuts added to the storage while solving.
  const int& total_number_of_cuts_found() const {
    return total_number_of_cuts_found_;
  }

 private:
  std::vector<Cut> cuts_;
  std::vector<int> active_cut_indices_;
  int total_number_of_cuts_found_;
};

}  // namespace minimip
#endif  // SRC_DATA_STRUCTURES_CUTS_DATA_H_
