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
class CutData {
 public:
  // Constructor that takes a CutData object
  explicit CutData(SparseRow row, double right_hand_side,
                   int number_of_non_zeros, int number_of_integer_variables,
                   double objective_parallelism, double efficacy,
                   std::string name, bool is_forced = false)
      : row_(std::move(row)),
        right_hand_side_(right_hand_side),
        number_of_non_zeros_(number_of_non_zeros),
        number_of_integer_variables_(number_of_integer_variables),
        objective_parallelism_(objective_parallelism),
        name_(std::move(name)),
        is_forced_(is_forced),
        efficacy_(efficacy) {
    // Check the validity of the parameters
    DCHECK(row_.IsClean());
    DCHECK_GT(right_hand_side_, -std::numeric_limits<double>::infinity());
    DCHECK_GT(number_of_non_zeros_, 0);
    DCHECK_GE(number_of_integer_variables_, 0);
    DCHECK_GE(objective_parallelism_, -1.0);
    DCHECK_LE(objective_parallelism_, 1.0);
    DCHECK_GE(efficacy_, 0.0);
    DCHECK(!name_.empty());
  }

  // Setter methods for the changing data
  void SetEfficacy(double efficacy) { efficacy_ = efficacy; }
  void SetScore(double score) { score_ = score; }
  void SetIndex(int index) { cut_index_ = index; }
  void SetActive(bool is_active) { is_active_ = is_active; }

  // Getter methods for the changing data
  double score() const {
    if (score_.has_value()) {
      return *score_;
    }
    throw std::runtime_error("Score has not been set.");
  }

  int index() const {
    if (cut_index_.has_value()) {
      return *cut_index_;
    }
    throw std::runtime_error("Index has not been set.");
  }

  bool is_active() const { return is_active_; }

  double efficacy() const { return efficacy_; }

  // Getter methods for the constant data
  const SparseRow& row() const { return row_; }
  double right_hand_side() const { return right_hand_side_; }
  int number_of_non_zeros() const { return number_of_non_zeros_; }
  int number_of_integer_variables() const {
    return number_of_integer_variables_;
  }
  double objective_parallelism() const { return objective_parallelism_; }
  const std::string& name() const { return name_; }
  bool is_forced() const { return is_forced_; }

 private:
  // The coefficients of the cutting plane are stored as a SparseRow
  SparseRow row_;
  double right_hand_side_;

  // Important characteristics of a cutting plane.
  int number_of_non_zeros_;
  int number_of_integer_variables_;

  // The relative parallelism to the objective function.
  double objective_parallelism_;

  // The cut name is set by the separator it originates from.
  std::string name_;

  // If the cut will always be selected and activated. This is necessary for
  // cut selection (e.g., bound changes returned from separators).
  bool is_forced_;

  // If the cut is currently applied to the problem, is_active is true.
  bool is_active_{false};

  // The efficacy of the cut is the orthogonal projection of the LP optimum onto
  // the cutting plane.
  double efficacy_;

  // The index of the cut in the storage, for easy reference.
  // The index is set to the current number of cuts in storage once the cut is
  // added to the storage (i.e. cut_index += cuts().size()).
  std::optional<int> cut_index_;

  // The original score is the score the cut was given when first selected.
  std::optional<double> score_;
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
  CutStorage(std::vector<CutData> cuts, std::vector<int> cut_indices);

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
  int AddCut(CutData&& cut_data);

  // Activate stored cut.
  void ActivateCut(int cut_index);

  // Remove a single cut from storage.
  void RemoveCut(const int& cut_index);

  // Getter for individual cuts.
  const CutData& GetCut(const int& cut_index) const;

  // Getter for all cuts.
  const std::vector<CutData>& cuts() const { return cuts_; }

  // Getter for active cut indices.
  const std::vector<int>& active_cuts() const { return active_cut_indices_; }

  // Getter for total number of cuts added to the storage while solving.
  const int& total_number_of_cuts_found() const {
    return total_number_of_cuts_found_;
  }

 private:
  std::vector<CutData> cuts_;
  std::vector<int> active_cut_indices_;
  int total_number_of_cuts_found_;
};

}  // namespace minimip
#endif  // SRC_DATA_STRUCTURES_CUTS_DATA_H_
