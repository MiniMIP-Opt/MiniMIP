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

class CutRegistry;  // Forward declaration

// ============================================================================
// The CutData struct stores all information of a cutting plane a^Tx <= b.
// ============================================================================
class CutData {
  // Allow CutRegistry to create CutData objects
  friend class CutRegistry;

 public:
  // The efficacy is a measure of how much the cut improves the LP solution.
  // It describes the length of orthogonal projection of the LP solution onto
  // the cutting plane.
  void SetEfficacy(double efficacy) { efficacy_ = efficacy; }
  double efficacy() const { return efficacy_; }

  // The score is a measure of preference that is set by the cut selector.
  void SetScore(double score) { score_ = score; }
  double score() const {
    if (score_.has_value()) {
      return *score_;
    }
    throw std::runtime_error("Score has not been set.");
  }

  // The index is the position of the cut in the cut registry.
  // It is set by the Cut registry, after the cut has been added.
  void SetIndex(int index) { cut_index_ = index; }
  int index() const {
    if (cut_index_.has_value()) {
      return *cut_index_;
    }
    throw std::runtime_error("Index has not been set.");
  }

  // The active flag is set by the cut selector. It is used to indicate whether
  // the cut is active in the LP. This is useful to avoid adding the same or
  // strictly dominated cuts multiple times.
  void SetActive(bool is_active) { is_active_ = is_active; }
  bool is_active() const { return is_active_; }

  // The row is the coefficient vector a of the cutting plane a^Tx <= b.
  const SparseRow& row() const { return row_; }

  // The right hand side is the constant b of the cutting plane a^Tx <= b.
  double right_hand_side() const { return right_hand_side_; }

  // Important characteristics of a cutting plane.
  int number_of_non_zeros() const { return number_of_non_zeros_; }
  int number_of_integer_variables() const {
    return number_of_integer_variables_;
  }
  double objective_parallelism() const { return objective_parallelism_; }

  // The cut name is set by the generator it originates from.
  const std::string& name() const { return name_; }

  // If the cut will always be selected and activated. This is useful for single
  // variable cuts, also called bound changes.
  bool is_forced() const { return is_forced_; }

 private:
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
  // The coefficients of the cutting plane are stored as a SparseRow
  SparseRow row_;
  double right_hand_side_;

  // Important characteristics of a cutting plane.
  int number_of_non_zeros_;
  int number_of_integer_variables_;

  // The relative parallelism to the objective function.
  double objective_parallelism_;

  // The cut name is set by the generator it originates from.
  std::string name_;

  // If the cut will always be selected and activated. This is necessary for
  // cut selection (e.g., bound changes returned from generators).
  bool is_forced_;

  // If the cut is currently applied to the problem, is_active is true.
  bool is_active_ = false;

  // The efficacy of the cut is the orthogonal projection of the LP optimum onto
  // the cutting plane.
  double efficacy_;

  // The index of the cut in the registry, for easy reference.
  // The index is set to the current number of cuts in registry once the cut is
  // added to the registry (i.e. cut_index += cuts().size()).
  std::optional<int> cut_index_;

  // The original score is the score the cut was given when first selected.
  std::optional<double> score_;
};

// ============================================================================
// CutRegistry initializes and contains the generated cutting planes and all
// relevant meta-data for cutting plane management (e.g., it tracks which cuts
// are currently active).
//
// NOTE: As of 2022/11/09, MiniMip uses only globally valid cuts
// (though, a cutting plane might have been generated in the inner nodes).
// ============================================================================

class CutRegistry {
 public:
  // ==========================================================================
  // Constructors
  // ==========================================================================

  CutRegistry();

  // Initialize CutRegistry from initial separation round.
  CutRegistry(std::vector<CutData> cuts, std::vector<int> cut_indices);

  // CutRegistry is not copyable to make sure a copy will not be
  // triggered by accident (copy constructor and assign operator are private).
  // CutRegistry is (no-throw) movable.
  CutRegistry(CutRegistry&&) noexcept = default;
  CutRegistry& operator=(CutRegistry&&) noexcept = default;

  // Use this to initialize by deep copy from another cut registry.
  // Under-the-hood we just use a private copy / move constructor and
  // assignment operator.
  void PopulateFromCutRegistry(CutRegistry cut_registry);

  // ==========================================================================
  // Methods for managing the cuts in the cut registry.
  // ==========================================================================

  CutData CreateCut(const MipData& mip_data, const SparseRow& lp_optimum,
                    SparseRow row, double right_hand_side, std::string name,
                    bool is_forced = false) const;

  // Activate cutting plane with given index in the registry
  void ActivateCut(int cut_index);

  // Remove a single cut from registry.
  void RemoveCut(const int& cut_index);

  // Getter for individual cuts.
  const CutData& GetCut(const int& cut_index) const;

  double ComputeEfficacy(const SparseRow& row, double right_hand_side,
                         const SparseRow& lp_optimum) const {
    return (row.DotProduct(lp_optimum) - right_hand_side) /
           sqrt(row.DotProduct(row));
  }

  // Given a new LP optimum, update the efficacy of all cuts before reselection.
  void UpdateCutEfficacy(const SparseRow& lp_optimum);

  // Getter for all cuts.
  const std::vector<CutData>& cuts() const { return cuts_; }

  // Getter for active cut indices.
  const std::vector<int>& active_cuts() const { return active_cut_indices_; }

  // Getter for total number of cuts added to the registry while solving.
  const int& total_number_of_cuts_found() const {
    return total_number_of_cuts_found_;
  }

  // Add a cut to registry and return its index.
  int AddCut(CutData&& cut_data);

 private:
  std::vector<CutData> cuts_;
  std::vector<int> active_cut_indices_;
  int total_number_of_cuts_found_;
};

}  // namespace minimip
#endif  // SRC_DATA_STRUCTURES_CUTS_DATA_H_
