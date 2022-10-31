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

#ifndef SRC_CUTTING_INTERFACE_SELECTOR_H_
#define SRC_CUTTING_INTERFACE_SELECTOR_H_

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "ortools/base/logging.h"
#include "ortools/base/status_macros.h"
#include "src/data_structures/cuts_data.h"
#include "src/data_structures/mip_data.h"
#include "src/lp_interface/lpi.h"

namespace minimip {

// NOTE: This file should include the default scoring function as the efficacy
//      of a given cutting plane in regard to the current LP solution as a
//      general baseline for scoring cuts. Additionally, a basic greedy cut
//      selection method adding all non-duplicate cuts given to the cut
//      selector should be implemented here, allowing the user to turn off all
//      cut selection strategies for testing or analytics.

class SelectorInterface {
 public:
  // Generate up to `max_num_cuts` cutting planes. Returns an empty vector if no
  // cuts could be generated given the current state. If called multiple times,
  // the generator is expected to continue where it left of, if appropriate.
  virtual absl::StatusOr<int> SelectCuttingPlanes(
      const LPInterface* lpi, const MipData& mip_data,
      std::vector<CuttingPlane>& cuts) = 0;
};

class Selector : SelectorInterface {
 public:
  virtual ~Selector() = default;

  absl::StatusOr<int> SelectCuttingPlanes(
      const LPInterface* lpi, const MipData& mip_data,
      std::vector<CuttingPlane>& cuts) final {
    RETURN_IF_ERROR(PrepareSelection(lpi, mip_data, cuts));

    while (number_of_cuts_selected_ < max_number_of_cuts_to_select_) {
      RETURN_IF_ERROR(PrepareIteration(cuts));

      RETURN_IF_ERROR(Filtering(cuts));

      RETURN_IF_ERROR(Scoring(cuts));
      number_of_cuts_selected_++;
    }

    return number_of_cuts_selected_;
  }

 protected:
  // Extract the relevant data needed for the following cut selection loop.
  // This must include some sort of pre-scoring notion, refreshing the current
  // score of a cut.
  // TODO: add "isCutFresh()" like function corresponding to its current_score.
  virtual absl::Status PrepareSelection(
      const LPInterface* lpi, const MipData& mip_data,
      const std::vector<CuttingPlane>& cuts) const = 0;

  // Prepare the next iteration of filtering cuts and any preprocessing needed
  virtual absl::Status PrepareIteration(std::vector<CuttingPlane>& cuts) = 0;

  // Computing the cutting plane from the current iterative and any additional
  // data needed.
  virtual absl::Status Filtering(std::vector<CuttingPlane>& cuts) = 0;

  // Compute new scores after filtering if necessary.
  virtual absl::Status Scoring(std::vector<CuttingPlane>& cuts) = 0;

 private:
  double max_parallelism_between_cuts_ = 1.0;
  double min_selection_threshold_ = 0.0;
  int number_of_cuts_selected_ = 0;
  int max_number_of_cuts_to_select_ = 0;
  int current_number_of_cuts_ = 0;
  char filtermode_ = 'd';
  char scoringmode_ = 'd';
};

}  // namespace minimip
#endif  // SRC_CUTTING_INTERFACE_SELECTOR_H_
