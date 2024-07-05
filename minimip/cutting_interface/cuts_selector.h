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

#ifndef MINIMIP_CUTTING_INTERFACE_CUTS_SELECTOR_H_
#define MINIMIP_CUTTING_INTERFACE_CUTS_SELECTOR_H_

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "minimip/data_structures/cuts_data.h"
#include "ortools/base/logging.h"
#include "ortools/base/status_macros.h"

namespace minimip {

// Forward declaration of Solver. This is required to break the circular
// dependency between the solver and the selector.
class SolverContextInterface;

// NOTE: This file should include the default scoring function as the efficacy
//      of a given cutting plane in regard to the current LP solution as a
//      general baseline for scoring cuts. Additionally, a basic greedy cut
//      selection method adding all non-duplicate cuts given to the cut
//      selector should be implemented here, allowing the user to turn off all
//      cut selection strategies for testing or analytics.

// The CutSelectorInterface provides the interface to implement cut selectors
// to select the most effective cuts to be included in the LP model. It allows
// evaluating and ranking of different cuts based on predefined criteria,
// such that the most impactful cuts are chosen with regard to this metric.
class CutSelectorInterface {
 public:
  virtual ~CutSelectorInterface() = default;

  // Select up to `max_num_cuts` cutting planes.
  // Returns the selected cuts given the current state.
  virtual absl::StatusOr<std::vector<CutData>> SelectCuttingPlanes(
      const SolverContextInterface& context, std::vector<CutData>& cuts) = 0;
};

}  // namespace minimip

#endif  // MINIMIP_CUTTING_INTERFACE_CUTS_SELECTOR_H_
