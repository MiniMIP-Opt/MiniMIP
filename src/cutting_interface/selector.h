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

namespace minimip {

class CutSelector {
  // TODO: Implement abstract cut selector class to allow the implementation of
  //      different cut selection strategies including scoring functions.

  // NOTE: This file should include the default scoring function as the efficacy
  //      of a given cutting plane in regard to the current LP solution as a
  //      general baseline for scoring cuts. Additionally, a basic greedy cut
  //      selection method adding all non-duplicate cuts given to the cut
  //      selector should be implemented here, allowing the user to turn off all
  //      cut selection strategies for testing or analytics.

 public:
  virtual ~CutSelector() = default;

  virtual absl::Status MyCutSelectorFunction() = 0;
};

}  // namespace minimip
#endif  // SRC_CUTTING_INTERFACE_SELECTOR_H_
