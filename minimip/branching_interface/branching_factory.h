// Copyright (2024) the MiniMIP Project
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

#ifndef MINIMIP_BRANCHING_INTERFACE_BRANCHING_FACTORY_H_
#define MINIMIP_BRANCHING_INTERFACE_BRANCHING_FACTORY_H_

#include "minimip/branching_interface/branching_interface.h"
#include "minimip/branching_interface/maxfractional_branching.h"
#include "minimip/branching_interface/random_branching.h"
#include "minimip/parameters.pb.h"

namespace minimip {

inline absl::StatusOr<std::unique_ptr<BranchingRuleInterface>>
CreateBranchingRule(const BranchingParameters& branching_parameters) {
  if (branching_parameters.has_random_branching_parameters()) {
    return std::make_unique<RandomBranching>(branching_parameters);
  }
  if (branching_parameters.has_max_fractional_branching_parameters()) {
    return std::make_unique<MaxFractionalBranching>(branching_parameters);
  }
  return absl::InvalidArgumentError("No branching-specific parameters set.");
}

}  // namespace minimip
#endif  // MINIMIP_BRANCHING_INTERFACE_BRANCHING_FACTORY_H_
