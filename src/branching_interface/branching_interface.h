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

#ifndef SRC_BRANCHING_INTERFACE_BRANCHING_INTERFACE_H_
#define SRC_BRANCHING_INTERFACE_BRANCHING_INTERFACE_H_

#include "absl/status/statusor.h"
#include "src/data_structures/mip_types.h"
#include "src/parameters.pb.h"

namespace minimip {

// Forward declaration of Solver. This is required to break the circular
// dependency between the solver and the cutting interface.
class SolverContextInterface;

class BranchingInterface {
 public:
  virtual ~BranchingInterface() = default;

  virtual const absl::StatusOr<ColIndex> NextBranchingVariable(
      const SolverContextInterface& context) const = 0;
};

}  // namespace minimip
#endif  // SRC_BRANCHING_INTERFACE_BRANCHING_INTERFACE_H_
