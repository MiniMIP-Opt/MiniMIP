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

#ifndef MINIMIP_CUTTING_INTERFACE_CUTS_GENERATOR_H_
#define MINIMIP_CUTTING_INTERFACE_CUTS_GENERATOR_H_

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "minimip/data_structures/cuts_data.h"
#include "ortools/base/logging.h"

namespace minimip {

// Forward declaration of Solver. This is required to break the circular
// dependency between the solver and the generator.
class SolverContextInterface;

// The CutGeneratorInterface provides an interface for generating valid
// inequalities (cuts) from the current Linear Programming (LP) solution.
// It formulates potential cuts based on the current state of the solver
// and the LP solution's characteristics.
class CutGeneratorInterface {
 public:
  virtual ~CutGeneratorInterface() = default;

  // Generate up to `max_num_cuts` cutting planes. Returns an empty vector if no
  // cuts could be generated given the current state. If called multiple times,
  // the generator is expected to continue where it left of, if appropriate.
  virtual absl::StatusOr<std::vector<CutData>> GenerateCuttingPlanes(
      const SolverContextInterface& context) = 0;
};

}  // namespace minimip

#endif  // MINIMIP_CUTTING_INTERFACE_CUTS_GENERATOR_H_
