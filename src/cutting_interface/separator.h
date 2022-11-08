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

#ifndef SRC_CUTTING_INTERFACE_SEPARATOR_H_
#define SRC_CUTTING_INTERFACE_SEPARATOR_H_

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "ortools/base/logging.h"
#include "src/parameters.pb.h"

namespace minimip {

// Forward declaration of MiniMipSolver. This is required to break the circular
// dependency between the solver and the separator.
class MiniMipSolver;

class Separator {
 public:
  // Generate up to `max_num_cuts` cutting planes. Returns an empty vector if no
  // cuts could be generated given the current state. If called multiple times,
  // the generator is expected to continue where it left of, if appropriate.
  virtual absl::StatusOr<std::vector<CuttingPlane>> GenerateCuttingPlanes(
      const MiniMipSolver& solver) = 0;
};

// TODO(happlegr): This function currently doesn't do anything, since we don't
// have any separators yet. Rewrite when the first separator is integrated.
inline absl::StatusOr<std::unique_ptr<Separator>> ConfigureSeparatorFromProto(
    const SeparatorParameters& separator_parameters) {
  if (separator_parameters.has_awesome_separator_parameters()) {
    // return std::make_unique<AwesomeSeparator>(separator_parameters);
  } else if (separator_parameters.has_parameterless_separator_parameters()) {
    // return std::make_unique<ParameterlessSeparator>(separator_parameters);
  }
  return absl::InvalidArgumentError("No separator-specific parameters set.");
}

}  // namespace minimip

#endif  // SRC_CUTTING_INTERFACE_SEPARATOR_H_
