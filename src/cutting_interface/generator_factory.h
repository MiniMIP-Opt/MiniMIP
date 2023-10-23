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

#ifndef SRC_CUTTING_INTERFACE_GENERATOR_FACTORY_H_
#define SRC_CUTTING_INTERFACE_GENERATOR_FACTORY_H_

#include "src/cutting_interface/aggregating_generator.h"
#include "src/cutting_interface/cuts_generator.h"

namespace minimip {

inline absl::StatusOr<std::unique_ptr<CutGenerator>>
ConfigureGeneratorFromProto(const GeneratorParameters& generator_parameters) {
  if (generator_parameters.has_tableau_rounding_generator_parameters()) {
    return std::make_unique<TableauRoundingGenerator>(generator_parameters);
  }
  return absl::InvalidArgumentError("No generator-specific parameters set.");
}

}  // namespace minimip

#endif SRC_CUTTING_INTERFACE_GENERATOR_FACTORY_H_
