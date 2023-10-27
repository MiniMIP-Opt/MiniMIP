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

#ifndef SRC_CUTTING_INTERFACE_RUNNER_FACTORY_H_
#define SRC_CUTTING_INTERFACE_RUNNER_FACTORY_H_

#include "src/cutting_interface/cuts_runner.h"
#include "src/cutting_interface/default_runner.h"
#include "src/parameters.pb.h"

namespace minimip {

inline absl::StatusOr<std::unique_ptr<CutRunnerInterface>>
ConfigureRunnerFromProto(const CutRunnerParameters& runner_parameters) {
  if (runner_parameters.has_default_runner_parameters()) {
    // Here we create and return the specific CutRunnerInterface.
    return std::make_unique<DefaultRunner>(runner_parameters);
  }
  return absl::InvalidArgumentError("No cut runner implementation was chosen.");
}

}  // namespace minimip

#endif  // SRC_CUTTING_INTERFACE_RUNNER_FACTORY_H_
