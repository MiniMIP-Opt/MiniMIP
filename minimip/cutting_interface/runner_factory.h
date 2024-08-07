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

#ifndef MINIMIP_CUTTING_INTERFACE_RUNNER_FACTORY_H_
#define MINIMIP_CUTTING_INTERFACE_RUNNER_FACTORY_H_

#include "minimip/cutting_interface/cuts_runner.h"
#include "minimip/cutting_interface/default_runner.h"
#include "minimip/parameters.pb.h"

namespace minimip {

inline absl::StatusOr<std::unique_ptr<CutRunnerInterface>> CreateCutRunner(
    const CutRunnerParameters& runner_parameters) {
  VLOG(10) << "calling CreateCutRunner().";
  if (runner_parameters.has_default_runner_parameters()) {
    return std::make_unique<DefaultRunner>(runner_parameters);
  }
  return absl::InvalidArgumentError("No cut runner implementation was chosen.");
}

}  // namespace minimip

#endif  // MINIMIP_CUTTING_INTERFACE_RUNNER_FACTORY_H_
