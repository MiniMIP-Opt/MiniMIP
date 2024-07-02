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

#ifndef MINIMIP_SRC_CUTTING_INTERFACE_DEFAULT_RUNNER_H_
#define MINIMIP_SRC_CUTTING_INTERFACE_DEFAULT_RUNNER_H_

#include <vector>

#include "ortools/base/status_macros.h"
#include "src/cutting_interface/cuts_runner.h"
#include "src/cutting_interface/generator_factory.h"
#include "src/cutting_interface/selector_factory.h"
#include "src/parameters.pb.h"

namespace minimip {

// Forward declaration of Solver. This is required to break the circular
// dependency between the solver and the selector.
class SolverContextInterface;

class DefaultRunner : public CutRunnerInterface {
 public:
  explicit DefaultRunner(CutRunnerParameters params)
      : params_(std::move(params.default_runner_parameters())) {
    VLOG(10) << "calling DefaultRunner().";
    // Check and initialize the cut generators.
    for (const auto& generator_params : params.generator_parameters()) {
      auto generator_or_status = CreateCutGenerator(generator_params);
      if (!generator_or_status.ok()) {
        LOG(ERROR) << "Could not configure generator: "
                   << generator_or_status.status();
        continue;
      }
      AddGenerator(std::move(generator_or_status.value()));
    }

    auto selector_or_status =
        CreateCutSelector(params.selector_parameters().Get(0));
    if (!selector_or_status.ok()) {
      LOG(ERROR) << "Could not configure selector: "
                 << selector_or_status.status();
    }
    AddSelector(std::move(selector_or_status.value()));
  }

  bool MayRunOneMoreSeperationRound(const SolverContextInterface& context) final;

  absl::Status SeparateCurrentLPSolution(SolverContextInterface& context) final;

 private:
  const DefaultRunnerParameters params_;
};

}  // namespace minimip

#endif  // MINIMIP_SRC_CUTTING_INTERFACE_DEFAULT_RUNNER_H_
