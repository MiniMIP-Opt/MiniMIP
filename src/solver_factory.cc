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

#include "src/solver_factory.h"

#include "src/cutting_interface/runner_factory.h"
#include "src/cutting_interface/selector_factory.h"
#include "src/cutting_interface/separator_factory.h"
#include "src/lp_interface/lpi_factory.h"

namespace minimip {

absl::StatusOr<std::unique_ptr<MiniMipSolver>> ConfigureMiniMipSolverFromProto(
    const MiniMipParameters& params, const MiniMipProblem& problem) {
  const std::string problem_error = FindErrorInMiniMipProblem(problem);
  if (!problem_error.empty()) {
    return util::InvalidArgumentErrorBuilder()
           << "Error found in problem: " << problem_error;
  }
  auto mip_data = std::make_unique<MipData>(problem);
  auto cut_storage = std::make_unique<CutStorage>();
  ASSIGN_OR_RETURN(std::unique_ptr<LPInterface> lpi,
                   ConfigureLPSolverFromProto(params.lp_parameters()));
  ASSIGN_OR_RETURN(std::unique_ptr<CuttingInterface> cut_runner,
                   ConfigureCutInterfaceFromProto(params.cut_runner()));
  for (const SeparatorParameters& separator_params : params.separators()) {
    ASSIGN_OR_RETURN(std::unique_ptr<Separator> separator,
                     ConfigureSeparatorFromProto(separator_params));
    cut_runner->AddSeparator(std::move(separator));
  }
  for (const SelectorParameters& selector_params : params.selectors()) {
    ASSIGN_OR_RETURN(std::unique_ptr<Selector> selector,
                     ConfigureSelectorFromProto(selector_params));
    cut_runner->AddSelector(std::move(selector));
  }
  // TODO(cgraczy): Configure the selector(s) in a similar way.
  auto solver = std::make_unique<MiniMipSolver>(
      params, std::move(mip_data), std::move(cut_storage), std::move(lpi),
      std::move(cut_runner));
  return solver;
}

}  // namespace minimip
