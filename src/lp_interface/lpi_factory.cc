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

#include "src/lp_interface/lpi_factory.h"

#include <memory>

#include "lpi_soplex.h"
#include "ortools/base/logging.h"
#include "src/lp_interface/lpi_glop.h"
#include "src/lp_interface/lpi_soplex.h"

namespace minimip {

absl::StatusOr<std::unique_ptr<LPInterface>> ConfigureLPSolverFromProto(
    const LPParameters& params) {
  std::unique_ptr<LPInterface> lpi;
  switch (params.lp_solver_type()) {
    case LPParameters::LP_GLOP:
      lpi = std::make_unique<LPGlopInterface>();
      break;
    case LPParameters::LP_SOPLEX:
      lpi = std::make_unique<LPSoplexInterface>();
      break;
    default:
      return util::InvalidArgumentErrorBuilder()
             << "Invalid lp solver type: " << params.lp_solver_type();
  }

  RETURN_IF_ERROR(
      lpi->SetIntegerParameter(LPParameter::kScaling, params.scaling()));
  RETURN_IF_ERROR(lpi->SetIntegerParameter(LPParameter::kPresolving,
                                           params.use_presolve()));
  RETURN_IF_ERROR(lpi->SetIntegerParameter(LPParameter::kPolishing,
                                           params.use_polishing()));
  RETURN_IF_ERROR(
      lpi->SetIntegerParameter(LPParameter::kPricing, params.pricing()));
  if (params.feasibility_tolerance() >= 0) {
    RETURN_IF_ERROR(lpi->SetRealParameter(LPParameter::kFeasibilityTolerance,
                                          params.feasibility_tolerance()));
  }
  if (params.dual_feasibility_tolerance() >= 0) {
    RETURN_IF_ERROR(
        lpi->SetRealParameter(LPParameter::kDualFeasibilityTolerance,
                              params.dual_feasibility_tolerance()));
  }
  if (params.markowitz_tolerance() >= 0) {
    RETURN_IF_ERROR(lpi->SetRealParameter(LPParameter::kMarkowitz,
                                          params.markowitz_tolerance()));
  }
  RETURN_IF_ERROR(lpi->SetIntegerParameter(LPParameter::kRefactor,
                                           params.refactorization_interval()));
  if (params.use_objective_limit()) {
    RETURN_IF_ERROR(lpi->SetRealParameter(LPParameter::kObjectiveLimit,
                                          params.objective_limit()));
  }
  if (params.time_limit() > 0.0) {
    RETURN_IF_ERROR(
        lpi->SetRealParameter(LPParameter::kLPTimeLimit, params.time_limit()));
  }
  if (params.iteration_limit() > 0) {
    RETURN_IF_ERROR(lpi->SetIntegerParameter(LPParameter::kLPIterationLimit,
                                             params.iteration_limit()));
  }
  RETURN_IF_ERROR(
      lpi->SetIntegerParameter(LPParameter::kTiming, params.timing()));
  RETURN_IF_ERROR(
      lpi->SetIntegerParameter(LPParameter::kThreads, params.num_threads()));
  if (params.random_seed() != 0) {
    RETURN_IF_ERROR(lpi->SetIntegerParameter(LPParameter::kRandomSeed,
                                             params.random_seed()));
  }
  RETURN_IF_ERROR(lpi->SetIntegerParameter(
      LPParameter::kLPInfo, params.enable_internal_solver_output()));

  return lpi;
}

}  // namespace minimip
