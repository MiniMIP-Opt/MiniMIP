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

#include "src/lp_interface/lpi_factory.h"

#include <memory>

#include "ortools/base/logging.h"
#include "src/lp_interface/lpi_glop.h"
#include "src/lp_interface/lpi_soplex.h"
#include "src/parameters.pb.h"

namespace minimip {

absl::StatusOr<std::unique_ptr<LpInterface>> CreateLpSolver(
    const LpParameters& params) {
  std::unique_ptr<LpInterface> lpi;
  switch (params.lp_solver_type()) {
    case LpParameters::LP_GLOP:
      lpi = std::make_unique<LpGlopInterface>();
      break;
    case LpParameters::LP_SOPLEX:
      lpi = std::make_unique<LpSoplexInterface>();
      break;
    default:
      return util::InvalidArgumentErrorBuilder()
             << "Invalid lp solver type: " << params.lp_solver_type();
  }

  RETURN_IF_ERROR(lpi->SetLpParameters(params));
  return lpi;
}

}  // namespace minimip
