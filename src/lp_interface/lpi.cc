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

#include "src/lp_interface/lpi.h"

#include "src/minimip/minimip_def.h"

namespace minimip {

absl::Status LPInterface::SetObjectiveCoefficients(
    const std::vector<int>&
        indices,  // column indices to change objective value for
    const std::vector<double>&
        objective_coefficients  // new objective values for columns
) {
  for (size_t idx = 0; idx < indices.size(); ++idx) {
    MINIMIP_CALL(this->SetObjectiveCoefficient(indices[idx],
                                               objective_coefficients[idx]));
  }
  return absl::OkStatus();
}

}  // namespace minimip
