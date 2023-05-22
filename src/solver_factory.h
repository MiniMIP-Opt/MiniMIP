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

#include "absl/status/statusor.h"
#include "src/data_structures/problem.h"
#include "src/solver.h"

namespace minimip {

// Factory function that configures
absl::StatusOr<std::unique_ptr<MiniMipSolver>> ConfigureMiniMipSolverFromProto(
    const MiniMipParameters& params, const MiniMipProblem& problem);

}  // namespace minimip
