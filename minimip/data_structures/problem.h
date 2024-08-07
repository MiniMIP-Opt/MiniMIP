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

#ifndef MINIMIP_DATA_STRUCTURES_PROBLEM_H_
#define MINIMIP_DATA_STRUCTURES_PROBLEM_H_

#include <cstdint>
#include <limits>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_format.h"
#include "ortools/base/logging.h"
#include "ortools/base/status_builder.h"
#include "ortools/base/status_macros.h"

namespace minimip {

constexpr double kInf = std::numeric_limits<double>::infinity();

// ==========================================================================
// API Input Datastructures
// ==========================================================================

struct MiniMipVariable;
struct MiniMipConstraint;
struct MiniMipSolutionHint;

struct MiniMipProblem {
  std::string name;
  std::vector<MiniMipVariable> variables;
  std::vector<MiniMipConstraint> constraints;
  std::vector<MiniMipSolutionHint> hints;
  bool is_maximization = false;
  double objective_offset = 0.0;
};

struct MiniMipVariable {
  std::string name;
  double objective_coefficient = 0.0;
  double lower_bound = -kInf;
  double upper_bound = kInf;
  bool is_integer = false;
};

struct MiniMipConstraint {
  std::string name;
  std::vector<int> var_indices;
  std::vector<double> coefficients;
  double left_hand_side = -kInf;
  double right_hand_side = kInf;
};

struct MiniMipSolutionHint {
  std::vector<int> var_indices;
  std::vector<double> values;
};

// ==========================================================================
// API Output Datastructures
// ==========================================================================

enum class MiniMipSolveStatus;
enum class MiniMipStoppingReason;

struct MiniMipSolution {
  // Dense vector with the solution.
  std::vector<double> variable_values;
  double objective_value = kInf;

  // TODO(lpawel):
  // Some extra stuff like max constraints/integrality violations wrt the
  // original model to allow client to quickly see whether the solution is
  // precise enough.
};

struct MiniMipResult {
  MiniMipSolveStatus solve_status;
  MiniMipStoppingReason stopping_reason;
  std::string error_message;

  double best_primal_bound;
  double best_dual_bound;
  MiniMipSolution best_solution = MiniMipSolution();

  std::vector<MiniMipSolution> additional_solutions;

  // TODO(lpawel): Extra stuff like solve stats, wallclock time, number of
  // nodes, etc.

  absl::Status AddSolution(const MiniMipSolution& solution) {
    VLOG(10) << "calling AddSolution().";
    if (solution.objective_value < best_solution.objective_value) {
      if (best_solution.objective_value != kInf) {
        additional_solutions.push_back(best_solution);
      }
      best_solution = solution;
    } else {
      // Store all primal solutions found
      additional_solutions.push_back(solution);
    }
    return absl::OkStatus();
  }
};

enum class MiniMipSolveStatus {
  // The provided input problem was invalid.
  kProblemInvalid = 0,

  // No solution, nor proof of infeasibility or unboundedness was found.
  kNotSolved = 1,

  // The solver found a solution and proven its optimality (within given gap
  // limits)
  kOptimal = 2,

  // The solver found a solution, but has not proven its optimality (whiting
  // given gap limits)
  kFeasible = 3,

  // The solver found a solution (perhaps even proved its "optimality"),
  // but the solution turned out to be imprecise (e.g., it violates primal
  // feasibility tolerance).
  kImprecise = 4,

  // The solver proved the problem is infeasible, i.e., no solution exists.
  kInfeasible = 5,

  // The solver proved the problem is unbounded, i.e., for each solution there
  // exists another solution with a better objective value (ad infinitum).
  // This should happen rarely for real-life problems and may indicate a "bug"
  // in the input problem itself.
  kUnbounded = 5,

  // The solver proved the problem is either infeasible or unbounded. Note, it's
  // sometimes very time-consuming to distinguish between infeasible and
  // unbounded
  // for a mixed-integer program (in contrast to pure linear programming).
  // In such cases, we stop and return this status.
  kInfeasibleOrUnbounded = 6,
};

enum class MiniMipStoppingReason {
  kError = 0,
  kWallclockTimeLimitReached = 1,
  kDeterministicTimeLimitReached = 2,
  kSolutionLimitReached = 3,
  kNodeLimitReached = 4,
  kSimplexIterationLimitReached = 5,
  kPrimalBoundLimitReached = 6,
  kDualBoundLimitReached = 7,
  kAbsoluteGapLimitReached = 8,
  kRelativeGapLimitReached = 9,
  kMemoryLimitReached = 10,
  kUserCallbackLimitReached = 11,
};

// ==========================================================================
// API Functionality
// ==========================================================================
std::string FindErrorInMiniMipProblem(const MiniMipProblem& problem);

}  // namespace minimip
#endif  // MINIMIP_DATA_STRUCTURES_PROBLEM_H_
