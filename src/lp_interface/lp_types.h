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

#ifndef SRC_LP_INTERFACE_LP_TYPES_H_
#define SRC_LP_INTERFACE_LP_TYPES_H_

#include <cstddef>
#include <string>
#include <vector>

namespace minimip {

// LP objective sense
enum class LPObjectiveSense {
  kMaximization = -1,  // maximize objective function
  kMinimization = +1,  // minimize objective function
};

// LP solver parameters
enum class LPParameter {
  kFromScratch          = 0,  // Solver should start from scratch at next call?
  kScaling              = 1,  // Should LP solver use scaling?
  kPresolving           = 2,  // Should LP solver use presolving?
  kPolishing            = 3,  // set solution polishing (0: disable, 1: enable)
  kPricing              = 4,  // pricing strategy
  kFeasibilityTolerance = 5,  // feasibility tolerance for primal variables and
                              // slacks, strictly positive
  kDualFeasibilityTolerance = 6,  // feasibility tolerance for dual variables
                                  // and reduced costs, strictly positive
  kMarkowitz      = 7,            // Markowitz tolerance
  kRefactor       = 8,            // set refactorization interval (0: automatic)
  kObjectiveLimit = 9,     // objective limit (stop if objective is known be
                           // larger/smaller than limit for min/max-imization)
  kLPIterationLimit = 10,  // LP iteration limit (> 0)
  kLPTimeLimit      = 11,  // LP time limit, positive
  kThreads    = 12,  // number of threads used to solve the LP (-1: automatic)
  kTiming     = 13,  // type of timer (1: cpu, 2: wallclock, 0: off)
  kRandomSeed = 14,  // inital random seed, e.g., for perturbations in the
                     // simplex (0: LP default)
  kLPInfo = 15,      // Should LP solver output information to the screen?
};

// LP pricing strategy
enum class LPPricing {
  kDefault = 0,  // the MiniMIP/LP interface should use its preferred strategy
  kAuto    = 1,  // the LP solver should use its preferred strategy
  kFull    = 2,  // full pricing
  kPartial = 3,  // partial pricing
  kSteep   = 4,  // steepest edge pricing
  kSteepQStart = 5,  // steepest edge pricing without initial dual norms
  kDevex       = 6,  // devex pricing
};

// LP basis status for columns and rows
enum class LPBasisStatus {
  kBasic        = 0,  // (slack) variable is basic
  kAtLowerBound = 1,  // (slack) variable is at its lower bound
  kAtUpperBound = 2,  // (slack) variable is at its upper bound
  kFixed        = 3,  // variable fixed to identical bounds
  kFree         = 4,  // free variable is non-basic and set to zero
};

}  // namespace minimip
#endif  // SRC_LP_INTERFACE_LP_TYPES_H_
