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

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "absl/status/status.h"
#include "ortools/base/status_macros.h"
#include "unit_tests/utils.h"

namespace minimip {

TEST(SolveLoopTest, RootNode) {
  // Run a single round of main loop
  // Creates a very simple model where all variables are non-negative.
  //
  // max: z = x2
  //  3*x1 + 2*x2 <= 6
  // -3*x1 + 2*x2 <= 0
  // x1, x2 >= 0
  // x1, x2 integer
  MiniMipProblem problem{
      .variables = {MiniMipVariable{.name = "x1",
                                    .objective_coefficient = 0.0,
                                    .lower_bound = 0,
                                    .upper_bound = kInf,
                                    .is_integer = true},
                    MiniMipVariable{.name = "x2",
                                    .objective_coefficient = 1.0,
                                    .lower_bound = 0,
                                    .upper_bound = kInf,
                                    .is_integer = true}},
      .constraints = {MiniMipConstraint{
                          .name = "ct1",
                          .var_indices = {0, 1},
                          .coefficients = {3.0, 2.0},
                          .left_hand_side = -kInf,
                          .right_hand_side = 6.0,
                      },
                      MiniMipConstraint{
                          .name = "ct2",
                          .var_indices = {0, 1},
                          .coefficients = {-3.0, 2.0},
                          .left_hand_side = -kInf,
                          .right_hand_side = 0.0,
                      }},
      .is_maximization = true};

  // Call the Create function to create a Solver object
  ASSERT_OK_AND_ASSIGN(std::unique_ptr<Solver> solver, Solver::Create(problem));

  ASSERT_OK(solver->Solve());
  ASSERT_EQ(solver->result().solve_status, MiniMipSolveStatus::kOptimal);
  ASSERT_FLOAT_EQ(solver->result().best_solution.objective_value, -1.0);
  ASSERT_FLOAT_EQ(solver->result().best_solution.variable_values[0], 1.0);
  ASSERT_FLOAT_EQ(solver->result().best_solution.variable_values[1], 1.0);
}
// TODO(CG): Add Test for multiple rounds of main loop, as cutting solves
//           the rootnode to optimality and branching is not tested anymore.

}  // namespace minimip