// Copyright (2024) the MiniMIP Project
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

#include "minimip/solver.h"

#include <gtest/gtest.h>

#include "absl/status/status.h"
#include "ortools/base/status_macros.h"
#include "unit_tests/utils.h"

namespace minimip {
namespace {
// max: z = x2
//  3*x1 + 2*x2 <= 6
// -3*x1 + 2*x2 <= 0
// x1, x2 >= 0
// x1, x2 integer
MiniMipProblem CreateProblem() {
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
      .is_maximization = false};
  return problem;
}

TEST(SolverTest, CreateWithValidInput) {
  MiniMipProblem problem = CreateProblem();

  auto solver_result = Solver::Create(problem);
  ASSERT_TRUE(solver_result.ok());
  EXPECT_NE(solver_result.value(), nullptr);
}

TEST(SolverTest, InitializeSolver) {
  MiniMipProblem problem = CreateProblem();

  ASSERT_OK_AND_ASSIGN(std::unique_ptr<Solver> solver, Solver::Create(problem));

  MipData& mip_data = solver->mutable_mip_data();
  MipData mip_data_compare = MipData(problem);

  ASSERT_TRUE(AreEqual(mip_data.objective(), mip_data_compare.objective()));

  LpInterface* lp = solver->mutable_lpi();

  std::unique_ptr<LpInterface> lp_compare =
      CreateLpSolver(lp->GetLpParameters()).value();
  ASSERT_OK(lp_compare->PopulateFromMipData(mip_data));

  for (ColIndex col : mip_data.objective().indices()) {
    ASSERT_EQ(lp->GetObjectiveCoefficient(col),
              lp_compare->GetObjectiveCoefficient(col));
  }

  CHECK_OK(lp->SolveLpWithDualSimplex());
  CHECK_OK(lp_compare->SolveLpWithDualSimplex());

  ASSERT_EQ(lp->IsOptimal(), lp_compare->IsOptimal());
  ASSERT_EQ(lp->GetObjectiveValue(), lp_compare->GetObjectiveValue());

  for (ColIndex col : mip_data.objective().indices()) {
    ASSERT_EQ(lp->GetObjectiveCoefficient(col),
              lp_compare->GetObjectiveCoefficient(col));
  }

  for (ColIndex col : mip_data.integer_variables()) {
    ASSERT_EQ(lp->GetLowerBound(col), lp_compare->GetLowerBound(col));
    ASSERT_EQ(lp->GetUpperBound(col), lp_compare->GetUpperBound(col));
  }
}
}  // namespace

}  // namespace minimip