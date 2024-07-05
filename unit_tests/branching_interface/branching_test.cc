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

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "minimip/branching_interface/branching_factory.h"
#include "minimip/branching_interface/branching_interface.h"
#include "minimip/branching_interface/maxfractional_branching.h"
#include "minimip/branching_interface/random_branching.h"
#include "minimip/parameters.pb.h"
#include "minimip/solver.h"
#include "unit_tests/utils.h"

namespace minimip {

// ==========================================================================
// Branching Factory Tests
// ==========================================================================

TEST(BranchingFactoryTest, CreatesRandomBranchingWithValidParams) {
  BranchingParameters params;
  RandomBranchingParameters* random_params =
      params.mutable_random_branching_parameters();
  random_params->set_random_seed(42);

  absl::StatusOr<std::unique_ptr<BranchingRuleInterface>> branching =
      CreateBranchingRule(params);
  EXPECT_TRUE(branching.ok());
  EXPECT_TRUE(dynamic_cast<RandomBranching*>(branching.value().get()) !=
              nullptr);
}

TEST(BranchingFactoryTest, FailsToCreateBranchingRuleWithInvalidParams) {
  BranchingParameters params;  // Do not set random_branching_parameters
  auto branching = minimip::CreateBranchingRule(params);
  EXPECT_FALSE(branching.ok());
  EXPECT_EQ(branching.status().message(),
            "No branching-specific parameters set.");
}

// ==========================================================================
// Random Branching Tests
// ==========================================================================

TEST(RandomBranchingTest, NextBranchingVariableWithoutIntegerVariables) {
  BranchingParameters params;
  RandomBranchingParameters* random_params =
      params.mutable_random_branching_parameters();
  random_params->set_random_seed(42);

  absl::StatusOr<std::unique_ptr<BranchingRuleInterface>> branching =
      CreateBranchingRule(params);

  auto* random_branching =
      dynamic_cast<RandomBranching*>(branching.value().get());

  ASSERT_OK_AND_ASSIGN(std::unique_ptr<Solver> solver,
                       Solver::Create(MiniMipProblem()));

  absl::StatusOr<BranchingVariable> branching_variable =
      random_branching->NextBranchingVariable(*solver);
  EXPECT_FALSE(branching_variable.ok());
  EXPECT_EQ(branching_variable.status().message(),
            "No integer variables to branch on.");
}

TEST(RandomBranchingTest, NextBranchingVariableReturnsValidVariable) {
  BranchingParameters params;
  RandomBranchingParameters* random_params =
      params.mutable_random_branching_parameters();
  random_params->set_random_seed(42);

  absl::StatusOr<std::unique_ptr<BranchingRuleInterface>> branching =
      CreateBranchingRule(params);

  auto* random_branching =
      dynamic_cast<RandomBranching*>(branching.value().get());

  MiniMipProblem problem;
  problem.variables = std::vector<MiniMipVariable>{
      MiniMipVariable{.name = "1", .is_integer = true},
      MiniMipVariable{.name = "2", .is_integer = true},
      MiniMipVariable{.name = "3", .is_integer = true},
      MiniMipVariable{.name = "4", .is_integer = false},
  };

  ASSERT_OK_AND_ASSIGN(std::unique_ptr<Solver> solver, Solver::Create(problem));

  absl::StatusOr<BranchingVariable> branching_variable =
      random_branching->NextBranchingVariable(*solver);
  EXPECT_TRUE(branching_variable.ok());
  EXPECT_LT(branching_variable.value().index, ColIndex(4));
  EXPECT_GE(branching_variable.value().index, ColIndex(0));
}

// ==========================================================================
// MaxFractional Branching Tests
// ==========================================================================

// TODO (CG): Use Mockup to set specific primal value arrays to test the
// branching
TEST(MaxFractionalBranchingTest, NextBranchingVariableReturnsValidVariable) {
  // Creates a very simple model where all variables are non-negative.
  //
  // max: z = x2
  //  3*x1 + 2*x2 <= 6
  //  -3*x1 + 2*x2 <= 0
  // x1, x2 >= 0
  // x1, x2 integer
  MiniMipProblem problem;
  problem.variables.push_back(MiniMipVariable{.name = "x1",
                                              .objective_coefficient = 0.0,
                                              .lower_bound = 0,
                                              .upper_bound = kInf,
                                              .is_integer = true});
  problem.variables.push_back(MiniMipVariable{.name = "x2",
                                              .objective_coefficient = 1.0,
                                              .lower_bound = 0,
                                              .upper_bound = kInf,
                                              .is_integer = true});
  problem.constraints.push_back(MiniMipConstraint{
      .name = "ct1",
      .var_indices = {0, 1},
      .coefficients = {3.0, 2.0},
      .left_hand_side = -kInf,
      .right_hand_side = 6.0,
  });
  problem.constraints.push_back(MiniMipConstraint{
      .name = "ct2",
      .var_indices = {0, 1},
      .coefficients = {-3.0, 2.0},
      .left_hand_side = -kInf,
      .right_hand_side = 0.0,
  });
  problem.is_maximization = true;

  ASSERT_OK_AND_ASSIGN(std::unique_ptr<Solver> solver, Solver::Create(problem));

  ASSERT_OK(solver->mutable_lpi()->SolveLpWithDualSimplex())
  ASSERT_TRUE(solver->lpi()->IsSolved());

  BranchingParameters params;
  MaxFractionalBranchingParameters* fractional_branching_parameters =
      params.mutable_max_fractional_branching_parameters();

  absl::StatusOr<std::unique_ptr<BranchingRuleInterface>> branching =
      CreateBranchingRule(params);

  auto* max_fractional_branching =
      dynamic_cast<MaxFractionalBranching*>(branching.value().get());

  absl::StatusOr<BranchingVariable> branching_variable =
      max_fractional_branching->NextBranchingVariable(*solver);

  EXPECT_TRUE(branching_variable.ok());
  EXPECT_LT(branching_variable.value().index, ColIndex(4));
  EXPECT_GE(branching_variable.value().index, ColIndex(0));
}
}  // namespace minimip