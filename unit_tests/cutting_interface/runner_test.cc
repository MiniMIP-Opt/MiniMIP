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

#include "gmock/gmock.h"
#include "google/protobuf/text_format.h"
#include "gtest/gtest.h"
#include "src/cutting_interface/cuts_runner.h"
#include "src/parameters.pb.h"
#include "src/solver.h"
#include "unit_tests/utils.h"

namespace minimip {

enum class CutGeneratorType { kGomoryMixedInteger, kGomoryStrongCG };

std::unique_ptr<CutGeneratorInterface> CreateTestCutGenerator(
    CutGeneratorType type) {
  switch (type) {
    case CutGeneratorType::kGomoryMixedInteger: {
      CutGeneratorParameters params;
      CHECK(google::protobuf::TextFormat::ParseFromString(
          R"pb(
            tableau_rounding_generator_parameters: {
              use_mixed_integer_rounding: true
            })pb",
          &params));
      return std::make_unique<TableauRoundingGenerator>(params);
    }
    case CutGeneratorType::kGomoryStrongCG: {
      CutGeneratorParameters params;
      CHECK(google::protobuf::TextFormat::ParseFromString(
          R"pb(
            tableau_rounding_generator_parameters: {
              use_strong_cg_rounding: true
            })pb",
          &params));
      return std::make_unique<TableauRoundingGenerator>(params);
    }
  }
}

TEST(CutRunnerTests, CreateCutRunner) {
  CutRunnerParameters default_runner_params = DefaultCutRunnerParameters();

  // checks the default number of available generators and selectors
  ASSERT_EQ(default_runner_params.generator_parameters().size(), 1);
  ASSERT_EQ(default_runner_params.selector_parameters().size(), 1);

  // At this point, the default_params object has explicitly instantiated the
  // nested messages, so they are "set" even though their fields are at default
  // values.

  ASSERT_OK_AND_ASSIGN(std::unique_ptr<CutRunnerInterface> runner,
                       CreateCutRunner(default_runner_params));

  ASSERT_TRUE(runner->generator(0) != nullptr);
  ASSERT_TRUE(runner->selector() != nullptr);
}

TEST(CutRunnerTests, SimpleSolve) {
  MiniMipProblem problem;

  // Another small integer problem that requires proper use of slack
  // variables.
  //
  // max: z = 3x_1 + 2x_2
  //  x_1 + 2x_2 <= 4
  // 2x_1 +  x_2 <= 6
  // x_1, x_2 >= 0
  // x_1 integer
  //
  // ==> max(z) = 9 at x_1 = 3, x_2 = 0

  problem.variables.push_back(MiniMipVariable{.name = "x1",
                                              .objective_coefficient = 3.0,
                                              .lower_bound = 0,
                                              .upper_bound = kInf,
                                              .is_integer = true});
  problem.variables.push_back(MiniMipVariable{.name = "x2",
                                              .objective_coefficient = 2.0,
                                              .lower_bound = 0,
                                              .upper_bound = kInf,
                                              .is_integer = false});
  problem.constraints.push_back(MiniMipConstraint{
      .name = "ct1",
      .var_indices = {0, 1},
      .coefficients = {1.0, 2.0},
      .left_hand_side = -kInf,
      .right_hand_side = 4.0,
  });
  problem.constraints.push_back(MiniMipConstraint{
      .name = "ct2",
      .var_indices = {0, 1},
      .coefficients = {2.0, 1.0},
      .left_hand_side = -kInf,
      .right_hand_side = 6,
  });

  problem.is_maximization = true;

  ASSERT_OK_AND_ASSIGN(std::unique_ptr<Solver> solver, Solver::Create(problem));

  CHECK_OK(solver->mutable_lpi()->PopulateFromMipData(solver->mip_data()));
  CHECK_OK(solver->mutable_lpi()->SolveLpWithDualSimplex());

  ASSERT_TRUE(solver->lpi()->IsSolved());
  ASSERT_TRUE(solver->lpi()->IsOptimal());

  ASSERT_TRUE(solver->IsEqualToWithinTolerance(
      solver->lpi()->GetObjectiveValue(), -28.0 / 3.0));

  ASSERT_OK(
      solver->mutable_cut_runner()->SeparateCurrentLPSolution(*solver.get()));

  ASSERT_TRUE(solver->IsEqualToWithinTolerance(
      solver->lpi()->GetObjectiveValue(), -9.0));

  ASSERT_TRUE(solver->lpi()->IsSolved());
  ASSERT_TRUE(solver->lpi()->IsOptimal());
}

}  // namespace minimip
