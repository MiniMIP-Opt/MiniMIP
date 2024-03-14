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

#include "src/cutting_interface/aggregating_generator.h"

#include <algorithm>

#include "gmock/gmock.h"
#include "google/protobuf/text_format.h"
#include "gtest/gtest.h"
#include "src/cutting_interface/cuts_generator.h"
#include "src/parameters.pb.h"
#include "src/solver.h"
#include "unit_tests/utils.h"

namespace minimip {
namespace {

using testing::Gt;
using testing::IsEmpty;
using testing::Le;
using testing::Not;

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

class SmallModelSmokeTest
    : public testing::TestWithParam<std::tuple<CutGeneratorType, int>> {
 protected:
  void SetUp() final {
    generator_ = CreateTestCutGenerator(std::get<0>(GetParam()));
    std::pair<MiniMipProblem, SparseRow> data =
        CreateSampleProblem(std::get<1>(GetParam()));

    // Call the Create function to create a Solver object
    ASSERT_OK_AND_ASSIGN(std::unique_ptr<Solver> solver,
                         Solver::Create(data.first));

    solver_ = std::move(solver);
    optimum_ = std::move(data.second);
    CHECK_OK(solver_->mutable_lpi()->PopulateFromMipData(solver_->mip_data()));
  }

  static std::pair<MiniMipProblem, SparseRow> CreateSampleProblem(
      int model_index) {
    MiniMipProblem problem;
    SparseRow optimum;
    switch (model_index) {
      case 0: {
        // Creates a very simple model where all variables are non-negative.
        //
        // max: z = x2
        //  3*x1 + 2*x2 <= 6
        // -3*x1 + 2*x2 <= 0
        // x1, x2 >= 0
        // x1, x2 integer
        problem.variables.push_back(
            MiniMipVariable{.name = "x1",
                            .objective_coefficient = 0.0,
                            .lower_bound = 0,
                            .upper_bound = kInf,
                            .is_integer = true});
        problem.variables.push_back(
            MiniMipVariable{.name = "x2",
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
        optimum = CreateSparseRow({{0, 1}, {1, 1}});
        break;
      }
      case 1: {
        // Another small integer problem that requires proper use of slack
        // variables.
        //
        // max: z = 0.2*x1 + 0.1*x2
        // 0.1*x1 + 0.1*x2 <= 1.0
        // 0.1*x1          <= 0.75
        // 0.8*x1 + 2.0*x2 <= 10
        // x1, x2 >= 0
        // x1, x2 integer
        problem.variables.push_back(
            MiniMipVariable{.name = "x1",
                            .objective_coefficient = 0.2,
                            .lower_bound = 0,
                            .upper_bound = kInf,
                            .is_integer = true});
        problem.variables.push_back(
            MiniMipVariable{.name = "x2",
                            .objective_coefficient = 0.1,
                            .lower_bound = 0,
                            .upper_bound = kInf,
                            .is_integer = true});
        problem.constraints.push_back(MiniMipConstraint{
            .name = "ct1",
            .var_indices = {0, 1},
            .coefficients = {0.1, 0.1},
            .left_hand_side = -kInf,
            .right_hand_side = 1.0,
        });
        problem.constraints.push_back(MiniMipConstraint{
            .name = "ct2",
            .var_indices = {0},
            .coefficients = {0.1},
            .left_hand_side = -kInf,
            .right_hand_side = 0.75,
        });
        problem.constraints.push_back(MiniMipConstraint{
            .name = "ct3",
            .var_indices = {0, 1},
            .coefficients = {0.8, 2.0},
            .left_hand_side = -kInf,
            .right_hand_side = 10.0,
        });
        problem.is_maximization = true;
        optimum = CreateSparseRow({{0, 7}, {1, 2}});
        break;
      }
      case 2: {
        // Problem where some constraints use the left hand side.
        //
        // min: z = x1 + 2.0*x2
        // 1.0 <= 1.5*x1 + 1.5*x2
        //            x1          <= 2
        //                     x2 <= 2
        // x1, x2 >= 0
        // x1, x2 integer
        problem.variables.push_back(
            MiniMipVariable{.name = "x1",
                            .objective_coefficient = 1.0,
                            .lower_bound = 0,
                            .upper_bound = kInf,
                            .is_integer = true});
        problem.variables.push_back(
            MiniMipVariable{.name = "x2",
                            .objective_coefficient = 2.0,
                            .lower_bound = 0,
                            .upper_bound = kInf,
                            .is_integer = true});
        problem.constraints.push_back(MiniMipConstraint{
            .name = "ct1",
            .var_indices = {0, 1},
            .coefficients = {1.5, 1.5},
            .left_hand_side = 1.0,
            .right_hand_side = kInf,
        });
        problem.constraints.push_back(MiniMipConstraint{
            .name = "ct2",
            .var_indices = {0},
            .coefficients = {1.0},
            .left_hand_side = -kInf,
            .right_hand_side = 2,
        });
        problem.constraints.push_back(MiniMipConstraint{
            .name = "ct3",
            .var_indices = {1},
            .coefficients = {1.0},
            .left_hand_side = -kInf,
            .right_hand_side = 2.0,
        });
        problem.is_maximization = false;
        optimum = CreateSparseRow({{0, 1}});
        break;
      }
      case 3: {
        // Problem where the variables have lower and upper bounds.
        //
        // min: z = x1 + 1.5*x2
        // -1.0 <= 1.5*x1 + 1.5*x2
        // -1.0 <=     x1          <= 3
        // -1.0 <=              x2 <= 7
        // x1, x2 integer
        problem.variables.push_back(
            MiniMipVariable{.name = "x1",
                            .objective_coefficient = 1.0,
                            .lower_bound = -1.0,
                            .upper_bound = 3,
                            .is_integer = true});
        problem.variables.push_back(
            MiniMipVariable{.name = "x2",
                            .objective_coefficient = 1.5,
                            .lower_bound = -1.0,
                            .upper_bound = 7,
                            .is_integer = true});
        problem.constraints.push_back(MiniMipConstraint{
            .name = "ct1",
            .var_indices = {0, 1},
            .coefficients = {1.5, 1.5},
            .left_hand_side = -1.0,
            .right_hand_side = kInf,
        });
        problem.is_maximization = false;
        optimum = CreateSparseRow({{0, 1}, {1, -1}});
        break;
      }
      default:
        LOG(FATAL) << "Invalid model: " << model_index;
    }
    return std::make_pair(problem, optimum);
  }

  absl::Status AddCutsAndResolveLp(const std::vector<CutData>& cuts) {
    for (const CutData& cut : cuts) {
      RETURN_IF_ERROR(solver_->mutable_lpi()->AddRow(
          cut.row(), -solver_->lpi()->Infinity(), cut.right_hand_side(), ""));
    }
    return solver_->mutable_lpi()->SolveLpWithDualSimplex();
  }

  absl::StatusOr<bool> SolutionIsMipFeasible() {
    ASSIGN_OR_RETURN((const absl::StrongVector<ColIndex, double> primal_values),
                     solver_->lpi()->GetPrimalValues());
    return std::all_of(
        solver_->mip_data().integer_variables().begin(),
        solver_->mip_data().integer_variables().end(),
        [&primal_values, solver = solver_.get()](ColIndex col) {
          return solver->IsIntegerWithinTolerance(primal_values[col]);
        });
  }

  std::unique_ptr<CutGeneratorInterface> generator_;
  std::unique_ptr<Solver> solver_;
  SparseRow optimum_;
};

INSTANTIATE_TEST_SUITE_P(
    GomoryMixedInteger, SmallModelSmokeTest,
    testing::Combine(testing::Values(CutGeneratorType::kGomoryMixedInteger),
                     testing::Range(0, 4)));

INSTANTIATE_TEST_SUITE_P(
    StrongCG, SmallModelSmokeTest,
    testing::Combine(testing::Values(CutGeneratorType::kGomoryStrongCG),
                     testing::Range(0, 4)));

TEST_P(SmallModelSmokeTest, SmokeTest) {
  ASSERT_OK(solver_->mutable_lpi()->SolveLpWithPrimalSimplex());
  ASSERT_TRUE(solver_->lpi()->IsSolved());
  ASSERT_TRUE(solver_->lpi()->IsOptimal());
  for (int i = 0; i < 10; ++i) {
    ASSERT_OK_AND_ASSIGN(const std::vector<CutData> cuts,
                         generator_->GenerateCuttingPlanes(*solver_));
    ASSERT_THAT(cuts, Not(IsEmpty()));

    ASSERT_OK_AND_ASSIGN(
        (const absl::StrongVector<ColIndex, double> primal_values),
        solver_->lpi()->GetPrimalValues());
    const SparseRow lp_optimum = SparseRow(primal_values);

    double tolerance = 1e-12;
    for (const CutData& cut : cuts) {
      // All cuts should remove the LP optimum.
      ASSERT_THAT(cut.row(), Activation(lp_optimum, Gt(cut.right_hand_side() - tolerance)));

      // No cuts should remove the MIP optimum.
      ASSERT_THAT(cut.row(), Activation(optimum_, Le(cut.right_hand_side() + tolerance)));
    }

    ASSERT_OK(AddCutsAndResolveLp(cuts));
    ASSERT_TRUE(solver_->lpi()->IsSolved());
    ASSERT_TRUE(solver_->lpi()->IsOptimal());
    ASSERT_OK_AND_ASSIGN(bool is_mip_feasible, SolutionIsMipFeasible());
    if (is_mip_feasible) break;
  }

  // The cutting plane algorithm is complete when using Gomory cuts, so we
  // expect that the MIP is solved. Note that this is not guaranteed for all
  // generator families.
  ASSERT_OK_AND_ASSIGN(bool is_mip_feasible, SolutionIsMipFeasible());
  ASSERT_TRUE(is_mip_feasible);
}
}  // namespace
}  // namespace minimip
