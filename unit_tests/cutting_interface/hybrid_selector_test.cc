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

#include "src/cutting_interface/hybrid_selector.h"

#include <algorithm>

#include "gmock/gmock.h"
#include "google/protobuf/text_format.h"
#include "gtest/gtest.h"
#include "src/parameters.pb.h"
#include "src/solver.h"
#include "unit_tests/utils.h"

namespace minimip {
namespace {

using testing::Gt;
using testing::IsEmpty;
using testing::Le;
using testing::Not;

// This is a helper function to create the selector parameters.
enum class CutSelectorType { kHybridSelectorUnsigned, kHybridSelectorSigned };

std::unique_ptr<CutSelector> CreateCutSelector(CutSelectorType type) {
  switch (type) {
    case CutSelectorType::kHybridSelectorUnsigned: {
      CutSelectorParameters params;
      CHECK(google::protobuf::TextFormat::ParseFromString(
          R"pb(
            hybrid_selector_parameters: { signed_orthogonality: true })pb",
          &params));
      return std::make_unique<HybridSelector>(params);
    }
    case CutSelectorType::kHybridSelectorSigned: {
      CutSelectorParameters params;
      CHECK(google::protobuf::TextFormat::ParseFromString(
          R"pb(
            hybrid_selector_parameters: { signed_orthogonality: false })pb",
          &params));
      return std::make_unique<HybridSelector>(params);
    }
  }
}

class MinimalCutSelectorTest
    : public testing::TestWithParam<std::tuple<CutSelectorType>> {
 protected:
  void SetUp() final {
    // Here we initialize the cut selector.
    CutSelectorParameters params;
    selector_ = CreateCutSelector(std::get<0>(GetParam()));
    MiniMipProblem problem = CreateSampleProblem();

    // Call the Create function to create a Solver object
    ASSERT_OK_AND_ASSIGN(std::unique_ptr<Solver> solver,
                         Solver::Create(MiniMipParameters{}, problem));

    solver_ = std::move(solver);
  }

  // Creates a very simple model where all variables are non-negative.
  //
  // max: z = x1 + x2
  //  3*x1 + 2*x2 <= 6
  // -3*x1 + 2*x2 <= 0
  // x1, x2 >= 0
  // x1 integer

  static MiniMipProblem CreateSampleProblem() {
    MiniMipProblem problem;
    SparseRow optimum;
    problem.variables.push_back(MiniMipVariable{.name = "x1",
                                                .objective_coefficient = 1.0,
                                                .lower_bound = 0,
                                                .upper_bound = kInf,
                                                .is_integer = true});
    problem.variables.push_back(MiniMipVariable{.name = "x2",
                                                .objective_coefficient = 0.0,
                                                .lower_bound = 0,
                                                .upper_bound = kInf,
                                                .is_integer = false});
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
    return problem;
  }

  std::unique_ptr<CutSelector> selector_;
  std::unique_ptr<Solver> solver_;
};

INSTANTIATE_TEST_SUITE_P(
    SimpleCutExamples, MinimalCutSelectorTest,
    testing::Values(CutSelectorType::kHybridSelectorUnsigned,
                    CutSelectorType::kHybridSelectorSigned));

TEST_P(MinimalCutSelectorTest, SelectsFirstCutIfAllAreIdentical) {
  ASSERT_OK(solver_->mutable_lpi()->SolveLPWithPrimalSimplex());
  ASSERT_TRUE(solver_->lpi()->IsSolved());
  ASSERT_TRUE(solver_->lpi()->IsOptimal());

  SparseRow lp_optimum = CreateSparseRow({{0, 1.0}, {1, 1.5}});

  std::vector<CutData> cuts = {solver_->cut_registry().CreateCut(
                                   solver_->mip_data(), lp_optimum,
                                   CreateSparseRow({{1, 1.0}}), 1.0, "cut_1"),
                               solver_->cut_registry().CreateCut(
                                   solver_->mip_data(), lp_optimum,
                                   CreateSparseRow({{1, 1.0}}), 1.0, "cut_2")};

  absl::StatusOr<std::vector<CutData>> selected_cuts =
      selector_->SelectCuttingPlanes(*solver_, cuts);
  ASSERT_EQ(selected_cuts.value().size(), 1);
  ASSERT_EQ(selected_cuts.value()[0].name(), "cut_1");
}

TEST_P(MinimalCutSelectorTest, SelectsIdenticalCutIfItIsForced) {
  ASSERT_OK_AND_ASSIGN(std::unique_ptr<Solver> solver,
                       Solver::Create(MiniMipParameters{}, MiniMipProblem{}));

  SparseRow lp_optimum = CreateSparseRow({{0, 1.0}, {1, 1.5}});

  std::vector<CutData> cuts = {
      solver_->cut_registry().CreateCut(solver_->mip_data(), lp_optimum,
                                        CreateSparseRow({{1, 1.0}}), 1.0,
                                        "cut_1"),
      solver_->cut_registry().CreateCut(solver_->mip_data(), lp_optimum,
                                        CreateSparseRow({{1, 1.0}}), 1.0,
                                        "cut_2", true)};

  absl::StatusOr<std::vector<CutData>> selected_cuts =
      selector_->SelectCuttingPlanes(*solver, cuts);
  ASSERT_EQ(selected_cuts.value().size(), 2);
  ASSERT_EQ(selected_cuts.value()[0].name(), "cut_1");
  ASSERT_EQ(selected_cuts.value()[1].name(), "cut_2");
}

TEST_P(MinimalCutSelectorTest, SelectsParallelCutWithHigherScore) {
  ASSERT_OK_AND_ASSIGN(std::unique_ptr<Solver> solver,
                       Solver::Create(MiniMipParameters{}, MiniMipProblem{}));

  SparseRow lp_optimum = CreateSparseRow({{0, 1.0}, {1, 1.5}});

  std::vector<CutData> cuts = {solver_->cut_registry().CreateCut(
                                   solver_->mip_data(), lp_optimum,
                                   CreateSparseRow({{1, 1.0}}), 1.0, "cut_1"),
                               solver_->cut_registry().CreateCut(
                                   solver_->mip_data(), lp_optimum,
                                   CreateSparseRow({{1, 1.0}}), 0.9, "cut_2")};

  absl::StatusOr<std::vector<CutData>> selected_cuts =
      selector_->SelectCuttingPlanes(*solver, cuts);
  ASSERT_EQ(selected_cuts.value().size(), 1);
  ASSERT_EQ(selected_cuts.value()[0].name(), "cut_2");
}

TEST_P(MinimalCutSelectorTest, SelectsForcedParallelCutWithLowerScoreLast) {
  ASSERT_OK_AND_ASSIGN(std::unique_ptr<Solver> solver,
                       Solver::Create(MiniMipParameters{}, MiniMipProblem{}));

  SparseRow lp_optimum = CreateSparseRow({{0, 1.0}, {1, 1.5}});

  std::vector<CutData> cuts = {
      solver_->cut_registry().CreateCut(solver_->mip_data(), lp_optimum,
                                        CreateSparseRow({{1, 1.0}}), 1.0,
                                        "cut_1", true),
      solver_->cut_registry().CreateCut(solver_->mip_data(), lp_optimum,
                                        CreateSparseRow({{1, 1.0}}), 0.9,
                                        "cut_2")};

  absl::StatusOr<std::vector<CutData>> selected_cuts =
      selector_->SelectCuttingPlanes(*solver, cuts);
  ASSERT_EQ(selected_cuts.value().size(), 2);
  ASSERT_EQ(selected_cuts.value()[0].name(), "cut_2");
  ASSERT_EQ(selected_cuts.value()[1].name(), "cut_1");
}

TEST_P(MinimalCutSelectorTest, SelectsAllOrthgonalCuts) {
  ASSERT_OK_AND_ASSIGN(std::unique_ptr<Solver> solver,
                       Solver::Create(MiniMipParameters{}, MiniMipProblem{}));

  SparseRow lp_optimum = CreateSparseRow({{0, 1.0}, {1, 1.5}});

  std::vector<CutData> cuts = {solver_->cut_registry().CreateCut(
                                   solver_->mip_data(), lp_optimum,
                                   CreateSparseRow({{1, 1.0}}), 1.0, "cut_1"),
                               solver_->cut_registry().CreateCut(
                                   solver_->mip_data(), lp_optimum,
                                   CreateSparseRow({{0, 1.0}}), 0.9, "cut_2")};

  absl::StatusOr<std::vector<CutData>> selected_cuts =
      selector_->SelectCuttingPlanes(*solver, cuts);

  ASSERT_EQ(selected_cuts.value().size(), 2);
  ASSERT_EQ(selected_cuts.value()[0].name(), "cut_1");
  ASSERT_EQ(selected_cuts.value()[1].name(), "cut_2");

  std::vector<CutData> reverse_ordered_cuts = {
      solver_->cut_registry().CreateCut(solver_->mip_data(), lp_optimum,
                                        CreateSparseRow({{0, 1.0}}), 1.0,
                                        "cut_2"),
      solver_->cut_registry().CreateCut(solver_->mip_data(), lp_optimum,
                                        CreateSparseRow({{1, 1.0}}), 0.9,
                                        "cut_1")};

  absl::StatusOr<std::vector<CutData>> selected_reverse_cuts =
      selector_->SelectCuttingPlanes(*solver, reverse_ordered_cuts);

  ASSERT_EQ(selected_reverse_cuts.value().size(), 2);
  ASSERT_EQ(selected_reverse_cuts.value()[0].name(), "cut_1");
  ASSERT_EQ(selected_reverse_cuts.value()[1].name(), "cut_2");
}

TEST_P(MinimalCutSelectorTest, FiltersOrthogonalCutBelowScoreThreshold) {
  ASSERT_OK_AND_ASSIGN(std::unique_ptr<Solver> solver,
                       Solver::Create(MiniMipParameters{}, MiniMipProblem{}));

  SparseRow lp_optimum = CreateSparseRow({{0, 1.5}, {1, 1.0}});

  std::vector<CutData> cuts = {solver_->cut_registry().CreateCut(
                                   solver_->mip_data(), lp_optimum,
                                   CreateSparseRow({{0, 1.0}}), 1.0, "cut_1"),
                               solver_->cut_registry().CreateCut(
                                   solver_->mip_data(), lp_optimum,
                                   CreateSparseRow({{1, 1.0}}), 1.0, "cut_2")};

  absl::StatusOr<std::vector<CutData>> selected_cuts =
      selector_->SelectCuttingPlanes(*solver, cuts);

  ASSERT_EQ(selected_cuts.value().size(), 1);
  ASSERT_EQ(selected_cuts.value()[0].name(), "cut_1");
}

TEST_P(MinimalCutSelectorTest, SelectsNoCutsIfNoneAreAboveScoreThreshold) {
  ASSERT_OK_AND_ASSIGN(std::unique_ptr<Solver> solver,
                       Solver::Create(MiniMipParameters{}, MiniMipProblem{}));

  SparseRow lp_optimum = CreateSparseRow({{0, 1.5}, {1, 1.0}});

  std::vector<CutData> cuts = {solver_->cut_registry().CreateCut(
                                   solver_->mip_data(), lp_optimum,
                                   CreateSparseRow({{1, 1.0}}), 1.0, "cut_1"),
                               solver_->cut_registry().CreateCut(
                                   solver_->mip_data(), lp_optimum,
                                   CreateSparseRow({{1, 1.0}}), 1.0, "cut_2")};

  absl::StatusOr<std::vector<CutData>> selected_cuts =
      selector_->SelectCuttingPlanes(*solver, cuts);
  ASSERT_EQ(selected_cuts.value().size(), 0);
}

TEST_P(MinimalCutSelectorTest, SelectBetterOrthogonalCut) {
  ASSERT_OK_AND_ASSIGN(std::unique_ptr<Solver> solver,
                       Solver::Create(MiniMipParameters{}, MiniMipProblem{}));

  SparseRow lp_optimum = CreateSparseRow({{0, 1.5}, {1, 1.5}});

  std::vector<CutData> cuts = {solver_->cut_registry().CreateCut(
                                   solver_->mip_data(), lp_optimum,
                                   CreateSparseRow({{0, 1.0}}), 0.4, "cut_1"),
                               solver_->cut_registry().CreateCut(
                                   solver_->mip_data(), lp_optimum,
                                   CreateSparseRow({{1, 1.0}}), 0.5, "cut_2"),
                               solver_->cut_registry().CreateCut(
                                   solver_->mip_data(), lp_optimum,
                                   CreateSparseRow({{1, 1.0}}), 1.0, "cut_3")};

  absl::StatusOr<std::vector<CutData>> selected_cuts =
      selector_->SelectCuttingPlanes(*solver, cuts);
  ASSERT_EQ(selected_cuts.value().size(), 2);
  ASSERT_EQ(selected_cuts.value()[0].name(), "cut_1");
  ASSERT_EQ(selected_cuts.value()[1].name(), "cut_2");
}

TEST_P(MinimalCutSelectorTest, SelectBestCutsFirstAndForcedCutsInPosition) {
  ASSERT_OK_AND_ASSIGN(std::unique_ptr<Solver> solver,
                       Solver::Create(MiniMipParameters{}, MiniMipProblem{}));

  SparseRow lp_optimum = CreateSparseRow({{0, 1.5}, {1, 1.5}});

  std::vector<CutData> cuts = {
      solver_->cut_registry().CreateCut(solver_->mip_data(), lp_optimum,
                                        CreateSparseRow({{0, 1.0}}), 1.0,
                                        "cut_1"),
      solver_->cut_registry().CreateCut(solver_->mip_data(), lp_optimum,
                                        CreateSparseRow({{1, 1.0}}), 0.5,
                                        "cut_2"),
      solver_->cut_registry().CreateCut(solver_->mip_data(), lp_optimum,
                                        CreateSparseRow({{1, 1.0}}), 0.9,
                                        "cut_3", true)};

  absl::StatusOr<std::vector<CutData>> selected_cuts =
      selector_->SelectCuttingPlanes(*solver, cuts);
  ASSERT_EQ(selected_cuts.value().size(), 3);
  ASSERT_EQ(selected_cuts.value()[0].name(), "cut_2");
  ASSERT_EQ(selected_cuts.value()[1].name(), "cut_1");
  ASSERT_EQ(selected_cuts.value()[2].name(), "cut_3");
}

TEST_P(MinimalCutSelectorTest, SelectInvertedCutIfSignedOrthogonalityIsUsed) {
  ASSERT_OK_AND_ASSIGN(std::unique_ptr<Solver> solver,
                       Solver::Create(MiniMipParameters{}, MiniMipProblem{}));

  SparseRow lp_optimum = CreateSparseRow({{0, 1.5}, {1, 1.0}});

  std::vector<CutData> cuts = {
      solver_->cut_registry().CreateCut(solver_->mip_data(), lp_optimum,
                                        CreateSparseRow({{1, 1.0}}), 0.4,
                                        "cut_1"),
      solver_->cut_registry().CreateCut(solver_->mip_data(), lp_optimum,
                                        CreateSparseRow({{0, 1.0}, {1, -0.5}}),
                                        0.5, "cut_2")};

  absl::StatusOr<std::vector<CutData>> selected_cuts =
      selector_->SelectCuttingPlanes(*solver, cuts);

  if (std::get<0>(GetParam()) == CutSelectorType::kHybridSelectorUnsigned) {
    ASSERT_EQ(selected_cuts.value().size(), 2);
    ASSERT_EQ(selected_cuts.value()[0].name(), "cut_1");
    ASSERT_EQ(selected_cuts.value()[1].name(), "cut_2");
  } else if (std::get<0>(GetParam()) ==
             CutSelectorType::kHybridSelectorSigned) {
    ASSERT_EQ(selected_cuts.value().size(), 1);
    ASSERT_EQ(selected_cuts.value()[0].name(), "cut_1");
  }
}

}  // namespace
}  // namespace minimip
