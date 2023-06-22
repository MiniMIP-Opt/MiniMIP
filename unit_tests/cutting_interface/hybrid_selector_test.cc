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
  }

  std::unique_ptr<CutSelector> selector_;
};

INSTANTIATE_TEST_SUITE_P(
    SimpleCutExamples, MinimalCutSelectorTest,
    testing::Values(CutSelectorType::kHybridSelectorUnsigned,
                    CutSelectorType::kHybridSelectorSigned));

TEST_P(MinimalCutSelectorTest, SelectsFirstCutIfAllAreIdentical) {
  ASSERT_OK_AND_ASSIGN(std::unique_ptr<Solver> solver,
                       Solver::Create(MiniMipParameters{}, MiniMipProblem{}));
  std::vector<CutData> cuts = {
      CutData{CreateSparseRow({{1, 1.0}}), 1.0, 1, 1, 1.0, 1.0, "cut_1", false},
      CutData{CreateSparseRow({{1, 1.0}}), 1.0, 1, 1, 1.0, 1.0, "cut_2",
              false}};

  absl::StatusOr<std::vector<CutData>> selected_cuts =
      selector_->SelectCuttingPlanes(*solver, cuts);
  ASSERT_EQ(selected_cuts.value().size(), 1);
  ASSERT_EQ(selected_cuts.value()[0].name(), "cut_1");
}

TEST_P(MinimalCutSelectorTest, SelectsIdenticalCutIfItIsForced) {
  ASSERT_OK_AND_ASSIGN(std::unique_ptr<Solver> solver,
                       Solver::Create(MiniMipParameters{}, MiniMipProblem{}));
  std::vector<CutData> cuts = {
      CutData{CreateSparseRow({{1, 1.0}}), 1.0, 1, 1, 1.0, 1.0, "cut_1", false},
      CutData{CreateSparseRow({{1, 1.0}}), 1.0, 1, 1, 1.0, 1.0, "cut_2", true}};

  absl::StatusOr<std::vector<CutData>> selected_cuts =
      selector_->SelectCuttingPlanes(*solver, cuts);
  ASSERT_EQ(selected_cuts.value().size(), 2);
  ASSERT_EQ(selected_cuts.value()[0].name(), "cut_1");
  ASSERT_EQ(selected_cuts.value()[1].name(), "cut_2");
}

TEST_P(MinimalCutSelectorTest, SelectsSecondIdenticalCutWithHigherScore) {
  ASSERT_OK_AND_ASSIGN(std::unique_ptr<Solver> solver,
                       Solver::Create(MiniMipParameters{}, MiniMipProblem{}));
  std::vector<CutData> cuts = {
      CutData{CreateSparseRow({{1, 1.0}}), 1.0, 1, 1, 1.0, 1.0, "cut_1", false},
      CutData{CreateSparseRow({{1, 1.0}}), 1.0, 1, 1, 1.0, 1.1, "cut_2",
              false}};

  absl::StatusOr<std::vector<CutData>> selected_cuts =
      selector_->SelectCuttingPlanes(*solver, cuts);
  ASSERT_EQ(selected_cuts.value().size(), 1);
  ASSERT_EQ(selected_cuts.value()[0].name(), "cut_2");
}

TEST_P(MinimalCutSelectorTest, SelectsForcedIdenticalCutWithLowerScoreLast) {
  ASSERT_OK_AND_ASSIGN(std::unique_ptr<Solver> solver,
                       Solver::Create(MiniMipParameters{}, MiniMipProblem{}));
  std::vector<CutData> cuts = {
      CutData{CreateSparseRow({{1, 1.0}}), 1.0, 1, 1, 1.0, 1.0, "cut_1", true},
      CutData{CreateSparseRow({{1, 1.0}}), 1.0, 1, 1, 1.0, 1.1, "cut_2",
              false}};

  absl::StatusOr<std::vector<CutData>> selected_cuts =
      selector_->SelectCuttingPlanes(*solver, cuts);
  ASSERT_EQ(selected_cuts.value().size(), 2);
  ASSERT_EQ(selected_cuts.value()[0].name(), "cut_2");
  ASSERT_EQ(selected_cuts.value()[1].name(), "cut_1");
}

TEST_P(MinimalCutSelectorTest, SelectsAllOrthgonalCuts) {
  ASSERT_OK_AND_ASSIGN(std::unique_ptr<Solver> solver,
                       Solver::Create(MiniMipParameters{}, MiniMipProblem{}));
  std::vector<CutData> cuts = {
      CutData{CreateSparseRow({{1, 1.0}}), 1.0, 1, 1, 1.0, 1.1, "cut_1", false},
      CutData{CreateSparseRow({{0, 1.0}}), 1.0, 1, 1, 1.0, 1.0, "cut_2",
              false}};

  absl::StatusOr<std::vector<CutData>> selected_cuts =
      selector_->SelectCuttingPlanes(*solver, cuts);
  ASSERT_EQ(selected_cuts.value().size(), 2);
  ASSERT_EQ(selected_cuts.value()[0].name(), "cut_1");
  ASSERT_EQ(selected_cuts.value()[1].name(), "cut_2");
}

TEST_P(MinimalCutSelectorTest, FiltersOrthogonalCutBelowScoreThreshold) {
  ASSERT_OK_AND_ASSIGN(std::unique_ptr<Solver> solver,
                       Solver::Create(MiniMipParameters{}, MiniMipProblem{}));
  std::vector<CutData> cuts = {
      CutData{CreateSparseRow({{1, 1.0}}), 1.0, 1, 1, 1.0, 1.0, "cut_1", false},
      CutData{CreateSparseRow({{0, 1.0}}), 1.0, 1, 0, 0.0, 0.0, "cut_2",
              false}};

  absl::StatusOr<std::vector<CutData>> selected_cuts =
      selector_->SelectCuttingPlanes(*solver, cuts);
  ASSERT_EQ(selected_cuts.value().size(), 1);
  ASSERT_EQ(selected_cuts.value()[0].name(), "cut_1");
}

TEST_P(MinimalCutSelectorTest, SelectsNoCutsIfNoneAreAboveScoreThreshold) {
  ASSERT_OK_AND_ASSIGN(std::unique_ptr<Solver> solver,
                       Solver::Create(MiniMipParameters{}, MiniMipProblem{}));
  std::vector<CutData> cuts = {
      CutData{CreateSparseRow({{1, 1.0}}), 1.0, 1, 0, 0.0, 0.0, "cut_1", false},
      CutData{CreateSparseRow({{0, 1.0}}), 1.0, 1, 0, 0.0, 0.0, "cut_2",
              false}};

  absl::StatusOr<std::vector<CutData>> selected_cuts =
      selector_->SelectCuttingPlanes(*solver, cuts);
  ASSERT_EQ(selected_cuts.value().size(), 0);
}

TEST_P(MinimalCutSelectorTest, SelectBetterOrthogonalCut) {
  ASSERT_OK_AND_ASSIGN(std::unique_ptr<Solver> solver,
                       Solver::Create(MiniMipParameters{}, MiniMipProblem{}));
  std::vector<CutData> cuts = {
      CutData{CreateSparseRow({{1, 1.0}}), 1.0, 1, 1, 1.0, 1.1, "cut_1", false},
      CutData{CreateSparseRow({{0, 1.0}}), 1.0, 1, 1, 1.0, 1.0, "cut_2", false},
      CutData{CreateSparseRow({{0, 2.0}}), 1.0, 1, 1, 1.0, 0.5, "cut_3",
              false}};

  absl::StatusOr<std::vector<CutData>> selected_cuts =
      selector_->SelectCuttingPlanes(*solver, cuts);
  ASSERT_EQ(selected_cuts.value().size(), 2);
  ASSERT_EQ(selected_cuts.value()[0].name(), "cut_1");
  ASSERT_EQ(selected_cuts.value()[1].name(), "cut_2");
}

TEST_P(MinimalCutSelectorTest, SelectBestCutsFirstAndForcedCutsInPosition) {
  ASSERT_OK_AND_ASSIGN(std::unique_ptr<Solver> solver,
                       Solver::Create(MiniMipParameters{}, MiniMipProblem{}));
  std::vector<CutData> cuts = {
      CutData{CreateSparseRow({{1, 1.0}}), 1.0, 1, 1, 1.0, 0.5, "cut_1", false},
      CutData{CreateSparseRow({{0, 1.0}}), 1.0, 1, 1, 1.0, 1.0, "cut_2", false},
      CutData{CreateSparseRow({{0, 2.0}}), 1.0, 1, 1, 1.0, 0.6, "cut_3", true}};

  absl::StatusOr<std::vector<CutData>> selected_cuts =
      selector_->SelectCuttingPlanes(*solver, cuts);
  ASSERT_EQ(selected_cuts.value().size(), 3);
  ASSERT_EQ(selected_cuts.value()[0].name(), "cut_2");
  ASSERT_EQ(selected_cuts.value()[1].name(), "cut_3");
  ASSERT_EQ(selected_cuts.value()[2].name(), "cut_1");
}

TEST_P(MinimalCutSelectorTest, SelectInvertedCutIfSignedOrthogonalityIsUsed) {
  ASSERT_OK_AND_ASSIGN(std::unique_ptr<Solver> solver,
                       Solver::Create(MiniMipParameters{}, MiniMipProblem{}));
  std::vector<CutData> cuts = {
      CutData{CreateSparseRow({{1, 1.0}}), 1.0, 1, 1, 1.0, 1.1, "cut_1", false},
      CutData{CreateSparseRow({{0, -1.0}}), 1.0, 1, 1, 1.0, 1.0, "cut_2",
              false}};

  absl::StatusOr<std::vector<CutData>> selected_cuts =
      selector_->SelectCuttingPlanes(*solver, cuts);

  if (std::get<0>(GetParam()) == CutSelectorType::kHybridSelectorUnsigned) {
    ASSERT_EQ(selected_cuts.value().size(), 1);
    ASSERT_EQ(selected_cuts.value()[0].name(), "cut_1");
  } else if (std::get<0>(GetParam()) ==
             CutSelectorType::kHybridSelectorSigned) {
    ASSERT_EQ(selected_cuts.value().size(), 2);
    ASSERT_EQ(selected_cuts.value()[0].name(), "cut_1");
    ASSERT_EQ(selected_cuts.value()[1].name(), "cut_2");
  }
}

}  // namespace
}  // namespace minimip
