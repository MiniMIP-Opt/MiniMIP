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

std::unique_ptr<Selector> CreateCutSelector(CutSelectorType type) {
  switch (type) {
    case CutSelectorType::kHybridSelectorUnsigned: {
      SelectorParameters params;
      CHECK(google::protobuf::TextFormat::ParseFromString(
          R"pb(
            hybrid_selector_parameters: { unsigned_orthogonality: true })pb",
          &params));
      return std::make_unique<HybridSelector>(params);
    }
    case CutSelectorType::kHybridSelectorSigned: {
      SelectorParameters params;
      CHECK(google::protobuf::TextFormat::ParseFromString(
          R"pb(
            hybrid_selector_parameters: { unsigned_orthogonality: false })pb",
          &params));
      return std::make_unique<HybridSelector>(params);
    }
  }
}

class MinimalCutSelectorTest
    : public testing::TestWithParam<std::tuple<CutSelectorType, int>> {
 protected:
  void SetUp() final {
    // Here we initialize the cut selector.
    SelectorParameters params;
    selector_ = CreateCutSelector(std::get<0>(GetParam()));

    std::vector<CutData> cuts = CreateCuttingplanes(std::get<1>(GetParam()));
  }

  static std::vector<CutData> CreateCuttingplanes(int switch_variable_index) {
    std::vector<CutData> cuts;

    switch (switch_variable_index) {
      case 0: {
        // cut_1 is identical to cut_2
        CutData cut1(CreateSparseRow({{1, 1.0}}), 1.0, 1, 1, 1.0, 1.0, "cut_1",
                     false);
        CutData cut2(CreateSparseRow({{1, 1.0}}), 1.0, 1, 1, 1.0, 1.0, "cut_2",
                     false);
        cuts.push_back(cut1);
        cuts.push_back(cut2);
        break;
      }
      case 1: {
        // cut_1 is identical to cut_2, but cut_2 is forced
        CutData cut1(CreateSparseRow({{1, 1.0}}), 1.0, 1, 1, 1.0, 1.0, "cut_1",
                     false);
        CutData cut2(CreateSparseRow({{1, 1.0}}), 1.0, 1, 1, 1.0, 1.0, "cut_2",
                     true);
        cuts.push_back(cut1);
        cuts.push_back(cut2);
        break;
      }
      case 2: {
        // cut_1 is scored lower than parallel cut_2
        CutData cut1(CreateSparseRow({{1, 1.0}}), 1.0, 1, 1, 1.0, 1.0, "cut_1",
                     false);
        CutData cut2(CreateSparseRow({{1, 1.0}}), 1.0, 1, 1, 1.0, 1.1, "cut_2",
                     false);
        cuts.push_back(cut1);
        cuts.push_back(cut2);
        break;
      }
      case 3: {
        // cut_1 is scored lower than cut_2, but cut_1 is forced
        CutData cut1(CreateSparseRow({{1, 1.0}}), 1.0, 1, 1, 1.0, 1.0, "cut_1",
                     true);
        CutData cut2(CreateSparseRow({{1, 1.0}}), 1.0, 1, 1, 1.0, 1.1, "cut_2",
                     false);
        cuts.push_back(cut1);
        cuts.push_back(cut2);
        break;
      }
      case 4: {
        // cut_1 is orthogonal to cut_2
        CutData cut1(CreateSparseRow({{1, 1.0}}), 1.0, 1, 1, 1.0, 1.1, "cut_1",
                     false);
        CutData cut2(CreateSparseRow({{0, 1.0}}), 1.0, 1, 1, 1.0, 1.0, "cut_2",
                     false);
        cuts.push_back(cut1);
        cuts.push_back(cut2);
        break;
      }
      case 5: {
        // cut_1 is orthogonal to cut_2, but cut_2 is scored below threshold
        CutData cut1(CreateSparseRow({{1, 1.0}}), 1.0, 1, 1, 1.0, 1.0, "cut_1",
                     false);
        CutData cut2(CreateSparseRow({{0, 1.0}}), 1.0, 1, 0, 0.0, 0.0, "cut_2",
                     false);
        cuts.push_back(cut1);
        cuts.push_back(cut2);
        break;
      }
      case 6: {
        // cut_1 is orthogonal to cut_2, but both are scored below threshold
        CutData cut1(CreateSparseRow({{1, 1.0}}), 1.0, 1, 0, 0.0, 0.0, "cut_1",
                     false);
        CutData cut2(CreateSparseRow({{0, 1.0}}), 1.0, 1, 0, 0.0, 0.0, "cut_2",
                     false);
        cuts.push_back(cut1);
        cuts.push_back(cut2);
        break;
      }
      case 7: {
        // cut_1 is orthogonal to 2 cut, one that is parallel to cut_2
        CutData cut1(CreateSparseRow({{1, 1.0}}), 1.0, 1, 1, 1.0, 1.1, "cut_1",
                     false);
        CutData cut2(CreateSparseRow({{0, 1.0}}), 1.0, 1, 1, 1.0, 1.0, "cut_2",
                     false);
        CutData cut3(CreateSparseRow({{0, 2.0}}), 1.0, 1, 1, 1.0, 0.5, "cut_3",
                     false);
        cuts.push_back(cut1);
        cuts.push_back(cut2);
        cuts.push_back(cut3);
        break;
      }
      case 8: {
        // cut_1 is orthogonal to 2 cut, one that is parallel to cut_2, but is
        // forced
        CutData cut1(CreateSparseRow({{1, 1.0}}), 1.0, 1, 1, 1.0, 0.6, "cut_1",
                     false);
        CutData cut2(CreateSparseRow({{0, 1.0}}), 1.0, 1, 1, 1.0, 1.0, "cut_2",
                     false);
        CutData cut3(CreateSparseRow({{0, 2.0}}), 1.0, 1, 1, 1.0, 0.5, "cut_3",
                     true);
        cuts.push_back(cut1);
        cuts.push_back(cut2);
        cuts.push_back(cut3);
        break;
      }
      case 9: {
        // cut_1 is parallel to cut_2 with cos = -1
        CutData cut1(CreateSparseRow({{1, 1.0}}), 1.0, 1, 1, 1.0, 1.1, "cut_1",
                     false);
        CutData cut2(CreateSparseRow({{1, -1.0}}), 1.0, 1, 1, 1.0, 1.0, "cut_2",
                     false);
        cuts.push_back(cut1);
        cuts.push_back(cut2);
        break;
      }
      default:
        LOG(FATAL) << "Invalid Cutting Plane Setup: " << switch_variable_index;
    }
    return cuts;
  }

  std::unique_ptr<Selector> selector_;
};

INSTANTIATE_TEST_SUITE_P(IdenticalCuts, MinimalCutSelectorTest,
                         testing::Combine(testing::Values(0, 1),
                                          testing::Range(0, 3)));

INSTANTIATE_TEST_SUITE_P(OrthogonalCuts, MinimalCutSelectorTest,
                         testing::Combine(testing::Values(0, 1),
                                          testing::Range(4, 8)));

INSTANTIATE_TEST_SUITE_P(SignedOrthogonalCuts, MinimalCutSelectorTest,
                         testing::Combine(testing::Values(0, 1),
                                          testing::Range(9, 10)));

TEST_P(MinimalCutSelectorTest, CutSelectorTest) {
  // Call the Create function to create a Solver object
  ASSERT_OK_AND_ASSIGN(std::unique_ptr<Solver> solver,
                       Solver::Create(MiniMipParameters{}, MiniMipProblem{}));

  std::vector<CutData> cuts = CreateCuttingplanes(std::get<1>(GetParam()));
  absl::StatusOr<std::vector<CutData>> selected_cuts =
      selector_->SelectCuttingPlanes(*solver, cuts);

  if (std::get<1>(GetParam()) == 0) {
    ASSERT_EQ(selected_cuts.value().size(), 1);
    ASSERT_EQ(selected_cuts.value()[0].name(), "cut_1");
  } else if (std::get<1>(GetParam()) == 1) {
    ASSERT_EQ(selected_cuts.value().size(), 2);
    ASSERT_EQ(selected_cuts.value()[0].name(), "cut_1");
    ASSERT_EQ(selected_cuts.value()[1].name(), "cut_2");
  } else if (std::get<1>(GetParam()) == 2) {
    ASSERT_EQ(selected_cuts.value().size(), 1);
    ASSERT_EQ(selected_cuts.value()[0].name(), "cut_2");
  } else if (std::get<1>(GetParam()) == 3) {
    ASSERT_EQ(selected_cuts.value().size(), 2);
    ASSERT_EQ(selected_cuts.value()[0].name(), "cut_2");
    ASSERT_EQ(selected_cuts.value()[1].name(), "cut_1");
  } else if (std::get<1>(GetParam()) == 4) {
    ASSERT_EQ(selected_cuts.value().size(), 2);
    ASSERT_EQ(selected_cuts.value()[0].name(), "cut_1");
    ASSERT_EQ(selected_cuts.value()[1].name(), "cut_2");
  } else if (std::get<1>(GetParam()) == 5) {
    ASSERT_EQ(selected_cuts.value().size(), 1);
    ASSERT_EQ(selected_cuts.value()[0].name(), "cut_1");
  } else if (std::get<1>(GetParam()) == 6) {
    ASSERT_EQ(selected_cuts.value().size(), 0);
  } else if (std::get<1>(GetParam()) == 7) {
    ASSERT_EQ(selected_cuts.value().size(), 2);
    ASSERT_EQ(selected_cuts.value()[0].name(), "cut_1");
    ASSERT_EQ(selected_cuts.value()[1].name(), "cut_2");
  } else if (std::get<1>(GetParam()) == 8) {
    ASSERT_EQ(selected_cuts.value().size(), 3);
    ASSERT_EQ(selected_cuts.value()[0].name(), "cut_2");
    ASSERT_EQ(selected_cuts.value()[1].name(), "cut_1");
    ASSERT_EQ(selected_cuts.value()[2].name(), "cut_3");
  } else if (std::get<1>(GetParam()) == 9) {
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
}

}  // namespace
}  // namespace minimip
