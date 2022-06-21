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

#include "src/lp_interface/lpi.h"

#include <gtest/gtest.h>

#include "absl/status/status.h"
#include "src/data_structures/strong_sparse_vector.h"
#include "src/lp_interface/lpi_factory.h"
#include "unit_tests/utils.h"

namespace minimip {

namespace {

class LPInterfaceImplementationTest
    : public ::testing::TestWithParam<LPInterfaceCode> {
 protected:
  void SetUp() override { lpi_ = CreateLPInterface(GetParam()); }

  void PopulateSmallLp() {
    // max 4 x        (objective)
    // 1 <= x <= 2  (linear constraint)
    // 0 <= x <= 3  (bounds)
    ASSERT_OK(lpi_->SetObjectiveSense(/*is_maximization=*/true));
    ASSERT_OK(lpi_->AddColumn({}, 0.0, 3.0, 4.0, "x"));
    SparseRow row({{0, 1.0}});
    ASSERT_OK(lpi_->AddRow(row, 1.0, 2.0, "r1"));
  }

  std::unique_ptr<LPInterface> lpi_;
};

INSTANTIATE_TEST_SUITE_P(All, LPInterfaceImplementationTest,
                         testing::ValuesIn({LPInterfaceCode::kGlop,
                                            LPInterfaceCode::kSoplex}));

TEST_P(LPInterfaceImplementationTest, CreateSimpleLp) {
  PopulateSmallLp();
  EXPECT_TRUE(lpi_->IsMaximization());

  EXPECT_EQ(lpi_->GetNumberOfColumns(), ColIndex(1));
  EXPECT_EQ(lpi_->GetObjectiveCoefficient(ColIndex(0)), 4.0);
  EXPECT_EQ(lpi_->GetLowerBound(ColIndex(0)), 0.0);
  EXPECT_EQ(lpi_->GetUpperBound(ColIndex(0)), 3.0);

  EXPECT_EQ(lpi_->GetNumberOfRows(), RowIndex(1));
  EXPECT_EQ(lpi_->GetLeftHandSide(RowIndex(0)), 1.0);
  EXPECT_EQ(lpi_->GetRightHandSide(RowIndex(0)), 2.0);

  EXPECT_EQ(lpi_->GetMatrixCoefficient(ColIndex(0), RowIndex(0)), 1.0);
}

// TODO(lpawel, cgraczy): Add missing tests from `bases.cc`, `change.cc`,
// `bounds_changes.cc`, `solve.cc` (or put them in seperate unit test file with
// a descriptive name).

}  // namespace

}  // namespace minimip
