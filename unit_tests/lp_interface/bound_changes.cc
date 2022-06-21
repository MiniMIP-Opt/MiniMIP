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

// Unit test for checking bound changes
//
// We perform two tests:
// - We change the bounds by a very small value and check whether this has an
// effect on the LP solver interface.
// - We fix a variable to infinity.
//
// In both cases it is unclear what happens. Also LP-solvers might react
// differently.
//
// These tests can be used for debugging or checking the behavior of LP-solvers.

#include <gtest/gtest.h>

#include "absl/status/status.h"
#include "src/lp_interface/lpi_factory.h"
#include "unit_tests/utils.h"

#define DEF_INTERFACE \
  0  // 0 = Glop Interface (Default),
     // 1 = SoPlex Interface,

namespace minimip {

static LPInterface* lp_interface_ = nullptr;

class BoundChanges : public ::testing::Test {
 protected:
  double objective_coefficients_, lower_bound_, upper_bound_, new_lower_bound_,
      new_upper_bound_;

  void SetUp() override {
    // build interface factory
    auto* interface_factory = new LPInterfaceFactory();
    InterfaceCode interface_code;
    switch (DEF_INTERFACE) {
      case 1:
        interface_code = InterfaceCode::kSoplex;
        break;
      default:
        interface_code = InterfaceCode::kGlop;
        break;
    }
    lp_interface_ = interface_factory->CreateLPInterface(interface_code);
    ASSERT_OK(
        lp_interface_->SetObjectiveSense(LPObjectiveSense::kMaximization));

    objective_coefficients_ = 1.0;
    lower_bound_            = 0.0;
    upper_bound_            = 1.0;

    SparseVector empty_coefficients = {{}, {}};

    // add one empty column
    ASSERT_OK(lp_interface_->AddColumn(empty_coefficients, lower_bound_,
                                       upper_bound_, objective_coefficients_,
                                       "x1"));
  }
};

// TESTS
TEST_F(BoundChanges, SimpleBoundTest) {
  lower_bound_ = 1.0;
  upper_bound_ = 2.0;

  // change bounds to some value
  ASSERT_OK(lp_interface_->SetColumnBounds(0, lower_bound_, upper_bound_));

  // get bounds and compare
  new_lower_bound_ = lp_interface_->GetLowerBound(0);
  new_upper_bound_ = lp_interface_->GetUpperBound(0);

  ASSERT_FLOAT_EQ(lower_bound_, new_lower_bound_);
  ASSERT_FLOAT_EQ(upper_bound_, new_upper_bound_);
}

TEST_F(BoundChanges, ChangeBoundBySmallValue) {
  // change bound by small value
  lower_bound_ = 1e-11;
  upper_bound_ = 1.0 - 1e-11;
  ASSERT_OK(lp_interface_->SetColumnBounds(0, lower_bound_, upper_bound_));

  // get bounds and compare
  new_lower_bound_ = lp_interface_->GetLowerBound(0);
  new_upper_bound_ = lp_interface_->GetUpperBound(0);

  ASSERT_FLOAT_EQ(lower_bound_, new_lower_bound_);
  ASSERT_FLOAT_EQ(upper_bound_, new_upper_bound_);
}

TEST_F(BoundChanges, FixToInfinity) {
  absl::Status absl_code;

  // try to fix variables to infinity
  lower_bound_ = lp_interface_->Infinity();
  upper_bound_ = lp_interface_->Infinity();

  // calling should return an LPERROR
  ASSERT_EQ(lp_interface_->SetColumnBounds(0, lower_bound_, upper_bound_),
            absl::Status(absl::StatusCode::kInternal, "LP Error"));
}

}  // namespace minimip
