
// @file   boundchg.c
// @brief  unit test for checking bound changes
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

#define DEF_INTERFACE \
  1  // 0 = Glop Interface (Default),
     // 1 = SoPlex Interface,

namespace minimip {
// TEST SUITE SIMPLE

static LPInterface* lp_interface_ = nullptr;

class BoundChanges : public ::testing::Test {
 protected:
  std::vector<double> obj, lb, ub, lbnew, ubnew, empty_vals;
  std::vector<int> ind, empty_indices;
  std::vector<std::string> empty_names;

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
    lp_interface_->ChangeObjectiveSense(LPObjectiveSense::kMaximize);

    obj.push_back(1.0);
    lb.push_back(0.0);
    ub.push_back(1.0);
    ind.push_back(0);
    lbnew.reserve(1);
    ubnew.reserve(1);

    // add one column
    ASSERT_EQ(
        lp_interface_->AddColumns(1, obj, lb, ub, empty_names, 0, empty_indices,
                                  empty_indices, empty_vals),
        absl::OkStatus());
  }
};

// TESTS
TEST_F(BoundChanges, SimpleBoundTest) {
  lb[0] = 1.0;
  ub[0] = 2.0;

  // change bounds to some value
  ASSERT_EQ(lp_interface_->SetColumnBounds(1, ind, lb, ub), absl::OkStatus());

  // get bounds and compare
  ASSERT_EQ(lp_interface_->GetBounds(0, 0, lbnew, ubnew), absl::OkStatus());

  ASSERT_FLOAT_EQ(lb[0], lbnew[0]);
  ASSERT_FLOAT_EQ(ub[0], ubnew[0]);
}

TEST_F(BoundChanges, ChangeBoundBySmallValue) {
  // change bound to small value
  lb[0] = 1e-11;
  ub[0] = 1.0 - 1e-11;
  ASSERT_EQ(lp_interface_->SetColumnBounds(1, ind, lb, ub), absl::OkStatus());

  // get bounds and compare
  ASSERT_EQ(lp_interface_->GetBounds(0, 0, lbnew, ubnew), absl::OkStatus());

  ASSERT_FLOAT_EQ(lb[0], lbnew[0]);
  ASSERT_FLOAT_EQ(ub[0], ubnew[0]);
}

TEST_F(BoundChanges, FixToInfinity) {
  absl::Status absl_code;

  // try to fix variables to infinity
  lb[0] = lp_interface_->Infinity();
  ub[0] = lp_interface_->Infinity();

  // calling should return an LPERROR
  absl_code = lp_interface_->SetColumnBounds(1, ind, lb, ub);

  ASSERT_EQ(absl_code, absl::Status(absl::StatusCode::kInternal, "LP Error"))
      << "Fixing variables to infinity does not return an error.";
}

}  // namespace minimip
