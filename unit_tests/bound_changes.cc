
/**@file   boundchg.c
 * @brief  unit test for checking bound changes
 *
 * We perform two tests:
 * - We change the bounds by a very small value and check whether this has an effect on the LP solver interface.
 * - We fix a variable to infinity.
 *
 * In both cases it is unclear what happens. Also LP-solvers might react differently.
 *
 * These tests can be used for debugging or checking the behavior of LP-solvers.
 */

#include "src/lp_interface/lpi_factory.h"
#include <gtest/gtest.h>

#define DEF_INTERFACE 1 /** 0 = Glop Interface (Default),
                          * 1 = SoPlex Interface, **/

namespace minimip {
/*** TEST SUITE SIMPLE ***/

static LPInterface* lp_interface_ = nullptr;

class BoundChanges : public ::testing::Test {
 protected:
  LPValueArray obj, lb, ub, lbnew, ubnew, empty_vals;
  LPIndexArray ind, empty_indices;
  StringArray empty_names;

  void SetUp() override {
    /* build interface factory */
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

    /* add one column */
    ASSERT_EQ(lp_interface_->AddColumns(1, obj, lb, ub, empty_names, 0, empty_indices, empty_indices, empty_vals), RetCode::kOkay);
  }
};

/** TESTS **/
TEST_F(BoundChanges, SimpleBoundTest) {
  lb[0] = 1.0;
  ub[0] = 2.0;

  /* change bounds to some value */
  ASSERT_EQ(lp_interface_->ChangeBounds(1, ind, lb, ub), RetCode::kOkay);

  /* get bounds and compare */
  ASSERT_EQ(lp_interface_->GetBounds(0, 0, lbnew, ubnew), RetCode::kOkay);

  ASSERT_FLOAT_EQ(lb[0], lbnew[0]);
  ASSERT_FLOAT_EQ(ub[0], ubnew[0]);
}

TEST_F(BoundChanges, ChangeBoundBySmallValue) {
  /* change bound to small value */
  lb[0] = 1e-11;
  ub[0] = 1.0 - 1e-11;
  ASSERT_EQ(lp_interface_->ChangeBounds(1, ind, lb, ub), RetCode::kOkay);

  /* get bounds and compare */
  ASSERT_EQ(lp_interface_->GetBounds(0, 0, lbnew, ubnew), RetCode::kOkay);

  ASSERT_FLOAT_EQ(lb[0], lbnew[0]);
  ASSERT_FLOAT_EQ(ub[0], ubnew[0]);
}

TEST_F(BoundChanges, FixToInfinity) {
  RetCode retcode;

  /* try to fix variables to infinity */
  lb[0] = lp_interface_->Infinity();
  ub[0] = lp_interface_->Infinity();

  /* calling should return an LPERROR */
  retcode = lp_interface_->ChangeBounds(1, ind, lb, ub);

  ASSERT_EQ(retcode, RetCode::kLPError) << "Fixing variables to infinity does not return an error.";
}

} /*namespace minimip */

