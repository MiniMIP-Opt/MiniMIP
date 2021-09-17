
/**@file   solve.c
 * @brief  unit test for checking lpi solution
 * @author Marc Pfetsch
 *
 * We perform tests with solving several examples. These are inspired by the unit tests of OSI in COIN-OR.
 */

#include <tuple>

#include "src/lp_interface/lpi_factory.h"
#include <gtest/gtest.h>

#define EPS 1e-6
#define DEF_INTERFACE 1 /** 0 = Glop Interface (Default),
                          * 1 = SoPlex Interface, **/

namespace minimip {

/** expected feasibility status for primal or dual problem */
enum LPFeasibilityStat {
  FEASIBLE = 0,  /**< the problem is feasible */
  UNBOUNDED = 1, /**< the problem is unbounded */
  INFEASIBLE = 2 /**< the problem is infeasible */
};
typedef enum LPFeasibilityStat LPFeasibilityStat;

/*** TEST SUITE SOLVE ***/

static LPInterface* lp_interface_ = nullptr;

class Solve : public ::testing::Test {
 protected:
  /* perform test */
  LPValueArray obj;          /**< objective function values of columns */
  LPValueArray lb;           /**< lower bounds of columns */
  LPValueArray ub;           /**< upper bounds of columns */
  LPValueArray lhs;          /**< left hand sides of rows */
  LPValueArray rhs;          /**< right hand sides of rows */
  LPIndexArray beg;          /**< start index of each column in ind- and val-array */
  LPIndexArray ind;          /**< row indices of constraint matrix entries */
  LPValueArray val;          /**< values of constraint matrix entries */
  LPValueArray exp_primsol;  /**< expected primal optimal solution or primal ray if primal is unbounded or NULL */
  LPValueArray exp_dualsol;  /**< expected dual optimal solution or dual ray if dual is unbounded or NULL */
  LPValueArray exp_activity; /**< expected activity of optimal solution or NULL */
  LPValueArray exp_redcost;  /**< expected reduced cost of optimal solution or NULL */

  void SetUp() override {
    /* build interface factory */
    auto* interface_factory = new LPInterfaceFactory();
    InterfaceCode interface_code;
    switch (DEF_INTERFACE) {
      case 1:
        interface_code = InterfaceCode::SOPLEX;
        break;
      default:
        interface_code = InterfaceCode::GLOP;
        break;
    }
    lp_interface_ = interface_factory->CreateLPInterface(interface_code);
    lp_interface_->ChangeObjectiveSense(LPObjectiveSense::OBJ_SENSE_MAXIMIZE);
  }
  /* local functions */

  /** solve problem */
  static void solveTest(
    bool solveprimal,                 /**< use primal simplex */
    LPNum ncols,                      /**< number of columns */
    LPNum nrows,                      /**< number of rows */
    LPFeasibilityStat exp_primalfeas, /**< expected primal feasibility status */
    LPFeasibilityStat exp_dualfeas,   /**< expected primal feasibility status */
    const LPValueArray& exp_primsol,  /**< expected primal optimal solution or primal ray if primal is unbounded or NULL */
    const LPValueArray& exp_dualsol,  /**< expected dual optimal solution or dual ray if dual is unbounded or NULL */
    const LPValueArray& exp_activity, /**< expected activity of optimal solution or NULL */
    const LPValueArray& exp_redcost   /**< expected reduced cost of optimal solution or NULL */
  ) {
    /* solution data */
    LPValue objval;
    LPValueArray primsol;
    LPValueArray dualsol;
    LPValueArray activity;
    LPValueArray redcost;

    /* auxiliary data */
    bool primalfeasible;
    bool dualfeasible;
    LPNum ntmprows;
    LPNum ntmpcols;
    LPNum i;
    LPNum j;

    /* check size */
    ntmprows = lp_interface_->GetNumberOfRows();
    ntmpcols = lp_interface_->GetNumberOfColumns();
    ASSERT_EQ(nrows, ntmprows);
    ASSERT_EQ(ncols, ntmpcols);

    /* solve problem */
    if (solveprimal) {
      ASSERT_EQ(lp_interface_->SolvePrimal(), RetCode::OKAY);
    } else {
      ASSERT_EQ(lp_interface_->SolveDual(), RetCode::OKAY);
    }

    /* check status */
    ASSERT_TRUE(lp_interface_->WasSolved());
    ASSERT_TRUE(!lp_interface_->IsObjectiveLimitExceeded());
    ASSERT_TRUE(!lp_interface_->IsIterationLimitExceeded());
    ASSERT_TRUE(!lp_interface_->IsTimeLimitExceeded());

    /* check feasibility status */
    primalfeasible = lp_interface_->IsPrimalFeasible();
    dualfeasible = lp_interface_->IsDualFeasible();

    /* if we are feasible, we should be optimal */
    if (exp_primalfeas == LPFeasibilityStat::FEASIBLE && exp_dualfeas == LPFeasibilityStat::FEASIBLE) {
      ASSERT_TRUE(lp_interface_->IsOptimal());
    }

    /* check more primal statuses */
    switch (exp_primalfeas) {
      case LPFeasibilityStat::FEASIBLE:
        ASSERT_TRUE(primalfeasible);
        ASSERT_TRUE(!lp_interface_->ExistsPrimalRay());
        ASSERT_TRUE(!lp_interface_->HasPrimalRay());
        ASSERT_TRUE(!lp_interface_->IsPrimalUnbounded());
        ASSERT_TRUE(!lp_interface_->IsPrimalInfeasible());
        ASSERT_TRUE(lp_interface_->IsPrimalFeasible());
        break;

      case LPFeasibilityStat::UNBOUNDED:
        /* Because of SoPlex, cannot always determine feasibility status here, even if we want to apply the primal
         * simplex. In any case, the results of primalfeasible and IsPrimalFeasible(lpi) should coincide. */
        ASSERT_EQ(primalfeasible, lp_interface_->IsPrimalFeasible());

        /* It seems that we cannot guarantee that the primal is shown to be unbounded. */
        /* cr_assert( IsPrimalUnbounded(lpi) ); */

        /* primal ray should exist if the primal simplex ran */
        ASSERT_TRUE(!solveprimal || lp_interface_->ExistsPrimalRay());
        ASSERT_TRUE(!lp_interface_->IsPrimalInfeasible());
        break;

      case LPFeasibilityStat::INFEASIBLE:
        ASSERT_TRUE(!primalfeasible);
        /* It seems that we cannot always prove that primal is infeasible. */
        /* cr_assert( IsPrimalInfeasible(lpi) ); */

        /* It seems that we cannot always prove that primal is not unbounded. */
        /* cr_assert( ! IsPrimalUnbounded(lpi) ); */
        ASSERT_TRUE(!lp_interface_->IsPrimalFeasible());
        break;

      default:
        abort();
    }

    /* check more dual statuses */
    switch (exp_dualfeas) {
      case LPFeasibilityStat::FEASIBLE:
        ASSERT_TRUE(dualfeasible);
        ASSERT_TRUE(!lp_interface_->ExistsDualRay());
        ASSERT_TRUE(!lp_interface_->HasDualRay());
        ASSERT_TRUE(!lp_interface_->IsDualUnbounded());
        ASSERT_TRUE(!lp_interface_->IsDualInfeasible());
        ASSERT_TRUE(lp_interface_->IsDualFeasible());
        break;

      case LPFeasibilityStat::UNBOUNDED:
        /* Because of SoPlex, cannot always determine feasibility status here, even if we want to apply the dual
         * simplex. In any case, the results of dualfeasible and IsDualFeasible(lpi) should coincide. */
        ASSERT_EQ(dualfeasible, lp_interface_->IsDualFeasible());

        /* It seems that we cannot guarantee that the dual is shown to be unbounded. */
        /* cr_assert( IsDualUnbounded(lpi) ); */

        /* dual ray should exist if the dual simplex ran */
        ASSERT_TRUE(solveprimal || lp_interface_->ExistsDualRay());
        ASSERT_TRUE(!lp_interface_->IsDualInfeasible());
        break;

      case LPFeasibilityStat::INFEASIBLE:
        ASSERT_TRUE(!dualfeasible);
        ASSERT_TRUE(!lp_interface_->IsDualUnbounded());
        /* It seems that we cannot always prove that dual is infeasible. */
        /* cr_assert( IsDualInfeasible(lpi) ); */
        ASSERT_TRUE(!lp_interface_->IsDualFeasible());
        break;

      default:
        abort();
    }

    primsol.reserve(ncols);
    dualsol.reserve(nrows);
    activity.reserve(nrows);
    redcost.reserve(ncols);

    /* check solution */
    if (exp_primalfeas == LPFeasibilityStat::FEASIBLE) {
      /* get solution */
      ASSERT_EQ(lp_interface_->GetSolution(objval, primsol, dualsol, activity, redcost), RetCode::OKAY);

      for (j = 0; j < ncols; ++j) {
        ASSERT_FLOAT_EQ(primsol[j], exp_primsol[j]);// EPS, "Violation of primal solution %d: %g != %g\n", j, primsol[j], exp_primsol[j]);
        ASSERT_FLOAT_EQ(redcost[j], exp_redcost[j]);// EPS, "Violation of reduced cost of solution %d: %g != %g\n", j, redcost[j], exp_redcost[j]);
      }
    } else if (exp_primalfeas == LPFeasibilityStat::UNBOUNDED) {

      if (lp_interface_->HasPrimalRay()) {
        LPValue scalingfactor = 1.0;

        ASSERT_EQ(lp_interface_->GetPrimalRay(primsol), RetCode::OKAY);

        /* loop until scaling factor can be determined */
        for (j = 0; j < ncols; ++j) {
          if (REALABS(exp_primsol[j]) < EPS)
            ASSERT_FLOAT_EQ(primsol[j], exp_primsol[j]);// EPS, "Violation of primal ray %d: %g != %g\n", j, primsol[j], exp_primsol[j]);
          else {
            scalingfactor = primsol[j] / exp_primsol[j];
            break;
          }
        }

        /* again loop over ray */
        for (j = 0; j < ncols; ++j) {
          ASSERT_FLOAT_EQ(primsol[j], scalingfactor * exp_primsol[j]);// EPS, "Violation of primal ray %d: %g != %g\n", j, primsol[j], scalingfactor * exp_primsol[j]);
        }
      }
    }

    if (exp_dualfeas == LPFeasibilityStat::FEASIBLE) {
      /* get solution */
      ASSERT_EQ(lp_interface_->GetSolution(objval, primsol, dualsol, activity, redcost), RetCode::OKAY);

      for (i = 0; i < nrows; ++i) {
        ASSERT_FLOAT_EQ(dualsol[i], exp_dualsol[i]);  // EPS, "Violation of dual solution %d: %g != %g\n", i, dualsol[i], exp_dualsol[i]);
        ASSERT_FLOAT_EQ(activity[i], exp_activity[i]);// EPS, "Violation of activity of solution %d: %g != %g\n", i, activity[i], exp_activity[i]);
      }
    } else if (exp_dualfeas == LPFeasibilityStat::UNBOUNDED) {

      if (lp_interface_->HasDualRay()) {
        LPValue scalingfactor = 1.0;
        LPValueArray lhs(nrows);
        LPValueArray rhs(nrows);

        /* get lhs/rhs for check of dual ray */
        ASSERT_EQ(lp_interface_->GetSides(0, nrows - 1, lhs, rhs), RetCode::OKAY);

        /* get dual ray */
        ASSERT_EQ(lp_interface_->GetDualFarkasMultiplier(dualsol), RetCode::OKAY);

        /* loop until scaling factor can be determined */
        for (i = 0; i < nrows; ++i) {
          if (REALABS(exp_dualsol[i]) < EPS)
            ASSERT_FLOAT_EQ(dualsol[i], exp_dualsol[i]);// EPS, "Violation of dual ray %d: %g != %g\n", i, dualsol[i], exp_dualsol[i]);
          else {
            scalingfactor = dualsol[i] / exp_dualsol[i];
            break;
          }
        }

        /* again loop over ray */
        for (i = 0; i < nrows; ++i) {
          ASSERT_FLOAT_EQ(dualsol[i], scalingfactor * exp_dualsol[i]);// EPS, "Violation of dual ray %d: %g != %g\n", i, dualsol[i], scalingfactor * exp_dualsol[i]);
          ASSERT_TRUE(!lp_interface_->IsInfinity(-lhs[i]) || dualsol[i] <= -EPS);
          ASSERT_TRUE(!lp_interface_->IsInfinity(rhs[i]) || dualsol[i] >= EPS);
        }
      }
    }
  }

  /** perform basic test for the given problem */
  static void performTest(
    bool solveprimal,                 /**< use primal simplex */
    LPObjectiveSense objsen,          /**< objective sense */
    LPNum ncols,                      /**< number of columns */
    const LPValueArray& obj,          /**< objective function values of columns */
    const LPValueArray& lb,           /**< lower bounds of columns */
    const LPValueArray& ub,           /**< upper bounds of columns */
    LPNum nrows,                      /**< number of rows */
    const LPValueArray& lhs,          /**< left hand sides of rows */
    const LPValueArray& rhs,          /**< right hand sides of rows */
    const LPIndexArray& beg,          /**< start index of each column in ind- and val-array */
    const LPIndexArray& ind,          /**< row indices of constraint matrix entries */
    const LPValueArray& val,          /**< values of constraint matrix entries */
    LPFeasibilityStat exp_primalfeas, /**< expected primal feasibility status */
    LPFeasibilityStat exp_dualfeas,   /**< expected primal feasibility status */
    const LPValueArray& exp_primsol,  /**< expected primal optimal solution or primal ray if primal is unbounded or NULL */
    const LPValueArray& exp_dualsol,  /**< expected dual optimal solution or dual ray if dual is unbounded or NULL */
    const LPValueArray& exp_activity, /**< expected activity of optimal solution or NULL */
    const LPValueArray& exp_redcost   /**< expected reduced cost of optimal solution or NULL */
  ) {
    StringArray empty_names;

    /* load problem */
    ASSERT_EQ(lp_interface_->LoadColumnLP(objsen, 2, obj, lb, ub, empty_names, 2, lhs, rhs, empty_names, 4, beg, ind, val), RetCode::OKAY);
    ASSERT_TRUE(!lp_interface_->WasSolved());

    /* solve problem */
    ASSERT_NO_FATAL_FAILURE(solveTest(solveprimal, ncols, nrows, exp_primalfeas, exp_dualfeas, exp_primsol, exp_dualsol, exp_activity, exp_redcost));
  }

  /** check whether data in LP solver aggrees with original data */
  static void checkData(
    LPObjectiveSense objsen, /**< objective sense */
    LPNum ncols,             /**< number of columns */
    const LPValueArray& obj, /**< objective function values of columns */
    const LPValueArray& lb,  /**< lower bounds of columns */
    const LPValueArray& ub,  /**< upper bounds of columns */
    LPNum nrows,             /**< number of rows */
    const LPValueArray& lhs, /**< left hand sides of rows */
    const LPValueArray& rhs, /**< right hand sides of rows */
    LPNum nnonz,             /**< number of nonzero elements in the constraint matrix */
    const LPIndexArray& beg, /**< start index of each column in ind- and val-array */
    const LPIndexArray& ind, /**< row indices of constraint matrix entries */
    const LPValueArray& val  /**< values of constraint matrix entries */
  ) {
    LPObjectiveSense check_objsen;
    LPValueArray check_val;
    LPValueArray check_lb;
    LPValueArray check_ub;
    LPValueArray check_obj;
    LPValueArray check_lhs;
    LPValueArray check_rhs;
    LPIndexArray check_beg;
    LPIndexArray check_ind;
    LPNum check_ncols;
    LPNum check_nrows;
    LPNum check_nnonz;
    LPNum check_nnonz2;
    LPNum i;
    LPNum j;

    /* check number of rows and columns */
    check_nrows = lp_interface_->GetNumberOfRows();
    check_ncols = lp_interface_->GetNumberOfColumns();
    ASSERT_EQ(check_nrows, nrows);
    ASSERT_EQ(check_ncols, ncols);

    /* check objective sense */
    check_objsen = lp_interface_->GetObjectiveSense();
    ASSERT_EQ(objsen, check_objsen);

    /* get number of nonzeros in matrix */
    check_nnonz = lp_interface_->GetNumberOfNonZeros();
    ASSERT_EQ(check_nnonz, nnonz);

    check_val.reserve(check_nnonz);
    check_lb.reserve(ncols);
    check_ub.reserve(ncols);
    check_obj.reserve(ncols);
    check_beg.reserve(ncols);
    check_ind.reserve(check_nnonz);

    /* get matrix data */
    ASSERT_EQ(lp_interface_->GetColumns(0, ncols - 1, check_lb, check_ub, check_nnonz2, check_beg, check_ind, check_val), RetCode::OKAY);
    ASSERT_EQ(lp_interface_->GetObjective(0, ncols - 1, check_obj), RetCode::OKAY);

    /* compare data */
    for (j = 0; j < ncols; ++j) {
      if ( fabs(check_lb[j]) < 1e30 && fabs(lb[j]) < 1e30 ){
        ASSERT_FLOAT_EQ(check_lb[j], lb[j]);
      }
      if ( fabs(check_ub[j]) < 1e30 && fabs(ub[j]) < 1e30 ){
        ASSERT_FLOAT_EQ(check_ub[j], ub[j]);
      }

      ASSERT_FLOAT_EQ(check_obj[j], obj[j]);// EPS, "Violation of objective coefficient %d: %g != %g\n", j, check_obj[j], obj[j]);

      ASSERT_EQ(check_beg[j], beg[j]);
    }

    /* compare matrix */
    for (j = 0; j < nnonz; ++j) {
      ASSERT_EQ(check_ind[j], ind[j]);
      ASSERT_FLOAT_EQ(check_val[j], val[j]);// EPS, "Violation of matrix entry (%d, %d): %g != %g\n", ind[j], j, check_val[j], val[j]);
    }

    check_lhs.reserve(nrows);
    check_rhs.reserve(nrows);

    ASSERT_EQ(lp_interface_->GetSides(0, nrows - 1, check_lhs, check_rhs), RetCode::OKAY);

    for (i = 0; i < nrows; ++i) {

      if ( fabs(check_lhs[i]) < 1e30 && fabs(lhs[i]) < 1e30 ){
        ASSERT_FLOAT_EQ(check_lhs[i], lhs[i]);
      }
      if ( fabs(check_lhs[i]) < 1e30 && fabs(lhs[i]) < 1e30 ){
        ASSERT_FLOAT_EQ(check_rhs[i], rhs[i]);
      }
    }
  }
};

/** TESTS **/

/** Test 1
 *
 * max 3 x1 +   x2
 *     2 x1 +   x2 <= 10
 *       x1 + 3 x2 <= 15
 *       x1,    x2 >= 0
 *
 * with primal optimal solution (5, 0), dual optimal solution (1.5, 0), activity (10, 5), and redcost (0, -0.5).
 *
 * Then use objective (1, 1) with primal optimal solution (3,4), dual optimal solution (0.4, 0.2), activity (10, 15), and redcost (0, 0).
 */
TEST_F(Solve, test1) {
  /* data to be filled */
  obj.reserve(2);
  lb.reserve(2);
  rhs.reserve(2);
  beg.reserve(2);
  ind.reserve(4);
  val.reserve(4);

  exp_primsol.reserve(2);
  exp_dualsol.reserve(2);
  exp_activity.reserve(2);
  exp_redcost.reserve(2);

  ub.reserve(2);
  lhs.reserve(2);

  /* data with fixed values: */
  obj = {3.0, 1.0};
  lb = {0.0, 0.0};
  rhs = {10.0, 15.0};
  beg = {0, 2};
  ind = {0, 1, 0, 1};
  val = {2.0, 1.0, 1.0, 3.0};

  /* expected solutions */
  exp_primsol = {5.0, 0.0};
  exp_dualsol = {1.5, 0.0};
  exp_activity = {10.0, 5.0};
  exp_redcost = {0.0, -0.5};

  /* fill variable data */
  ub[0] = lp_interface_->Infinity();
  ub[1] = lp_interface_->Infinity();
  for(int j = 0; j < 2; j++){
  ASSERT_TRUE(lp_interface_->IsInfinity(ub[j]));
  }
  lhs[0] = -lp_interface_->Infinity();
  lhs[1] = -lp_interface_->Infinity();
  for(int j = 0; j < 2; j++){
    ASSERT_EQ(lhs[j], -lp_interface_->Infinity());
  }
  /* solve problem with primal simplex */
  ASSERT_NO_FATAL_FAILURE(performTest(true, LPObjectiveSense::OBJ_SENSE_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, beg, ind, val,
              LPFeasibilityStat::FEASIBLE, LPFeasibilityStat::FEASIBLE, exp_primsol, exp_dualsol, exp_activity, exp_redcost));

  /* check that data stored in lpi is still the same */
  ASSERT_NO_FATAL_FAILURE(checkData(LPObjectiveSense::OBJ_SENSE_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val));

  /* clear basis status */
  ASSERT_EQ(lp_interface_->ClearState(), RetCode::OKAY);

  /* solve problem with dual simplex */
  ASSERT_NO_FATAL_FAILURE(solveTest(false, 2, 2, LPFeasibilityStat::FEASIBLE, LPFeasibilityStat::FEASIBLE, exp_primsol, exp_dualsol, exp_activity, exp_redcost));

  /* check that data stored in lpi is still the same */
  ASSERT_NO_FATAL_FAILURE(checkData(LPObjectiveSense::OBJ_SENSE_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val));

  /* clear basis status */
  ASSERT_EQ(lp_interface_->ClearState(), RetCode::OKAY);

  /* change objective */
  obj[0] = 1.0;
  ASSERT_EQ(lp_interface_->ChangeObjective(1, ind, obj), RetCode::OKAY);

  /* change expected solution */
  exp_primsol[0] = 3;
  exp_primsol[1] = 4;
  exp_dualsol[0] = 0.4;
  exp_dualsol[1] = 0.2;
  exp_activity[0] = 10;
  exp_activity[1] = 15;
  exp_redcost[1] = 0;

  /* check changed problem with primal simplex */
  ASSERT_NO_FATAL_FAILURE(performTest(true, LPObjectiveSense::OBJ_SENSE_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, beg, ind, val,
              LPFeasibilityStat::FEASIBLE, LPFeasibilityStat::FEASIBLE, exp_primsol, exp_dualsol, exp_activity, exp_redcost));

  /* check that data stored in lpi is still the same */
  ASSERT_NO_FATAL_FAILURE(checkData(LPObjectiveSense::OBJ_SENSE_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val));
}

/** Test 2
 *
 * max 3 x1 +   x2
 *     2 x1 +   x2 <= 10
 *       x1 + 3 x2 <= 15
 *       x1, x2 free
 *
 * which is unbounded (the only difference to Test 1 is that the variables are free).
 *
 * Then use objective (1, 1) with primal optimal solution (3,4), dual optimal solution (0.4, 0.2), activity (10, 15), and redcost (0, 0).
 */
TEST_F(Solve, test2) {
  /* data to be filled */
  obj.reserve(2);
  rhs.reserve(2);
  beg.reserve(2);
  ind.reserve(4);
  val.reserve(4);

  LPValueArray exp_primray(2);
  exp_primsol.reserve(2);
  exp_dualsol.reserve(2);
  exp_activity.reserve(2);
  exp_redcost.reserve(2);

  lb.reserve(2);
  ub.reserve(2);
  lhs.reserve(2);

  /* data with fixed values: */
  obj = {3, 1};
  rhs = {10, 15};
  beg = {0, 2};
  ind = {0, 1, 0, 1};
  val = {2, 1, 1, 3};

  /* expected ray for first run */
  exp_primray = {0.5, -1};

  /* expected solutions */
  exp_primsol = {3, 4};
  exp_dualsol = {0.4, 0.2};
  exp_activity = {10, 15};
  exp_redcost = {0, 0};

  /* fill variable data */
  lb[0] = -lp_interface_->Infinity();
  lb[1] = -lp_interface_->Infinity();
  ub[0] = lp_interface_->Infinity();
  ub[1] = lp_interface_->Infinity();
  lhs[0] = -lp_interface_->Infinity();
  lhs[1] = -lp_interface_->Infinity();

  /* empty_placeholders */
  LPValueArray empty_vals;

  /* solve problem with primal simplex */
  ASSERT_NO_FATAL_FAILURE(performTest(true, LPObjectiveSense::OBJ_SENSE_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, beg, ind, val,
              LPFeasibilityStat::UNBOUNDED, LPFeasibilityStat::INFEASIBLE, exp_primray, empty_vals, empty_vals, empty_vals));

  /* check that data stored in lpi is still the same */
  ASSERT_NO_FATAL_FAILURE(checkData(LPObjectiveSense::OBJ_SENSE_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val));

  /* clear basis status */
  ASSERT_EQ(lp_interface_->ClearState(), RetCode::OKAY);

  /* solve problem with dual simplex */
  ASSERT_NO_FATAL_FAILURE(solveTest(false, 2, 2, LPFeasibilityStat::UNBOUNDED, LPFeasibilityStat::INFEASIBLE, exp_primray, empty_vals, empty_vals, empty_vals));

  /* check that data stored in lpi is still the same */
  ASSERT_NO_FATAL_FAILURE(checkData(LPObjectiveSense::OBJ_SENSE_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val));

  /* clear basis status */
  ASSERT_EQ(lp_interface_->ClearState(), RetCode::OKAY);

  /* change objective */
  obj[0] = 1.0;
  ASSERT_EQ(lp_interface_->ChangeObjective(1, ind, obj), RetCode::OKAY);

  /* solve with primal simplex */
  ASSERT_NO_FATAL_FAILURE(solveTest(true, 2, 2, LPFeasibilityStat::FEASIBLE, LPFeasibilityStat::FEASIBLE, exp_primsol, exp_dualsol, exp_activity, exp_redcost));

  /* check that data stored in lpi is still the same */
  ASSERT_NO_FATAL_FAILURE(checkData(LPObjectiveSense::OBJ_SENSE_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val));
}

/** Test 3
 *
 * min 10 y1 + 15 y2
 *      2 y1 +   y2 == 3
 *        y1 + 3 y2 == 1
 *        y1,    y2 >= 0
 *
 * which is dual unbounded (this is the dual of the problem in Test 2).
 *
 * Then use rhs (1, 1) with primal optimal solution (0.4,0.2), dual optimal solution (3, 4), activity (0, 0), and redcost (0, 0).
 */
TEST_F(Solve, test3) {
  /* data to be filled */
  obj.reserve(2);
  rhs.reserve(2);
  beg.reserve(2);
  ind.reserve(4);
  val.reserve(4);
  LPValueArray exp_dualray(2);
  exp_primsol.reserve(2);
  exp_dualsol.reserve(2);
  exp_activity.reserve(2);
  exp_redcost.reserve(2);

  lb.reserve(2);
  ub.reserve(2);
  lhs.reserve(2);

  /* data with fixed values: */
  obj = {10, 15};
  rhs = {3, 1};
  lhs = {3, 1};
  beg = {0, 2};
  ind = {0, 1, 0, 1};
  val = {2, 1, 1, 3};
  lb = {0, 0};

  /* expected ray */
  exp_dualray = {0.5, -1};

  /* expected solutions */
  exp_primsol = {0.4, 0.2};
  exp_dualsol = {3, 4};
  exp_activity = {1, 1};
  exp_redcost = {0, 0};

  /* fill variable data */
  ub[0] = lp_interface_->Infinity();
  ub[1] = lp_interface_->Infinity();

  /* empty_placeholders */
  LPValueArray empty_vals;

  /* check problem with primal simplex */
  ASSERT_NO_FATAL_FAILURE(performTest(true, LPObjectiveSense::OBJ_SENSE_MINIMIZE, 2, obj, lb, ub, 2, lhs, rhs, beg, ind, val,
              LPFeasibilityStat::INFEASIBLE, LPFeasibilityStat::UNBOUNDED, empty_vals, exp_dualray, empty_vals, empty_vals));

  /* check that data stored in lpi is still the same */
  ASSERT_NO_FATAL_FAILURE(checkData(LPObjectiveSense::OBJ_SENSE_MINIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val));

  /* clear basis status */
  ASSERT_EQ(lp_interface_->ClearState(), RetCode::OKAY);

  /* check problem with dual simplex */
  ASSERT_NO_FATAL_FAILURE(solveTest(false, 2, 2, LPFeasibilityStat::INFEASIBLE, LPFeasibilityStat::UNBOUNDED, empty_vals, exp_dualray, empty_vals, empty_vals));

  /* check that data stored in lpi is still the same */
  ASSERT_NO_FATAL_FAILURE(checkData(LPObjectiveSense::OBJ_SENSE_MINIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val));

  /* clear basis status */
  ASSERT_EQ(lp_interface_->ClearState(), RetCode::OKAY);

  /* change lhs/rhs */
  lhs[0] = 1.0;
  rhs[0] = 1.0;
  ASSERT_EQ(lp_interface_->ChangeSides(1, ind, lhs, rhs), RetCode::OKAY);
  ASSERT_NO_FATAL_FAILURE(solveTest(true, 2, 2, LPFeasibilityStat::FEASIBLE, LPFeasibilityStat::FEASIBLE, exp_primsol, exp_dualsol, exp_activity, exp_redcost));

  /* check that data stored in lpi is still the same */
  ASSERT_NO_FATAL_FAILURE(checkData(LPObjectiveSense::OBJ_SENSE_MINIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val));
}

/** Test 4
 *
 * max x1 + x2
 *     x1 - x2 <= 0
 *   - x1 + x2 <= -1
 *     x1,  x2 free
 *
 * which primal and dual infeasible.
 */
TEST_F(Solve, test4) {
  /* data to be filled */
  obj.reserve(2);
  rhs.reserve(2);
  beg.reserve(2);
  ind.reserve(4);
  val.reserve(4);
  lb.reserve(2);
  ub.reserve(2);
  lhs.reserve(2);

  /* data with fixed values: */
  obj = {1, 1};
  rhs = {0, -1};
  beg = {0, 2};
  ind = {0, 1, 0, 1};
  val = {1, -1, -1, 1};

  /* fill variable data */
  lb[0] = -lp_interface_->Infinity();
  lb[1] = -lp_interface_->Infinity();
  ub[0] = lp_interface_->Infinity();
  ub[1] = lp_interface_->Infinity();
  lhs[0] = -lp_interface_->Infinity();
  lhs[1] = -lp_interface_->Infinity();

  /* empty_placeholders */
  LPValueArray empty_vals;

  /* check problem with primal simplex */
  ASSERT_NO_FATAL_FAILURE(  performTest(true, LPObjectiveSense::OBJ_SENSE_MINIMIZE, 2, obj, lb, ub, 2, lhs, rhs, beg, ind, val,
              LPFeasibilityStat::INFEASIBLE, LPFeasibilityStat::INFEASIBLE, empty_vals, empty_vals, empty_vals, empty_vals));

  /* check that data stored in lpi is still the same */
  ASSERT_NO_FATAL_FAILURE(checkData(LPObjectiveSense::OBJ_SENSE_MINIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val));

  /* check problem with dual simplex */
  ASSERT_NO_FATAL_FAILURE(solveTest(false, 2, 2, LPFeasibilityStat::INFEASIBLE, LPFeasibilityStat::INFEASIBLE, empty_vals, empty_vals, empty_vals, empty_vals));

  /* check that data stored in lpi is still the same */
  ASSERT_NO_FATAL_FAILURE(checkData(LPObjectiveSense::OBJ_SENSE_MINIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val));
}

/** Test 5: Test objective limit
 *
 * Use second problem from Test 1 and set objective limit.
 *
 * This is a quite weak test. For instance SoPlex directly finishes with the optimal solution.
 */
TEST_F(Solve, test5) {
  /* data to be filled */
  obj.reserve(2);
  rhs.reserve(2);
  beg.reserve(2);
  ind.reserve(4);
  val.reserve(4);

  exp_primsol.reserve(2);
  exp_dualsol.reserve(2);
  exp_activity.reserve(2);
  exp_redcost.reserve(2);

  lb.reserve(2);
  ub.reserve(2);
  lhs.reserve(2);

  /* data with fixed values: */
  obj = {1, 1};
  lb = {0, 0};
  rhs = {10, 15};
  beg = {0, 2};
  ind = {0, 1, 0, 1};
  val = {2, 1, 1, 3};

  LPValue objval;
  LPBaseStatArray cstat(2);
  LPBaseStatArray rstat(2);
  cstat = {LPBaseStat::BASESTAT_LOWER, LPBaseStat::BASESTAT_LOWER};
  rstat = {LPBaseStat::BASESTAT_BASIC, LPBaseStat::BASESTAT_BASIC};
  LPValue exp_objval = 5.0;

  /* expected solutions */
  exp_primsol = {0.4, 0.2};
  exp_dualsol = {3, 4};
  exp_activity = {1, 1};
  exp_redcost = {0, 0};

  /* fill variable data */
  ub[0] = lp_interface_->Infinity();
  ub[1] = lp_interface_->Infinity();
  lhs[0] = -lp_interface_->Infinity();
  lhs[1] = -lp_interface_->Infinity();

  /* empty_placeholders */
  StringArray empty_names;

  /* load problem */
  ASSERT_EQ(lp_interface_->LoadColumnLP(LPObjectiveSense::OBJ_SENSE_MAXIMIZE, 2, obj, lb, ub, empty_names, 2, lhs, rhs, empty_names, 4, beg, ind, val), RetCode::OKAY);

  /* set objective limit */
  ASSERT_TRUE((lp_interface_->SetIntegerParameter(LPParameter::FROM_SCRATCH, 1) == RetCode::OKAY) || (lp_interface_->SetIntegerParameter(LPParameter::FROM_SCRATCH, 1) == RetCode::PARAMETER_UNKNOWN));
  ASSERT_TRUE((lp_interface_->SetIntegerParameter(LPParameter::PRESOLVING, 0) == RetCode::OKAY) || (lp_interface_->SetIntegerParameter(LPParameter::PRESOLVING, 0) == RetCode::PARAMETER_UNKNOWN));
  ASSERT_EQ(lp_interface_->SetRealParameter(LPParameter::OBJECTIVE_LIMIT, 0.0), RetCode::OKAY);

  /* set basis */
  ASSERT_EQ(lp_interface_->SetBase(cstat, rstat), RetCode::OKAY);

  /* solve problem */
  ASSERT_EQ(lp_interface_->SolveDual(), RetCode::OKAY);

  /* check status */
  ASSERT_TRUE(lp_interface_->WasSolved());
  ASSERT_TRUE(lp_interface_->IsObjectiveLimitExceeded() || lp_interface_->IsOptimal());
  ASSERT_TRUE(!lp_interface_->IsIterationLimitExceeded());
  ASSERT_TRUE(!lp_interface_->IsTimeLimitExceeded());

  /* the objective should be equal to the objective limit */
  ASSERT_EQ(lp_interface_->GetObjectiveValue(objval), RetCode::OKAY);
  ASSERT_GE(objval, exp_objval);// << "Objective value not equal to objective limit: %g != %g\n", objval, exp_objval);

  /* check that data stored in lpi is still the same */
  ASSERT_NO_FATAL_FAILURE(checkData(LPObjectiveSense::OBJ_SENSE_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val));
}

/** Test 6: More complex example
 *
 * The original problem was the following (arising from the qpkktref unit test), which displays a bug in CPLEX 12.7.0
 * w.r.t. scaling:
 *   Minimize t_objvar
 *   Subject To
 *     KKTBinary1_y:                 - t_dual_y_bin1 + t_dual_y_bin2 + t_dual_y_slackbin1 = 0
 *     KKTlin_lower_1:               - t_x - t_y + t_slack_lhs_lower + t_slack_ub_z       = 0.75
 *     KKTBinary1_x:                 - t_dual_x_bin1 + t_dual_x_bin2 + t_dual_x_slackbin1 = 0
 *     KKTlin_lower_0:               - t_x - t_y + t_slack_ub_z - t_slack_rhs_lower       = 0.25
 *     quadratic_side1_estimation_0: 2.75 t_x - 3.75 t_y + t_objvar + 2.28 t_slack_ub_z  <= 5.0496
 *     quadratic_side0_estimation_0: 1.25 t_x - 0.25 t_y + t_objvar + 2 t_slack_ub_z     >= 2.6875
 *     quadratic_side1_estimation_0: 0.75 t_x - 0.25 t_y + t_objvar + 0.68 t_slack_ub_z  <= 4.2056
 *     quadratic_side0_estimation_0: 2.75 t_x - 0.25 t_y + t_objvar + 3 t_slack_ub_z     >= 4.4375
 *   Bounds
 *     t_x = 1
 *     t_y = 0
 *     -2.562500001 <= t_objvar <= -0.0624999989999999
 *     0 <= t_slack_lhs_lower <= 0.5
 *     t_dual_x_bin1 Free
 *     t_dual_x_bin2 Free
 *     t_dual_x_slackbin1 = 0
 *     t_dual_y_bin1 = 0
 *     t_dual_y_bin2 Free
 *     t_dual_y_slackbin1 Free
 *     1.25 <= t_slack_ub_z <= 1.75
 *     0 <= t_slack_rhs_lower <= 0.5
 *   End
 *
 *  We use the following mapping between variables and indices:
 *  0:  t_x = 1
 *  1:  t_y = 0
 *  2:  t_objvar
 *  3:  t_slack_lhs_lower
 *  4:  t_dual_x_bin1
 *  5:  t_dual_x_bin2
 *  6:  t_dual_x_slackbin1
 *  7:  t_dual_y_bin1
 *  8:  t_dual_y_bin2
 *  9:  t_dual_y_slackbin1
 *  10: t_slack_ub_z
 *  11: t_slack_rhs_lower
 */
TEST_F(Solve, test6) {
  /* data to be filled */
  obj.reserve(12);
  lb.reserve(12);
  ub.reserve(12);
  lhs.reserve(8);
  rhs.reserve(8);
  beg.reserve(8);
  ind.reserve(30);
  val.reserve(30);

  exp_primsol.reserve(2);
  exp_dualsol.reserve(2);
  exp_activity.reserve(2);
  exp_redcost.reserve(2);

  /* data with fixed values: */

  LPValue exp_objval = -2.0625;
  LPValue objval;

  LPIndexArray varind(12);
  varind = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

  /* LP data: */
  obj = {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  lb = {1, 0, -2.5625, 0, -1e20, -1e20, 0.0, 0.0, -1e20, -1e20, 1.25, 0};
  ub = {1, 0, -0.0625, 0.5, 1e20, 1e20, 0.0, 0.0, 1e20, 1e20, 1.75, 0.5};

  lhs = {0, 0.75, 0, 0.25, -1e20, 2.6875, -1e20, 4.4375};
  rhs = {0, 0.75, 0, 0.25, 5.0496, 1e20, 4.2056, 1e20};

  /* matrix */
  beg = {0, 6, 12, 16, 17, 18, 19, 20, 21, 22, 23, 29};
  /*             x0                x1                x2          x3 x4 x5 x6 x7 x8 x9 x10               x11 */
  ind = {1, 3, 4, 5, 6, 7, 1, 3, 4, 5, 6, 7, 4, 5, 6, 7, 1, 2, 2, 2, 0, 0, 0, 1, 3, 4, 5, 6, 7, 3};
  /*                   x0                              x1                                  x2          x3 x4  x5 x6 x7  x8 x9 x10                     x11 */
  val = {-1, -1, 2.75, 1.25, 0.75, 2.75, -1, -1, -3.75, -0.25, -0.25, -0.25, 1, 1, 1, 1, 1, -1, 1, 1, -1, 1, 1, 1, 1, 2.28, 2, 0.68, 3, -1.0};
  int j;

  /* possibly convert |1e20| to infinity of LPI */
  for (j = 0; j < 12; ++j) {
    if (lb[j] == -1e20)
      lb[j] = -lp_interface_->Infinity();
    if (ub[j] == 1e20)
      ub[j] = lp_interface_->Infinity();
  }
  for (j = 0; j < 8; ++j) {
    if (lhs[j] == -1e20)
      lhs[j] = -lp_interface_->Infinity();
    if (rhs[j] == 1e20)
      rhs[j] = lp_interface_->Infinity();
  }
  /* empty_placeholders */
  StringArray empty_names;

  /* load problem */
  ASSERT_EQ(lp_interface_->LoadColumnLP(LPObjectiveSense::OBJ_SENSE_MINIMIZE, 12, obj, lb, ub, empty_names, 8, lhs, rhs, empty_names, 30, beg, ind, val), RetCode::OKAY);

  /* set some parameters - simulate settings in MiniMIP */
  ASSERT_TRUE((lp_interface_->SetIntegerParameter(LPParameter::FROM_SCRATCH, 0) == RetCode::OKAY) || (lp_interface_->SetIntegerParameter(LPParameter::FROM_SCRATCH, 0) == RetCode::PARAMETER_UNKNOWN));
  ASSERT_TRUE((lp_interface_->SetIntegerParameter(LPParameter::SCALING, 1) == RetCode::OKAY) || (lp_interface_->SetIntegerParameter(LPParameter::SCALING, 1) == RetCode::PARAMETER_UNKNOWN));
  ASSERT_TRUE((lp_interface_->SetIntegerParameter(LPParameter::PRESOLVING, 1) == RetCode::OKAY) || (lp_interface_->SetIntegerParameter(LPParameter::PRESOLVING, 1) == RetCode::PARAMETER_UNKNOWN));
  ASSERT_TRUE((lp_interface_->SetIntegerParameter(LPParameter::PRICING, 0) == RetCode::OKAY) || (lp_interface_->SetIntegerParameter(LPParameter::PRICING, 0) == RetCode::PARAMETER_UNKNOWN));

  ASSERT_TRUE((lp_interface_->SetRealParameter(LPParameter::FEASIBLITY_TOLERANCE, 1e-06) == RetCode::OKAY) || (lp_interface_->SetRealParameter(LPParameter::FEASIBLITY_TOLERANCE, 1e-06) == RetCode::PARAMETER_UNKNOWN));
  ASSERT_TRUE((lp_interface_->SetRealParameter(LPParameter::DUAL_FEASIBILITY_TOLERANCE, 1e-07) == RetCode::OKAY) || (lp_interface_->SetRealParameter(LPParameter::DUAL_FEASIBILITY_TOLERANCE, 1e-07) == RetCode::PARAMETER_UNKNOWN));

  ASSERT_EQ(lp_interface_->ClearState(), RetCode::OKAY);

  /* set objlimit */
  ASSERT_EQ(lp_interface_->SetRealParameter(LPParameter::OBJECTIVE_LIMIT, 4.320412501), RetCode::OKAY);

  /* solve problem */
  ASSERT_EQ(lp_interface_->SolveDual(), RetCode::OKAY);

  /* check status */
  ASSERT_TRUE(lp_interface_->WasSolved());

  /* the objective should be equal to the objective limit */
  ASSERT_EQ(lp_interface_->GetObjectiveValue(objval), RetCode::OKAY);
  ASSERT_FLOAT_EQ(objval, exp_objval);

  /* check that data stored in lpi is still the same */
  ASSERT_NO_FATAL_FAILURE(checkData(LPObjectiveSense::OBJ_SENSE_MINIMIZE, 12, obj, lb, ub, 8, lhs, rhs, 30, beg, ind, val));

  /* change some bounds */
  lb[0] = 1;
  ub[0] = 1;
  lb[1] = 0;
  ub[1] = 0;
  lb[2] = -2.06255;
  ub[2] = -2.0625;
  lb[3] = 0;
  ub[3] = 4.94694e-05;
  lb[6] = 0;
  ub[6] = 0;
  lb[7] = 0;
  ub[7] = 0;
  lb[10] = 1.74995;
  ub[10] = 1.750;
  lb[11] = 0.499951;
  ub[11] = 0.5;
  ASSERT_EQ(lp_interface_->ChangeBounds(12, varind, lb, ub), RetCode::OKAY);

  /* set objlimit */
  ASSERT_EQ(lp_interface_->SetRealParameter(LPParameter::OBJECTIVE_LIMIT, -2.0625), RetCode::OKAY);

  ASSERT_EQ(lp_interface_->ClearState(), RetCode::OKAY);

  /* solve problem */
  ASSERT_EQ(lp_interface_->SolveDual(), RetCode::OKAY);

  /* check that data stored in lpi is still the same */
  ASSERT_NO_FATAL_FAILURE(checkData(LPObjectiveSense::OBJ_SENSE_MINIMIZE, 12, obj, lb, ub, 8, lhs, rhs, 30, beg, ind, val));
}

/** Test 7
 *
 *  min 10 x1 + 15 x2
 *       2 x1 +   x2 >= 3
 *         x1 + 3 x2 <= 1
 *         x1,    x2 >= 0
 *
 *  which is dual unbounded (this is a variant of Test 3 in which the equations have been replaced by inequalities).
 *
 *  The dual is:
 *  max  3 y1 +   y2
 *       2 y1 +   y2 <= 10
 *         y1 + 3 y2 <= 15
 *         y1 >= 0, y2 <= 0
 */
TEST_F(Solve, test7) {
  /* data to be filled */
  obj.reserve(2);
  lb.reserve(2);
  ub.reserve(2);
  lhs.reserve(2);
  rhs.reserve(2);
  beg.reserve(2);
  ind.reserve(4);
  val.reserve(4);

  exp_primsol.reserve(2);
  exp_dualsol.reserve(2);
  exp_activity.reserve(2);
  exp_redcost.reserve(2);

  /* data with fixed values: */
  obj = {10, 15};
  lb = {0, 0};
  lhs = {3, 1};
  rhs = {3, 1};
  beg = {0, 2};
  ind = {0, 1, 0, 1};
  val = {2, 1, 1, 3};

  /* expected ray */
  LPValueArray exp_dualray(2);
  exp_dualray= {0.5, -1};

  /* fill data */
  rhs[0] = lp_interface_->Infinity();
  lhs[1] = -lp_interface_->Infinity();
  ub[0] = lp_interface_->Infinity();
  ub[1] = lp_interface_->Infinity();

  /* empty_placeholders */
  LPValueArray empty_vals;

  /* check problem with primal simplex */
  ASSERT_NO_FATAL_FAILURE(  performTest(true, LPObjectiveSense::OBJ_SENSE_MINIMIZE, 2, obj, lb, ub, 2, lhs, rhs, beg, ind, val,
              LPFeasibilityStat::INFEASIBLE, LPFeasibilityStat::UNBOUNDED, empty_vals, exp_dualray, empty_vals, empty_vals));

  /* check that data stored in lpi is still the same */
  ASSERT_NO_FATAL_FAILURE(checkData(LPObjectiveSense::OBJ_SENSE_MINIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val));

  /* clear basis status */
  ASSERT_EQ(lp_interface_->ClearState(), RetCode::OKAY);

  /* check problem with dual simplex */
  ASSERT_NO_FATAL_FAILURE(solveTest(false, 2, 2, LPFeasibilityStat::INFEASIBLE, LPFeasibilityStat::UNBOUNDED, empty_vals, exp_dualray, empty_vals, empty_vals));

  /* check that data stored in lpi is still the same */
  ASSERT_NO_FATAL_FAILURE(checkData(LPObjectiveSense::OBJ_SENSE_MINIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val));
}
} /* namespace minimip */
