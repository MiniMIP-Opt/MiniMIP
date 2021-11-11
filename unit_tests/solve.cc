
// Unit tests for checking LP solutions.
// We perform tests with solving several examples. These are inspired by the
// unit tests of OSI in COIN-OR.

#include <gtest/gtest.h>

#include <tuple>

#include "absl/status/status.h"
#include "src/lp_interface/lpi_factory.h"
#include "unit_tests/utils.h"

#define EPS 1e-6
#define DEF_INTERFACE \
  1  // 0 = Glop Interface (Default),
     // 1 = SoPlex Interface,

namespace minimip {

// expected feasibility status for primal or dual problem
enum LPFeasibilityStat {
  FEASIBLE   = 0,  // the problem is feasible
  UNBOUNDED  = 1,  // the problem is unbounded
  INFEASIBLE = 2   // the problem is infeasible
};
typedef enum LPFeasibilityStat LPFeasibilityStat;

// TEST SUITE SOLVE

static LPInterface* lp_interface_ = nullptr;

class Solve : public ::testing::Test {
 protected:
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
    ASSERT_OK(lp_interface_->SetObjectiveSense(LPObjectiveSense::kMaximization));
  }
  // local functions

  // solve problem
  static void solveTest(bool solve_primal, int num_columns, int num_rows,
                        LPFeasibilityStat expected_primal_feasibility_status,
                        LPFeasibilityStat expected_dual_feasibility_status,
                        const std::vector<double>& expected_primal_solution,
                        const std::vector<double>& expected_dual_solution,
                        const std::vector<double>& expected_row_activities,
                        const std::vector<double>& expected_reduced_costs) {
    // solution data
    std::vector<double> primal_solution;
    std::vector<double> dual_solution;
    std::vector<double> row_activities;
    std::vector<double> reduced_costs;

    // check size
    ASSERT_EQ(num_rows, lp_interface_->GetNumberOfRows());
    ASSERT_EQ(num_columns, lp_interface_->GetNumberOfColumns());

    // solve problem
    if (solve_primal) {
      ASSERT_OK(lp_interface_->SolveLPWithPrimalSimplex());
    } else {
      ASSERT_OK(lp_interface_->SolveLPWithDualSimplex());
    }

    // check status
    ASSERT_TRUE(lp_interface_->IsSolved());
    ASSERT_TRUE(!lp_interface_->ObjectiveLimitIsExceeded());
    ASSERT_TRUE(!lp_interface_->IterationLimitIsExceeded());
    ASSERT_TRUE(!lp_interface_->TimeLimitIsExceeded());

    // check feasibility status
    auto primal_feasible = lp_interface_->IsPrimalFeasible();
    auto dual_feasible   = lp_interface_->IsDualFeasible();

    // if we are feasible, we should be optimal
    if (expected_primal_feasibility_status == LPFeasibilityStat::FEASIBLE &&
        expected_dual_feasibility_status == LPFeasibilityStat::FEASIBLE) {
      ASSERT_TRUE(lp_interface_->IsOptimal());
    }

    // check more primal statuses
    switch (expected_primal_feasibility_status) {
      case LPFeasibilityStat::FEASIBLE:
        ASSERT_TRUE(primal_feasible);
        ASSERT_TRUE(!lp_interface_->ExistsPrimalRay());
        ASSERT_TRUE(!lp_interface_->HasPrimalRay());
        ASSERT_TRUE(!lp_interface_->IsPrimalUnbounded());
        ASSERT_TRUE(!lp_interface_->IsPrimalInfeasible());
        ASSERT_TRUE(lp_interface_->IsPrimalFeasible());
        break;

      case LPFeasibilityStat::UNBOUNDED:
        // Because of SoPlex, cannot always determine feasibility status here,
        // even if we want to apply the primal simplex. In any case, the results
        // of primal_feasible and IsPrimalFeasible() should coincide.
        ASSERT_EQ(primal_feasible, lp_interface_->IsPrimalFeasible());

        // primal ray should exist if the primal simplex ran
        ASSERT_TRUE(!solve_primal || lp_interface_->ExistsPrimalRay());
        ASSERT_TRUE(!lp_interface_->IsPrimalInfeasible());
        break;

      case LPFeasibilityStat::INFEASIBLE:
        ASSERT_TRUE(!primal_feasible);
        ASSERT_TRUE(!lp_interface_->IsPrimalFeasible());
        break;

      default:
        abort();
    }

    // check more dual statuses
    switch (expected_dual_feasibility_status) {
      case LPFeasibilityStat::FEASIBLE:
        ASSERT_TRUE(dual_feasible);
        ASSERT_TRUE(!lp_interface_->ExistsDualRay());
        ASSERT_TRUE(!lp_interface_->HasDualRay());
        ASSERT_TRUE(!lp_interface_->IsDualUnbounded());
        ASSERT_TRUE(!lp_interface_->IsDualInfeasible());
        ASSERT_TRUE(lp_interface_->IsDualFeasible());
        break;

      case LPFeasibilityStat::UNBOUNDED:
        // Because of SoPlex, cannot always determine feasibility status here,
        // even if we want to apply the dual simplex. In any case, the results
        // of dual_feasible and IsDualFeasible() should coincide.
        ASSERT_EQ(dual_feasible, lp_interface_->IsDualFeasible());

        // dual ray should exist if the dual simplex ran
        ASSERT_TRUE(solve_primal || lp_interface_->ExistsDualRay());
        ASSERT_TRUE(!lp_interface_->IsDualInfeasible());
        break;

      case LPFeasibilityStat::INFEASIBLE:
        ASSERT_TRUE(!dual_feasible);
        ASSERT_TRUE(!lp_interface_->IsDualUnbounded());
        ASSERT_TRUE(!lp_interface_->IsDualFeasible());
        break;

      default:
        abort();
    }

    primal_solution.reserve(num_columns);
    dual_solution.reserve(num_rows);
    row_activities.reserve(num_rows);
    reduced_costs.reserve(num_columns);

    // check solution
    if (expected_primal_feasibility_status == LPFeasibilityStat::FEASIBLE) {
      // get solution

      absl::StatusOr<std::vector<double>> absl_tmp;

      absl_tmp = lp_interface_->GetPrimalSolution();
      ASSERT_OK(absl_tmp.status());
      primal_solution = *absl_tmp;

      absl_tmp = lp_interface_->GetRowActivity();
      ASSERT_OK(absl_tmp.status());
      row_activities = *absl_tmp;

      absl_tmp = lp_interface_->GetDualSolution();
      ASSERT_OK(absl_tmp.status());
      dual_solution = *absl_tmp;

      absl_tmp = lp_interface_->GetReducedCost();
      ASSERT_OK(absl_tmp.status());
      reduced_costs = *absl_tmp;

      for (int j = 0; j < num_columns; ++j) {
        ASSERT_FLOAT_EQ(primal_solution[j], expected_primal_solution[j]);
        ASSERT_FLOAT_EQ(reduced_costs[j], expected_reduced_costs[j]);
      }
    } else if (expected_primal_feasibility_status ==
               LPFeasibilityStat::UNBOUNDED) {
      if (lp_interface_->HasPrimalRay()) {
        double scaling_factor = 1.0;

        absl::StatusOr<std::vector<double>> absl_tmp;

        absl_tmp = lp_interface_->GetPrimalRay();
        ASSERT_OK(absl_tmp.status());
        primal_solution = *absl_tmp;

        // loop until scaling factor can be determined
        for (int j = 0; j < num_columns; ++j) {
          if (REALABS(expected_primal_solution[j]) < EPS) {
            ASSERT_FLOAT_EQ(primal_solution[j], expected_primal_solution[j]);
          } else {
            scaling_factor = primal_solution[j] / expected_primal_solution[j];
            break;
          }
        }

        // again loop over ray
        for (int j = 0; j < num_columns; ++j) {
          ASSERT_FLOAT_EQ(primal_solution[j],
                          scaling_factor * expected_primal_solution[j]);
        }
      }
    }

    if (expected_dual_feasibility_status == LPFeasibilityStat::FEASIBLE) {
      // get solution
      absl::StatusOr<std::vector<double>> absl_tmp;

      absl_tmp = lp_interface_->GetPrimalSolution();
      ASSERT_OK(absl_tmp.status());
      primal_solution = *absl_tmp;

      absl_tmp = lp_interface_->GetRowActivity();
      ASSERT_OK(absl_tmp.status());
      row_activities = *absl_tmp;

      absl_tmp = lp_interface_->GetDualSolution();
      ASSERT_OK(absl_tmp.status());
      dual_solution = *absl_tmp;

      absl_tmp = lp_interface_->GetReducedCost();
      ASSERT_OK(absl_tmp.status());
      reduced_costs = *absl_tmp;

      for (int i = 0; i < num_rows; ++i) {
        ASSERT_FLOAT_EQ(dual_solution[i], expected_dual_solution[i]);
        ASSERT_FLOAT_EQ(row_activities[i], expected_row_activities[i]);
      }
    } else if (expected_dual_feasibility_status ==
               LPFeasibilityStat::UNBOUNDED) {
      if (lp_interface_->HasDualRay()) {
        double scaling_factor = 1.0;
        std::vector<double> left_hand_sides(num_rows);
        std::vector<double> right_hand_sides(num_rows);

        // get left_hand_sides/right_hand_sides for check of dual ray
        for (int i = 0; i < num_rows; i++) {
          left_hand_sides.push_back(lp_interface_->GetLeftHandSide(i));
          right_hand_sides.push_back(lp_interface_->GetRightHandSide(i));
        }

        // get dual ray
        absl::StatusOr<std::vector<double>> absl_tmp;

        absl_tmp = lp_interface_->GetDualFarkasMultiplier();
        ASSERT_OK(absl_tmp.status());
        dual_solution = *absl_tmp;

        // loop until scaling factor can be determined
        for (int i = 0; i < num_rows; ++i) {
          if (REALABS(expected_dual_solution[i]) < EPS) {
            ASSERT_FLOAT_EQ(dual_solution[i], expected_dual_solution[i]);
          } else {
            scaling_factor = dual_solution[i] / expected_dual_solution[i];
            break;
          }
        }

        // again loop over ray
        for (int i = 0; i < num_rows; ++i) {
          ASSERT_FLOAT_EQ(dual_solution[i],
                          scaling_factor * expected_dual_solution[i]);
          ASSERT_TRUE(!lp_interface_->IsInfinity(-left_hand_sides[i]) ||
                      dual_solution[i] <= -EPS);
          ASSERT_TRUE(!lp_interface_->IsInfinity(right_hand_sides[i]) ||
                      dual_solution[i] >= EPS);
        }
      }
    }
  }

  // perform basic test for the given problem
  static void performTest(bool solve_primal, LPObjectiveSense objective_sense,
                          int num_columns,
                          const std::vector<double>& objective_coefficients,
                          const std::vector<double>& lower_bounds,
                          const std::vector<double>& upper_bounds, int num_rows,
                          const std::vector<double>& left_hand_sides,
                          const std::vector<double>& right_hand_sides,
                          std::vector<SparseVector> sparse_columns,
                          LPFeasibilityStat expected_primal_feasibility_status,
                          LPFeasibilityStat expected_dual_feasibility_status,
                          const std::vector<double>& expected_primal_solution,
                          const std::vector<double>& expected_dual_solution,
                          const std::vector<double>& expected_row_activities,
                          const std::vector<double>& expected_reduced_costs) {
    std::vector<std::string> empty_names;

    ASSERT_OK(lp_interface_->LoadSparseColumnLP(
        objective_sense, objective_coefficients, lower_bounds, upper_bounds,
        empty_names, left_hand_sides, right_hand_sides, empty_names,
        sparse_columns));

    // load problem
    ASSERT_TRUE(!lp_interface_->IsSolved());

    // solve problem
    ASSERT_NO_FATAL_FAILURE(solveTest(
        solve_primal, num_columns, num_rows, expected_primal_feasibility_status,
        expected_dual_feasibility_status, expected_primal_solution,
        expected_dual_solution, expected_row_activities,
        expected_reduced_costs));
  }

  // check whether data in LP solver aggrees with original data
  static void checkData(LPObjectiveSense objective_sense, int num_columns,
                        const std::vector<double>& objective_coefficients,
                        const std::vector<double>& lower_bounds,
                        const std::vector<double>& upper_bounds, int num_rows,
                        const std::vector<double>& left_hand_sides,
                        const std::vector<double>& right_hand_sides,
                        int num_nonzeros,
                        const std::vector<SparseVector>& sparse_columns) {
    LPObjectiveSense expected_objective_sense;
    std::vector<double> check_val;
    std::vector<double> expected_lower_bounds;
    std::vector<double> expected_upper_bounds;
    std::vector<double> expected_objective_coefficients;
    std::vector<double> expected_left_hand_sides;
    std::vector<double> expected_right_hand_sides;
    std::vector<SparseVector> expected_sparse_columns;
    int expected_num_columns;
    int expected_num_rows;
    int expected_num_nonzeros;

    // check number of rows and columns
    expected_num_rows    = lp_interface_->GetNumberOfRows();
    expected_num_columns = lp_interface_->GetNumberOfColumns();
    ASSERT_EQ(expected_num_rows, num_rows);
    ASSERT_EQ(expected_num_columns, num_columns);

    // check objective sense
    expected_objective_sense = lp_interface_->GetObjectiveSense();
    ASSERT_EQ(objective_sense, expected_objective_sense);

    // get number of nonzeros in matrix
    expected_num_nonzeros = lp_interface_->GetNumberOfNonZeros();
    ASSERT_EQ(expected_num_nonzeros, num_nonzeros);

    check_val.reserve(expected_num_nonzeros);
    expected_lower_bounds.reserve(num_columns);
    expected_upper_bounds.reserve(num_columns);
    expected_objective_coefficients.reserve(num_columns);

    // get matrix data
    expected_sparse_columns.reserve(num_columns);

    for (int i = 0; i < num_columns; i++) {
      expected_lower_bounds.push_back(lp_interface_->GetLowerBound(i));
      expected_upper_bounds.push_back(lp_interface_->GetUpperBound(i));

      expected_sparse_columns.push_back(
          lp_interface_->GetSparseColumnCoefficients(i));
    }

    for (int i = 0; i < num_columns; i++) {
      expected_objective_coefficients.push_back(
          lp_interface_->GetObjectiveCoefficient(i));
    }

    // compare data
    for (int j = 0; j < num_columns; ++j) {
      if (fabs(expected_lower_bounds[j]) < 1e30 &&
          fabs(lower_bounds[j]) < 1e30) {
        ASSERT_FLOAT_EQ(expected_lower_bounds[j], lower_bounds[j]);
      }
      if (fabs(expected_upper_bounds[j]) < 1e30 &&
          fabs(upper_bounds[j]) < 1e30) {
        ASSERT_FLOAT_EQ(expected_upper_bounds[j], upper_bounds[j]);
      }

      ASSERT_FLOAT_EQ(expected_objective_coefficients[j],
                      objective_coefficients[j]);
    }

    // compare sparse rows
    ASSERT_EQ(expected_sparse_columns.size(), sparse_columns.size());

    for (size_t j = 0; j < expected_sparse_columns.size(); ++j) {
      ASSERT_EQ(expected_sparse_columns[j].indices.size(),
                sparse_columns[j].indices.size());
      for (size_t i = 0; i < expected_sparse_columns[j].indices.size(); ++i) {
        expected_sparse_columns[j].values[i];
        ASSERT_EQ(expected_sparse_columns[j].indices[i],
                  sparse_columns[j].indices[i]);
        ASSERT_EQ(expected_sparse_columns[j].values[i],
                  sparse_columns[j].values[i]);
      }
    }

    expected_left_hand_sides.reserve(num_rows);
    expected_right_hand_sides.reserve(num_rows);

    for (int i = 0; i < num_rows; i++) {
      expected_left_hand_sides.push_back(lp_interface_->GetLeftHandSide(i));
      expected_right_hand_sides.push_back(lp_interface_->GetRightHandSide(i));
    }

    for (int i = 0; i < num_rows; ++i) {
      if (fabs(expected_left_hand_sides[i]) < 1e30 &&
          fabs(left_hand_sides[i]) < 1e30) {
        ASSERT_FLOAT_EQ(expected_left_hand_sides[i], left_hand_sides[i]);
      }
      if (fabs(expected_left_hand_sides[i]) < 1e30 &&
          fabs(left_hand_sides[i]) < 1e30) {
        ASSERT_FLOAT_EQ(expected_right_hand_sides[i], right_hand_sides[i]);
      }
    }
  }
};

// TESTS

// Test 1

// max 3 x1 +   x2
//    2 x1 +   x2 <= 10
//      x1 + 3 x2 <= 15
//      x1,    x2 >= 0

// with primal optimal solution (5, 0), dual optimal solution (1.5, 0), activity
// (10, 5), and reduced_costs (0, -0.5).

// Then use objective (1, 1) with primal optimal solution (3,4), dual optimal
// solution (0.4, 0.2), activity (10, 15), and reduced_costs (0, 0).
TEST_F(Solve, test1) {
  // LP data
  std::vector<double> objective_coefficients = {3.0, 1.0};
  std::vector<double> lower_bounds           = {0.0, 0.0};
  std::vector<double> upper_bounds           = {lp_interface_->Infinity(),
                                      lp_interface_->Infinity()};
  std::vector<double> right_hand_sides       = {10.0, 15.0};
  std::vector<double> left_hand_sides        = {-lp_interface_->Infinity(),
                                         -lp_interface_->Infinity()};

  SparseVector sparse_col_1                = {{0, 1}, {2, 1}};
  SparseVector sparse_col_2                = {{0, 1}, {1, 3}};
  std::vector<SparseVector> sparse_columns = {sparse_col_1, sparse_col_2};

  // expected solutions
  std::vector<double> expected_primal_solution = {5.0, 0.0};
  std::vector<double> expected_dual_solution   = {1.5, 0.0};
  std::vector<double> expected_row_activities  = {10.0, 5.0};
  std::vector<double> expected_reduced_costs   = {0.0, -0.5};

  // assert infinity values
  for (int j = 0; j < 2; j++) {
    ASSERT_TRUE(lp_interface_->IsInfinity(upper_bounds[j]));
    ASSERT_EQ(left_hand_sides[j], -lp_interface_->Infinity());
  }

  // solve problem with primal simplex
  ASSERT_NO_FATAL_FAILURE(performTest(
      true, LPObjectiveSense::kMaximization, 2, objective_coefficients,
      lower_bounds, upper_bounds, 2, left_hand_sides, right_hand_sides,
      sparse_columns, LPFeasibilityStat::FEASIBLE, LPFeasibilityStat::FEASIBLE,
      expected_primal_solution, expected_dual_solution, expected_row_activities,
      expected_reduced_costs));

  // check that data stored in lpi is still the same
  ASSERT_NO_FATAL_FAILURE(checkData(
      LPObjectiveSense::kMaximization, 2, objective_coefficients, lower_bounds,
      upper_bounds, 2, left_hand_sides, right_hand_sides, 4, sparse_columns));

  // clear basis status
  ASSERT_OK(lp_interface_->ClearState());

  // solve problem with dual simplex
  ASSERT_NO_FATAL_FAILURE(solveTest(
      false, 2, 2, LPFeasibilityStat::FEASIBLE, LPFeasibilityStat::FEASIBLE,
      expected_primal_solution, expected_dual_solution, expected_row_activities,
      expected_reduced_costs));

  // check that data stored in lpi is still the same
  ASSERT_NO_FATAL_FAILURE(checkData(
      LPObjectiveSense::kMaximization, 2, objective_coefficients, lower_bounds,
      upper_bounds, 2, left_hand_sides, right_hand_sides, 4, sparse_columns));

  // clear basis status
  ASSERT_OK(lp_interface_->ClearState());

  // change objective
  objective_coefficients[0] = 1.0;
  ASSERT_OK(
      lp_interface_->SetObjectiveCoefficient(0, objective_coefficients[0]));

  // change expected solution
  expected_primal_solution[0] = 3;
  expected_primal_solution[1] = 4;
  expected_dual_solution[0]   = 0.4;
  expected_dual_solution[1]   = 0.2;
  expected_row_activities[0]  = 10;
  expected_row_activities[1]  = 15;
  expected_reduced_costs[1]   = 0;

  // check changed problem with primal simplex
  ASSERT_NO_FATAL_FAILURE(performTest(
      true, LPObjectiveSense::kMaximization, 2, objective_coefficients,
      lower_bounds, upper_bounds, 2, left_hand_sides, right_hand_sides,
      sparse_columns, LPFeasibilityStat::FEASIBLE, LPFeasibilityStat::FEASIBLE,
      expected_primal_solution, expected_dual_solution, expected_row_activities,
      expected_reduced_costs));

  // check that data stored in lpi is still the same
  ASSERT_NO_FATAL_FAILURE(checkData(
      LPObjectiveSense::kMaximization, 2, objective_coefficients, lower_bounds,
      upper_bounds, 2, left_hand_sides, right_hand_sides, 4, sparse_columns));
}

// Test 2
//
// max 3 x1 +   x2
//    2 x1 +   x2 <= 10
//      x1 + 3 x2 <= 15
//      x1, x2 free
//
// which is unbounded,
// the only difference to Test 1 is that the variables are free.
//
// Then use objective (1, 1) with primal optimal solution (3,4), dual optimal
// solution (0.4, 0.2), activity (10, 15), and reduced_costs (0, 0).
TEST_F(Solve, test2) {
  // data to be filled
  std::vector<double> objective_coefficients = {3.0, 1.0};
  std::vector<double> lower_bounds           = {-lp_interface_->Infinity(),
                                      -lp_interface_->Infinity()};
  std::vector<double> upper_bounds           = {lp_interface_->Infinity(),
                                      lp_interface_->Infinity()};
  std::vector<double> right_hand_sides       = {10.0, 15.0};
  std::vector<double> left_hand_sides        = {-lp_interface_->Infinity(),
                                         -lp_interface_->Infinity()};

  SparseVector sparse_col_1                = {{0, 1}, {2, 1}};
  SparseVector sparse_col_2                = {{0, 1}, {1, 3}};
  std::vector<SparseVector> sparse_columns = {sparse_col_1, sparse_col_2};

  // expected ray for first run
  std::vector<double> expected_primal_ray = {0.5, -1};

  // expected solutions
  std::vector<double> expected_primal_solution = {3, 4};
  std::vector<double> expected_dual_solution   = {0.4, 0.2};
  std::vector<double> expected_row_activities  = {10, 15};
  std::vector<double> expected_reduced_costs   = {0, 0};

  // fill variable data
  lower_bounds[0]    = -lp_interface_->Infinity();
  lower_bounds[1]    = -lp_interface_->Infinity();
  upper_bounds[0]    = lp_interface_->Infinity();
  upper_bounds[1]    = lp_interface_->Infinity();
  left_hand_sides[0] = -lp_interface_->Infinity();
  left_hand_sides[1] = -lp_interface_->Infinity();

  // empty placeholders
  std::vector<double> empty_vals;

  // solve problem with primal simplex
  ASSERT_NO_FATAL_FAILURE(
      performTest(true, LPObjectiveSense::kMaximization, 2,
                  objective_coefficients, lower_bounds, upper_bounds, 2,
                  left_hand_sides, right_hand_sides, sparse_columns,
                  LPFeasibilityStat::UNBOUNDED, LPFeasibilityStat::INFEASIBLE,
                  expected_primal_ray, empty_vals, empty_vals, empty_vals));

  // check that data stored in lpi is still the same
  ASSERT_NO_FATAL_FAILURE(checkData(
      LPObjectiveSense::kMaximization, 2, objective_coefficients, lower_bounds,
      upper_bounds, 2, left_hand_sides, right_hand_sides, 4, sparse_columns));

  // clear basis status
  ASSERT_OK(lp_interface_->ClearState());

  // solve problem with dual simplex
  ASSERT_NO_FATAL_FAILURE(solveTest(
      false, 2, 2, LPFeasibilityStat::UNBOUNDED, LPFeasibilityStat::INFEASIBLE,
      expected_primal_ray, empty_vals, empty_vals, empty_vals));

  // check that data stored in lpi is still the same
  ASSERT_NO_FATAL_FAILURE(checkData(
      LPObjectiveSense::kMaximization, 2, objective_coefficients, lower_bounds,
      upper_bounds, 2, left_hand_sides, right_hand_sides, 4, sparse_columns));

  // clear basis status
  ASSERT_OK(lp_interface_->ClearState());

  // change objective
  objective_coefficients[0] = 1.0;
  ASSERT_OK(
      lp_interface_->SetObjectiveCoefficient(0, objective_coefficients[0]));

  // solve with primal simplex
  ASSERT_NO_FATAL_FAILURE(solveTest(
      true, 2, 2, LPFeasibilityStat::FEASIBLE, LPFeasibilityStat::FEASIBLE,
      expected_primal_solution, expected_dual_solution, expected_row_activities,
      expected_reduced_costs));

  // check that data stored in lpi is still the same
  ASSERT_NO_FATAL_FAILURE(checkData(
      LPObjectiveSense::kMaximization, 2, objective_coefficients, lower_bounds,
      upper_bounds, 2, left_hand_sides, right_hand_sides, 4, sparse_columns));
}

// Test 3
//
// min 10 y1 + 15 y2
//     2 y1 +   y2 == 3
//       y1 + 3 y2 == 1
//       y1,    y2 >= 0
//
// which is dual unbounded (this is the dual of the problem in Test 2).
//
// Then use right_hand_sides (1, 1) with primal optimal solution (0.4,0.2), dual
// optimal solution (3, 4), activity (0, 0), and reduced_costs (0, 0).
TEST_F(Solve, test3) {
  // LP data
  std::vector<double> objective_coefficients = {10, 15};
  std::vector<double> lower_bounds           = {0, 0};
  std::vector<double> upper_bounds           = {lp_interface_->Infinity(),
                                      lp_interface_->Infinity()};
  std::vector<double> right_hand_sides       = {3.0, 1.0};
  std::vector<double> left_hand_sides        = {3.0, 1.0};

  SparseVector sparse_col_1                = {{0, 1}, {2, 1}};
  SparseVector sparse_col_2                = {{0, 1}, {1, 3}};
  std::vector<SparseVector> sparse_columns = {sparse_col_1, sparse_col_2};

  // expected ray
  std::vector<double> expected_dual_ray = {0.5, -1};

  // expected solutions
  std::vector<double> expected_primal_solution = {0.4, 0.2};
  std::vector<double> expected_dual_solution   = {3, 4};
  std::vector<double> expected_row_activities  = {1, 1};
  std::vector<double> expected_reduced_costs   = {0, 0};

  // empty placeholders
  std::vector<double> empty_vals;

  // check problem with primal simplex
  ASSERT_NO_FATAL_FAILURE(
      performTest(true, LPObjectiveSense::kMinimization, 2,
                  objective_coefficients, lower_bounds, upper_bounds, 2,
                  left_hand_sides, right_hand_sides, sparse_columns,
                  LPFeasibilityStat::INFEASIBLE, LPFeasibilityStat::UNBOUNDED,
                  empty_vals, expected_dual_ray, empty_vals, empty_vals));

  // check that data stored in lpi is still the same
  ASSERT_NO_FATAL_FAILURE(checkData(
      LPObjectiveSense::kMinimization, 2, objective_coefficients, lower_bounds,
      upper_bounds, 2, left_hand_sides, right_hand_sides, 4, sparse_columns));

  // clear basis status
  ASSERT_OK(lp_interface_->ClearState());

  // check problem with dual simplex
  ASSERT_NO_FATAL_FAILURE(solveTest(false, 2, 2, LPFeasibilityStat::INFEASIBLE,
                                    LPFeasibilityStat::UNBOUNDED, empty_vals,
                                    expected_dual_ray, empty_vals, empty_vals));

  // check that data stored in lpi is still the same
  ASSERT_NO_FATAL_FAILURE(checkData(
      LPObjectiveSense::kMinimization, 2, objective_coefficients, lower_bounds,
      upper_bounds, 2, left_hand_sides, right_hand_sides, 4, sparse_columns));

  // clear basis status
  ASSERT_OK(lp_interface_->ClearState());

  // change left_hand_sides/right_hand_sides
  left_hand_sides[0]  = 1.0;
  right_hand_sides[0] = 1.0;
  ASSERT_OK(
      lp_interface_->SetRowSides(0, left_hand_sides[0], right_hand_sides[0]));
  ASSERT_NO_FATAL_FAILURE(solveTest(
      true, 2, 2, LPFeasibilityStat::FEASIBLE, LPFeasibilityStat::FEASIBLE,
      expected_primal_solution, expected_dual_solution, expected_row_activities,
      expected_reduced_costs));

  // check that data stored in lpi is still the same
  ASSERT_NO_FATAL_FAILURE(checkData(
      LPObjectiveSense::kMinimization, 2, objective_coefficients, lower_bounds,
      upper_bounds, 2, left_hand_sides, right_hand_sides, 4, sparse_columns));
}

// Test 4
//
// max x1 + x2
//    x1 - x2 <= 0
//  - x1 + x2 <= -1
//    x1,  x2 free
//
// which primal and dual infeasible.
TEST_F(Solve, test4) {
  // LP data
  std::vector<double> objective_coefficients{1, 1};
  std::vector<double> lower_bounds     = {-lp_interface_->Infinity(),
                                      -lp_interface_->Infinity()};
  std::vector<double> upper_bounds     = {lp_interface_->Infinity(),
                                      lp_interface_->Infinity()};
  std::vector<double> right_hand_sides = {0, -1};
  std::vector<double> left_hand_sides  = {-lp_interface_->Infinity(),
                                         -lp_interface_->Infinity()};

  SparseVector sparse_col_1                = {{0, 1}, {1, -1}};
  SparseVector sparse_col_2                = {{0, 1}, {-1, 1}};
  std::vector<SparseVector> sparse_columns = {sparse_col_1, sparse_col_2};

  // empty placeholders
  std::vector<double> empty_vals;

  // check problem with primal simplex
  ASSERT_NO_FATAL_FAILURE(
      performTest(true, LPObjectiveSense::kMinimization, 2,
                  objective_coefficients, lower_bounds, upper_bounds, 2,
                  left_hand_sides, right_hand_sides, sparse_columns,
                  LPFeasibilityStat::INFEASIBLE, LPFeasibilityStat::INFEASIBLE,
                  empty_vals, empty_vals, empty_vals, empty_vals));

  // check that data stored in lpi is still the same
  ASSERT_NO_FATAL_FAILURE(checkData(
      LPObjectiveSense::kMinimization, 2, objective_coefficients, lower_bounds,
      upper_bounds, 2, left_hand_sides, right_hand_sides, 4, sparse_columns));

  // check problem with dual simplex
  ASSERT_NO_FATAL_FAILURE(solveTest(false, 2, 2, LPFeasibilityStat::INFEASIBLE,
                                    LPFeasibilityStat::INFEASIBLE, empty_vals,
                                    empty_vals, empty_vals, empty_vals));

  // check that data stored in lpi is still the same
  ASSERT_NO_FATAL_FAILURE(checkData(
      LPObjectiveSense::kMinimization, 2, objective_coefficients, lower_bounds,
      upper_bounds, 2, left_hand_sides, right_hand_sides, 4, sparse_columns));
}

// Test 5: Test objective limit
//
// Use second problem from Test 1 and set objective limit.
//
// This is a quite weak test. For instance SoPlex directly finishes with the
// optimal solution.
TEST_F(Solve, test5) {
  // LP data
  std::vector<double> objective_coefficients = {1, 1};
  std::vector<double> lower_bounds           = {0, 0};
  std::vector<double> upper_bounds           = {lp_interface_->Infinity(),
                                      lp_interface_->Infinity()};
  std::vector<double> right_hand_sides       = {10, 15};
  std::vector<double> left_hand_sides        = {-lp_interface_->Infinity(),
                                         -lp_interface_->Infinity()};

  double objective_value;
  std::vector<LPBasisStatus> column_basis_status(2);
  std::vector<LPBasisStatus> row_basis_status(2);
  column_basis_status = {LPBasisStatus::kAtLowerBound,
                         LPBasisStatus::kAtLowerBound};
  row_basis_status    = {LPBasisStatus::kBasic, LPBasisStatus::kBasic};
  double expected_objective_value = 5.0;

  // expected solutions
  std::vector<double> expected_primal_solution = {0.4, 0.2};
  std::vector<double> expected_dual_solution   = {3, 4};
  std::vector<double> expected_row_activities  = {1, 1};
  std::vector<double> expected_reduced_costs   = {0, 0};

  // empty placeholders
  std::vector<std::string> empty_names;

  SparseVector sparse_col_1                = {{0, 1}, {2, 1}};
  SparseVector sparse_col_2                = {{0, 1}, {1, 3}};
  std::vector<SparseVector> sparse_columns = {sparse_col_1, sparse_col_2};

  // load problem
  ASSERT_OK(lp_interface_->LoadSparseColumnLP(
      LPObjectiveSense::kMaximization, objective_coefficients, lower_bounds,
      upper_bounds, empty_names, left_hand_sides, right_hand_sides, empty_names,
      sparse_columns));

  // set objective limit
  ASSERT_TRUE(
      (lp_interface_->SetIntegerParameter(LPParameter::kFromScratch, 1) ==
       absl::OkStatus()) ||
      (lp_interface_->SetIntegerParameter(LPParameter::kFromScratch, 1) ==
       absl::Status(absl::StatusCode::kInvalidArgument, "Parameter Unknown")));
  ASSERT_TRUE(
      (lp_interface_->SetIntegerParameter(LPParameter::kPresolving, 0) ==
       absl::OkStatus()) ||
      (lp_interface_->SetIntegerParameter(LPParameter::kPresolving, 0) ==
       absl::Status(absl::StatusCode::kInvalidArgument, "Parameter Unknown")));
  ASSERT_OK(lp_interface_->SetRealParameter(LPParameter::kObjectiveLimit, 0.0));

  // set basis
  ASSERT_OK(
      lp_interface_->SetBasisStatus(column_basis_status, row_basis_status));

  // solve problem
  ASSERT_OK(lp_interface_->SolveLPWithDualSimplex());

  // check status
  ASSERT_TRUE(lp_interface_->IsSolved());
  ASSERT_TRUE(lp_interface_->ObjectiveLimitIsExceeded() ||
              lp_interface_->IsOptimal());
  ASSERT_TRUE(!lp_interface_->IterationLimitIsExceeded());
  ASSERT_TRUE(!lp_interface_->TimeLimitIsExceeded());

  // the objective should be equal to the objective limit
  objective_value = lp_interface_->GetObjectiveValue();
  ASSERT_GE(
      objective_value,
      expected_objective_value);  // << "Objective value not equal to objective
                                  // limit: %g != %g\n", objective_value,
                                  // expected_objective_value);

  // check that data stored in lpi is still the same
  ASSERT_NO_FATAL_FAILURE(checkData(
      LPObjectiveSense::kMaximization, 2, objective_coefficients, lower_bounds,
      upper_bounds, 2, left_hand_sides, right_hand_sides, 4, sparse_columns));
}

// Test 6: More complex example
//
// The original problem was the following (arising from the qpkktref unit test),
// which displays a bug in CPLEX 12.7.0 w.r.t. scaling:
//  Minimize t_objvar
//  Subject To
//    KKTBinary1_y:                 - t_dual_y_bin1 + t_dual_y_bin2 +
//    t_dual_y_slackbin1 = 0 KKTlin_lower_1:               - t_x - t_y +
//    t_slack_lhs_lower + t_slack_ub_z       = 0.75 KKTBinary1_x: -
//    t_dual_x_bin1 + t_dual_x_bin2 + t_dual_x_slackbin1 = 0 KKTlin_lower_0: -
//    t_x - t_y + t_slack_ub_z - t_slack_rhs_lower       = 0.25
//    quadratic_side1_estimation_0: 2.75 t_x - 3.75 t_y + t_objvar + 2.28
//    t_slack_ub_z  <= 5.0496 quadratic_side0_estimation_0: 1.25 t_x - 0.25 t_y
//    + t_objvar + 2 t_slack_ub_z     >= 2.6875 quadratic_side1_estimation_0:
//    0.75 t_x - 0.25 t_y + t_objvar + 0.68 t_slack_ub_z  <= 4.2056
//    quadratic_side0_estimation_0: 2.75 t_x - 0.25 t_y + t_objvar + 3
//    t_slack_ub_z     >= 4.4375
//  Bounds
//    t_x = 1
//    t_y = 0
//    -2.562500001 <= t_objvar <= -0.0624999989999999
//    0 <= t_slack_lhs_lower <= 0.5
//    t_dual_x_bin1 Free
//    t_dual_x_bin2 Free
//    t_dual_x_slackbin1 = 0
//    t_dual_y_bin1 = 0
//    t_dual_y_bin2 Free
//    t_dual_y_slackbin1 Free
//    1.25 <= t_slack_ub_z <= 1.75
//    0 <= t_slack_rhs_lower <= 0.5
//  End
//
// We use the following mapping between variables and indices:
// 0:  t_x = 1
// 1:  t_y = 0
// 2:  t_objvar
// 3:  t_slack_lhs_lower
// 4:  t_dual_x_bin1
// 5:  t_dual_x_bin2
// 6:  t_dual_x_slackbin1
// 7:  t_dual_y_bin1
// 8:  t_dual_y_bin2
// 9:  t_dual_y_slackbin1
// 10: t_slack_ub_z
// 11: t_slack_rhs_lower
TEST_F(Solve, test6) {
  // LP data
  std::vector<double> objective_coefficients = {0, 0, 1, 0, 0, 0,
                                                0, 0, 0, 0, 0, 0};
  std::vector<double> lower_bounds = {1,   0,   -2.5625, 0,     -1e20, -1e20,
                                      0.0, 0.0, -1e20,   -1e20, 1.25,  0};
  std::vector<double> upper_bounds = {1,   0,   -0.0625, 0.5,  1e20, 1e20,
                                      0.0, 0.0, 1e20,    1e20, 1.75, 0.5};

  std::vector<double> left_hand_sides  = {0,     0.75,   0,     0.25,
                                         -1e20, 2.6875, -1e20, 4.4375};
  std::vector<double> right_hand_sides = {0,      0.75, 0,      0.25,
                                          5.0496, 1e20, 4.2056, 1e20};

  SparseVector sparse_col_1  = {{1, 3, 4, 5, 6, 7},
                               {-1, -1, 2.75, 1.25, 0.75, 2.75}};
  SparseVector sparse_col_2  = {{1, 3, 4, 5, 6, 7},
                               {-1, -1, -3.75, -0.25, -0.25, -0.25}};
  SparseVector sparse_col_3  = {{4, 5, 6, 7}, {1, 1, 1, 1}};
  SparseVector sparse_col_4  = {{1}, {1}};
  SparseVector sparse_col_5  = {{2}, {-1}};
  SparseVector sparse_col_6  = {{2}, {1}};
  SparseVector sparse_col_7  = {{2}, {1}};
  SparseVector sparse_col_8  = {{0}, {-1}};
  SparseVector sparse_col_9  = {{0}, {1}};
  SparseVector sparse_col_10 = {{0}, {1}};
  SparseVector sparse_col_11 = {{1, 3, 4, 5, 6, 7}, {1, 1, 2.28, 2, 0.68, 3}};
  SparseVector sparse_col_12 = {{3}, {-1.0}};
  std::vector<SparseVector> sparse_columns = {
      sparse_col_1, sparse_col_2,  sparse_col_3,  sparse_col_4,
      sparse_col_5, sparse_col_6,  sparse_col_7,  sparse_col_8,
      sparse_col_9, sparse_col_10, sparse_col_11, sparse_col_12};

  // possibly convert |1e20| to infinity of LPI
  for (int j = 0; j < 12; ++j) {
    if (lower_bounds[j] == -1e20) lower_bounds[j] = -lp_interface_->Infinity();
    if (upper_bounds[j] == 1e20) upper_bounds[j] = lp_interface_->Infinity();
  }
  for (int j = 0; j < 8; ++j) {
    if (left_hand_sides[j] == -1e20)
      left_hand_sides[j] = -lp_interface_->Infinity();
    if (right_hand_sides[j] == 1e20)
      right_hand_sides[j] = lp_interface_->Infinity();
  }
  // empty placeholders
  std::vector<std::string> empty_names;

  // load problem
  ASSERT_OK(lp_interface_->LoadSparseColumnLP(
      LPObjectiveSense::kMinimization, objective_coefficients, lower_bounds,
      upper_bounds, empty_names, left_hand_sides, right_hand_sides, empty_names,
      sparse_columns));

  // set some parameters - simulate settings in MiniMIP
  ASSERT_TRUE(
      (lp_interface_->SetIntegerParameter(LPParameter::kFromScratch, 0) ==
       absl::OkStatus()) ||
      (lp_interface_->SetIntegerParameter(LPParameter::kFromScratch, 0) ==
       absl::Status(absl::StatusCode::kInvalidArgument, "Parameter Unknown")));
  ASSERT_TRUE(
      (lp_interface_->SetIntegerParameter(LPParameter::kScaling, 1) ==
       absl::OkStatus()) ||
      (lp_interface_->SetIntegerParameter(LPParameter::kScaling, 1) ==
       absl::Status(absl::StatusCode::kInvalidArgument, "Parameter Unknown")));
  ASSERT_TRUE(
      (lp_interface_->SetIntegerParameter(LPParameter::kPresolving, 1) ==
       absl::OkStatus()) ||
      (lp_interface_->SetIntegerParameter(LPParameter::kPresolving, 1) ==
       absl::Status(absl::StatusCode::kInvalidArgument, "Parameter Unknown")));
  ASSERT_TRUE(
      (lp_interface_->SetIntegerParameter(LPParameter::kPricing, 0) ==
       absl::OkStatus()) ||
      (lp_interface_->SetIntegerParameter(LPParameter::kPricing, 0) ==
       absl::Status(absl::StatusCode::kInvalidArgument, "Parameter Unknown")));

  ASSERT_TRUE(
      (lp_interface_->SetRealParameter(LPParameter::kFeasibilityTolerance,
                                       1e-06) == absl::OkStatus()) ||
      (lp_interface_->SetRealParameter(LPParameter::kFeasibilityTolerance,
                                       1e-06) ==
       absl::Status(absl::StatusCode::kInvalidArgument, "Parameter Unknown")));
  ASSERT_TRUE(
      (lp_interface_->SetRealParameter(LPParameter::kDualFeasibilityTolerance,
                                       1e-07) == absl::OkStatus()) ||
      (lp_interface_->SetRealParameter(LPParameter::kDualFeasibilityTolerance,
                                       1e-07) ==
       absl::Status(absl::StatusCode::kInvalidArgument, "Parameter Unknown")));

  ASSERT_OK(lp_interface_->ClearState());

  // set objlimit
  ASSERT_OK(lp_interface_->SetRealParameter(LPParameter::kObjectiveLimit,
                                            4.320412501));

  // solve problem
  ASSERT_OK(lp_interface_->SolveLPWithDualSimplex());

  // check status
  ASSERT_TRUE(lp_interface_->IsSolved());

  // expected objective
  double expected_objective_value = -2.0625;

  // the objective should be equal to the objective limit
  ASSERT_FLOAT_EQ(lp_interface_->GetObjectiveValue(), expected_objective_value);

  // check that data stored in lpi is still the same
  ASSERT_NO_FATAL_FAILURE(checkData(
      LPObjectiveSense::kMinimization, 12, objective_coefficients, lower_bounds,
      upper_bounds, 8, left_hand_sides, right_hand_sides, 30, sparse_columns));

  // change some bounds
  lower_bounds[0]  = 1;
  upper_bounds[0]  = 1;
  lower_bounds[1]  = 0;
  upper_bounds[1]  = 0;
  lower_bounds[2]  = -2.06255;
  upper_bounds[2]  = -2.0625;
  lower_bounds[3]  = 0;
  upper_bounds[3]  = 4.94694e-05;
  lower_bounds[6]  = 0;
  upper_bounds[6]  = 0;
  lower_bounds[7]  = 0;
  upper_bounds[7]  = 0;
  lower_bounds[10] = 1.74995;
  upper_bounds[10] = 1.750;
  lower_bounds[11] = 0.499951;
  upper_bounds[11] = 0.5;

  for (size_t j = 0; j < sparse_columns.size(); ++j) {
    ASSERT_OK(
        lp_interface_->SetColumnBounds(j, lower_bounds[j], upper_bounds[j]))
  }
  // set objlimit
  ASSERT_OK(
      lp_interface_->SetRealParameter(LPParameter::kObjectiveLimit, -2.0625));

  ASSERT_OK(lp_interface_->ClearState());

  // solve problem
  ASSERT_OK(lp_interface_->SolveLPWithDualSimplex());

  // check that data stored in lpi is still the same
  ASSERT_NO_FATAL_FAILURE(checkData(
      LPObjectiveSense::kMinimization, 12, objective_coefficients, lower_bounds,
      upper_bounds, 8, left_hand_sides, right_hand_sides, 30, sparse_columns));
}

// Test 7
//
// min 10 x1 + 15 x2
//      2 x1 +   x2 >= 3
//        x1 + 3 x2 <= 1
//        x1,    x2 >= 0
//
// which is dual unbounded (this is a variant of Test 3 in which the equations
// have been replaced by inequalities).
//
// The dual is:
// max  3 y1 +   y2
//      2 y1 +   y2 <= 10
//        y1 + 3 y2 <= 15
//        y1 >= 0, y2 <= 0
TEST_F(Solve, test7) {
  // LP data
  std::vector<double> objective_coefficients{10, 15};
  std::vector<double> lower_bounds     = {0, 0};
  std::vector<double> upper_bounds     = {lp_interface_->Infinity(),
                                      lp_interface_->Infinity()};
  std::vector<double> right_hand_sides = {lp_interface_->Infinity(), 1};
  std::vector<double> left_hand_sides  = {3, -lp_interface_->Infinity()};

  SparseVector sparse_col_1                = {{0, 1}, {2, 1}};
  SparseVector sparse_col_2                = {{0, 1}, {1, 3}};
  std::vector<SparseVector> sparse_columns = {sparse_col_1, sparse_col_2};

  // expected ray
  std::vector<double> expected_dual_ray(2);
  expected_dual_ray = {0.5, -1};

  // empty placeholders
  std::vector<double> empty_vals;

  // check problem with primal simplex
  ASSERT_NO_FATAL_FAILURE(
      performTest(true, LPObjectiveSense::kMinimization, 2,
                  objective_coefficients, lower_bounds, upper_bounds, 2,
                  left_hand_sides, right_hand_sides, sparse_columns,
                  LPFeasibilityStat::INFEASIBLE, LPFeasibilityStat::UNBOUNDED,
                  empty_vals, expected_dual_ray, empty_vals, empty_vals));

  // check that data stored in lpi is still the same
  ASSERT_NO_FATAL_FAILURE(checkData(
      LPObjectiveSense::kMinimization, 2, objective_coefficients, lower_bounds,
      upper_bounds, 2, left_hand_sides, right_hand_sides, 4, sparse_columns));

  // clear basis status
  ASSERT_OK(lp_interface_->ClearState());

  // check problem with dual simplex
  ASSERT_NO_FATAL_FAILURE(solveTest(false, 2, 2, LPFeasibilityStat::INFEASIBLE,
                                    LPFeasibilityStat::UNBOUNDED, empty_vals,
                                    expected_dual_ray, empty_vals, empty_vals));

  // check that data stored in lpi is still the same
  ASSERT_NO_FATAL_FAILURE(checkData(
      LPObjectiveSense::kMinimization, 2, objective_coefficients, lower_bounds,
      upper_bounds, 2, left_hand_sides, right_hand_sides, 4, sparse_columns));
}
}  // namespace minimip
