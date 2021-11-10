
// Unit tests for testing the methods that change and get coefficients
// The tested methods are:
// GetCoefficient,
// ChangeObjective, GetObjective,
// ChangeBounds, GetBounds,
// ChangeSides, GetSides,
// ChangeObjectiveSense, GetObjectiveSense,
// GetNumberOfColumns, GetNumberOfRows, GetNumberOfNonZeros,
// GetSparseColumnCoefficients, GetSparseRowCoefficients,
// ClearState,
// WriteLP, ReadLP, Clear

#include <gtest/gtest.h>

#include "absl/status/status.h"
#include "src/lp_interface/lpi_factory.h"
#include "unit_tests/utils.h"

#define TEST_ERRORS 0  // if 0 then skip tests expected to fail.
#define DEF_INTERFACE \
  1  // 0 = Glop Interface (Default),
     // 1 = SoPlex Interface,

namespace minimip {

// TEST SUITE SOLVE

static LPInterface* lp_interface_ = nullptr;

class Change : public ::testing::Test {
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
    lp_interface_->SetObjectiveSense(LPObjectiveSense::kMaximization);
  }

  // write num_columns, num_rows and objective_sense into variables to check later
  static void initProb(int pos, int& num_columns, int& num_rows, int& num_nonzeros,
                       LPObjectiveSense& objective_sense) {
    std::vector<double> objective_coefficients = {1.0, 1.0};
    std::vector<double> lower_bounds  = {0.0, 0.0};
    std::vector<double> upper_bounds  = {lp_interface_->Infinity(),
                              lp_interface_->Infinity()};
    std::vector<double> left_hand_sides = {-lp_interface_->Infinity(),
                               -lp_interface_->Infinity()};
    std::vector<double> right_hand_sides = {1.0, 1.0};
    std::vector<SparseVector> sparse_rows;

    // maximization problems, num_columns is 1, num_rows is 1
    switch (pos) {
      case 0:
        // unbounded - infeasible
        // (P):  max x
        // -x <= 1 (constr)
        //  0 <= x (bound)
        //
        // (D):  min y
        // 1 <= -y (constr)
        // 0 <= y (bound)
        num_columns  = 1;
        num_rows  = 1;
        num_nonzeros  = 1;
        objective_sense = LPObjectiveSense::kMaximization;
        sparse_rows = {{{0},{-1}}};
        break;

      case 1:
        // optimal - optimal
        // (P):  max x
        //  x <= 0 (constr)
        //  0 <= x (bound)
        //
        // (D):  min 0
        // 1 <= y (constr)
        // 0 <= y (bound)
        num_columns  = 1;
        num_rows  = 1;
        num_nonzeros  = 1;
        objective_sense = LPObjectiveSense::kMaximization;
        sparse_rows = {{{0},{1}}};
        right_hand_sides[0] = 0.0;
        break;

      case 2:
        // minimization problems (duals of the above)
        num_columns  = 1;
        num_rows  = 1;
        num_nonzeros  = 1;
        objective_sense = LPObjectiveSense::kMinimization;
        right_hand_sides[0] = lp_interface_->Infinity();
        left_hand_sides[0] = 1;
        sparse_rows = {{{0},{-1}}};
        break;

      case 3:
        num_columns  = 1;
        num_rows  = 1;
        num_nonzeros  = 1;
        objective_sense = LPObjectiveSense::kMinimization;
        right_hand_sides[0] = lp_interface_->Infinity();
        left_hand_sides[0] = 1;
        objective_coefficients[0] = 0.0;
        sparse_rows = {{{0},{1}}};
        break;

      case 4:
        // maximization problems, num_columns is 2, *num_rows is 2
        // unbounded - infeasible
        // (P):  max x+y
        // -x    <= 1 (constr)
        //    -y <= 1 (constr)
        //
        //  0 <= x (bound)
        //  0 <= y (bound)
        //
        // (D):  min x+y
        // 1 <= -x   (constr)
        // 1 <=   -y (constr)
        //
        // 0 <= x (bound)
        // 0 <= y (bound)
        num_columns  = 2;
        num_rows  = 2;
        num_nonzeros  = 2;
        objective_sense = LPObjectiveSense::kMaximization;
        sparse_rows = {{{0},{-1}}, {{1},{-1}}};
        break;

      case 5:
        // optimal - optimal
        // (P):  max x+y
        // x     <= 1 (constr)
        //     y <= 1 (constr)
        //
        //  0 <= x (bound)
        //  0 <= y (bound)
        //
        // (D):  min x+y
        // 1 <= x    (constr)
        // 1 <=    y (constr)
        //
        // 0 <= x (bound)
        // 0 <= y (bound)
        num_columns  = 2;
        num_rows  = 2;
        num_nonzeros  = 2;
        objective_sense = LPObjectiveSense::kMaximization;
        sparse_rows = {{{0},{1}}, {{1},{1}}};
        break;

      case 6:
        // infeasible - infeasible
        // (P):  max x+y
        // -x    <= -1 (constr)
        //     y <= -1 (constr)
        //
        //  0 <= x (bound)
        //  0 <= y (bound)
        //
        // (D):  min -x-y
        // 1 <= -x    (constr)
        // 1 <=     y (constr)
        //
        // 0 <= x (bound)
        // 0 <= y (bound)
        num_columns  = 2;
        num_rows  = 2;
        num_nonzeros  = 2;
        objective_sense = LPObjectiveSense::kMaximization;
        right_hand_sides[0] = -1.0;
        right_hand_sides[1] = -1.0;
        sparse_rows = {{{0},{-1}}, {{1},{1}}};
        break;

      case 7:
        // minimization problems (duals of the above)
        num_columns  = 2;
        num_rows  = 2;
        num_nonzeros  = 2;
        objective_sense = LPObjectiveSense::kMinimization;
        right_hand_sides[0] = lp_interface_->Infinity();
        right_hand_sides[1] = lp_interface_->Infinity();
        left_hand_sides[0] = 1.0;
        left_hand_sides[1] = 1.0;
        sparse_rows = {{{0},{-1}}, {{1},{-1}}};
        break;

      case 8:
        num_columns  = 2;
        num_rows  = 2;
        num_nonzeros  = 2;
        objective_sense = LPObjectiveSense::kMinimization;
        right_hand_sides[0] = lp_interface_->Infinity();
        right_hand_sides[1] = lp_interface_->Infinity();
        left_hand_sides[0] = 1.0;
        left_hand_sides[1] = 1.0;
        sparse_rows = {{{0},{1}}, {{1},{1}}};
        break;

      case 9:
        num_columns  = 2;
        num_rows  = 2;
        num_nonzeros  = 2;
        objective_sense = LPObjectiveSense::kMinimization;
        right_hand_sides[0] = lp_interface_->Infinity();
        right_hand_sides[1] = lp_interface_->Infinity();
        left_hand_sides[0] = 1.0;
        left_hand_sides[1] = 1.0;
        objective_coefficients[0] = -1.0;
        objective_coefficients[1] = -1.0;
        sparse_rows = {{{0},{-1}}, {{1},{1}}};
        break;

      default:
        return;
    }
    // empty placeholders
    std::vector<std::string> empty_names;
    std::vector<double> empty_values;
    std::vector<int> empty_indices;

    ASSERT_OK(lp_interface_->SetObjectiveSense(objective_sense));
    
    std::vector<SparseVector> empty_columns;
    for (int i = 0; i < num_columns; ++i){
      empty_columns.push_back({{},{}});
    }

    ASSERT_OK(lp_interface_->AddColumns(empty_columns, lower_bounds, upper_bounds, objective_coefficients, empty_names));

    ASSERT_OK(lp_interface_->AddRows(sparse_rows, left_hand_sides, right_hand_sides, empty_names));
    ASSERT_TRUE(!lp_interface_->IsSolved());
    ASSERT_OK(lp_interface_->SolveLPWithPrimalSimplex());
    ASSERT_TRUE(lp_interface_->IsSolved());
  }

  // We treat -2 and 2 as -infinity and infinity resp.
  static double substituteInfinity(double inf) {
    if (inf == 2)
      return lp_interface_->Infinity();
    else if (inf == -2)
      return -lp_interface_->Infinity();
    else
      return inf;
  }
};

// TESTS

// Test ChangeObjectives
class ChangeObjective
    : public Change,
      public ::testing::WithParamInterface<std::tuple<double, double, int>> {
 protected:
  static void checkChangeObjective(int& last_col, std::vector<int>& ind,
                          std::vector<double>& set_objective_coefficients) {
    std::vector<double> objective_coefficients(2);
    ASSERT_LE(last_col, 2);

    std::vector<int> current_indices;
    current_indices.reserve(last_col);
    for (int i = 0; i < last_col; i++) {
      ASSERT_OK(lp_interface_->SetObjectiveCoefficient(ind[i], set_objective_coefficients[i]));
    }
    ASSERT_TRUE(!lp_interface_->IsSolved());
    for (int i = 0; i < last_col; i++) {
      objective_coefficients[i] = lp_interface_->GetObjectiveCoefficient(i);
    }

    for (size_t i; i < set_objective_coefficients.size(); i++) {
      ASSERT_EQ(objective_coefficients[i], set_objective_coefficients[i]);
    }
  }
};

TEST_P(ChangeObjective, checkChangeObjective) {
  int prob = std::get<2>(GetParam());
  int num_rows, num_columns, num_nonzeros;
  LPObjectiveSense objective_sense;
  std::vector<int> ind = {0, 1};
  std::vector<double> set_objective_coefficients(2);
  bool deathflag = false;

  set_objective_coefficients[0] = substituteInfinity(std::get<0>(GetParam()));
  set_objective_coefficients[1] = substituteInfinity(std::get<1>(GetParam()));

  ASSERT_NO_FATAL_FAILURE(initProb(prob, num_columns, num_rows, num_nonzeros, objective_sense));

  if (DEF_INTERFACE == 0) {
    if (TEST_ERRORS) {
      for (int i = 0; i < num_columns; i++) {
        if (lp_interface_->IsInfinity(fabs(set_objective_coefficients[i]))) deathflag = true;
      }
    } else {
      for (int i = 0; i < num_columns; i++) {
        if (lp_interface_->IsInfinity(fabs(set_objective_coefficients[i]))) GTEST_SKIP();
      }
    }
    if (deathflag)
      ASSERT_DEATH(checkChangeObjective(num_columns, ind, set_objective_coefficients), "");
    else
      ASSERT_NO_FATAL_FAILURE(checkChangeObjective(num_columns, ind, set_objective_coefficients));
  } else {
    ASSERT_NO_FATAL_FAILURE(checkChangeObjective(num_columns, ind, set_objective_coefficients));
  }
}

INSTANTIATE_TEST_SUITE_P(ChangeObjectivesCombinations, ChangeObjective,
                         ::testing::Combine(::testing::Values(0, 1, -1, 2, -2),
                                            ::testing::Values(0, 1, -1, 2, -2),
                                            ::testing::Values(0, 1, 2, 3, 4, 5,
                                                              6, 7, 8, 9)));

// Test ChangeBounds
class ChangeBounds : public Change,
                     public ::testing::WithParamInterface<
                         std::tuple<double, double, double, double, int>> {
 protected:
  static void CheckChangeBounds(int& last_col, std::vector<int>& ind,
                             std::vector<double>& set_lower_bounds,
                             std::vector<double>& set_upper_bounds,
                             absl::Status error = absl::OkStatus()) {
    std::vector<double> upper_bounds(2);
    std::vector<double> lower_bounds(2);

    ASSERT_LE(last_col, 2);
    if (error != absl::OkStatus()) {
      for (int i = 0; i < last_col; ++i) {
      ASSERT_EQ(lp_interface_->SetColumnBounds(i, set_lower_bounds[i], set_upper_bounds[i]),
                error);
      }
      abort();
    } else {
      for (int i = 0; i < last_col; ++i) {
      ASSERT_OK(lp_interface_->SetColumnBounds(i, set_lower_bounds[i], set_upper_bounds[i]));
      }
    }
    ASSERT_TRUE(!lp_interface_->IsSolved());
    for (int i = 0; i < last_col; i++) {
      lower_bounds[i] = lp_interface_->GetLowerBound(i);
      upper_bounds[i] = lp_interface_->GetUpperBound(i);
    }

    for (size_t i; i < last_col; i++) {
      ASSERT_EQ(upper_bounds[i], set_upper_bounds[i]);
    }
    for (size_t i; i < last_col; i++) {
      ASSERT_EQ(lower_bounds[i], set_lower_bounds[i]);
    }
  }
};

TEST_P(ChangeBounds, CheckChangeBounds) {
  int prob = std::get<4>(GetParam());

  int num_rows, num_columns, num_nonzeros;
  std::vector<int> ind = {0, 1};
  LPObjectiveSense objective_sense;
  std::vector<double> set_upper_bounds(2);
  std::vector<double> set_lower_bounds(2);
  bool deathflag = false;

  set_lower_bounds[0] = substituteInfinity(std::get<0>(GetParam()));
  set_lower_bounds[1] = substituteInfinity(std::get<1>(GetParam()));
  set_upper_bounds[0] = substituteInfinity(std::get<2>(GetParam()));
  set_upper_bounds[1] = substituteInfinity(std::get<3>(GetParam()));

  ASSERT_NO_FATAL_FAILURE(initProb(prob, num_columns, num_rows, num_nonzeros, objective_sense));

  if (TEST_ERRORS) {
    for (int i = 0; i < num_columns; i++) {
      if (set_upper_bounds[i] < set_lower_bounds[i]) deathflag = true;
      if (lp_interface_->IsInfinity(set_lower_bounds[i]) ||
          lp_interface_->IsInfinity(-set_upper_bounds[i]))
        deathflag = true;
    }
  } else {
    for (int i = 0; i < num_columns; i++) {
      if (set_upper_bounds[i] < set_lower_bounds[i]) GTEST_SKIP();
      if (lp_interface_->IsInfinity(set_lower_bounds[i]) ||
          lp_interface_->IsInfinity(-set_upper_bounds[i]))
        GTEST_SKIP();
    }
  }

  if (deathflag)
    ASSERT_DEATH(
        CheckChangeBounds(num_columns, ind, set_lower_bounds, set_upper_bounds,
                       absl::Status(absl::StatusCode::kInternal, "LP Error")),
        "");
  else
    ASSERT_NO_FATAL_FAILURE(CheckChangeBounds(num_columns, ind, set_lower_bounds, set_upper_bounds));
}

INSTANTIATE_TEST_SUITE_P(ChangeBoundsCombinations, ChangeBounds,
                         ::testing::Combine(::testing::Values(0, 1, -1, 2, -2),
                                            ::testing::Values(0, 1, -1, 2, -2),
                                            ::testing::Values(0, 1, -1, 2, -2),
                                            ::testing::Values(0, 1, -1, 2, -2),
                                            ::testing::Values(0, 1, 2, 3, 4, 5,
                                                              6, 7, 8, 9)));

// Test ChangeSides
class ChangeSides : public Change,
                    public ::testing::WithParamInterface<
                        std::tuple<double, double, double, double, int>> {
 protected:
  static void checkChangedSides(int& lastrow, std::vector<int>& ind,
                            std::vector<double>& set_left_hand_sides,
                            std::vector<double>& set_right_hand_sides) {
    std::vector<double> ls(2);
    std::vector<double> rs(2);

    ASSERT_LE(lastrow, 2);
    for (int i = 0; i < lastrow; i++) {
      ASSERT_OK(lp_interface_->SetRowSides(ind[i], set_left_hand_sides[i], set_right_hand_sides[i]));
      // ls.push_back(lp_interface_->GetLeftHandSide(i));
      // rs.push_back(lp_interface_->GetRightHandSide(i));
    }

    // ASSERT_OK(lp_interface_->SetRowSides(lastrow, ind, set_left_hand_sides, set_right_hand_sides));
    ASSERT_TRUE(!lp_interface_->IsSolved());

    for (int i = 0; i < lastrow; i++) {
      ls[i] = lp_interface_->GetLeftHandSide(i);
      rs[i] = lp_interface_->GetRightHandSide(i);
    }

    for (int i; i < lastrow; i++) {
      ASSERT_EQ(ls[i], set_left_hand_sides[i]);
    }
    for (int i; i < lastrow; i++) {
      ASSERT_EQ(rs[i], set_right_hand_sides[i]);
    }
  }
};

TEST_P(ChangeSides, checkChangedSides) {
  double left1  = substituteInfinity(std::get<0>(GetParam()));
  double left2  = substituteInfinity(std::get<1>(GetParam()));
  double right1 = substituteInfinity(std::get<2>(GetParam()));
  double right2 = substituteInfinity(std::get<3>(GetParam()));
  int prob      = std::get<4>(GetParam());

  int num_rows, num_columns, num_nonzeros;
  std::vector<int> ind = {0, 1};
  LPObjectiveSense objective_sense;
  std::vector<double> setlhs = {left1, left2};
  std::vector<double> setrhs = {right1, right2};
  bool deathflag             = false;

  ASSERT_NO_FATAL_FAILURE(initProb(prob, num_columns, num_rows, num_nonzeros, objective_sense));

  if (TEST_ERRORS) {
    for (int i = 0; i < num_rows; i++) {
      if (setrhs[i] < setlhs[i]) deathflag = true;
      if (lp_interface_->IsInfinity(fabs(setlhs[i])) &&
          lp_interface_->IsInfinity(fabs(setrhs[i])))
        deathflag = true;
    }
  } else {
    for (int i = 0; i < num_rows; i++) {
      if (setrhs[i] < setlhs[i]) GTEST_SKIP();
      if (lp_interface_->IsInfinity(fabs(setlhs[i])) &&
          lp_interface_->IsInfinity(fabs(setrhs[i])))
        GTEST_SKIP();
    }
  }

  if (deathflag)
    ASSERT_DEATH(checkChangedSides(num_rows, ind, setlhs, setrhs), "");
  else
    ASSERT_NO_FATAL_FAILURE(checkChangedSides(num_rows, ind, setlhs, setrhs));
}

INSTANTIATE_TEST_SUITE_P(ChangeSidesCombinations, ChangeSides,
                         ::testing::Combine(::testing::Values(0, 1, -1, 2, -2),
                                            ::testing::Values(0, 1, -1, 2, -2),
                                            ::testing::Values(0, 1, -1, 2, -2),
                                            ::testing::Values(0, 1, -1, 2, -2),
                                            ::testing::Values(0, 1, 2, 3, 4, 5,
                                                              6, 7, 8, 9)));

// Test ChangeObjectiveSense
class ChangeObjectiveSense
    : public Change,
      public ::testing::WithParamInterface<std::tuple<LPObjectiveSense, int>> {
};

TEST_P(ChangeObjectiveSense, ChangeObjectiveSense) {
  LPObjectiveSense newsense = std::get<0>(GetParam());
  int prob                  = std::get<1>(GetParam());

  int num_rows, num_columns, num_nonzeros;
  LPObjectiveSense objective_sense;
  LPObjectiveSense probsense;

  ASSERT_NO_FATAL_FAILURE(initProb(prob, num_columns, num_rows, num_nonzeros, objective_sense));

  ASSERT_OK(lp_interface_->SetObjectiveSense(newsense));
  ASSERT_TRUE(!lp_interface_->IsSolved());

  probsense = lp_interface_->GetObjectiveSense();

  ASSERT_EQ(newsense, probsense);
}

INSTANTIATE_TEST_SUITE_P(
    ChangeObjectiveSenseCombinations, ChangeObjectiveSense,
    ::testing::Combine(::testing::Values(LPObjectiveSense::kMaximization,
                                         LPObjectiveSense::kMinimization),
                       ::testing::Values(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)));

// Test for AddRows, DeleteRowSet, DeleteRows
// TEST_F(Change, TestRowMethods) {
//   // problem data
//   std::vector<double> objective_coefficients = {1.0, 1.0, 1.0, 1.0, 1.0};
//   std::vector<double> lower_bounds  = {-1.0, -lp_interface_->Infinity(), 0.0,
//                             -lp_interface_->Infinity(), 0.0};
//   std::vector<double> upper_bounds  = {10.0, lp_interface_->Infinity(),
//                             lp_interface_->Infinity(), 29.0, 0.0};
//   int ncolsbefore, ncolsafter;
//   int nrowsbefore, nrowsafter;
//   std::vector<double> lhsvals = {
//       -lp_interface_->Infinity(), -1.0, -3e-10, 0.0, 1.0, 3e10};
//   std::vector<double> rhsvals = {-1.0, -3e-10, 0.0,
//                                  1.0,  3e10,   lp_interface_->Infinity()};
//   std::vector<int> nnonzs = {1, 10, -1, 6, -1};
//   std::vector<int> begvals = {0, 2, 3, 5, 8, 9};
//   std::vector<int> indvals = {0, 1, 3, 2, 1, 1, 2, 4, 0, 3};
//   std::vector<double> vals = {1.0, 5.0, -1.0, 3e5,  2.0,
//                               1.0, 20,  10,   -1.9, 1e-2};

//   std::vector<SparseVector> sparse_rows = {{{0,1},{1.0,5.0}}, {{3},{-1.0}}, {{2,1},{3e5,2.0}}, 
//                                           {{1,2,4},{1.0,20,10}}, {{0},{-1.9}}, {{3},{1e-2}}};

//   int iterations              = 5;
//   std::vector<int> k          = {1, 6, -1, 4, -2};
//   std::vector<int> nnonzsdiff = {1, 10, -1, 6, -3};
//   int i;
//   int j;

//   // empty placeholders
//   std::vector<std::string> empty_names;
//   std::vector<double> empty_values;
//   std::vector<int> empty_indices;

//   // empty columns
//     std::vector<SparseVector> empty_columns;
//     // empty_columns.reserve(num_columns);
//     for (int i = 0; i < objective_coefficients.size(); ++i){
//       empty_columns.push_back({{},{}});
//     }

//   // create original lp
//   ASSERT_OK(lp_interface_->AddColumns(empty_columns, lower_bounds, upper_bounds, objective_coefficients, empty_names));
//   ncolsbefore = lp_interface_->GetNumberOfColumns();

//   for (i = 0; i < iterations; i++) {
//     // setup row values
//     int num_rows;
//     int nnonzsbefore;
//     int nnonzsafter;

//     num_rows = k[i];
//     // get data before modification
//     nnonzsbefore = lp_interface_->GetNumberOfNonZeros();
//     nrowsbefore  = lp_interface_->GetNumberOfRows();

//     if (num_rows < 0) {
//       ASSERT_OK(lp_interface_->DeleteRows(0, -(1 + num_rows)));
//     } else {  // num_rows >= 0
//       std::vector<double> left_hand_sides(100);
//       std::vector<double> right_hand_sides(100);
//       std::vector<int> beg(100);

//       int num_nonzeros = nnonzs[i];
//       std::vector<int> ind(100);
//       std::vector<double> val(100);

//       std::vector<double> newlhs(100);
//       std::vector<double> newval(100);
//       std::vector<double> newrhs(100);
//       std::vector<int> newbeg(100);
//       std::vector<int> newind(100);
//       int newnnonz;
//       int indold;
//       int indnew;

//       ASSERT_LE(num_rows, 100);
//       for (j = 0; j < num_rows; j++) {
//         left_hand_sides[j] = lhsvals[j];
//         right_hand_sides[j] = rhsvals[j];
//         beg[j] = begvals[j];
//       }
//       ASSERT_LE(num_nonzeros, 100);
//       for (j = 0; j < num_nonzeros; j++) {
//         ind[j] = indvals[j];
//         val[j] = vals[j];
//       }
//       ASSERT_OK(lp_interface_->AddRows(sparse_rows, left_hand_sides, right_hand_sides, empty_names));

//       // checks
//       std::vector<SparseVector> check_sparse_rows;
//       newnnonz = 0;

//       for (int i = nrowsbefore; i < nrowsbefore + num_rows -1; i++) {
//         // ASSERT_EQ(left_hand_sides[i], lp_interface_->GetLeftHandSide(i));
//         // ASSERT_EQ(right_hand_sides[i], lp_interface_->GetRightHandSide(i));
//         newlhs[i] = lp_interface_->GetLeftHandSide(i);
//         newrhs[i] = lp_interface_->GetRightHandSide(i);

//         check_sparse_rows.push_back(lp_interface_->GetSparseRowCoefficients(i));
//         int entries = static_cast<int>(check_sparse_rows[i].indices.size());
//         // for (int j = 0; j < entries; j++) {
//         //   newind[i] = check_sparse_rows[i].indices[j];
//         //   newval[i] = check_sparse_rows[i].values[j];
//         // }
//         newbeg.push_back(entries);
//         newnnonz += entries;
//       }
//       ASSERT_EQ(num_nonzeros, newnnonz);


//       for (j = 0; j < num_rows; j++) {
//         ASSERT_EQ(beg[j], newbeg[j]);
//       }
//       beg[num_rows]    = num_nonzeros;
//       newbeg[num_rows] = newnnonz;

//       // check each row seperately
//       for (j = 0; j < num_rows; j++) {
//         if (fabs(left_hand_sides[j]) < 1e30 && fabs(newlhs[j]) < 1e30) {
//           ASSERT_DOUBLE_EQ(left_hand_sides[j], newlhs[j]);
//         }
//         if (fabs(right_hand_sides[j]) < 1e30 && fabs(newrhs[j]) < 1e30) {
//           ASSERT_DOUBLE_EQ(right_hand_sides[j], newrhs[j]);
//         }

//     // compare sparse rows

//     ASSERT_EQ(expected_sparse_columns.size(), sparse_columns.size());

//     for (int j = 0; j < expected_sparse_columns.size(); ++j) {
//       ASSERT_EQ(expected_sparse_columns[j].indices.size(),
//                 sparse_columns[j].indices.size());
//       for (int i = 0; i < expected_sparse_columns[j].indices.size(); ++i) {
//         expected_sparse_columns[j].values[i];
//         ASSERT_EQ(expected_sparse_columns[j].indices[i],
//                   sparse_columns[j].indices[i]);
//         ASSERT_EQ(expected_sparse_columns[j].values[i],
//                   sparse_columns[j].values[i]);
//       }
//     }

//         // // We add a row where the indices are not sorted, some lp solvers give
//         // // them back sorted (e.g. soplex), some others don't (e.g. cplex).
//         // // Therefore we cannot simply assert the ind and val arrays to be equal,
//         // // but have to search for and check each value individually.
//         // for (indold = beg[j]; indold < beg[j + 1]; indold++) {
//         //   int occurrences = 0;

//         //   // for each value ind associated to the current row search for it in
//         //   // the newind array
//         //   for (indnew = beg[j]; indnew < beg[j + 1]; indnew++) {
//         //     if (ind[indold] == newind[indnew]) {
//         //       occurrences = occurrences + 1;
//         //       ASSERT_DOUBLE_EQ(val[indold], newval[indnew]);
//         //     }
//         //   }
//         //   // assert that we found only one occurrence in the current row
//         //   ASSERT_EQ(occurrences, 1);
//         // }

//     }

//     }

//     // checks
//     nrowsafter = lp_interface_->GetNumberOfRows();
//     ASSERT_EQ(nrowsbefore + num_rows, nrowsafter);

//     nnonzsafter = lp_interface_->GetNumberOfNonZeros();
//     ASSERT_EQ(nnonzsbefore + nnonzsdiff[i], nnonzsafter);

//     ncolsafter = lp_interface_->GetNumberOfColumns();
//     ASSERT_EQ(ncolsbefore, ncolsafter);
//   }

//   // delete rowsets
//   // should have 8 rows now
//   nrowsbefore = lp_interface_->GetNumberOfRows();
//   ASSERT_EQ(8, nrowsbefore);
//   for (i = 3; i > 0; i--) {
//     std::vector<bool> rows = {false, false, false, false,
//                               false, false, false, false};

//     for (j = 0; j < i; j++) rows[(2 * j) + 1] = true;

//     nrowsbefore = lp_interface_->GetNumberOfRows();
//     ASSERT_OK(lp_interface_->DeleteRowSet(rows));

//     nrowsafter = lp_interface_->GetNumberOfRows();
//     ASSERT_EQ(nrowsbefore - i, nrowsafter);
//     // assert that the rows that are left are the ones I intended
//   }
// }

// // Test for AddColumns, DeleteColumns
// TEST_F(Change, TestColumnMethods) {
//   // problem data
//   std::vector<double> objective_coefficients = {1.0, 1.0, 1.0, 1.0, 1.0};
//   std::vector<double> left_hand_sides = {-1.0, -lp_interface_->Infinity(), 0.0,
//                              -lp_interface_->Infinity(), 0.0};
//   std::vector<double> right_hand_sides = {10.0, lp_interface_->Infinity(),
//                              lp_interface_->Infinity(), 29.0, 0.0};
//   int ncolsbefore, ncolsafter;
//   int nrowsbefore, nrowsafter;
//   std::vector<double> lbvals = {
//       -lp_interface_->Infinity(), -1.0, -3e-10, 0.0, 1.0, 3e10};
//   std::vector<double> ubvals = {-1.0, -3e-10, 0.0,
//                                 1.0,  3e10,   lp_interface_->Infinity()};
//   std::vector<double> vals   = {1.0, 5.0, -1.0, 3e5,  2.0,
//                               1.0, 20,  10,   -1.9, 1e-2};
//   std::vector<int> nnonzs    = {1, 10, -1, 6, -1};
//   std::vector<int> begvals   = {0, 2, 3, 5, 8, 9};
//   std::vector<int> indvals   = {0, 1, 3, 2, 1, 1, 2, 4, 0, 3};

//   int iterations              = 5;
//   std::vector<int> k          = {1, 6, -1, 4, -2};
//   std::vector<int> nnonzsdiff = {1, 10, -1, 6, -3};

//   // empty placeholders
//   std::vector<std::string> empty_names;
//   std::vector<double> empty_values;
//   std::vector<int> empty_indices;

//   // create original lp
//   ASSERT_OK(lp_interface_->AddRows(5, left_hand_sides, right_hand_sides, empty_names, 0, empty_indices,
//                                    empty_indices, empty_values));
//   nrowsbefore = lp_interface_->GetNumberOfRows();

//   for (int i = 0; i < iterations; i++) {
//     // setup col values
//     int num_columns;
//     int nnonzsbefore;
//     int nnonzsafter;

//     num_columns = k[i];

//     // get data before modification
//     nnonzsbefore = lp_interface_->GetNumberOfNonZeros();
//     ncolsbefore  = lp_interface_->GetNumberOfColumns();

//     if (k[i] < 0) {
//       ASSERT_EQ(lp_interface_->DeleteColumns(0, -(1 + num_columns)),
//                 absl::OkStatus());
//     } else {  // num_columns >= 0
//       std::vector<double> lower_bounds(100);
//       std::vector<double> upper_bounds(100);
//       std::vector<int> beg(100);

//       assert(nnonzs[i] >= 0);
//       int num_nonzeros = nnonzs[i];
//       std::vector<int> ind(100);
//       std::vector<double> val(100);

//       std::vector<double> newlb(100);
//       std::vector<double> newval(100);
//       std::vector<double> newub(100);
//       std::vector<int> newbeg(100);
//       std::vector<int> newind(100);
//       int newnnonz;

//       ASSERT_LE(num_columns, 100);
//       for (int j = 0; j < num_columns; j++) {
//         lower_bounds[j]  = lbvals[j];
//         upper_bounds[j]  = ubvals[j];
//         beg[j] = begvals[j];
//       }

//       ASSERT_LE(num_nonzeros, 100);
//       for (int j = 0; j < num_nonzeros; j++) {
//         ind[j] = indvals[j];
//         val[j] = vals[j];
//       }
//       ASSERT_OK(lp_interface_->AddColumns(num_columns, objective_coefficients, lower_bounds, upper_bounds, empty_names,
//                                           num_nonzeros, beg, ind, val));

//       // checks
//       std::vector<SparseVector> sparse_rows;
//       newnnonz = 0;

//       for (int i = 0; i < ncolsbefore; i++) {
//         newlb.push_back(lp_interface_->GetLowerBound(i));
//         newub.push_back(lp_interface_->GetUpperBound(i));

//         sparse_rows.push_back(lp_interface_->GetSparseColumnCoefficients(i));
//         int entries = static_cast<int>(sparse_rows[i].indices.size());
//         for (int j = 0; j < entries; j++) {
//           newind.push_back(sparse_rows[i].indices[j]);
//           newval.push_back(sparse_rows[i].values[j]);
//         }
//         newbeg.push_back(entries);
//         newnnonz += entries;
//       }

//       ASSERT_EQ(num_nonzeros, newnnonz);

//       for (int j = 0; j < num_columns; j++) {
//         ASSERT_EQ(lower_bounds[j], newlb[j]);
//         ASSERT_EQ(upper_bounds[j], newub[j]);
//         ASSERT_EQ(beg[j], newbeg[j]);
//       }
//       for (int j = 0; j < num_nonzeros; j++) {
//         ASSERT_EQ(ind[j], newind[j]);
//         ASSERT_EQ(val[j], newval[j]);
//       }
//     }

//     // checks
//     nrowsafter = lp_interface_->GetNumberOfRows();
//     ASSERT_EQ(nrowsbefore, nrowsafter);

//     nnonzsafter = lp_interface_->GetNumberOfNonZeros();
//     ASSERT_EQ(nnonzsbefore + nnonzsdiff[i], nnonzsafter);

//     ncolsafter = lp_interface_->GetNumberOfColumns();
//     ASSERT_EQ(ncolsbefore + num_columns, ncolsafter);
//   }
// }

// // Test adding zero coeffs cols
// TEST_F(Change, testzerosincols) {
//   int num_columns;
//   int num_rows;
//   int num_nonzeros = 2;
//   LPObjectiveSense objective_sense;
//   std::vector<double> lower_bounds  = {0};
//   std::vector<double> upper_bounds  = {20};
//   std::vector<int> beg    = {0};
//   std::vector<int> ind    = {0, 1};
//   std::vector<double> val = {0, 3};
//   std::vector<double> objective_coefficients = {1};

//   // empty placeholders
//   std::vector<std::string> empty_names;

//   // 2x2 problem
//   initProb(4, num_columns, num_rows, num_nonzeros, objective_sense);
//   ASSERT_EQ(2, num_rows);
//   ASSERT_EQ(2, num_columns);

// #ifndef NDEBUG
//   ASSERT_DEATH(lp_interface_->AddColumns(1, objective_coefficients, lower_bounds, upper_bounds, empty_names, num_nonzeros,
//                                          beg, ind, val),
//                "");
// #endif
//   // this test can only work in debug mode, so we make it pass in opt mode
// #ifdef NDEBUG
//   ASSERT_OK(lp_interface_->AddColumns(1, objective_coefficients, lower_bounds, upper_bounds, empty_names, num_nonzeros, beg,
//                                       ind, val));
//   SUCCEED();  // return SIGABORT
// #endif
// }

// // Test adding zero coeffs in rows, expecting an assert in debug mode
// //
// // This test should fail with an assert from the which causes SIGABRT to be
// // issued. Thus, this test should pass.
// TEST_F(Change, testzerosinrows) {
//   int num_rows;
//   int num_columns;
//   int num_nonzeros = 2;
//   LPObjectiveSense objective_sense;
//   std::vector<double> left_hand_sides = {0};
//   std::vector<double> right_hand_sides = {20};
//   std::vector<int> beg    = {0};
//   std::vector<int> ind    = {0, 1};
//   std::vector<double> val = {0, 3};
//   // empty placeholders
//   std::vector<std::string> empty_names;

//   // 2x2 problem
//   initProb(4, num_columns, num_rows, num_nonzeros, objective_sense);
//   ASSERT_EQ(2, num_rows);
//   ASSERT_EQ(2, num_columns);

// #ifndef NDEBUG
//   ASSERT_DEATH(
//       lp_interface_->AddRows(1, left_hand_sides, right_hand_sides, empty_names, num_nonzeros, beg, ind, val),
//       "");
// #else
//   // this test can only work in debug mode, so we make it pass in opt mode
//   ASSERT_OK(
//       lp_interface_->AddRows(1, left_hand_sides, right_hand_sides, empty_names, num_nonzeros, beg, ind, val));
//   SUCCEED();
// #endif
// }

// // test WriteLP, ReadLP, Clear
// TEST_F(Change, testlpiwritereadlpmethods) {
//   int num_rows, num_columns, num_nonzeros;
//   double objval;
//   std::vector<double> primsol(2);
//   std::vector<double> dualsol(2);
//   std::vector<double> activity(2);
//   std::vector<double> redcost(2);
//   double objval2;
//   std::vector<double> primsol2(2);
//   std::vector<double> dualsol2(2);
//   std::vector<double> activity2(2);
//   std::vector<double> redcost2(2);
//   LPObjectiveSense objective_sense;

//   // 2x2 problem
//   ASSERT_NO_FATAL_FAILURE(initProb(5, num_columns, num_rows, num_nonzeros, objective_sense));

//   ASSERT_OK(lp_interface_->SolveLPWithPrimalSimplex());
//   objval = lp_interface_->GetObjectiveValue();
//   absl::StatusOr<std::vector<double>> absl_tmp;

//   absl_tmp = lp_interface_->GetPrimalSolution();
//   if (absl_tmp.ok())
//     primsol = *absl_tmp;
//   else
//     std::cout << absl_tmp.status();

//   absl_tmp = lp_interface_->GetRowActivity();
//   if (absl_tmp.ok())
//     activity = *absl_tmp;
//   else
//     std::cout << absl_tmp.status();

//   absl_tmp = lp_interface_->GetDualSolution();
//   if (absl_tmp.ok())
//     dualsol = *absl_tmp;
//   else
//     std::cout << absl_tmp.status();

//   absl_tmp = lp_interface_->GetReducedCost();
//   if (absl_tmp.ok())
//     redcost = *absl_tmp;
//   else
//     std::cout << absl_tmp.status();

//   ASSERT_EQ(lp_interface_->WriteLP("lpi_change_test_problem.lp"),
//             absl::OkStatus());
//   ASSERT_OK(lp_interface_->Clear());

//   if (DEF_INTERFACE == 0) ASSERT_OK(lp_interface_->ReadLP("lpi_change_test_problem.lp.gz"));
//   else ASSERT_OK(lp_interface_->ReadLP("lpi_change_test_problem.lp"));

//   ASSERT_OK(lp_interface_->SolveLPWithPrimalSimplex());
//   objval2  = lp_interface_->GetObjectiveValue();
//   absl_tmp = lp_interface_->GetPrimalSolution();
//   if (absl_tmp.ok())
//     primsol2 = *absl_tmp;
//   else
//     std::cout << absl_tmp.status();

//   absl_tmp = lp_interface_->GetRowActivity();
//   if (absl_tmp.ok())
//     activity2 = *absl_tmp;
//   else
//     std::cout << absl_tmp.status();

//   absl_tmp = lp_interface_->GetDualSolution();
//   if (absl_tmp.ok())
//     dualsol2 = *absl_tmp;
//   else
//     std::cout << absl_tmp.status();

//   absl_tmp = lp_interface_->GetReducedCost();
//   if (absl_tmp.ok())
//     redcost2 = *absl_tmp;
//   else {
//     std::cout << absl_tmp.status();  //@TODO: should break for all
//                                      //absl_tmp.status() return
//   }

//   ASSERT_FLOAT_EQ(objval, objval2);

//   ASSERT_EQ(primsol.size(), primsol2.size());
//   for (size_t j = 0; j < primsol.size(); j++) {
//     ASSERT_EQ(primsol[j], primsol2[j]);
//   }

//   ASSERT_EQ(dualsol.size(), dualsol2.size());
//   for (size_t j = 0; j < dualsol.size(); j++) {
//     ASSERT_EQ(dualsol[j], dualsol2[j]);
//   }

//   ASSERT_EQ(activity.size(), activity2.size());
//   for (size_t j = 0; j < activity.size(); j++) {
//     ASSERT_EQ(activity[j], activity2[j]);
//   }

//   ASSERT_EQ(redcost.size(), redcost2.size());
//   for (size_t j = 0; j < redcost.size(); j++) {
//     ASSERT_EQ(redcost[j], redcost2[j]);
//   }

//   if (DEF_INTERFACE == 0) {
//     remove("lpi_change_test_problem.lp.gz");
//     SUCCEED();
//   } else {
//     ASSERT_OK(lp_interface_->WriteLP("lpi_change_test_problem2.lp"));
//     ASSERT_OK(lp_interface_->Clear());

//     std::ifstream t1("lpi_change_test_problem.lp");
//     std::stringstream buffer1;
//     buffer1 << t1.rdbuf();

//     std::ifstream t2("lpi_change_test_problem2.lp");
//     std::stringstream buffer2;
//     buffer2 << t2.rdbuf();
//     ASSERT_EQ(buffer1.str(), buffer2.str());

//     remove("lpi_change_test_problem.lp");
//     remove("lpi_change_test_problem2.lp");
//   }
//  }
}  // namespace minimip
