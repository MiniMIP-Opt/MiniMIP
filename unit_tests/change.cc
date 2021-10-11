// @file   change.c
// @brief  unit tests for testing the methods that change and get coefficients
// in the .
//
// The tested methods are:
// GetCoefficient,
// ChangeObjective, GetObjective,
// ChangeBounds, GetBounds,
// ChangeSides, GetSides,
// ChangeObjectiveSense, GetObjectiveSense,
// GetNumberOfColumns, GetNumberOfRows, GetNumberOfNonZeros,
// GetColumns, GetRows,
// ClearState,
// WriteLP, ReadLP, Clear

#include <gtest/gtest.h>

#include "absl/status/status.h"
#include "src/lp_interface/lpi_factory.h"

#define TEST_ERRORS 0  // if 0 then skip tests expected to fail.
#define DEF_INTERFACE \
  0  // 0 = Glop Interface (Default),
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

  // write ncols, nrows and objsen into variables to check later
  static void initProb(int pos, int& ncols, int& nrows, int& nnonz,
                       LPObjectiveSense& objsen) {
    std::vector<double> obj = {1.0, 1.0};
    std::vector<double> lb  = {0.0, 0.0};
    std::vector<double> ub  = {lp_interface_->Infinity(),
                              lp_interface_->Infinity()};
    std::vector<double> lhs = {-lp_interface_->Infinity(),
                               -lp_interface_->Infinity()};
    std::vector<double> rhs = {1.0, 1.0};
    std::vector<double> val = {1.0, 1.0};
    std::vector<int> beg    = {0, 1};
    std::vector<int> ind    = {0, 1};

    // maximization problems, ncols is 1, nrows is 1
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
        ncols  = 1;
        nrows  = 1;
        nnonz  = 1;
        objsen = LPObjectiveSense::kMaximization;
        val[0] = -1.0;
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
        ncols  = 1;
        nrows  = 1;
        nnonz  = 1;
        objsen = LPObjectiveSense::kMaximization;
        rhs[0] = 0.0;
        break;

      case 2:
        // minimization problems (duals of the above)
        ncols  = 1;
        nrows  = 1;
        nnonz  = 1;
        objsen = LPObjectiveSense::kMinimization;
        rhs[0] = lp_interface_->Infinity();
        lhs[0] = 1;
        val[0] = -1.0;
        break;

      case 3:
        ncols  = 1;
        nrows  = 1;
        nnonz  = 1;
        objsen = LPObjectiveSense::kMinimization;
        rhs[0] = lp_interface_->Infinity();
        lhs[0] = 1;
        obj[0] = 0.0;
        break;

      case 4:
        // maximization problems, ncols is 2, *nrows is 2
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
        ncols  = 2;
        nrows  = 2;
        nnonz  = 2;
        objsen = LPObjectiveSense::kMaximization;
        val[0] = -1.0;
        val[1] = -1.0;
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
        ncols  = 2;
        nrows  = 2;
        nnonz  = 2;
        objsen = LPObjectiveSense::kMaximization;
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
        ncols  = 2;
        nrows  = 2;
        nnonz  = 2;
        objsen = LPObjectiveSense::kMaximization;
        rhs[0] = -1.0;
        rhs[1] = -1.0;
        val[0] = -1.0;
        break;

      case 7:
        // minimization problems (duals of the above)
        ncols  = 2;
        nrows  = 2;
        nnonz  = 2;
        objsen = LPObjectiveSense::kMinimization;
        rhs[0] = lp_interface_->Infinity();
        rhs[1] = lp_interface_->Infinity();
        lhs[0] = 1.0;
        lhs[1] = 1.0;
        val[0] = -1.0;
        val[1] = -1.0;
        break;

      case 8:
        ncols  = 2;
        nrows  = 2;
        nnonz  = 2;
        objsen = LPObjectiveSense::kMinimization;
        rhs[0] = lp_interface_->Infinity();
        rhs[1] = lp_interface_->Infinity();
        lhs[0] = 1.0;
        lhs[1] = 1.0;
        break;

      case 9:
        ncols  = 2;
        nrows  = 2;
        nnonz  = 2;
        objsen = LPObjectiveSense::kMinimization;
        rhs[0] = lp_interface_->Infinity();
        rhs[1] = lp_interface_->Infinity();
        lhs[0] = 1.0;
        lhs[1] = 1.0;
        obj[0] = -1.0;
        obj[1] = -1.0;
        val[0] = -1.0;
        break;

      default:
        return;
    }
    // empty placeholders
    std::vector<std::string> empty_names;
    std::vector<double> empty_vals;
    std::vector<int> empty_indices;

    ASSERT_EQ(lp_interface_->SetObjectiveSense(objsen), absl::OkStatus());
    ASSERT_EQ(
        lp_interface_->AddColumns(ncols, obj, lb, ub, empty_names, 0,
                                  empty_indices, empty_indices, empty_vals),
        absl::OkStatus());
    ASSERT_EQ(lp_interface_->AddRows(nrows, lhs, rhs, empty_names, nnonz, beg,
                                     ind, val),
              absl::OkStatus());
    ASSERT_TRUE(!lp_interface_->IsSolved());
    ASSERT_EQ(lp_interface_->SolvePrimal(), absl::OkStatus());
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
  static void checkChgObj(int& lastcol, std::vector<int>& ind,
                          std::vector<double>& setobj) {
    std::vector<double> obj(2);
    ASSERT_LE(lastcol, 2);
    ASSERT_EQ(lp_interface_->SetObjectiveCoefficients(lastcol, ind, setobj),
              absl::OkStatus());
    ASSERT_TRUE(!lp_interface_->IsSolved());
    ASSERT_EQ(lp_interface_->GetObjective(0, lastcol - 1, obj),
              absl::OkStatus());

    for (size_t i; i < setobj.size(); i++) {
      ASSERT_EQ(obj[i], setobj[i]);
    }
  }
};

TEST_P(ChangeObjective, checkChgObj) {
  int prob = std::get<2>(GetParam());
  int nrows, ncols, nnonz;
  LPObjectiveSense sense;
  std::vector<int> ind = {0, 1};
  std::vector<double> setobj(2);
  bool deathflag = false;

  setobj[0] = substituteInfinity(std::get<0>(GetParam()));
  setobj[1] = substituteInfinity(std::get<1>(GetParam()));

  ASSERT_NO_FATAL_FAILURE(initProb(prob, ncols, nrows, nnonz, sense));

  if (DEF_INTERFACE == 0) {
    if (TEST_ERRORS) {
      for (int i = 0; i < ncols; i++) {
        if (lp_interface_->IsInfinity(fabs(setobj[i]))) deathflag = true;
      }
    } else {
      for (int i = 0; i < ncols; i++) {
        if (lp_interface_->IsInfinity(fabs(setobj[i]))) GTEST_SKIP();
      }
    }
    if (deathflag)
      ASSERT_DEATH(checkChgObj(ncols, ind, setobj), "");
    else
      ASSERT_NO_FATAL_FAILURE(checkChgObj(ncols, ind, setobj));
  } else {
    ASSERT_NO_FATAL_FAILURE(checkChgObj(ncols, ind, setobj));
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
  static void checkChgBounds(int& lastcol, std::vector<int>& ind,
                             std::vector<double>& setlb,
                             std::vector<double>& setub,
                             absl::Status error = absl::OkStatus()) {
    std::vector<double> ub(2);
    std::vector<double> lb(2);

    ASSERT_LE(lastcol, 2);
    if (error != absl::OkStatus()) {
      ASSERT_EQ(lp_interface_->SetColumnBounds(lastcol, ind, setlb, setub), error);
      abort();
    } else {
      ASSERT_EQ(lp_interface_->SetColumnBounds(lastcol, ind, setlb, setub),
                absl::OkStatus());
    }
    ASSERT_TRUE(!lp_interface_->IsSolved());
    ASSERT_EQ(lp_interface_->GetBounds(0, lastcol - 1, lb, ub),
              absl::OkStatus());

    for (size_t i; i < setub.size(); i++) {
      ASSERT_EQ(ub[i], setub[i]);
    }
    for (size_t i; i < setlb.size(); i++) {
      ASSERT_EQ(lb[i], setlb[i]);
    }
  }
};

TEST_P(ChangeBounds, checkChgBounds) {
  int prob = std::get<4>(GetParam());

  int nrows, ncols, nnonz;
  std::vector<int> ind = {0, 1};
  LPObjectiveSense sense;
  std::vector<double> setub(2);
  std::vector<double> setlb(2);
  bool deathflag = false;

  setlb[0] = substituteInfinity(std::get<0>(GetParam()));
  setlb[1] = substituteInfinity(std::get<1>(GetParam()));
  setub[0] = substituteInfinity(std::get<2>(GetParam()));
  setub[1] = substituteInfinity(std::get<3>(GetParam()));

  ASSERT_NO_FATAL_FAILURE(initProb(prob, ncols, nrows, nnonz, sense));

  if (TEST_ERRORS) {
    for (int i = 0; i < ncols; i++) {
      if (setub[i] < setlb[i]) deathflag = true;
      if (lp_interface_->IsInfinity(setlb[i]) ||
          lp_interface_->IsInfinity(-setub[i]))
        deathflag = true;
    }
  } else {
    for (int i = 0; i < ncols; i++) {
      if (setub[i] < setlb[i]) GTEST_SKIP();
      if (lp_interface_->IsInfinity(setlb[i]) ||
          lp_interface_->IsInfinity(-setub[i]))
        GTEST_SKIP();
    }
  }

  if (deathflag)
    ASSERT_DEATH(
        checkChgBounds(ncols, ind, setlb, setub,
                       absl::Status(absl::StatusCode::kInternal, "LP Error")),
        "");
  else
    ASSERT_NO_FATAL_FAILURE(checkChgBounds(ncols, ind, setlb, setub));
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
  static void checkChgSides(int& lastcol, std::vector<int>& ind,
                            std::vector<double>& setls,
                            std::vector<double>& setrs) {
    std::vector<double> ls(2);
    std::vector<double> rs(2);

    ASSERT_LE(lastcol, 2);
    ASSERT_EQ(lp_interface_->SetRowSides(lastcol, ind, setls, setrs),
              absl::OkStatus());
    ASSERT_TRUE(!lp_interface_->IsSolved());
    ASSERT_EQ(lp_interface_->GetSides(0, lastcol - 1, ls, rs),
              absl::OkStatus());

    for (int i; i < setls.size(); i++) {
      ASSERT_EQ(ls[i], setls[i]);
    }
    for (int i; i < setrs.size(); i++) {
      ASSERT_EQ(rs[i], setrs[i]);
    }
  }
};

TEST_P(ChangeSides, checkChgSides) {
  double left1  = substituteInfinity(std::get<0>(GetParam()));
  double left2  = substituteInfinity(std::get<1>(GetParam()));
  double right1 = substituteInfinity(std::get<2>(GetParam()));
  double right2 = substituteInfinity(std::get<3>(GetParam()));
  int prob      = std::get<4>(GetParam());

  int nrows, ncols, nnonz;
  std::vector<int> ind = {0, 1};
  LPObjectiveSense sense;
  std::vector<double> setlhs = {left1, left2};
  std::vector<double> setrhs = {right1, right2};
  bool deathflag             = false;

  ASSERT_NO_FATAL_FAILURE(initProb(prob, ncols, nrows, nnonz, sense));

  if (TEST_ERRORS) {
    for (int i = 0; i < nrows; i++) {
      if (setrhs[i] < setlhs[i]) deathflag = true;
      if (lp_interface_->IsInfinity(fabs(setlhs[i])) &&
          lp_interface_->IsInfinity(fabs(setrhs[i])))
        deathflag = true;
    }
  } else {
    for (int i = 0; i < nrows; i++) {
      if (setrhs[i] < setlhs[i]) GTEST_SKIP();
      if (lp_interface_->IsInfinity(fabs(setlhs[i])) &&
          lp_interface_->IsInfinity(fabs(setrhs[i])))
        GTEST_SKIP();
    }
  }

  if (deathflag)
    ASSERT_DEATH(checkChgSides(nrows, ind, setlhs, setrhs), "");
  else
    ASSERT_NO_FATAL_FAILURE(checkChgSides(nrows, ind, setlhs, setrhs));
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

TEST_P(ChangeObjectiveSense, ChgObjSen) {
  LPObjectiveSense newsense = std::get<0>(GetParam());
  int prob                  = std::get<1>(GetParam());

  int nrows, ncols, nnonz;
  LPObjectiveSense sense;
  LPObjectiveSense probsense;

  ASSERT_NO_FATAL_FAILURE(initProb(prob, ncols, nrows, nnonz, sense));

  ASSERT_EQ(lp_interface_->SetObjectiveSense(newsense), absl::OkStatus());
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
TEST_F(Change, testrowmethods) {
  // problem data
  std::vector<double> obj = {1.0, 1.0, 1.0, 1.0, 1.0};
  std::vector<double> lb  = {-1.0, -lp_interface_->Infinity(), 0.0,
                            -lp_interface_->Infinity(), 0.0};
  std::vector<double> ub  = {10.0, lp_interface_->Infinity(),
                            lp_interface_->Infinity(), 29.0, 0.0};
  int ncolsbefore, ncolsafter;
  int nrowsbefore, nrowsafter;
  std::vector<double> lhsvals = {
      -lp_interface_->Infinity(), -1.0, -3e-10, 0.0, 1.0, 3e10};
  std::vector<double> rhsvals = {-1.0, -3e-10, 0.0,
                                 1.0,  3e10,   lp_interface_->Infinity()};
  std::vector<int> nnonzs     = {1, 10, -1, 6, -1};
  std::vector<int> begvals    = {0, 2, 3, 5, 8, 9};
  std::vector<int> indvals    = {0, 1, 3, 2, 1, 1, 2, 4, 0, 3};
  std::vector<double> vals    = {1.0, 5.0, -1.0, 3e5,  2.0,
                              1.0, 20,  10,   -1.9, 1e-2};

  int iterations              = 5;
  std::vector<int> k          = {1, 6, -1, 4, -2};
  std::vector<int> nnonzsdiff = {1, 10, -1, 6, -3};
  int i;
  int j;

  // empty placeholders
  std::vector<std::string> empty_names;
  std::vector<double> empty_vals;
  std::vector<int> empty_indices;

  // create original lp
  ASSERT_EQ(lp_interface_->AddColumns(5, obj, lb, ub, empty_names, 0,
                                      empty_indices, empty_indices, empty_vals),
            absl::OkStatus());
  ncolsbefore = lp_interface_->GetNumberOfColumns();

  for (i = 0; i < iterations; i++) {
    // setup row values
    int nrows;
    int nnonzsbefore;
    int nnonzsafter;

    nrows = k[i];
    // get data before modification
    nnonzsbefore = lp_interface_->GetNumberOfNonZeros();
    nrowsbefore  = lp_interface_->GetNumberOfRows();

    if (nrows < 0) {
      ASSERT_EQ(lp_interface_->DeleteRows(0, -(1 + nrows)), absl::OkStatus());
    } else {  // nrows >= 0
      std::vector<double> lhs(100);
      std::vector<double> rhs(100);
      std::vector<int> beg(100);

      int nnonz = nnonzs[i];
      std::vector<int> ind(100);
      std::vector<double> val(100);

      std::vector<double> newlhs(100);
      std::vector<double> newval(100);
      std::vector<double> newrhs(100);
      std::vector<int> newbeg(100);
      std::vector<int> newind(100);
      int newnnonz;
      int indold;
      int indnew;

      ASSERT_LE(nrows, 100);
      for (j = 0; j < nrows; j++) {
        lhs[j] = lhsvals[j];
        rhs[j] = rhsvals[j];
        beg[j] = begvals[j];
      }

      ASSERT_LE(nnonz, 100);
      for (j = 0; j < nnonz; j++) {
        ind[j] = indvals[j];
        val[j] = vals[j];
      }
      ASSERT_EQ(lp_interface_->AddRows(nrows, lhs, rhs, empty_names, nnonz, beg,
                                       ind, val),
                absl::OkStatus());

      // checks
      ASSERT_EQ(
          lp_interface_->GetRows(nrowsbefore, nrowsbefore - 1 + nrows, newlhs,
                                 newrhs, newnnonz, newbeg, newind, newval),
          absl::OkStatus());
      ASSERT_EQ(nnonz, newnnonz);

      for (j = 0; j < nrows; j++) {
        ASSERT_EQ(beg[j], newbeg[j]);
      }
      beg[nrows]    = nnonz;
      newbeg[nrows] = newnnonz;

      // check each row seperately
      for (j = 0; j < nrows; j++) {
        if (fabs(lhs[j]) < 1e30 && fabs(newlhs[j]) < 1e30) {
          ASSERT_DOUBLE_EQ(lhs[j], newlhs[j]);
        }
        if (fabs(rhs[j]) < 1e30 && fabs(newrhs[j]) < 1e30) {
          ASSERT_DOUBLE_EQ(rhs[j], newrhs[j]);
        }

        // We add a row where the indices are not sorted, some lp solvers give
        // them back sorted (e.g. soplex), some others don't (e.g. cplex).
        // Therefore we cannot simply assert the ind and val arrays to be equal,
        // but have to search for and check each value individually.
        for (indold = beg[j]; indold < beg[j + 1]; indold++) {
          int occurrences = 0;

          // for each value ind associated to the current row search for it in
          // the newind array
          for (indnew = beg[j]; indnew < beg[j + 1]; indnew++) {
            if (ind[indold] == newind[indnew]) {
              occurrences = occurrences + 1;
              ASSERT_DOUBLE_EQ(val[indold], newval[indnew]);
            }
          }
          // assert that we found only one occurrence in the current row
          ASSERT_EQ(occurrences, 1);
        }
      }
    }

    // checks
    nrowsafter = lp_interface_->GetNumberOfRows();
    ASSERT_EQ(nrowsbefore + nrows, nrowsafter);

    nnonzsafter = lp_interface_->GetNumberOfNonZeros();
    ASSERT_EQ(nnonzsbefore + nnonzsdiff[i], nnonzsafter);

    ncolsafter = lp_interface_->GetNumberOfColumns();
    ASSERT_EQ(ncolsbefore, ncolsafter);
  }

  // delete rowsets
  // should have 8 rows now
  nrowsbefore = lp_interface_->GetNumberOfRows();
  ASSERT_EQ(8, nrowsbefore);
  for (i = 3; i > 0; i--) {
    std::vector<bool> rows = {false, false, false, false,
                              false, false, false, false};

    for (j = 0; j < i; j++) rows[(2 * j) + 1] = true;

    nrowsbefore = lp_interface_->GetNumberOfRows();
    ASSERT_EQ(lp_interface_->DeleteRowSet(rows), absl::OkStatus());

    nrowsafter = lp_interface_->GetNumberOfRows();
    ASSERT_EQ(nrowsbefore - i, nrowsafter);
    // assert that the rows that are left are the ones I intended
  }
}

// Test for AddColumns, DeleteColumns
TEST_F(Change, testcolmethods) {
  // problem data
  std::vector<double> obj = {1.0, 1.0, 1.0, 1.0, 1.0};
  std::vector<double> lhs = {-1.0, -lp_interface_->Infinity(), 0.0,
                             -lp_interface_->Infinity(), 0.0};
  std::vector<double> rhs = {10.0, lp_interface_->Infinity(),
                             lp_interface_->Infinity(), 29.0, 0.0};
  int ncolsbefore, ncolsafter;
  int nrowsbefore, nrowsafter;
  std::vector<double> lbvals = {
      -lp_interface_->Infinity(), -1.0, -3e-10, 0.0, 1.0, 3e10};
  std::vector<double> ubvals = {-1.0, -3e-10, 0.0,
                                1.0,  3e10,   lp_interface_->Infinity()};
  std::vector<double> vals   = {1.0, 5.0, -1.0, 3e5,  2.0,
                              1.0, 20,  10,   -1.9, 1e-2};
  std::vector<int> nnonzs    = {1, 10, -1, 6, -1};
  std::vector<int> begvals   = {0, 2, 3, 5, 8, 9};
  std::vector<int> indvals   = {0, 1, 3, 2, 1, 1, 2, 4, 0, 3};

  int iterations              = 5;
  std::vector<int> k          = {1, 6, -1, 4, -2};
  std::vector<int> nnonzsdiff = {1, 10, -1, 6, -3};

  // empty placeholders
  std::vector<std::string> empty_names;
  std::vector<double> empty_vals;
  std::vector<int> empty_indices;

  // create original lp
  ASSERT_EQ(lp_interface_->AddRows(5, lhs, rhs, empty_names, 0, empty_indices,
                                   empty_indices, empty_vals),
            absl::OkStatus());
  nrowsbefore = lp_interface_->GetNumberOfRows();

  for (int i = 0; i < iterations; i++) {
    // setup col values
    int ncols;
    int nnonzsbefore;
    int nnonzsafter;

    ncols = k[i];

    // get data before modification
    nnonzsbefore = lp_interface_->GetNumberOfNonZeros();
    ncolsbefore  = lp_interface_->GetNumberOfColumns();

    if (k[i] < 0) {
      ASSERT_EQ(lp_interface_->DeleteColumns(0, -(1 + ncols)),
                absl::OkStatus());
    } else {  // ncols >= 0
      std::vector<double> lb(100);
      std::vector<double> ub(100);
      std::vector<int> beg(100);

      assert(nnonzs[i] >= 0);
      int nnonz = nnonzs[i];
      std::vector<int> ind(100);
      std::vector<double> val(100);

      std::vector<double> newlb(100);
      std::vector<double> newval(100);
      std::vector<double> newub(100);
      std::vector<int> newbeg(100);
      std::vector<int> newind(100);
      int newnnonz;

      ASSERT_LE(ncols, 100);
      for (int j = 0; j < ncols; j++) {
        lb[j]  = lbvals[j];
        ub[j]  = ubvals[j];
        beg[j] = begvals[j];
      }

      ASSERT_LE(nnonz, 100);
      for (int j = 0; j < nnonz; j++) {
        ind[j] = indvals[j];
        val[j] = vals[j];
      }
      ASSERT_EQ(lp_interface_->AddColumns(ncols, obj, lb, ub, empty_names,
                                          nnonz, beg, ind, val),
                absl::OkStatus());

      // checks
      ASSERT_EQ(
          lp_interface_->GetColumns(ncolsbefore, ncolsbefore - 1 + ncols, newlb,
                                    newub, newnnonz, newbeg, newind, newval),
          absl::OkStatus());
      ASSERT_EQ(nnonz, newnnonz);

      for (int j = 0; j < ncols; j++) {
        ASSERT_EQ(lb[j], newlb[j]);
        ASSERT_EQ(ub[j], newub[j]);
        ASSERT_EQ(beg[j], newbeg[j]);
      }
      for (int j = 0; j < nnonz; j++) {
        ASSERT_EQ(ind[j], newind[j]);
        ASSERT_EQ(val[j], newval[j]);
      }
    }

    // checks
    nrowsafter = lp_interface_->GetNumberOfRows();
    ASSERT_EQ(nrowsbefore, nrowsafter);

    nnonzsafter = lp_interface_->GetNumberOfNonZeros();
    ASSERT_EQ(nnonzsbefore + nnonzsdiff[i], nnonzsafter);

    ncolsafter = lp_interface_->GetNumberOfColumns();
    ASSERT_EQ(ncolsbefore + ncols, ncolsafter);
  }

}

// Test adding zero coeffs cols
TEST_F(Change, testzerosincols) {
  int ncols;
  int nrows;
  int nnonz = 2;
  LPObjectiveSense sense;
  std::vector<double> lb  = {0};
  std::vector<double> ub  = {20};
  std::vector<int> beg    = {0};
  std::vector<int> ind    = {0, 1};
  std::vector<double> val = {0, 3};
  std::vector<double> obj = {1};

  // empty placeholders
  std::vector<std::string> empty_names;

  // 2x2 problem
  initProb(4, ncols, nrows, nnonz, sense);
  ASSERT_EQ(2, nrows);
  ASSERT_EQ(2, ncols);

#ifndef NDEBUG
  ASSERT_DEATH(lp_interface_->AddColumns(1, obj, lb, ub, empty_names, nnonz,
                                         beg, ind, val),
               "");
#endif
  // this test can only work in debug mode, so we make it pass in opt mode
#ifdef NDEBUG
  ASSERT_EQ(lp_interface_->AddColumns(1, obj, lb, ub, empty_names, nnonz, beg,
                                      ind, val),
            absl::OkStatus());
  SUCCEED();  // return SIGABORT
#endif
}

// Test adding zero coeffs in rows, expecting an assert in debug mode
//
// This test should fail with an assert from the which causes SIGABRT to be
// issued. Thus, this test should pass.
TEST_F(Change, testzerosinrows) {
  int nrows;
  int ncols;
  int nnonz = 2;
  LPObjectiveSense sense;
  std::vector<double> lhs = {0};
  std::vector<double> rhs = {20};
  std::vector<int> beg    = {0};
  std::vector<int> ind    = {0, 1};
  std::vector<double> val = {0, 3};
  // empty placeholders
  std::vector<std::string> empty_names;

  // 2x2 problem
  initProb(4, ncols, nrows, nnonz, sense);
  ASSERT_EQ(2, nrows);
  ASSERT_EQ(2, ncols);

#ifndef NDEBUG
  ASSERT_DEATH(
      lp_interface_->AddRows(1, lhs, rhs, empty_names, nnonz, beg, ind, val),
      "");
#else
  // this test can only work in debug mode, so we make it pass in opt mode
  ASSERT_EQ(
      lp_interface_->AddRows(1, lhs, rhs, empty_names, nnonz, beg, ind, val),
      absl::OkStatus());
  SUCCEED();
#endif
}

// test WriteLP, ReadLP, Clear
TEST_F(Change, testlpiwritereadlpmethods) {
  int nrows, ncols, nnonz;
  double objval;
  std::vector<double> primsol(2);
  std::vector<double> dualsol(2);
  std::vector<double> activity(2);
  std::vector<double> redcost(2);
  double objval2;
  std::vector<double> primsol2(2);
  std::vector<double> dualsol2(2);
  std::vector<double> activity2(2);
  std::vector<double> redcost2(2);
  LPObjectiveSense sense;

  // 2x2 problem
  ASSERT_NO_FATAL_FAILURE(initProb(5, ncols, nrows, nnonz, sense));

  ASSERT_EQ(lp_interface_->SolvePrimal(), absl::OkStatus());
  ASSERT_EQ(
      lp_interface_->GetSolution(objval, primsol, dualsol, activity, redcost),
      absl::OkStatus());

  ASSERT_EQ(lp_interface_->WriteLP("lpi_change_test_problem.lp"),
            absl::OkStatus());
  ASSERT_EQ(lp_interface_->Clear(), absl::OkStatus());

  if (DEF_INTERFACE == 0)
    ASSERT_EQ(lp_interface_->ReadLP("lpi_change_test_problem.lp.gz"),
              absl::OkStatus());
  else
    ASSERT_EQ(lp_interface_->ReadLP("lpi_change_test_problem.lp"),
              absl::OkStatus());

  ASSERT_EQ(lp_interface_->SolvePrimal(), absl::OkStatus());
  ASSERT_EQ(lp_interface_->GetSolution(objval2, primsol2, dualsol2, activity2,
                                       redcost2),
            absl::OkStatus());
  ASSERT_FLOAT_EQ(objval, objval2);

  ASSERT_EQ(primsol.size(), primsol2.size());
  for (size_t j = 0; j < primsol.size(); j++) {
    ASSERT_EQ(primsol[j], primsol2[j]);
  }

  ASSERT_EQ(dualsol.size(), dualsol2.size());
  for (size_t j = 0; j < dualsol.size(); j++) {
    ASSERT_EQ(dualsol[j], dualsol2[j]);
  }

  ASSERT_EQ(activity.size(), activity2.size());
  for (size_t j = 0; j < activity.size(); j++) {
    ASSERT_EQ(activity[j], activity2[j]);
  }

  ASSERT_EQ(redcost.size(), redcost2.size());
  for (size_t j = 0; j < redcost.size(); j++) {
    ASSERT_EQ(redcost[j], redcost2[j]);
  }

  if (DEF_INTERFACE == 0) {
    remove("lpi_change_test_problem.lp.gz");
    SUCCEED();
  } else {
    ASSERT_EQ(lp_interface_->WriteLP("lpi_change_test_problem2.lp"),
              absl::OkStatus());
    ASSERT_EQ(lp_interface_->Clear(), absl::OkStatus());

    std::ifstream t1("lpi_change_test_problem.lp");
    std::stringstream buffer1;
    buffer1 << t1.rdbuf();

    std::ifstream t2("lpi_change_test_problem2.lp");
    std::stringstream buffer2;
    buffer2 << t2.rdbuf();
    ASSERT_EQ(buffer1.str(), buffer2.str());

    remove("lpi_change_test_problem.lp");
    remove("lpi_change_test_problem2.lp");
  }
}
}  // namespace minimip
