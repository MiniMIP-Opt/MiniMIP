/**@file   change.c
 * @brief  unit tests for testing the methods that change and get coefficients in the .
 *
 * The tested methods are:
 * GetCoefficient,
 * ChangeObjective, GetObjective,
 * ChangeBounds, GetBounds,
 * ChangeSides, GetSides,
 * ChangeObjectiveSense, GetObjectiveSense,
 * GetNumberOfColumns, GetNumberOfRows, GetNumberOfNonZeros,
 * GetColumns, GetRows,
 * ClearState,
 * WriteLP, ReadLP, Clear
 */

#include "src/lp_interface/lpi_factory.h"

#include <gtest/gtest.h>

#define TEST_ERRORS 0   /* if 0 then skip tests expected to fail. */
#define DEF_INTERFACE 1 /* 0 = Glop Interface (Default), \
                          * 1 = SoPlex Interface, */

namespace minimip {

/*** TEST SUITE SOLVE ***/

static LPInterface* lp_interface_ = nullptr;

class Change : public ::testing::Test {
 protected:
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
  }

  /** write ncols, nrows and objsen into variables to check later */
  static void initProb(LPNum pos, LPNum& ncols, LPNum& nrows, LPNum& nnonz, LPObjectiveSense& objsen) {

    LPValueArray obj = {1.0, 1.0};
    LPValueArray lb = {0.0, 0.0};
    LPValueArray ub = {lp_interface_->Infinity(), lp_interface_->Infinity()};
    LPValueArray lhs = {-lp_interface_->Infinity(), -lp_interface_->Infinity()};
    LPValueArray rhs = {1.0, 1.0};
    LPValueArray val = {1.0, 1.0};
    LPIndexArray beg = {0, 1};
    LPIndexArray ind = {0, 1};

    /* maximization problems, ncols is 1, nrows is 1*/
    switch (pos) {
      case 0:
        /* unbounded - infeasible
         * (P):  max x
         * -x <= 1 (constr)
         *  0 <= x (bound)
         *
         * (D):  min y
         * 1 <= -y (constr)
         * 0 <= y (bound)
         */
        ncols = 1;
        nrows = 1;
        nnonz = 1;
        objsen = LPObjectiveSense::kMaximize;
        val[0] = -1.0;
        break;

      case 1:
        /* optimal - optimal
         * (P):  max x
         *  x <= 0 (constr)
         *  0 <= x (bound)
         *
         * (D):  min 0
         * 1 <= y (constr)
         * 0 <= y (bound)
         */
        ncols = 1;
        nrows = 1;
        nnonz = 1;
        objsen = LPObjectiveSense::kMaximize;
        rhs[0] = 0.0;
        break;

      case 2:
        /* minimization problems (duals of the above) */
        ncols = 1;
        nrows = 1;
        nnonz = 1;
        objsen = LPObjectiveSense::kMinimize;
        rhs[0] = lp_interface_->Infinity();
        lhs[0] = 1;
        val[0] = -1.0;
        break;

      case 3:
        ncols = 1;
        nrows = 1;
        nnonz = 1;
        objsen = LPObjectiveSense::kMinimize;
        rhs[0] = lp_interface_->Infinity();
        lhs[0] = 1;
        obj[0] = 0.0;
        break;

      case 4:
        /* maximization problems, ncols is 2, *nrows is 2 */
        /* unbounded - infeasible
         * (P):  max x+y
         * -x    <= 1 (constr)
         *    -y <= 1 (constr)
         *
         *  0 <= x (bound)
         *  0 <= y (bound)
         *
         * (D):  min x+y
         * 1 <= -x   (constr)
         * 1 <=   -y (constr)
         *
         * 0 <= x (bound)
         * 0 <= y (bound)
         */
        ncols = 2;
        nrows = 2;
        nnonz = 2;
        objsen = LPObjectiveSense::kMaximize;
        val[0] = -1.0;
        val[1] = -1.0;
        break;

      case 5:
        /* optimal - optimal
         * (P):  max x+y
         * x     <= 1 (constr)
         *     y <= 1 (constr)
         *
         *  0 <= x (bound)
         *  0 <= y (bound)
         *
         * (D):  min x+y
         * 1 <= x    (constr)
         * 1 <=    y (constr)
         *
         * 0 <= x (bound)
         * 0 <= y (bound)
         */
        ncols = 2;
        nrows = 2;
        nnonz = 2;
        objsen = LPObjectiveSense::kMaximize;
        break;

      case 6:
        /* infeasible - infeasible
         * (P):  max x+y
         * -x    <= -1 (constr)
         *     y <= -1 (constr)
         *
         *  0 <= x (bound)
         *  0 <= y (bound)
         *
         * (D):  min -x-y
         * 1 <= -x    (constr)
         * 1 <=     y (constr)
         *
         * 0 <= x (bound)
         * 0 <= y (bound)
         */
        ncols = 2;
        nrows = 2;
        nnonz = 2;
        objsen = LPObjectiveSense::kMaximize;
        rhs[0] = -1.0;
        rhs[1] = -1.0;
        val[0] = -1.0;
        break;

      case 7:
        /* minimization problems (duals of the above) */
        ncols = 2;
        nrows = 2;
        nnonz = 2;
        objsen = LPObjectiveSense::kMinimize;
        rhs[0] = lp_interface_->Infinity();
        rhs[1] = lp_interface_->Infinity();
        lhs[0] = 1.0;
        lhs[1] = 1.0;
        val[0] = -1.0;
        val[1] = -1.0;
        break;

      case 8:
        ncols = 2;
        nrows = 2;
        nnonz = 2;
        objsen = LPObjectiveSense::kMinimize;
        rhs[0] = lp_interface_->Infinity();
        rhs[1] = lp_interface_->Infinity();
        lhs[0] = 1.0;
        lhs[1] = 1.0;
        break;

      case 9:
        ncols = 2;
        nrows = 2;
        nnonz = 2;
        objsen = LPObjectiveSense::kMinimize;
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
    /* empty placeholders */
    StringArray empty_names;
    LPValueArray empty_vals;
    LPIndexArray empty_indices;

    ASSERT_EQ(lp_interface_->ChangeObjectiveSense(objsen), RetCode::kOkay);
    ASSERT_EQ(lp_interface_->AddColumns(ncols, obj, lb, ub, empty_names, 0, empty_indices, empty_indices, empty_vals), RetCode::kOkay);
    ASSERT_EQ(lp_interface_->AddRows(nrows, lhs, rhs, empty_names, nnonz, beg, ind, val), RetCode::kOkay);
    ASSERT_TRUE(!lp_interface_->WasSolved());
    ASSERT_EQ(lp_interface_->SolvePrimal(), RetCode::kOkay);
    ASSERT_TRUE(lp_interface_->WasSolved());
  }

   /* We treat -2 and 2 as -infinity and infinity resp. */
  static LPValue substituteInfinity(LPValue inf) {
    if (inf == 2)
      return lp_interface_->Infinity();
    else if (inf == -2)
      return -lp_interface_->Infinity();
    else
      return inf;
  }
};

/* TESTS **/

/** Test ChangeObjectives */
class ChangeObjective : public Change,
                        public ::testing::WithParamInterface<std::tuple<LPValue, LPValue, LPIndex>> {
 protected:
  static void checkChgObj(LPIndex& lastcol, LPIndexArray& ind, LPValueArray& setobj) {
    LPValueArray obj(2);
    ASSERT_LE(lastcol, 2);
    ASSERT_EQ(lp_interface_->ChangeObjective(lastcol, ind, setobj), RetCode::kOkay);
    ASSERT_TRUE(!lp_interface_->WasSolved());
    ASSERT_EQ(lp_interface_->GetObjective(0, lastcol - 1, obj), RetCode::kOkay);

    for (size_t i; i < setobj.size(); i++) {
      ASSERT_EQ(obj[i], setobj[i]);
    }
  }
};

TEST_P(ChangeObjective, checkChgObj) {
  LPNum prob = std::get<2>(GetParam());
  LPNum nrows, ncols, nnonz;
  LPObjectiveSense sense;
  LPIndexArray ind = {0, 1};
  LPValueArray setobj(2);
  bool deathflag = false;

  setobj[0] = substituteInfinity(std::get<0>(GetParam()));
  setobj[1] = substituteInfinity(std::get<1>(GetParam()));

  ASSERT_NO_FATAL_FAILURE(initProb(prob, ncols, nrows, nnonz, sense));

  if (DEF_INTERFACE == 0) {
    if (TEST_ERRORS) {
      for (LPIndex i = 0; i < ncols; i++) {
        if (lp_interface_->IsInfinity(fabs(setobj[i])))
          deathflag = true;
      }
    } else {
      for (LPIndex i = 0; i < ncols; i++) {
        if (lp_interface_->IsInfinity(fabs(setobj[i])))
          GTEST_SKIP();
      }
    }
    if (deathflag)
      ASSERT_DEATH(checkChgObj(ncols, ind, setobj), "");
    else
      ASSERT_NO_FATAL_FAILURE(checkChgObj(ncols, ind, setobj));

  } else
    ASSERT_NO_FATAL_FAILURE(checkChgObj(ncols, ind, setobj));
}

INSTANTIATE_TEST_SUITE_P(
  ChangeObjectivesCombinations,
  ChangeObjective,
  ::testing::Combine(
    ::testing::Values(0, 1, -1, 2, -2),
    ::testing::Values(0, 1, -1, 2, -2),
    ::testing::Values(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)));

/** Test ChangeBounds */
class ChangeBounds : public Change,
                     public ::testing::WithParamInterface<std::tuple<LPValue, LPValue, LPValue, LPValue, LPIndex>> {
 protected:
  static void checkChgBounds(LPIndex& lastcol, LPIndexArray& ind, LPValueArray& setlb, LPValueArray& setub, RetCode error = RetCode::kOkay) {
    LPValueArray ub(2);
    LPValueArray lb(2);

    ASSERT_LE(lastcol, 2);
    if (error != RetCode::kOkay) {
      ASSERT_EQ(lp_interface_->ChangeBounds(lastcol, ind, setlb, setub), error);
      abort();
    } else
      ASSERT_EQ(lp_interface_->ChangeBounds(lastcol, ind, setlb, setub), RetCode::kOkay);
    ASSERT_TRUE(!lp_interface_->WasSolved());
    ASSERT_EQ(lp_interface_->GetBounds(0, lastcol - 1, lb, ub), RetCode::kOkay);

    for (size_t i; i < setub.size(); i++) {
      ASSERT_EQ(ub[i], setub[i]);
    }
    for (size_t i; i < setlb.size(); i++) {
      ASSERT_EQ(lb[i], setlb[i]);
    }
  }
};

TEST_P(ChangeBounds, checkChgBounds) {
  LPNum prob = std::get<4>(GetParam());

  LPNum nrows, ncols, nnonz;
  LPIndexArray ind = {0, 1};
  LPObjectiveSense sense;
  LPValueArray setub(2);
  LPValueArray setlb(2);
  bool deathflag = false;

  setlb[0] = substituteInfinity(std::get<0>(GetParam()));
  setlb[1] = substituteInfinity(std::get<1>(GetParam()));
  setub[0] = substituteInfinity(std::get<2>(GetParam()));
  setub[1] = substituteInfinity(std::get<3>(GetParam()));

  ASSERT_NO_FATAL_FAILURE(initProb(prob, ncols, nrows, nnonz, sense));

  if (TEST_ERRORS) {
    for (LPIndex i = 0; i < ncols; i++) {
      if (setub[i] < setlb[i])
        deathflag = true;
      if (lp_interface_->IsInfinity(setlb[i]) || lp_interface_->IsInfinity(-setub[i]))
        deathflag = true;
    }
  } else {
    for (LPIndex i = 0; i < ncols; i++) {
      if (setub[i] < setlb[i])
        GTEST_SKIP();
      if (lp_interface_->IsInfinity(setlb[i]) || lp_interface_->IsInfinity(-setub[i]))
        GTEST_SKIP();
    }
  }

  if (deathflag)
    ASSERT_DEATH(checkChgBounds(ncols, ind, setlb, setub, RetCode::kLPError), "");
  else
    ASSERT_NO_FATAL_FAILURE(checkChgBounds(ncols, ind, setlb, setub));
}

INSTANTIATE_TEST_SUITE_P(
  ChangeBoundsCombinations,
  ChangeBounds,
  ::testing::Combine(
    ::testing::Values(0, 1, -1, 2, -2),
    ::testing::Values(0, 1, -1, 2, -2),
    ::testing::Values(0, 1, -1, 2, -2),
    ::testing::Values(0, 1, -1, 2, -2),
    ::testing::Values(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)));

/** Test ChangeSides */
class ChangeSides : public Change,
                    public ::testing::WithParamInterface<std::tuple<LPValue, LPValue, LPValue, LPValue, LPIndex>> {
 protected:
  static void checkChgSides(LPIndex& lastcol, LPIndexArray& ind, LPValueArray& setls, LPValueArray& setrs) {
    LPValueArray ls(2);
    LPValueArray rs(2);

    ASSERT_LE(lastcol, 2);
    ASSERT_EQ(lp_interface_->ChangeSides(lastcol, ind, setls, setrs), RetCode::kOkay);
    ASSERT_TRUE(!lp_interface_->WasSolved());
    ASSERT_EQ(lp_interface_->GetSides(0, lastcol - 1, ls, rs), RetCode::kOkay);

    for (LPIndex i; i < setls.size(); i++) {
      ASSERT_EQ(ls[i], setls[i]);
    }
    for (LPIndex i; i < setrs.size(); i++) {
      ASSERT_EQ(rs[i], setrs[i]);
    }
  }
};

TEST_P(ChangeSides, checkChgSides) {
  LPValue left1 = substituteInfinity(std::get<0>(GetParam()));
  LPValue left2 = substituteInfinity(std::get<1>(GetParam()));
  LPValue right1 = substituteInfinity(std::get<2>(GetParam()));
  LPValue right2 = substituteInfinity(std::get<3>(GetParam()));
  LPNum prob = std::get<4>(GetParam());

  LPNum nrows, ncols, nnonz;
  LPIndexArray ind = {0, 1};
  LPObjectiveSense sense;
  LPValueArray setlhs = {left1, left2};
  LPValueArray setrhs = {right1, right2};
  bool deathflag = false;

  ASSERT_NO_FATAL_FAILURE(initProb(prob, ncols, nrows, nnonz, sense));

  if (TEST_ERRORS) {
    for (LPIndex i = 0; i < nrows; i++) {
      if (setrhs[i] < setlhs[i])
        deathflag = true;
      if (lp_interface_->IsInfinity(fabs(setlhs[i])) && lp_interface_->IsInfinity(fabs(setrhs[i])))
        deathflag = true;
    }
  } else {
    for (LPIndex i = 0; i < nrows; i++) {
      if (setrhs[i] < setlhs[i])
        GTEST_SKIP();
      if (lp_interface_->IsInfinity(fabs(setlhs[i])) && lp_interface_->IsInfinity(fabs(setrhs[i])))
        GTEST_SKIP();
    }
  }

  if (deathflag)
    ASSERT_DEATH(checkChgSides(nrows, ind, setlhs, setrhs), "");
  else
    ASSERT_NO_FATAL_FAILURE(checkChgSides(nrows, ind, setlhs, setrhs));
}

INSTANTIATE_TEST_SUITE_P(
  ChangeSidesCombinations,
  ChangeSides,
  ::testing::Combine(
    ::testing::Values(0, 1, -1, 2, -2),
    ::testing::Values(0, 1, -1, 2, -2),
    ::testing::Values(0, 1, -1, 2, -2),
    ::testing::Values(0, 1, -1, 2, -2),
    ::testing::Values(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)));

/** Test ChangeObjectiveSense */
class ChangeObjectiveSense : public Change,
                             public ::testing::WithParamInterface<std::tuple<LPObjectiveSense, LPIndex>> {
};

TEST_P(ChangeObjectiveSense, ChgObjSen) {

  LPObjectiveSense newsense = std::get<0>(GetParam());
  LPNum prob = std::get<1>(GetParam());

  LPNum nrows, ncols, nnonz;
  LPObjectiveSense sense;
  LPObjectiveSense probsense;

  ASSERT_NO_FATAL_FAILURE(initProb(prob, ncols, nrows, nnonz, sense));

  ASSERT_EQ(lp_interface_->ChangeObjectiveSense(newsense), RetCode::kOkay);
  ASSERT_TRUE(!lp_interface_->WasSolved());

  probsense = lp_interface_->GetObjectiveSense();

  ASSERT_EQ(newsense, probsense);
}

INSTANTIATE_TEST_SUITE_P(
  ChangeObjectiveSenseCombinations,
  ChangeObjectiveSense,
  ::testing::Combine(
    ::testing::Values(LPObjectiveSense::kMaximize, LPObjectiveSense::kMinimize),
    ::testing::Values(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)));

/** Test for AddRows, DeleteRowSet, DeleteRows */
TEST_F(Change, testrowmethods) {
  /* problem data */
  LPValueArray obj = {1.0, 1.0, 1.0, 1.0, 1.0};
  LPValueArray lb = {-1.0, -lp_interface_->Infinity(), 0.0, -lp_interface_->Infinity(), 0.0};
  LPValueArray ub = {10.0, lp_interface_->Infinity(), lp_interface_->Infinity(), 29.0, 0.0};
  LPNum ncolsbefore, ncolsafter;
  LPNum nrowsbefore, nrowsafter;
  LPValueArray lhsvals = {-lp_interface_->Infinity(), -1.0, -3e-10, 0.0, 1.0, 3e10};
  LPValueArray rhsvals = {-1.0, -3e-10, 0.0, 1.0, 3e10, lp_interface_->Infinity()};
  IntArray nnonzs = {1, 10, -1, 6, -1};
  LPIndexArray begvals = {0, 2, 3, 5, 8, 9};
  LPIndexArray indvals = {0, 1, 3, 2, 1, 1, 2, 4, 0, 3};
  LPValueArray vals = {1.0, 5.0, -1.0, 3e5, 2.0, 1.0, 20, 10, -1.9, 1e-2};

  LPNum iterations = 5;
  IntArray k = {1, 6, -1, 4, -2};
  IntArray nnonzsdiff = {1, 10, -1, 6, -3};
  LPIndex i;
  LPIndex j;

  /* empty placeholders */
  StringArray empty_names;
  LPValueArray empty_vals;
  LPIndexArray empty_indices;

  /* create original lp */
  ASSERT_EQ(lp_interface_->AddColumns(5, obj, lb, ub, empty_names, 0, empty_indices, empty_indices, empty_vals), RetCode::kOkay);
  ncolsbefore = lp_interface_->GetNumberOfColumns();

  for (i = 0; i < iterations; i++) {
    /* setup row values */
    int nrows;
    LPNum nnonzsbefore;
    LPNum nnonzsafter;

    nrows = k[i];
    /* get data before modification */
    nnonzsbefore = lp_interface_->GetNumberOfNonZeros();
    nrowsbefore = lp_interface_->GetNumberOfRows();

    if (nrows < 0) {
      ASSERT_EQ(lp_interface_->DeleteRows(0, -(1 + nrows)), RetCode::kOkay);
    } else { /* nrows >= 0 */
      LPValueArray lhs(100);
      LPValueArray rhs(100);
      LPIndexArray beg(100);

      LPNum nnonz = nnonzs[i];
      LPIndexArray ind(100);
      LPValueArray val(100);

      LPValueArray newlhs(100);
      LPValueArray newval(100);
      LPValueArray newrhs(100);
      LPIndexArray newbeg(100);
      LPIndexArray newind(100);
      LPNum newnnonz;
      LPIndex indold;
      LPIndex indnew;

      ASSERT_LE(nrows, 100);
      for (j = 0; j < (unsigned) nrows; j++) {
        lhs[j] = lhsvals[j];
        rhs[j] = rhsvals[j];
        beg[j] = begvals[j];
      }

      ASSERT_LE(nnonz, 100);
      for (j = 0; j < nnonz; j++) {
        ind[j] = indvals[j];
        val[j] = vals[j];
      }
      ASSERT_EQ(lp_interface_->AddRows(nrows, lhs, rhs, empty_names, nnonz, beg, ind, val), RetCode::kOkay);

      /* checks */
      ASSERT_EQ(lp_interface_->GetRows(nrowsbefore, nrowsbefore - 1 + nrows, newlhs, newrhs, newnnonz, newbeg, newind, newval), RetCode::kOkay);
      ASSERT_EQ(nnonz, newnnonz);

      for (j = 0; j < (unsigned) nrows; j++) {
        ASSERT_EQ(beg[j], newbeg[j]);
      }
      beg[nrows] = nnonz;
      newbeg[nrows] = newnnonz;

      /* check each row seperately */
      for (j = 0; j < (unsigned) nrows; j++) {
        if (fabs(lhs[j]) < 1e30 && fabs(newlhs[j]) < 1e30) {
          ASSERT_DOUBLE_EQ(lhs[j], newlhs[j]);
        }
        if (fabs(rhs[j]) < 1e30 && fabs(newrhs[j]) < 1e30) {
          ASSERT_DOUBLE_EQ(rhs[j], newrhs[j]);
        }

        /* We add a row where the indices are not sorted, some lp solvers give them back sorted (e.g. soplex), some others don't (e.g. cplex).
         * Therefore we cannot simply assert the ind and val arrays to be equal, but have to search for and check each value individually. */
        for (indold = beg[j]; indold < beg[j + 1]; indold++) {
          int occurrences = 0;

          /* for each value ind associated to the current row search for it in the newind array */
          for (indnew = beg[j]; indnew < beg[j + 1]; indnew++) {
            if (ind[indold] == newind[indnew]) {
              occurrences = occurrences + 1;
              ASSERT_DOUBLE_EQ(val[indold], newval[indnew]);
            }
          }
          /* assert that we found only one occurrence in the current row */
          ASSERT_EQ(occurrences, 1);
        }
      }
    }

    /* checks */
    nrowsafter = lp_interface_->GetNumberOfRows();
    ASSERT_EQ(nrowsbefore + nrows, nrowsafter);

    nnonzsafter = lp_interface_->GetNumberOfNonZeros();
    ASSERT_EQ(nnonzsbefore + nnonzsdiff[i], nnonzsafter);

    ncolsafter = lp_interface_->GetNumberOfColumns();
    ASSERT_EQ(ncolsbefore, ncolsafter);
  }

  /* delete rowsets */
  /* should have 8 rows now */
  nrowsbefore = lp_interface_->GetNumberOfRows();
  ASSERT_EQ(8, nrowsbefore);
  for (i = 3; i > 0; i--) {
    BoolArray rows = {false, false, false, false, false, false, false, false};

    for (j = 0; j < i; j++)
      rows[(2 * j) + 1] = true;

    nrowsbefore = lp_interface_->GetNumberOfRows();
    ASSERT_EQ(lp_interface_->DeleteRowSet(rows), RetCode::kOkay);

    nrowsafter = lp_interface_->GetNumberOfRows();
    ASSERT_EQ(nrowsbefore - i, nrowsafter);
    /* assert that the rows that are left are the ones I intended */
  }
}

/** Test for AddColumns, DeleteColumnSet, DeleteColumns */
TEST_F(Change, testcolmethods) {
  /* problem data */
  LPValueArray obj = {1.0, 1.0, 1.0, 1.0, 1.0};
  LPValueArray lhs = {-1.0, -lp_interface_->Infinity(), 0.0, -lp_interface_->Infinity(), 0.0};
  LPValueArray rhs = {10.0, lp_interface_->Infinity(), lp_interface_->Infinity(), 29.0, 0.0};
  LPNum ncolsbefore, ncolsafter;
  LPNum nrowsbefore, nrowsafter;
  LPValueArray lbvals = {-lp_interface_->Infinity(), -1.0, -3e-10, 0.0, 1.0, 3e10};
  LPValueArray ubvals = {-1.0, -3e-10, 0.0, 1.0, 3e10, lp_interface_->Infinity()};
  LPValueArray vals = {1.0, 5.0, -1.0, 3e5, 2.0, 1.0, 20, 10, -1.9, 1e-2};
  IntArray nnonzs = {1, 10, -1, 6, -1};
  IntArray begvals = {0, 2, 3, 5, 8, 9};
  LPIndexArray indvals = {0, 1, 3, 2, 1, 1, 2, 4, 0, 3};

  LPNum iterations = 5;
  IntArray k = {1, 6, -1, 4, -2};
  IntArray nnonzsdiff = {1, 10, -1, 6, -3};

  /* empty placeholders */
  StringArray empty_names;
  LPValueArray empty_vals;
  LPIndexArray empty_indices;

  /* create original lp */
  ASSERT_EQ(lp_interface_->AddRows(5, lhs, rhs, empty_names, 0, empty_indices, empty_indices, empty_vals), RetCode::kOkay);
  nrowsbefore = lp_interface_->GetNumberOfRows();

  for (LPIndex i = 0; i < iterations; i++) {
    /* setup col values */
    LPNum ncols;
    LPNum nnonzsbefore;
    LPNum nnonzsafter;

    ncols = k[i];

    /* get data before modification */
    nnonzsbefore = lp_interface_->GetNumberOfNonZeros();
    ncolsbefore = lp_interface_->GetNumberOfColumns();

    if (k[i] < 0) {
      ASSERT_EQ(lp_interface_->DeleteColumns(0, -(1 + ncols)), RetCode::kOkay);
    } else { /* ncols >= 0 */
      LPValueArray lb(100);
      LPValueArray ub(100);
      LPIndexArray beg(100);

      assert( nnonzs[i] >= 0 );
      LPNum nnonz = nnonzs[i];
      LPIndexArray ind(100);
      LPValueArray val(100);

      LPValueArray newlb(100);
      LPValueArray newval(100);
      LPValueArray newub(100);
      LPIndexArray newbeg(100);
      LPIndexArray newind(100);
      LPNum newnnonz;

      ASSERT_LE(ncols, 100);
      for (LPIndex j = 0; j < ncols; j++) {
        lb[j] = lbvals[j];
        ub[j] = ubvals[j];
        beg[j] = begvals[j];
      }

      ASSERT_LE(nnonz, 100);
      for (LPIndex j = 0; j < nnonz; j++) {
        ind[j] = indvals[j];
        val[j] = vals[j];
      }
      ASSERT_EQ(lp_interface_->AddColumns(ncols, obj, lb, ub, empty_names, nnonz, beg, ind, val), RetCode::kOkay);

      /* checks */
      ASSERT_EQ(lp_interface_->GetColumns(ncolsbefore, ncolsbefore - 1 + ncols, newlb, newub, newnnonz, newbeg, newind, newval), RetCode::kOkay);
      ASSERT_EQ(nnonz, newnnonz);

      for (LPIndex j = 0; j < ncols; j++) {
        ASSERT_EQ(lb[j], newlb[j]);
        ASSERT_EQ(ub[j], newub[j]);
        ASSERT_EQ(beg[j], newbeg[j]);
      }
      for (LPIndex j = 0; j < nnonz; j++) {
        ASSERT_EQ(ind[j], newind[j]);
        ASSERT_EQ(val[j], newval[j]);
      }
    }

    /* checks */
    nrowsafter = lp_interface_->GetNumberOfRows();
    ASSERT_EQ(nrowsbefore, nrowsafter);

    nnonzsafter = lp_interface_->GetNumberOfNonZeros();
    ASSERT_EQ(nnonzsbefore + nnonzsdiff[i], nnonzsafter);

    ncolsafter = lp_interface_->GetNumberOfColumns();
    ASSERT_EQ(ncolsbefore + ncols, ncolsafter);
  }

  /* delete rowsets */
  /* should have 8 rows now */
  ncolsbefore = lp_interface_->GetNumberOfColumns();
  ASSERT_EQ(8, ncolsbefore);
  for (LPIndex i = 3; i > 0; i--) {
    BoolArray cols = {false, false, false, false, false, false, false, false};

    for (LPIndex j = 0; j < i; j++)
      cols[(2 * j) + 1] = true;

    ncolsbefore = lp_interface_->GetNumberOfColumns();
    ASSERT_EQ(lp_interface_->DeleteColumnSet(cols), RetCode::kOkay);
    ncolsafter = lp_interface_->GetNumberOfColumns();
    ASSERT_EQ(ncolsbefore - i, ncolsafter);
    /* assert that the rows that are left are the ones I intended */
  }
}

/** Test adding zero coeffs cols */
TEST_F(Change, testzerosincols) {
  LPNum ncols;
  LPNum nrows;
  LPNum nnonz = 2;
  LPObjectiveSense sense;
  LPValueArray lb = {0};
  LPValueArray ub = {20};
  LPIndexArray beg = {0};
  LPIndexArray ind = {0, 1};
  LPValueArray val = {0, 3};
  LPValueArray obj = {1};

  /* empty placeholders */
  StringArray empty_names;

  /* 2x2 problem */
  initProb(4, ncols, nrows, nnonz, sense);
  ASSERT_EQ(2, nrows);
  ASSERT_EQ(2, ncols);

#ifndef NDEBUG
  ASSERT_DEATH(lp_interface_->AddColumns(1, obj, lb, ub, empty_names, nnonz, beg, ind, val), "");
#endif
  /* this test can only work in debug mode, so we make it pass in opt mode */
#ifdef NDEBUG
  ASSERT_EQ(lp_interface_->AddColumns(1, obj, lb, ub, empty_names, nnonz, beg, ind, val), RetCode::kOkay);
  SUCCEED(); /* return SIGABORT */
#endif
}

/** Test adding zero coeffs in rows, expecting an assert in debug mode
 *
 *  This test should fail with an assert from the which causes SIGABRT to be issued. Thus, this test should pass.
 */
TEST_F(Change, testzerosinrows) {
  LPNum nrows;
  LPNum ncols;
  LPNum nnonz = 2;
  LPObjectiveSense sense;
  LPValueArray lhs = {0};
  LPValueArray rhs = {20};
  LPIndexArray beg = {0};
  LPIndexArray ind = {0, 1};
  LPValueArray val = {0, 3};
  /* empty placeholders */
  StringArray empty_names;

  /* 2x2 problem */
  initProb(4, ncols, nrows, nnonz, sense);
  ASSERT_EQ(2, nrows);
  ASSERT_EQ(2, ncols);

#ifndef NDEBUG
  ASSERT_DEATH(lp_interface_->AddRows(1, lhs, rhs, empty_names, nnonz, beg, ind, val), "");
#else
  /* this test can only work in debug mode, so we make it pass in opt mode */
  ASSERT_EQ(lp_interface_->AddRows(1, lhs, rhs, empty_names, nnonz, beg, ind, val), RetCode::kOkay);
  SUCCEED();
#endif
}

/** test WriteLP, ReadLP, Clear */
TEST_F(Change, testlpiwritereadlpmethods) {
  LPNum nrows, ncols, nnonz;
  LPValue objval;
  LPValueArray primsol(2);
  LPValueArray dualsol(2);
  LPValueArray activity(2);
  LPValueArray redcost(2);
  LPValue objval2;
  LPValueArray primsol2(2);
  LPValueArray dualsol2(2);
  LPValueArray activity2(2);
  LPValueArray redcost2(2);
  LPObjectiveSense sense;

  /* 2x2 problem */
  ASSERT_NO_FATAL_FAILURE(initProb(5, ncols, nrows, nnonz, sense));

  ASSERT_EQ(lp_interface_->SolvePrimal(), RetCode::kOkay);
  ASSERT_EQ(lp_interface_->GetSolution(objval, primsol, dualsol, activity, redcost), RetCode::kOkay);

  ASSERT_EQ(lp_interface_->WriteLP("lpi_change_test_problem.lp"), RetCode::kOkay);
  ASSERT_EQ(lp_interface_->Clear(), RetCode::kOkay);

  if (DEF_INTERFACE == 0)
    ASSERT_EQ(lp_interface_->ReadLP("lpi_change_test_problem.lp.gz"), RetCode::kOkay);
  else
    ASSERT_EQ(lp_interface_->ReadLP("lpi_change_test_problem.lp"), RetCode::kOkay);

  ASSERT_EQ(lp_interface_->SolvePrimal(), RetCode::kOkay);
  ASSERT_EQ(lp_interface_->GetSolution(objval2, primsol2, dualsol2, activity2, redcost2), RetCode::kOkay);
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
    ASSERT_EQ(lp_interface_->WriteLP("lpi_change_test_problem2.lp"), RetCode::kOkay);
    ASSERT_EQ(lp_interface_->Clear(), RetCode::kOkay);

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
} /* namespace minimip */
