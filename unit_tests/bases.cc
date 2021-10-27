#include <gtest/gtest.h>

#include <iostream>

#include "absl/status/status.h"
#include "src/lp_interface/lpi_factory.h"

#define DEF_INTERFACE \
  1  // 0 = Glop Interface (Default),
     // 1 = SoPlex Interface,

namespace minimip {

static LPInterface* lp_interface_ = nullptr;

// TEST SUITE SIMPLE
class SimpleTest : public ::testing::Test {
 protected:
  // setup for test
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

    // initialize program
    int num_rows;
    int num_cols;
    std::vector<int> begin_rows{0};
    std::vector<double> lower_bounds{0.0};
    std::vector<double> upper_bounds{3.0};
    std::vector<double> left_hand_sides{1.0};
    std::vector<double> right_hand_sides{2.0};
    std::vector<double> objective_values{1.0};
    std::vector<double> vals{1.0};
    std::vector<int> indices{0};

    // empty vectors
    std::vector<std::string> empty_names;
    std::vector<int> empty_indices;
    std::vector<double> empty_vals;

    // use the following LP as base:
    //   max x
    //       1 <= x <= 2  (linear constraint)
    //       0 <= x <= 3  (bounds)
    // add one column
    ASSERT_EQ(lp_interface_->AddColumns(
                  1, objective_values, lower_bounds, upper_bounds, empty_names,
                  0, empty_indices, empty_indices, empty_vals),
              absl::OkStatus());

    // add one row
    ASSERT_EQ(lp_interface_->AddRows(1, left_hand_sides, right_hand_sides,
                                     empty_names, 1, begin_rows, indices, vals),
              absl::OkStatus());

    // check size
    num_rows = lp_interface_->GetNumberOfRows();
    num_cols = lp_interface_->GetNumberOfColumns();
    ASSERT_EQ(num_rows, 1);
    ASSERT_EQ(num_cols, 1);
  }
};

// TESTS
TEST_F(SimpleTest, BasicAssertions) {
  std::vector<LPBasisStatus> column_basis_status(1);
  std::vector<LPBasisStatus> row_basis_status(1);

  // use LP from setup:
  //   min x
  //       1 <= x <= 2  (linear constraint)
  //       0 <= x <= 3  (bounds)
  // solve problem
  ASSERT_EQ(lp_interface_->SolveLPWithPrimalSimplex(), absl::OkStatus());

  // get basis
  absl::StatusOr<std::vector<LPBasisStatus>> absl_tmp;

  absl_tmp = lp_interface_->GetColumnBasisStatus();
  if(absl_tmp.ok())
    column_basis_status = *absl_tmp;
  else {
    std::cout<<absl_tmp.status();
  }
  absl_tmp = lp_interface_->GetRowBasisStatus();
  if(absl_tmp.ok())
    row_basis_status = *absl_tmp;
  else {
    std::cout<<absl_tmp.status();
  }

  // the variable should be basic and the slack variable at the upper bound
  ASSERT_EQ(column_basis_status[0], LPBasisStatus::kBasic);
  ASSERT_EQ(row_basis_status[0], LPBasisStatus::kAtUpperBound);
}

TEST_F(SimpleTest, test2) {
  std::vector<LPBasisStatus> column_basis_status(1);
  std::vector<LPBasisStatus> row_basis_status(1);

  // modify LP to:
  //   min x
  //       1 <= x <= 2  (linear constraint)
  //       0 <= x <= 3  (bounds)
  // change sense
  ASSERT_EQ(lp_interface_->SetObjectiveSense(LPObjectiveSense::kMinimization),
            absl::OkStatus());

  // solve problem
  ASSERT_EQ(lp_interface_->SolveLPWithPrimalSimplex(), absl::OkStatus());

  // get basis
  absl::StatusOr<std::vector<LPBasisStatus>> absl_tmp;

  absl_tmp = lp_interface_->GetColumnBasisStatus();
  if(absl_tmp.ok())
    column_basis_status = *absl_tmp;
  else {
    std::cout<<absl_tmp.status();
  }
  absl_tmp = lp_interface_->GetRowBasisStatus();
  if(absl_tmp.ok())
    row_basis_status = *absl_tmp;
  else {
    std::cout<<absl_tmp.status();
  }

  // the variable should be basic and the slack variable at the lower bound
  ASSERT_EQ(column_basis_status[0], LPBasisStatus::kBasic);
  ASSERT_EQ(row_basis_status[0], LPBasisStatus::kAtLowerBound);
}

TEST_F(SimpleTest, test3) {
  std::vector<LPBasisStatus> column_basis_status(1);
  std::vector<LPBasisStatus> row_basis_status(1);
  std::vector<double> left_hand_sides{1.0};
  std::vector<double> right_hand_sides(1.0);
  std::vector<int> indices{0};

  // modify LP to:
  //   min x
  //       1 <= x       (linear constraint)
  //       0 <= x <= 3  (bounds)
  // change sense
  ASSERT_EQ(lp_interface_->SetObjectiveSense(LPObjectiveSense::kMinimization),
            absl::OkStatus());

  // change row side
  right_hand_sides[0] = lp_interface_->Infinity();
  ASSERT_EQ(
      lp_interface_->SetRowSides(1, indices, left_hand_sides, right_hand_sides),
      absl::OkStatus());

  // solve problem
  ASSERT_EQ(lp_interface_->SolveLPWithPrimalSimplex(), absl::OkStatus());

  // get basis
  absl::StatusOr<std::vector<LPBasisStatus>> absl_tmp;

  absl_tmp = lp_interface_->GetColumnBasisStatus();
  if(absl_tmp.ok())
    column_basis_status = *absl_tmp;
  else {
    std::cout<<absl_tmp.status();
  }
  absl_tmp = lp_interface_->GetRowBasisStatus();
  if(absl_tmp.ok())
    row_basis_status = *absl_tmp;
  else {
    std::cout<<absl_tmp.status();
  }

  // the variable should be basic and the slack variable at the lower bound
  ASSERT_EQ(column_basis_status[0], LPBasisStatus::kBasic);
  ASSERT_EQ(row_basis_status[0], LPBasisStatus::kAtLowerBound);
}

TEST_F(SimpleTest, test4) {
  std::vector<LPBasisStatus> column_basis_status(1);
  std::vector<LPBasisStatus> row_basis_status(1);
  std::vector<double> left_hand_sides(1.0);
  std::vector<double> right_hand_sides{1.0};
  std::vector<int> indices{0};

  // modify LP to:
  //   max x
  //       x <= 1       (linear constraint)
  //       0 <= x <= 3  (bounds)

  // change row sides
  left_hand_sides[0] = -(lp_interface_->Infinity());
  ASSERT_EQ(
      lp_interface_->SetRowSides(1, indices, left_hand_sides, right_hand_sides),
      absl::OkStatus());

  // solve problem
  ASSERT_EQ(lp_interface_->SolveLPWithPrimalSimplex(), absl::OkStatus());

  // get basis
  absl::StatusOr<std::vector<LPBasisStatus>> absl_tmp;

  absl_tmp = lp_interface_->GetColumnBasisStatus();
  if(absl_tmp.ok())
    column_basis_status = *absl_tmp;
  else {
    std::cout<<absl_tmp.status();
  }
  absl_tmp = lp_interface_->GetRowBasisStatus();
  if(absl_tmp.ok())
    row_basis_status = *absl_tmp;
  else {
    std::cout<<absl_tmp.status();
  }

  // the variable should be basic and the slack variable at the upper bound
  ASSERT_EQ(column_basis_status[0], LPBasisStatus::kBasic);
  ASSERT_EQ(row_basis_status[0], LPBasisStatus::kAtUpperBound);
}

// TEST SUITE COMPLEX
class Complex : public ::testing::Test {
 protected:
  // empty vectors
  std::vector<std::string> empty_names;
  std::vector<int> empty_indices;
  std::vector<double> empty_vals;
  int null_int = 0;

  // setup for test
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

    // initialize program
    int ncols;
    std::vector<int> beg = {0};
    std::vector<int> inds(2);
    int nrows;
    std::vector<double> vals(2);
    std::vector<double> lb(1);
    std::vector<double> ub(1);
    std::vector<double> obj(1);
    std::vector<double> lhs(1);
    std::vector<double> rhs(1);

    // use the following LP:
    // max 1 x1 + 1 x2 + 1 x3
    //      -8 <= -x1 -          x3 <= -1
    //      -7 <= -x1 -   x2        <= -1
    //             x1 + 2 x2        <= 12
    //             x1,    x2,    x3 >= 0
    // add columns
    lb[0] = 0.0;
    ub[0] = lp_interface_->Infinity();
    obj[0] = 1.0;

    ASSERT_EQ(
        lp_interface_->AddColumns(1, obj, lb, ub, empty_names, 0, empty_indices,
                                  empty_indices, empty_vals),
        absl::OkStatus());
    ASSERT_EQ(
        lp_interface_->AddColumns(1, obj, lb, ub, empty_names, 0, empty_indices,
                                  empty_indices, empty_vals),
        absl::OkStatus());
    ASSERT_EQ(
        lp_interface_->AddColumns(1, obj, lb, ub, empty_names, 0, empty_indices,
                                  empty_indices, empty_vals),
        absl::OkStatus());

    // add rows
    lhs[0] = -8.0;
    rhs[0] = -1.0;
    inds[0] = 0;
    inds[1] = 2;
    vals[0] = -1.0;
    vals[1] = -1.0;
    ASSERT_EQ(
        lp_interface_->AddRows(1, lhs, rhs, empty_names, 2, beg, inds, vals),
        absl::OkStatus());

    lhs[0] = -7.0;
    rhs[0] = -1.0;
    inds[0] = 0;
    inds[1] = 1;
    vals[0] = -1.0;
    vals[1] = -1.0;
    ASSERT_EQ(
        lp_interface_->AddRows(1, lhs, rhs, empty_names, 2, beg, inds, vals),
        absl::OkStatus());

    lhs[0] = -lp_interface_->Infinity();
    rhs[0] = 12.0;
    inds[0] = 0;
    inds[1] = 1;
    vals[0] = 1.0;
    vals[1] = 2.0;
    ASSERT_EQ(
        lp_interface_->AddRows(1, lhs, rhs, empty_names, 2, beg, inds, vals),
        absl::OkStatus());

    // check size
    nrows = lp_interface_->GetNumberOfRows();
    ncols = lp_interface_->GetNumberOfColumns();
    ASSERT_EQ(nrows, 3);
    ASSERT_EQ(ncols, 3);
  }
};

// TESTS
TEST_F(Complex, test1) {
  std::vector<double> binvrow(3);
  std::vector<double> binvcol(3);
  std::vector<double> coef(3);
  std::vector<double> coeftwo(3);
  double objval;
  std::vector<LPBasisStatus> cstats(3);
  int nrows;
  std::vector<LPBasisStatus> rstats(3);
  std::vector<int> basinds(3);
  std::vector<int> inds(3);
  int ninds;
  int idx;
  int entry;
  int i;

  // expected values for the first column of BInv with corresponding variables
  std::vector<double> exp_vars = {-2, 1, 2};
  std::vector<double> exp_vals = {0.0, 0.0, -1.0};

  // expected values for the first column of BAInv with corresponding variables
  std::vector<double> exp_avals = {-0.5, 0.5, 1.0};

  // -------------------------------------
  // first solve problem
  ASSERT_EQ(lp_interface_->SolveLPWithPrimalSimplex(), absl::OkStatus());

  ASSERT_EQ(lp_interface_->GetObjectiveValue(objval), absl::OkStatus());
  ASSERT_FLOAT_EQ(objval, 14.0);

  // the optimal basis should be: {x2, x3, slack for second row}
  // get basis
  absl::StatusOr<std::vector<LPBasisStatus>> absl_tmp;

  absl_tmp = lp_interface_->GetColumnBasisStatus();
  if(absl_tmp.ok())
    cstats = *absl_tmp;
  else {
    std::cout<<absl_tmp.status();
  }
  absl_tmp = lp_interface_->GetRowBasisStatus();
  if(absl_tmp.ok())
    rstats = *absl_tmp;
  else {
    std::cout<<absl_tmp.status();
  }

  ASSERT_TRUE(cstats[0] == LPBasisStatus::kAtLowerBound);
  ASSERT_TRUE(cstats[1] == LPBasisStatus::kBasic);
  ASSERT_TRUE(cstats[2] == LPBasisStatus::kBasic);

  ASSERT_TRUE(rstats[0] == LPBasisStatus::kAtLowerBound);
  ASSERT_TRUE(rstats[1] == LPBasisStatus::kBasic);
  ASSERT_TRUE(rstats[2] == LPBasisStatus::kAtUpperBound);

  // get basis indices
  basinds = lp_interface_->GetBasisIndices();

  // search for slack variable in basis
  nrows = lp_interface_->GetNumberOfRows();
  for (i = 0; i < nrows; ++i) {
    if (basinds[i] < 0) break;
  }
  // assert that we found the slack variable in the basis
  ASSERT_LT(i, nrows);

  // check basis inverse for the row corresponding to the basic slack variable
  ASSERT_EQ(
      lp_interface_->GetRowOfBInverted(i, binvrow, empty_indices, null_int),
      absl::OkStatus());

  // row of basis inverse should be (0, 1, 0.5)
  ASSERT_FLOAT_EQ(binvrow[0], 0.0);
  ASSERT_FLOAT_EQ(binvrow[1], 1.0);
  ASSERT_FLOAT_EQ(binvrow[2], 0.5);

  // check whether sparse version is available and the same
  ASSERT_EQ(lp_interface_->GetRowOfBInverted(i, coef, inds, ninds),
            absl::OkStatus());
  if (ninds >= 0) {
    ASSERT_EQ(ninds, 2);
    for (entry = 0; entry < ninds; ++entry) {
      idx = inds[entry];
      ASSERT_TRUE(0 <= idx && idx < 3);
      ASSERT_FLOAT_EQ(coef[idx], binvrow[idx]);
    }
  }

  // check first column of basis inverse
  ASSERT_EQ(
      lp_interface_->GetColumnOfBInverted(0, binvcol, empty_indices, null_int),
      absl::OkStatus());

  // The columns will be in the same order, however, the rows might be permuted.
  // For each row/entry we check that it corresponds to the value of the
  // corresponding variable.
  // The correspondance variable to row/entry is given by basinds.
  for (entry = 0; entry < nrows; entry++) {
    // for the given entry try each variable in exp_vars
    for (idx = 0; idx < nrows; idx++) {
      // Check that the value is the expected one if the column corresponds to
      // the current variable given in exp_vars.
      if (exp_vars[idx] == basinds[entry]) {
        ASSERT_FLOAT_EQ(binvcol[entry], exp_vals[idx]);
      }
    }
  }

  // check whether number of nonzeros fits
  ASSERT_EQ(lp_interface_->GetColumnOfBInverted(0, coef, inds, ninds),
            absl::OkStatus());
  ASSERT_TRUE(ninds < 0 || ninds == 1);

  // check basis inverse times nonbasic matrix for row corresponding to the
  // basic slack variable
  ASSERT_TRUE(0 <= i && i < nrows);
  ASSERT_EQ(lp_interface_->GetRowOfBInvertedTimesA(i, coef,
                                                   empty_indices, null_int),
            absl::OkStatus());

  // row of basis inverse times nonbasic matrix should be (-0.5, 0, 0)
  ASSERT_FLOAT_EQ(coef[0], -0.5);
  ASSERT_FLOAT_EQ(coef[1], 0.0);
  ASSERT_FLOAT_EQ(coef[2], 0.0);

  // check nonzeros
  ASSERT_EQ(lp_interface_->GetRowOfBInvertedTimesA(i, coeftwo, inds,
                                                   ninds),
            absl::OkStatus());
  if (ninds >= 0) {
    ASSERT_EQ(ninds, 1);
    for (entry = 0; entry < ninds; ++entry) {
      idx = inds[entry];
      ASSERT_TRUE(0 <= idx && idx < 3);
      ASSERT_FLOAT_EQ(coeftwo[idx], coef[idx]);
    }
  }

  // check first column of basis inverse times nonbasic matrix
  ASSERT_EQ(lp_interface_->GetColumnOfBInvertedTimesA(0, coef, empty_indices,
                                                      null_int),
            absl::OkStatus());

  // The columns will be in the same order, however, the rows will be permuted.
  // For each row/entry we check that it corresponds to the value of the
  // corresponding variable.
  // The correspondance variable to row/entry is given by basinds.
  for (entry = 0; entry < nrows; entry++) {
    // for the given entry try each variable in exp_vars
    for (idx = 0; idx < nrows; idx++) {
      // Check that the value is the expected one if the column corresponds to
      // the current variable given in exp_vars.
      if (exp_vars[idx] == basinds[entry]) {
        ASSERT_FLOAT_EQ(coef[entry], exp_avals[idx]);
      }
    }
  }

  // check nonzeros
  ASSERT_EQ(lp_interface_->GetColumnOfBInvertedTimesA(0, coef, inds, ninds),
            absl::OkStatus());
  ASSERT_TRUE(ninds < 0 || ninds == 3);
}

// TEST SUITE MORE VARS THAN ROWS
class MoreVarsThanRows : public ::testing::Test {
 protected:
  // empty vectors
  std::vector<std::string> empty_names;
  std::vector<int> empty_indices;
  std::vector<double> empty_vals;
  int null_int = 0;

  // SetUp for Test
  void SetUp() override {
    // Build Interface Factory
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

    // initialize program
    int ncols;
    std::vector<int> beg = {0};
    std::vector<int> inds(4);
    int nrows;
    std::vector<double> vals(4);
    std::vector<double> lb(1);
    std::vector<double> ub(1);
    std::vector<double> lhs(1);
    std::vector<double> rhs(1);
    std::vector<double> obj(1);

    // use the following LP:
    // max 1 x1 + 1 x2 + 1 x3 + x4
    //      +x1          + x3 +   x4 >=   1
    //      +x1 +   x2   - x3 -   x4 >=   1
    //      -x1 - 2 x2        - 3 x4 >= -12
    //       x1,    x2,    x3,    x4 >=   0
    // add columns
    lb[0] = 0.0;
    ub[0] = lp_interface_->Infinity();
    obj[0] = 1.0;

    ASSERT_EQ(
        lp_interface_->AddColumns(1, obj, lb, ub, empty_names, 0, empty_indices,
                                  empty_indices, empty_vals),
        absl::OkStatus());
    ASSERT_EQ(
        lp_interface_->AddColumns(1, obj, lb, ub, empty_names, 0, empty_indices,
                                  empty_indices, empty_vals),
        absl::OkStatus());
    ASSERT_EQ(
        lp_interface_->AddColumns(1, obj, lb, ub, empty_names, 0, empty_indices,
                                  empty_indices, empty_vals),
        absl::OkStatus());
    ASSERT_EQ(
        lp_interface_->AddColumns(1, obj, lb, ub, empty_names, 0, empty_indices,
                                  empty_indices, empty_vals),
        absl::OkStatus());

    // add rows
    lhs[0] = 1.0;
    rhs[0] = lp_interface_->Infinity();
    inds[0] = 0;
    inds[1] = 2;
    inds[2] = 3;
    vals[0] = 1.0;
    vals[1] = 1.0;
    vals[2] = 1.0;
    ASSERT_EQ(
        lp_interface_->AddRows(1, lhs, rhs, empty_names, 3, beg, inds, vals),
        absl::OkStatus());

    lhs[0] = 1.0;
    rhs[0] = lp_interface_->Infinity();
    inds[0] = 0;
    inds[1] = 1;
    inds[2] = 2;
    inds[3] = 3;
    vals[0] = 1.0;
    vals[1] = 1.0;
    vals[2] = -1.0;
    vals[3] = -1.0;
    ASSERT_EQ(
        lp_interface_->AddRows(1, lhs, rhs, empty_names, 4, beg, inds, vals),
        absl::OkStatus());

    lhs[0] = -12;
    rhs[0] = lp_interface_->Infinity();
    inds[0] = 0;
    inds[1] = 1;
    inds[2] = 3;
    vals[0] = -1.0;
    vals[1] = -2.0;
    vals[2] = -3.0;
    ASSERT_EQ(
        lp_interface_->AddRows(1, lhs, rhs, empty_names, 3, beg, inds, vals),
        absl::OkStatus());

    // check size
    nrows = lp_interface_->GetNumberOfRows();
    ncols = lp_interface_->GetNumberOfColumns();
    ASSERT_EQ(nrows, 3);
    ASSERT_EQ(ncols, 4);
  }
};

// TESTS
TEST_F(MoreVarsThanRows, test1) {
  std::vector<double> binvarow(4);
  double objval;
  std::vector<LPBasisStatus> cstats(4);
  std::vector<LPBasisStatus> rstats(3);
  std::vector<int> basinds(3);
  int basicvarpos;

  // -------------------------------------
  // first solve problem
  ASSERT_EQ(lp_interface_->SolveLPWithPrimalSimplex(), absl::OkStatus());

  ASSERT_EQ(lp_interface_->GetObjectiveValue(objval), absl::OkStatus());
  ASSERT_FLOAT_EQ(objval, 23.0);

  // the optimal basis should be: {x1, x3, s1 = slack for first row}
  // get basis
  absl::StatusOr<std::vector<LPBasisStatus>> absl_tmp;

  absl_tmp = lp_interface_->GetColumnBasisStatus();
  if(absl_tmp.ok())
    cstats = *absl_tmp;
  else {
    std::cout<<absl_tmp.status();
  }
  absl_tmp = lp_interface_->GetRowBasisStatus();
  if(absl_tmp.ok())
    rstats = *absl_tmp;
  else {
    std::cout<<absl_tmp.status();
  }
  ASSERT_TRUE(cstats[0] == LPBasisStatus::kBasic);
  ASSERT_TRUE(cstats[1] == LPBasisStatus::kAtLowerBound);
  ASSERT_TRUE(cstats[2] == LPBasisStatus::kBasic);
  ASSERT_TRUE(cstats[3] == LPBasisStatus::kAtLowerBound);

  ASSERT_TRUE(rstats[0] == LPBasisStatus::kBasic);
  ASSERT_TRUE(rstats[1] == LPBasisStatus::kAtLowerBound);
  ASSERT_TRUE(rstats[2] == LPBasisStatus::kAtLowerBound);

  // binvarow should be
  //* 1.0   2.0  0.0   3.0  <- basic var x1
  //* 0.0   1.0  1.0   4.0  <- basic var x3
  //* 0.0  -3.0  0.0  -6.0  <- basic var s1

  // get basis indices
  basinds = lp_interface_->GetBasisIndices();

  // find position of x1 in basis indices; check binvarow of row where x1 is
  // basic
  for (basicvarpos = 0; basicvarpos < 3; ++basicvarpos) {
    if (basinds[basicvarpos] == 0) break;
  }
  ASSERT_LT(basicvarpos, 3);  // assert that we found the variable

  ASSERT_EQ(lp_interface_->GetRowOfBInvertedTimesA(
                basicvarpos, binvarow, empty_indices, null_int),
            absl::OkStatus());
  ASSERT_FLOAT_EQ(binvarow[0], 1.0);
  ASSERT_FLOAT_EQ(binvarow[1], 2.0);
  ASSERT_FLOAT_EQ(binvarow[2], 0.0);
  ASSERT_FLOAT_EQ(binvarow[3], 3.0);

  // find position of x3 in basis indices; check binvarow of row where x3 is
  // basic
  for (basicvarpos = 0; basicvarpos < 3; ++basicvarpos) {
    if (basinds[basicvarpos] == 2) break;
  }
  ASSERT_LT(basicvarpos, 3);  // assert that we found the variable

  ASSERT_EQ(lp_interface_->GetRowOfBInvertedTimesA(
                basicvarpos, binvarow, empty_indices, null_int),
            absl::OkStatus());
  ASSERT_FLOAT_EQ(binvarow[0], 0.0);
  ASSERT_FLOAT_EQ(binvarow[1], 1.0);
  ASSERT_FLOAT_EQ(binvarow[2], 1.0);
  ASSERT_FLOAT_EQ(binvarow[3], 4.0);

  // find position of s1 in basis indices; check binvarow of row where s1 is
  // basic
  for (basicvarpos = 0; basicvarpos < 3; ++basicvarpos) {
    if (basinds[basicvarpos] == -1) break;
  }
  ASSERT_LT(basicvarpos, 3);  // assert that we found the variable

  ASSERT_EQ(lp_interface_->GetRowOfBInvertedTimesA(
                basicvarpos, binvarow, empty_indices, null_int),
            absl::OkStatus());
  ASSERT_FLOAT_EQ(binvarow[0], 0.0);
  ASSERT_FLOAT_EQ(binvarow[1], -3.0);
  ASSERT_FLOAT_EQ(binvarow[2], 0.0);
  ASSERT_FLOAT_EQ(binvarow[3], -6.0);
}

}  // namespace minimip
