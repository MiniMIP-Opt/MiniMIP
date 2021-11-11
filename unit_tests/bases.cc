#include <gtest/gtest.h>

#include <iostream>

#include "absl/status/status.h"
#include "src/lp_interface/lpi_factory.h"
#include "unit_tests/utils.h"

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
    ASSERT_OK(
        lp_interface_->SetObjectiveSense(LPObjectiveSense::kMaximization));

    // use the following LP as base:
    //   max x
    //       1 <= x <= 2  (linear constraint)
    //       0 <= x <= 3  (bounds)

    minimip::SparseVector col_coefficients = {{}, {}};

    // add one column
    ASSERT_OK(lp_interface_->AddColumn(col_coefficients, 0.0, 3.0, 1.0, "x"));

    minimip::SparseVector row_coefficients = {{0}, {1.0}};
    // add one row
    ASSERT_OK(lp_interface_->AddRow(row_coefficients, 1.0, 2.0, "r1"));

    // check size
    auto num_rows = lp_interface_->GetNumberOfRows();
    auto num_cols = lp_interface_->GetNumberOfColumns();
    ASSERT_EQ(num_rows, 1);
    ASSERT_EQ(num_cols, 1);
  }
};

// TESTS
TEST_F(SimpleTest, BasicAssertions) {
  // use LP from setup:
  //   max x
  //       1 <= x <= 2  (linear constraint)
  //       0 <= x <= 3  (bounds)
  // solve problem
  ASSERT_OK(lp_interface_->SolveLPWithPrimalSimplex());

  ASSERT_FLOAT_EQ(lp_interface_->GetObjectiveValue(), 2.0);

  std::vector<LPBasisStatus> column_basis_status(1);
  std::vector<LPBasisStatus> row_basis_status(1);

  // get basis
  absl::StatusOr<std::vector<LPBasisStatus>> absl_tmp;

  absl_tmp = lp_interface_->GetColumnBasisStatus();
  ASSERT_OK(absl_tmp.status());
  column_basis_status = *absl_tmp;

  absl_tmp = lp_interface_->GetRowBasisStatus();
  ASSERT_OK(absl_tmp.status());
  row_basis_status = *absl_tmp;

  // the variable should be basic and the slack variable at the upper bound
  ASSERT_EQ(column_basis_status[0], LPBasisStatus::kBasic);
  ASSERT_EQ(row_basis_status[0], LPBasisStatus::kAtUpperBound);
}

TEST_F(SimpleTest, test2) {
  // modify LP to:
  //   min x
  //       1 <= x <= 2  (linear constraint)
  //       0 <= x <= 3  (bounds)
  // change sense
  ASSERT_OK(lp_interface_->SetObjectiveSense(LPObjectiveSense::kMinimization))

  // solve problem
  ASSERT_OK(lp_interface_->SolveLPWithPrimalSimplex());

  ASSERT_FLOAT_EQ(lp_interface_->GetObjectiveValue(), 1.0);

  std::vector<LPBasisStatus> column_basis_status(1);
  std::vector<LPBasisStatus> row_basis_status(1);

  // get basis
  absl::StatusOr<std::vector<LPBasisStatus>> absl_tmp;

  absl_tmp = lp_interface_->GetColumnBasisStatus();
  ASSERT_OK(absl_tmp.status());
  column_basis_status = *absl_tmp;

  absl_tmp = lp_interface_->GetRowBasisStatus();
  ASSERT_OK(absl_tmp.status());
  row_basis_status = *absl_tmp;

  // the variable should be basic and the slack variable at the lower bound
  ASSERT_EQ(column_basis_status[0], LPBasisStatus::kBasic);
  ASSERT_EQ(row_basis_status[0], LPBasisStatus::kAtLowerBound);
}

TEST_F(SimpleTest, test3) {
  // modify LP to:
  //   min x
  //       1 <= x       (linear constraint)
  //       0 <= x <= 3  (bounds)
  // change sense
  ASSERT_EQ(lp_interface_->SetObjectiveSense(LPObjectiveSense::kMinimization),
            absl::OkStatus());

  ASSERT_OK(lp_interface_->SetRowSides(0, 1.0, lp_interface_->Infinity()));

  // solve problem
  ASSERT_OK(lp_interface_->SolveLPWithPrimalSimplex());

  ASSERT_FLOAT_EQ(lp_interface_->GetObjectiveValue(), 1.0);

  std::vector<LPBasisStatus> column_basis_status(1);
  std::vector<LPBasisStatus> row_basis_status(1);

  // get basis
  absl::StatusOr<std::vector<LPBasisStatus>> absl_tmp;

  absl_tmp = lp_interface_->GetColumnBasisStatus();
  ASSERT_OK(absl_tmp.status());
  column_basis_status = *absl_tmp;
  absl_tmp            = lp_interface_->GetRowBasisStatus();
  ASSERT_OK(absl_tmp.status());
  row_basis_status = *absl_tmp;

  // the variable should be basic and the slack variable at the lower bound
  ASSERT_EQ(column_basis_status[0], LPBasisStatus::kBasic);
  ASSERT_EQ(row_basis_status[0], LPBasisStatus::kAtLowerBound);
}

TEST_F(SimpleTest, test4) {
  // modify LP to:
  //   max x
  //       x <= 1       (linear constraint)
  //       0 <= x <= 3  (bounds)

  // change row sides
  ASSERT_OK(lp_interface_->SetRowSides(0, -(lp_interface_->Infinity()), 1.0));

  // solve problem
  ASSERT_OK(lp_interface_->SolveLPWithPrimalSimplex());

  ASSERT_FLOAT_EQ(lp_interface_->GetObjectiveValue(), 1.0);

  std::vector<LPBasisStatus> column_basis_status(1);
  std::vector<LPBasisStatus> row_basis_status(1);

  // get basis
  absl::StatusOr<std::vector<LPBasisStatus>> absl_tmp;

  absl_tmp = lp_interface_->GetColumnBasisStatus();
  ASSERT_OK(absl_tmp.status());
  column_basis_status = *absl_tmp;

  absl_tmp = lp_interface_->GetRowBasisStatus();
  ASSERT_OK(absl_tmp.status());
  row_basis_status = *absl_tmp;

  // the variable should be basic and the slack variable at the upper bound
  ASSERT_EQ(column_basis_status[0], LPBasisStatus::kBasic);
  ASSERT_EQ(row_basis_status[0], LPBasisStatus::kAtUpperBound);
}

// TEST SUITE COMPLEX
class Complex : public ::testing::Test {
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
    ASSERT_OK(
        lp_interface_->SetObjectiveSense(LPObjectiveSense::kMaximization));

    // initialize program
    int num_columns;
    std::vector<int> beg = {0};
    std::vector<int> indices(2);
    int num_rows;
    std::vector<double> objective_coefficients(1);

    // use the following LP:
    // max 1 x1 + 1 x2 + 1 x3
    //      -8 <= -x1           -x3 <= -1
    //      -7 <= -x1 -   x2        <= -1
    //             x1 + 2 x2        <= 12
    //             x1,    x2,    x3 >= 0

    double inf = lp_interface_->Infinity();
    ASSERT_OK(lp_interface_->AddColumn({{}, {}}, 0.0, inf, 1.0, "x1"));
    std::vector<std::string> names = {"x2", "x3"};
    ASSERT_OK(lp_interface_->AddColumns({{{}, {}}, {{}, {}}}, {0.0, 0.0},
                                        {inf, inf}, {1.0, 1.0}, names));

    minimip::SparseVector row_coefficients = {{0, 2}, {-1.0, -1.0}};
    ASSERT_OK(lp_interface_->AddRow(row_coefficients, -8.0, -1.0, "r1"));

    std::vector<std::string> row_names          = {"r2", "r3"};
    std::vector<SparseVector> rows_coefficients = {{{0, 1}, {-1.0, -1.0}},
                                                   {{0, 1}, {1.0, 2.0}}};
    ASSERT_OK(lp_interface_->AddRows(rows_coefficients, {-7.0, -inf},
                                     {-1.0, 12.0}, row_names));

    // check size
    num_rows    = lp_interface_->GetNumberOfRows();
    num_columns = lp_interface_->GetNumberOfColumns();
    ASSERT_EQ(num_rows, 3);
    ASSERT_EQ(num_columns, 3);
  }
};

// TESTS
TEST_F(Complex, test1) {
  std::vector<double> b_inverted_row(3);
  std::vector<double> b_inverted_column(3);
  std::vector<double> b_inverted(3);
  std::vector<double> coeftwo(3);
  std::vector<LPBasisStatus> column_basis_status(3);
  int num_rows;
  std::vector<LPBasisStatus> row_basis_status(3);
  std::vector<int> basis_indices(3);
  std::vector<int> indices(3);
  int num_indices;
  int idx;
  int entry;
  int i;

  // expected values for the first column of BInv with corresponding variables
  std::vector<double> expected_b_invert_variables = {-2, 1, 2};
  std::vector<double> expected_b_invert_values    = {0.0, 0.0, -1.0};

  // expected values for the first column of BAInv with corresponding variables
  std::vector<double> expected_b_invert_times_a_values = {-0.5, 0.5, 1.0};

  // -------------------------------------
  // first solve problem
  ASSERT_OK(lp_interface_->SolveLPWithPrimalSimplex());

  ASSERT_FLOAT_EQ(lp_interface_->GetObjectiveValue(), 14.0);

  // the optimal basis should be: {x2, x3, slack for second row}
  // get basis
  absl::StatusOr<std::vector<LPBasisStatus>> absl_tmp;

  absl_tmp = lp_interface_->GetColumnBasisStatus();
  ASSERT_OK(absl_tmp.status());
  column_basis_status = *absl_tmp;

  absl_tmp = lp_interface_->GetRowBasisStatus();
  ASSERT_OK(absl_tmp.status());
  row_basis_status = *absl_tmp;

  ASSERT_TRUE(column_basis_status[0] == LPBasisStatus::kAtLowerBound);
  ASSERT_TRUE(column_basis_status[1] == LPBasisStatus::kBasic);
  ASSERT_TRUE(column_basis_status[2] == LPBasisStatus::kBasic);

  ASSERT_TRUE(row_basis_status[0] == LPBasisStatus::kAtLowerBound);
  ASSERT_TRUE(row_basis_status[1] == LPBasisStatus::kBasic);
  ASSERT_TRUE(row_basis_status[2] == LPBasisStatus::kAtUpperBound);

  // get basis indices
  basis_indices = lp_interface_->GetBasisIndices();

  // search for slack variable in basis
  num_rows = lp_interface_->GetNumberOfRows();
  for (i = 0; i < num_rows; ++i) {
    if (basis_indices[i] < 0) break;
  }
  // assert that we found the slack variable in the basis
  ASSERT_LT(i, num_rows);

  // check basis inverse for the row corresponding to the basic slack variable

  absl::StatusOr<SparseVector> absl_tmp_sparse =
      lp_interface_->GetSparseRowOfBInverted(i);
  ASSERT_OK(absl_tmp_sparse.status());
  b_inverted  = absl_tmp_sparse->values;
  indices     = absl_tmp_sparse->indices;
  num_indices = indices.size();

  // row of basis inverse should be (0, 1, 0.5)
  ASSERT_EQ(indices.size(), 2);
  ASSERT_EQ(indices[0], 1);
  ASSERT_EQ(indices[1], 2);
  ASSERT_FLOAT_EQ(b_inverted[0], 1.0);
  ASSERT_FLOAT_EQ(b_inverted[1], 0.5);

  // check first column of basis inverse
  absl_tmp_sparse = lp_interface_->GetSparseColumnOfBInverted(0);
  ASSERT_OK(absl_tmp_sparse.status());
  b_inverted  = absl_tmp_sparse->values;
  indices     = absl_tmp_sparse->indices;
  num_indices = indices.size();

  // column of basis inverse should be (0, 0, -1.0)
  ASSERT_EQ(indices.size(), 1);

  std::vector<double> b_inverted_dense(num_rows);
  for (size_t j = 0; j < indices.size(); ++j) {
    b_inverted_dense[indices[j]] = b_inverted[j];
  }

  // The columns will be in the same order, however, the rows might be
  // permuted.
  // For each row/entry we check that it corresponds to the value of the
  // corresponding variable.
  // The correspondance variable to row/entry is given by basis_indices.
  for (entry = 0; entry < num_rows; entry++) {
    // for the given entry try each variable in expected_b_invert_variables
    for (idx = 0; idx < num_rows; idx++) {
      // Check that the value is the expected one if the column corresponds
      // to
      // the current variable given in expected_b_invert_variables.
      if (expected_b_invert_variables[idx] == basis_indices[entry]) {
        ASSERT_FLOAT_EQ(b_inverted_dense[entry], expected_b_invert_values[idx]);
      }
    }
  }

  // check whether number of nonzeros fits
  ASSERT_TRUE(num_indices < 0 || num_indices == 1);

  // check basis inverse times nonbasic matrix for row corresponding to the
  // basic slack variable
  ASSERT_TRUE(0 <= i && i < num_rows);
  absl_tmp_sparse = lp_interface_->GetSparseRowOfBInvertedTimesA(i);

  ASSERT_OK(absl_tmp_sparse.status());
  indices    = absl_tmp_sparse->indices;
  b_inverted = absl_tmp_sparse->values;

  // row of basis inverse times nonbasic matrix should be (-0.5, 0, 0)
  ASSERT_EQ(indices.size(), 1);
  ASSERT_EQ(indices[0], 0);
  ASSERT_FLOAT_EQ(b_inverted[0], -0.5);

  // check first column of basis inverse times nonbasic matrix
  absl_tmp_sparse = lp_interface_->GetSparseColumnOfBInvertedTimesA(0);
  ASSERT_OK(absl_tmp_sparse.status());
  b_inverted  = absl_tmp_sparse->values;
  indices     = absl_tmp_sparse->indices;
  num_indices = absl_tmp_sparse->indices.size();

  b_inverted_dense.resize(num_rows, 0.0);
  for (size_t j = 0; j < indices.size(); ++j) {
    b_inverted_dense[indices[j]] = b_inverted[j];
  }

  // The columns will be in the same order, however, the rows will be permuted.
  // For each row/entry we check that it corresponds to the value of the
  // corresponding variable.
  // The correspondance variable to row/entry is given by basis_indices.
  for (entry = 0; entry < num_rows; entry++) {
    // for the given entry try each variable in expected_b_invert_variables
    for (idx = 0; idx < num_rows; idx++) {
      // Check that the value is the expected one if the column corresponds to
      // the current variable given in expected_b_invert_variables.
      if (expected_b_invert_variables[idx] == basis_indices[entry]) {
        ASSERT_FLOAT_EQ(b_inverted_dense[entry],
                        expected_b_invert_times_a_values[idx]);
      }
    }
  }

  // check nonzeros
  ASSERT_TRUE(num_indices == 3);
}

// TEST SUITE MORE VARS THAN ROWS
class MoreVarsThanRows : public ::testing::Test {
 protected:
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
    ASSERT_OK(
        lp_interface_->SetObjectiveSense(LPObjectiveSense::kMaximization));

    // use the following LP:
    // max 1 x1 + 1 x2 + 1 x3 + x4
    //      +x1          + x3 +   x4 >=   1
    //      +x1 +   x2   - x3 -   x4 >=   1
    //      -x1 - 2 x2        - 3 x4 >= -12
    //       x1,    x2,    x3,    x4 >=   0

    // initialize program
    std::vector<double> objective_coefficients = {1.0, 1.0, 1.0, 1.0};
    std::vector<double> lower_bounds           = {0.0, 0.0, 0.0, 0.0};
    std::vector<double> upper_bounds           = {
        lp_interface_->Infinity(), lp_interface_->Infinity(),
        lp_interface_->Infinity(), lp_interface_->Infinity()};

    // add columns
    ASSERT_OK(lp_interface_->AddColumn({{}, {}}, 0.0, lp_interface_->Infinity(),
                                       1.0, "x1"));
    ASSERT_OK(lp_interface_->AddColumn({{}, {}}, 0.0, lp_interface_->Infinity(),
                                       1.0, "x2"));
    ASSERT_OK(lp_interface_->AddColumn({{}, {}}, 0.0, lp_interface_->Infinity(),
                                       1.0, "x3"));
    ASSERT_OK(lp_interface_->AddColumn({{}, {}}, 0.0, lp_interface_->Infinity(),
                                       1.0, "x4"));

    // add one row
    ASSERT_OK(lp_interface_->AddRow({{0, 2, 3}, {1, 1, 1}}, 1,
                                    lp_interface_->Infinity(), "r1"));

    std::vector<std::string> row_names          = {"r2", "r3"};
    std::vector<SparseVector> rows_coefficients = {
        {{0, 1, 2, 3}, {1.0, 1.0, -1.0, -1.0}},
        {{0, 1, 3}, {-1.0, -2.0, -3.0}}};
    ASSERT_OK(lp_interface_->AddRows(
        rows_coefficients, {1.0, -12.0},
        {lp_interface_->Infinity(), lp_interface_->Infinity()}, row_names));

    int num_columns;
    int num_rows;

    // check size
    num_rows    = lp_interface_->GetNumberOfRows();
    num_columns = lp_interface_->GetNumberOfColumns();
    ASSERT_EQ(num_rows, 3);
    ASSERT_EQ(num_columns, 4);
  }
};

// TESTS
TEST_F(MoreVarsThanRows, test1) {
  // -------------------------------------
  // first solve problem
  ASSERT_EQ(lp_interface_->SolveLPWithPrimalSimplex(), absl::OkStatus());

  double objective_value;
  objective_value = lp_interface_->GetObjectiveValue();
  ASSERT_FLOAT_EQ(objective_value, 23.0);

  std::vector<LPBasisStatus> column_basis_status(4);
  std::vector<LPBasisStatus> row_basis_status(3);

  // the optimal basis should be: {x1, x3, s1 = slack for first row}
  // get basis
  absl::StatusOr<std::vector<LPBasisStatus>> absl_tmp;

  absl_tmp = lp_interface_->GetColumnBasisStatus();
  ASSERT_OK(absl_tmp.status());
  column_basis_status = *absl_tmp;

  absl_tmp = lp_interface_->GetRowBasisStatus();
  ASSERT_OK(absl_tmp.status());
  row_basis_status = *absl_tmp;

  ASSERT_TRUE(column_basis_status[0] == LPBasisStatus::kBasic);
  ASSERT_TRUE(column_basis_status[1] == LPBasisStatus::kAtLowerBound);
  ASSERT_TRUE(column_basis_status[2] == LPBasisStatus::kBasic);
  ASSERT_TRUE(column_basis_status[3] == LPBasisStatus::kAtLowerBound);

  ASSERT_TRUE(row_basis_status[0] == LPBasisStatus::kBasic);
  ASSERT_TRUE(row_basis_status[1] == LPBasisStatus::kAtLowerBound);
  ASSERT_TRUE(row_basis_status[2] == LPBasisStatus::kAtLowerBound);

  // // b_inverted_times_a_row should be
  // //* 1.0   2.0  0.0   3.0  <- basic var x1
  // //* 0.0   1.0  1.0   4.0  <- basic var x3
  // //* 0.0  -3.0  0.0  -6.0  <- basic var s1

  int basic_variable_positions;
  std::vector<int> basis_indices(3);

  // get basis indices
  basis_indices = lp_interface_->GetBasisIndices();

  // find position of x1 in basis indices; check b_inverted_times_a_row of row
  // where x1 is basic
  for (basic_variable_positions = 0; basic_variable_positions < 3;
       ++basic_variable_positions) {
    if (basis_indices[basic_variable_positions] == 0) break;
  }
  ASSERT_LT(basic_variable_positions, 3);  // assert that we found the variable

  std::vector<int> b_inverted_times_a_indices;
  std::vector<double> b_inverted_times_a_row;

  absl::StatusOr<SparseVector> absl_tmp_sparse =
      lp_interface_->GetSparseRowOfBInvertedTimesA(basic_variable_positions);
  ASSERT_OK(absl_tmp_sparse.status());

  b_inverted_times_a_indices = absl_tmp_sparse->indices;
  b_inverted_times_a_row     = absl_tmp_sparse->values;

  ASSERT_EQ(b_inverted_times_a_indices[0], 0);
  ASSERT_EQ(b_inverted_times_a_indices[1], 1);
  ASSERT_EQ(b_inverted_times_a_indices[2], 3);

  ASSERT_FLOAT_EQ(b_inverted_times_a_row[0], 1.0);
  ASSERT_FLOAT_EQ(b_inverted_times_a_row[1], 2.0);
  ASSERT_FLOAT_EQ(b_inverted_times_a_row[2], 3.0);

  // find position of x3 in basis indices; check b_inverted_times_a_row of row
  // where x3 is basic
  for (basic_variable_positions = 0; basic_variable_positions < 3;
       ++basic_variable_positions) {
    if (basis_indices[basic_variable_positions] == 2) break;
  }
  ASSERT_LT(basic_variable_positions, 3);  // assert that we found the variable

  absl_tmp_sparse =
      lp_interface_->GetSparseRowOfBInvertedTimesA(basic_variable_positions);
  ASSERT_OK(absl_tmp_sparse.status());
  b_inverted_times_a_indices = absl_tmp_sparse->indices;
  b_inverted_times_a_row     = absl_tmp_sparse->values;

  ASSERT_EQ(b_inverted_times_a_indices[0], 1);
  ASSERT_EQ(b_inverted_times_a_indices[1], 2);
  ASSERT_EQ(b_inverted_times_a_indices[2], 3);

  ASSERT_FLOAT_EQ(b_inverted_times_a_row[0], 1.0);
  ASSERT_FLOAT_EQ(b_inverted_times_a_row[1], 1.0);
  ASSERT_FLOAT_EQ(b_inverted_times_a_row[2], 4.0);

  // find position of s1 in basis indices; check b_inverted_times_a_row of row
  // where s1 is basic
  for (basic_variable_positions = 0; basic_variable_positions < 3;
       ++basic_variable_positions) {
    if (basis_indices[basic_variable_positions] == -1) break;
  }
  ASSERT_LT(basic_variable_positions, 3);  // assert that we found the variable

  absl_tmp_sparse =
      lp_interface_->GetSparseRowOfBInvertedTimesA(basic_variable_positions);
  ASSERT_OK(absl_tmp_sparse.status());
  b_inverted_times_a_indices = absl_tmp_sparse->indices;
  b_inverted_times_a_row     = absl_tmp_sparse->values;

  ASSERT_EQ(b_inverted_times_a_indices[0], 1);
  ASSERT_EQ(b_inverted_times_a_indices[1], 3);

  ASSERT_FLOAT_EQ(b_inverted_times_a_row[0], -3.0);
  ASSERT_FLOAT_EQ(b_inverted_times_a_row[1], -6.0);
}

}  // namespace minimip
