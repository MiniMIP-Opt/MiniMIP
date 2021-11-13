#include <gtest/gtest.h>

#include "absl/status/status.h"
#include "src/lp_interface/lpi_factory.h"
#include "unit_tests/utils.h"

#define DEF_INTERFACE \
  1  // 0 = Glop Interface (Default),
     // 1 = SoPlex Interface,

namespace minimip {
// TEST SUITE SIMPLE

static LPInterface* lp_interface_ = nullptr;

class matrix : public ::testing::Test {
  // simple problem
 protected:
  double lower_bound_            = 0.0;
  double upper_bound_            = 1.0;
  double left_hand_side_         = 1.0;
  double right_hand_side_        = 2.0;
  double objective_coefficients_ = 0.0;

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
  }
};

TEST_F(matrix, create_matrix) {
  SparseVector empty_coefficients = {{}, {}};
  // add one column
  ASSERT_OK(lp_interface_->AddColumn(empty_coefficients, lower_bound_,
                                     upper_bound_, objective_coefficients_,
                                     "x1"));

  // add additional identical column
  ASSERT_OK(lp_interface_->AddColumn(empty_coefficients, lower_bound_,
                                     upper_bound_, objective_coefficients_,
                                     "x2"));

  SparseVector row_coefficients = {{0}, {1.0}};
  // add one row
  ASSERT_OK(lp_interface_->AddRow(row_coefficients, left_hand_side_,
                                  right_hand_side_, "r1"));

  row_coefficients.indices[0] = 1;
  // add second row using a new variable
  ASSERT_OK(lp_interface_->AddRow(row_coefficients, left_hand_side_,
                                  right_hand_side_, "r2"));

  // ------------------------------------------------------------

  // check size
  auto num_rows    = lp_interface_->GetNumberOfRows();
  auto num_columns = lp_interface_->GetNumberOfColumns();
  ASSERT_EQ(num_rows, 2);
  ASSERT_EQ(num_columns, 2);

  // get rows
  std::vector<SparseVector> sparse_rows;
  int num_nonzeros = 0;
  std::vector<double> matrix_left_hand_sides;
  std::vector<double> matrix_right_hand_sides;

  for (int i = 0; i < num_rows; ++i) {
    matrix_left_hand_sides.push_back(lp_interface_->GetLeftHandSide(i));
    matrix_right_hand_sides.push_back(lp_interface_->GetRightHandSide(i));

    sparse_rows.push_back(lp_interface_->GetSparseRowCoefficients(i));
    num_nonzeros += static_cast<int>(sparse_rows[i].indices.size());
  }

  // check sparse matrix shape
  ASSERT_EQ(num_nonzeros, 2);
  ASSERT_EQ(sparse_rows.size(), 2);
  ASSERT_EQ(sparse_rows[0].indices[0], 0);
  ASSERT_EQ(sparse_rows[1].indices[0], 1);
  ASSERT_EQ(sparse_rows[0].indices.size(), 1);
  ASSERT_EQ(sparse_rows[1].indices.size(), 1);

  ASSERT_EQ(sparse_rows[0].indices[0], 0);
  ASSERT_EQ(sparse_rows[1].indices[0], 1);

  ASSERT_EQ(sparse_rows[0].values[0], 1.0);
  ASSERT_EQ(sparse_rows[1].values[0], 1.0);

  // check sparse matrix values
  ASSERT_FLOAT_EQ(matrix_left_hand_sides[0], 1.0);
  ASSERT_FLOAT_EQ(matrix_left_hand_sides[1], 1.0);

  ASSERT_FLOAT_EQ(matrix_right_hand_sides[0], 2.0);
  ASSERT_FLOAT_EQ(matrix_right_hand_sides[1], 2.0);
}
}  // namespace minimip
