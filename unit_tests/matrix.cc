#include <gtest/gtest.h>

#include "absl/status/status.h"
#include "src/lp_interface/lpi_factory.h"
#include "src/lp_interface/sparse_types.h"
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
  SparseVector empty_coefficients = {};
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

  row_coefficients.ReplaceEntry(0, 1, 1.0);
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
    num_nonzeros += static_cast<int>(sparse_rows[i].NumNonZeros());
  }

  // check sparse matrix shape
  ASSERT_EQ(num_nonzeros, 2);
  ASSERT_EQ(sparse_rows.size(), 2);
  ASSERT_EQ(sparse_rows[0].IndicesData()[0], 0);
  ASSERT_EQ(sparse_rows[1].IndicesData()[0], 1);
  ASSERT_EQ(sparse_rows[0].NumNonZeros(), 1);
  ASSERT_EQ(sparse_rows[1].NumNonZeros(), 1);

  ASSERT_EQ(sparse_rows[0].IndicesData()[0], 0);
  ASSERT_EQ(sparse_rows[1].IndicesData()[0], 1);

  ASSERT_EQ(sparse_rows[0].ValuesData()[0], 1.0);
  ASSERT_EQ(sparse_rows[1].ValuesData()[0], 1.0);

  // check sparse matrix values
  ASSERT_FLOAT_EQ(matrix_left_hand_sides[0], 1.0);
  ASSERT_FLOAT_EQ(matrix_left_hand_sides[1], 1.0);

  ASSERT_FLOAT_EQ(matrix_right_hand_sides[0], 2.0);
  ASSERT_FLOAT_EQ(matrix_right_hand_sides[1], 2.0);
}

TEST(SparseVectorTests, ViewFromSparse) {
  // SparseViewVector holds the same data as original array
  int indices[]   = {1, 5, 7};
  double values[] = {1.0, 5.0, 7.0};
  auto v          = SparseViewVector{indices, values, 3};
  ASSERT_EQ(v.NumNonZeros(), 3);
  ASSERT_FLOAT_EQ(v.ValueAt(1), 1.0);
  ASSERT_FLOAT_EQ(v.ValueAt(5), 5.0);
  ASSERT_FLOAT_EQ(v.ValueAt(7), 7.0);
  ASSERT_FLOAT_EQ(v.ValueAt(4), 0.0);
  values[2] = 42.0;
  ASSERT_FLOAT_EQ(v.ValueAt(7), 42.0);
  auto v_alloc = v.CopyToOwned();
  values[2]    = 1.0;
  ASSERT_FLOAT_EQ(v.ValueAt(7), 1.0);
  ASSERT_FLOAT_EQ(v_alloc.ValueAt(7), 42.0);
}

TEST(SparseVectorTests, InsertSorted) {
  auto v = SparseVector();
  ASSERT_EQ(v.NumNonZeros(), 0);
  v.InsertSorted(5, 4.5);
  v.InsertSorted(2, 2.5);
  // inserting should maintain order
  ASSERT_EQ(v.NumNonZeros(), 2);
  ASSERT_EQ(v.IndicesData()[0], 2);
  ASSERT_EQ(v.IndicesData()[1], 5);
  ASSERT_FLOAT_EQ(v.ValuesData()[1], 4.5);
  // inserting at existing index should replace existing value
  v.InsertSorted(5, 0.5);
  ASSERT_EQ(v.NumNonZeros(), 2);
  ASSERT_FLOAT_EQ(v.ValuesData()[1], 0.5);
  // inserting at intermediate index should add an entry
  v.InsertSorted(3, 1.5);
  ASSERT_EQ(v.NumNonZeros(), 3);
  ASSERT_EQ(v.IndicesData()[0], 2);
  ASSERT_EQ(v.IndicesData()[1], 3);
  ASSERT_EQ(v.IndicesData()[2], 5);
  ASSERT_FLOAT_EQ(v.ValuesData()[0], 2.5);
  ASSERT_FLOAT_EQ(v.ValuesData()[1], 1.5);
  ASSERT_FLOAT_EQ(v.ValuesData()[2], 0.5);
  // inserting at the end of vector
  v.InsertSorted(7, 3.5);
  ASSERT_EQ(v.NumNonZeros(), 4);
  ASSERT_FLOAT_EQ(v.ValuesData()[3], 3.5);
  ASSERT_EQ(v.IndicesData()[3], 7);
  // inserting 0 at new index should be a no-op
  v.InsertSorted(8, 0.0);
  ASSERT_FLOAT_EQ(v.ValueAt(8), 0.0);
  ASSERT_EQ(v.NumNonZeros(), 4);
  // inserting 0 at existing index should decrease number of entries
  v.InsertSorted(7, 0.0);
  ASSERT_FLOAT_EQ(v.ValueAt(7), 0.0);
  ASSERT_EQ(v.NumNonZeros(), 3);
}

TEST(SparseVectorTests, CreateFromEmpty) {
  auto v = SparseVector({}, {});
  ASSERT_TRUE(v.empty());
  ASSERT_EQ(v.NumNonZeros(), 0);
}

}  // namespace minimip
