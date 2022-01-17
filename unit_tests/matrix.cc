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

TEST(SparseMatrixTests, ConstructionAndAccessors) {
  auto m = SparseCompressedMatrix(3, 3);
  ASSERT_EQ(m.num_nonzeros(), 0);
  ASSERT_EQ(m.values.size(), 0);
  ASSERT_EQ(m.column_indices.size(), 4);
  ASSERT_EQ(m.row_indices.size(), 0);

  // test matrix with entries
  // 1.0   ⋅    ⋅
  // 3.0   ⋅    ⋅
  // 5.0  2.0  7.0
  std::vector<int> row_indices    = {0, 1, 2, 2, 2};
  std::vector<int> column_indices = {0, 3, 4, 5};
  std::vector<double> vals        = {1.0, 3.0, 5.0, 2.0, 7.0};

  auto m_filled =
      SparseCompressedMatrix(3, 3, column_indices, row_indices, vals);

  // elements match with a (row, col, val) ordering
  std::vector<int> col_ordering = {0, 0, 0, 1, 2};
  for (auto elem_idx = 0; elem_idx < 5; ++elem_idx) {
    ASSERT_FLOAT_EQ(m_filled.at(row_indices[elem_idx], col_ordering[elem_idx]),
                    vals[elem_idx]);
  }

  // verify structural zeros
  ASSERT_FLOAT_EQ(m_filled.at(0, 1), 0.0);
  ASSERT_FLOAT_EQ(m_filled.at(0, 2), 0.0);
  ASSERT_FLOAT_EQ(m_filled.at(1, 1), 0.0);
  ASSERT_FLOAT_EQ(m_filled.at(1, 2), 0.0);

  // test on non-square matrix
  // 1.0   ⋅    ⋅
  // 3.0   ⋅    ⋅
  // 5.0  2.0  7.0
  // 5.0   ⋅   7.0
  std::vector<int> row_indices2 = {0, 1, 2, 3, 2, 2, 3};
  std::vector<int> col_indices2 = {0, 4, 5, 7};
  std::vector<double> nzvals2   = {1.0, 3.0, 5.0, 5.0, 2.0, 7.0, 7.0};

  auto m_nonsquare =
      SparseCompressedMatrix(4, 3, col_indices2, row_indices2, nzvals2);
  std::vector<int> col_entries = {0, 0, 0, 0, 1, 2, 2};
  for (size_t elem_idx = 0; elem_idx < nzvals2.size(); ++elem_idx) {
    ASSERT_FLOAT_EQ(
        m_nonsquare.at(row_indices2[elem_idx], col_entries[elem_idx]),
        nzvals2[elem_idx]);
  }
  ASSERT_FLOAT_EQ(m_nonsquare.at(0, 1), 0.0);
  ASSERT_FLOAT_EQ(m_nonsquare.at(0, 2), 0.0);
  ASSERT_FLOAT_EQ(m_nonsquare.at(1, 1), 0.0);
  ASSERT_FLOAT_EQ(m_nonsquare.at(1, 2), 0.0);
  ASSERT_FLOAT_EQ(m_nonsquare.at(3, 1), 0.0);

  // copy constructor
  SparseCompressedMatrix m_ident = m_nonsquare;
  for (size_t elem_idx = 0; elem_idx < nzvals2.size(); ++elem_idx) {
    ASSERT_FLOAT_EQ(m_ident.at(row_indices2[elem_idx], col_entries[elem_idx]),
                    nzvals2[elem_idx]);
  }
  ASSERT_FLOAT_EQ(m_ident.at(0, 1), 0.0);
  ASSERT_FLOAT_EQ(m_ident.at(0, 2), 0.0);
  ASSERT_FLOAT_EQ(m_ident.at(1, 1), 0.0);
  ASSERT_FLOAT_EQ(m_ident.at(1, 2), 0.0);
  ASSERT_FLOAT_EQ(m_ident.at(3, 1), 0.0);
  // the copy is independent from the original
  m_ident.insert(1, 1, 5.0);
  ASSERT_EQ(m_ident.num_nonzeros(), nzvals2.size() + 1);
  ASSERT_EQ(m_nonsquare.num_nonzeros(), nzvals2.size());
  ASSERT_FLOAT_EQ(m_nonsquare.at(1, 1), 0.0);
  ASSERT_FLOAT_EQ(m_ident.at(1, 1), 5.0);
}

TEST(SparseMatrixTests, MatrixMutation) {
  auto nrows = 3;
  auto ncols = 4;
  auto mnew  = SparseCompressedMatrix(nrows, ncols);
  ASSERT_EQ(mnew.num_nonzeros(), 0);
  ASSERT_EQ(mnew.num_rows_, nrows);
  ASSERT_EQ(mnew.num_cols_, ncols);
  ASSERT_EQ(mnew.values.size(), 0);
  ASSERT_EQ(mnew.column_indices.size(), ncols + 1);
  ASSERT_EQ(mnew.row_indices.size(), 0);
  mnew.insert(0, 0, 1.0);
  ASSERT_EQ(mnew.num_nonzeros(), 1);
  ASSERT_EQ(mnew.row_indices[0], 0);
  ASSERT_FLOAT_EQ(mnew.values[0], 1.0);
  ASSERT_EQ(mnew.column_indices[0], 0);
  ASSERT_EQ(mnew.column_indices[1], 1);
  ASSERT_EQ(mnew.column_indices[2], 1);
  ASSERT_FLOAT_EQ(mnew.at(0, 0), 1.0);
  // add on same column after
  mnew.insert(2, 0, 3.0);
  ASSERT_EQ(mnew.num_nonzeros(), 2);
  ASSERT_EQ(mnew.row_indices.size(), 2);
  ASSERT_EQ(mnew.values.size(), 2);
  ASSERT_FLOAT_EQ(mnew.at(2, 0), 3.0);
  // inserting 0 over existing entry
  mnew.insert(2, 0, 0.0);
  ASSERT_FLOAT_EQ(mnew.at(2, 0), 0.0);
  ASSERT_EQ(mnew.num_nonzeros(), 2);
  // inserting 0 at new entry is a no-op
  ASSERT_FLOAT_EQ(mnew.at(2, 1), 0.0);
  ASSERT_EQ(mnew.num_nonzeros(), 2);
  ASSERT_FLOAT_EQ(mnew.at(2, 1), 0.0);

  // inserting at new column
  mnew.insert(2, 1, 4.0).insert(1, 1, 2.0);

  ASSERT_EQ(mnew.num_nonzeros(), 4);
  ASSERT_FLOAT_EQ(mnew.at(2, 1), 4.0);
  ASSERT_FLOAT_EQ(mnew.at(1, 1), 2.0);
}

}  // namespace minimip
