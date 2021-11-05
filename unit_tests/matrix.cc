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
  // 
 protected:
  double lower_bound = 0.0;
  double upper_bound = 1.0;
  double lhs = 1.0;
  double rhs = 2.0;
  double objective_coeff;

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
};

TEST_F(matrix, create_matrix) {
  minimip::SparseVector empty_coefficients = {{},{}};
  // add one column
  ASSERT_OK(lp_interface_->AddColumn(empty_coefficients, objective_coeff, lower_bound, upper_bound, "x1"));

  // add additional identical column
  ASSERT_OK(lp_interface_->AddColumn(empty_coefficients, objective_coeff, lower_bound, upper_bound, "x2"));

  minimip::SparseVector row_coefficients = {{0}, {1.0}};
  // add one row
  ASSERT_OK(lp_interface_->AddRow(row_coefficients, lhs, rhs, "r1"));

  row_coefficients.indices[0] = 1;
  // add second row using a new variable
  ASSERT_OK(lp_interface_->AddRow(row_coefficients, lhs, rhs, "r2"));

  // ------------------------------------------------------------

  // check size
  auto nrows = lp_interface_->GetNumberOfRows();
  auto ncols = lp_interface_->GetNumberOfColumns();
  ASSERT_EQ(nrows, 2);
  ASSERT_EQ(ncols, 2);

  // get rows
  std::vector<SparseVector> sparse_rows;
  int nnonz = 0;
  std::vector<double> matlhs;
  std::vector<double> matrhs;

  for (int i = 0; i < nrows; ++i) {
    matlhs.push_back(lp_interface_->GetLeftHandSide(i));
    matrhs.push_back(lp_interface_->GetRightHandSide(i));

    sparse_rows.push_back(lp_interface_->GetSparseRowCoefficients(i));
    int entries = static_cast<int>(sparse_rows[i].indices.size());
    nnonz += entries;
  }

  // check sparse matrix shape
  ASSERT_EQ(nnonz, 2);
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
  ASSERT_FLOAT_EQ(matlhs[0], 1.0);
  ASSERT_FLOAT_EQ(matlhs[1], 1.0);

  ASSERT_FLOAT_EQ(matrhs[0], 2.0);
  ASSERT_FLOAT_EQ(matrhs[1], 2.0);

}
}  // namespace minimip
