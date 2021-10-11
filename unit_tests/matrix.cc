#include <gtest/gtest.h>

#include "absl/status/status.h"
#include "src/lp_interface/lpi_factory.h"

#define DEF_INTERFACE \
  1  // 0 = Glop Interface (Default),
     // 1 = SoPlex Interface,

namespace minimip {
// TEST SUITE SIMPLE

static LPInterface* lp_interface_ = nullptr;

class matrix : public ::testing::Test {
 protected:
  std::vector<double> obj, lb, ub, lhs, rhs, val, matval, matlhs, matrhs, row1,
      row2, empty_vals;
  std::vector<int> matbeg, matind, beg, ind, empty_indices;
  std::vector<std::string> empty_names;

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

    obj.push_back(0.0);
    lb.push_back(0.0);
    ub.push_back(1.0);
    lhs.push_back(1.0);
    rhs.push_back(2.0);
    val.push_back(1.0);
    matval.reserve(2);
    row1.reserve(2);
    row2.reserve(2);
    matbeg.reserve(2);
    matind.reserve(2);
    matlhs.reserve(2);
    matrhs.reserve(2);
    beg.push_back(0);
    ind.push_back(0);
  }
};

TEST_F(matrix, create_matrix) {
  int nnonz, nrows, ncols;

  // add one column
  ASSERT_EQ(lp_interface_->AddColumns(1, obj, lb, ub, empty_names, 0,
                                      empty_indices, empty_indices, empty_vals),
            absl::OkStatus());

  // add additional column
  ASSERT_EQ(lp_interface_->AddColumns(1, obj, lb, ub, empty_names, 0,
                                      empty_indices, empty_indices, empty_vals),
            absl::OkStatus());

  // add one row
  ASSERT_EQ(lp_interface_->AddRows(1, lhs, rhs, empty_names, 1, beg, ind, val),
            absl::OkStatus());

  // add one more row using a new variable
  ind[0] = 1;
  ASSERT_EQ(lp_interface_->AddRows(1, lhs, rhs, empty_names, 1, beg, ind, val),
            absl::OkStatus());

  // ------------------------------------------------------------

  // check size
  nrows = lp_interface_->GetNumberOfRows();
  ncols = lp_interface_->GetNumberOfColumns();
  ASSERT_EQ(nrows, 2);
  ASSERT_EQ(ncols, 2);

  // get rows
  ASSERT_EQ(lp_interface_->GetRows(0, 1, matlhs, matrhs, nnonz, matbeg, matind,
                                   matval),
            absl::OkStatus());
  ASSERT_EQ(nnonz, 2);

  // equal, to within 4 ULPs ( Unit in the last place)
  ASSERT_FLOAT_EQ(matlhs[0], 1.0);
  ASSERT_FLOAT_EQ(matlhs[1], 1.0);

  ASSERT_FLOAT_EQ(matrhs[0], 2.0);
  ASSERT_FLOAT_EQ(matrhs[1], 2.0);

  ASSERT_EQ(matbeg[0], 0);
  ASSERT_EQ(matbeg[1], 1);

  ASSERT_EQ(matind[0], 0);
  ASSERT_EQ(matind[1], 1);

  ASSERT_FLOAT_EQ(matval[0], 1.0);
  ASSERT_FLOAT_EQ(matval[1], 1.0);
}
}  // namespace minimip
