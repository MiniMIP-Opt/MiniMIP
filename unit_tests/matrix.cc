#include "src/lp_interface/lpi_factory.h"

#include <gtest/gtest.h>
#define DEF_INTERFACE 1 /** 0 = Glop Interface (Default),
                          * 1 = SoPlex Interface, **/

namespace minimip {
/*** TEST SUITE SIMPLE ***/

static LPInterface* lp_interface_ = nullptr;

class matrix : public ::testing::Test {
 protected:
  LPValueArray obj, lb, ub, lhs, rhs, val, matval, matlhs, matrhs, row1, row2, empty_vals;
  LPIndexArray matbeg, matind, beg, ind, empty_indices;
  StringArray empty_names;

  void SetUp() override {
    /* build interface factory */
    auto* interface_factory = new LPInterfaceFactory();
    InterfaceCode interface_code;
    switch (DEF_INTERFACE) {
      case 1:
        interface_code = InterfaceCode::SOPLEX;
        break;
      default:
        interface_code = InterfaceCode::GLOP;
        break;
    }
    lp_interface_ = interface_factory->CreateLPInterface(interface_code);
    lp_interface_->ChangeObjectiveSense(LPObjectiveSense::OBJ_SENSE_MAXIMIZE);

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
  LPNum nnonz, nrows, ncols;

  /* add one column */
  ASSERT_EQ(lp_interface_->AddColumns(1, obj, lb, ub, empty_names, 0, empty_indices, empty_indices, empty_vals), RetCode::OKAY) << "hello";

  /* add additional column */
  ASSERT_EQ(lp_interface_->AddColumns(1, obj, lb, ub, empty_names, 0, empty_indices, empty_indices, empty_vals), RetCode::OKAY);

  /* add one row */
  ASSERT_EQ(lp_interface_->AddRows(1, lhs, rhs, empty_names, 1, beg, ind, val), RetCode::OKAY);

  /* add one more row using a new variable */
  ind[0] = 1;
  ASSERT_EQ(lp_interface_->AddRows(1, lhs, rhs, empty_names, 1, beg, ind, val), RetCode::OKAY);

  /* ------------------------------------------------------------ */

  /* check size */
  nrows = lp_interface_->GetNumberOfRows();
  ncols = lp_interface_->GetNumberOfColumns();
  ASSERT_EQ(nrows, 2);
  ASSERT_EQ(ncols, 2);

  /* get rows */
  ASSERT_EQ(lp_interface_->GetRows(0, 1, matlhs, matrhs, nnonz, matbeg, matind, matval), RetCode::OKAY);
  ASSERT_EQ(nnonz, 2);

  /* equal, to within 4 ULPs ( Unit in the last place) */
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
} /* namespace minimip */