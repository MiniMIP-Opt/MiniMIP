// Copyright 2023 the MiniMIP Project
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "src/lp_interface/lpi.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <cstdio>
#include <limits>
#include <numeric>
#include <sstream>
#include <tuple>

#include "absl/status/status.h"
#include "ortools/base/status_macros.h"
#include "src/data_structures/strong_sparse_vector.h"
#include "src/lp_interface/lpi_factory.h"
#include "src/parameters.pb.h"
#include "unit_tests/utils.h"

using testing::DoubleEq;
using testing::ElementsAre;
using testing::ExplainMatchResult;
using testing::Gt;
using testing::IsEmpty;
using testing::Le;
using testing::Lt;
using testing::Pointwise;
using testing::StartsWith;
using testing::UnorderedElementsAre;

namespace minimip {
namespace {

class LPInterfaceImplementationTest
    : public ::testing::TestWithParam<LPParameters::SolverType> {
 protected:
  void SetUp() override {
    LPParameters lp_params;
    lp_params.set_lp_solver_type(GetParam());
    lpi_ = ConfigureLPSolverFromProto(lp_params).value();
    solver_inf_ = lpi_->Infinity();
  }

  absl::Status UploadProblem(
      bool is_maximization,
      const absl::StrongVector<ColIndex, double>& lower_bounds,
      const absl::StrongVector<ColIndex, double>& upper_bounds,
      const absl::StrongVector<ColIndex, double>& objective_coefficients,
      const absl::StrongVector<RowIndex, double>& left_hand_sides,
      const absl::StrongVector<RowIndex, double>& right_hand_sides,
      const absl::StrongVector<ColIndex, absl::StrongVector<RowIndex, double>>&
          matrix) {
    // We use CHECK rather than ASSERT since they check if the test is written
    // properly, not that the lpi behaves correctly.
    CHECK_EQ(lower_bounds.size(), upper_bounds.size());
    CHECK_EQ(lower_bounds.size(), objective_coefficients.size());
    CHECK_EQ(left_hand_sides.size(), right_hand_sides.size());
    CHECK_EQ(lower_bounds.size(), matrix.size());
    for (const auto& row : matrix) CHECK_EQ(row.size(), left_hand_sides.size());

    RETURN_IF_ERROR(lpi_->Clear());
    RETURN_IF_ERROR(lpi_->SetObjectiveSense(is_maximization));
    for (ColIndex col(0); col < lower_bounds.size(); ++col) {
      RETURN_IF_ERROR(lpi_->AddColumn(SparseCol(), lower_bounds[col],
                                      upper_bounds[col],
                                      objective_coefficients[col], ""));
    }
    for (RowIndex row(0); row < left_hand_sides.size(); ++row) {
      SparseRow row_coefficients;
      for (ColIndex col(0); col < lower_bounds.size(); ++col) {
        if (matrix[col][row] != 0) {
          row_coefficients.AddEntry(col, matrix[col][row]);
        }
      }
      RETURN_IF_ERROR(lpi_->AddRow(row_coefficients, left_hand_sides[row],
                                   right_hand_sides[row], ""));
    }
    return absl::OkStatus();
  }

  absl::StrongVector<ColIndex, double> ExtractObjectiveCoefficients() {
    absl::StrongVector<ColIndex, double> objective_coefficients(
        lpi_->GetNumberOfColumns().value());
    for (ColIndex col(0); col < lpi_->GetNumberOfColumns(); ++col) {
      objective_coefficients[col] = lpi_->GetObjectiveCoefficient(col);
    }
    return objective_coefficients;
  }

  absl::StrongVector<ColIndex, double> ExtractLowerBounds() {
    absl::StrongVector<ColIndex, double> lower_bounds(
        lpi_->GetNumberOfColumns().value());
    for (ColIndex col(0); col < lpi_->GetNumberOfColumns(); ++col) {
      lower_bounds[col] = lpi_->GetLowerBound(col);
    }
    return lower_bounds;
  }

  absl::StrongVector<ColIndex, double> ExtractUpperBounds() {
    absl::StrongVector<ColIndex, double> upper_bounds(
        lpi_->GetNumberOfColumns().value());
    for (ColIndex col(0); col < lpi_->GetNumberOfColumns(); ++col) {
      upper_bounds[col] = lpi_->GetUpperBound(col);
    }
    return upper_bounds;
  }

  absl::StrongVector<RowIndex, double> ExtractLeftHandSides() {
    absl::StrongVector<RowIndex, double> left_hand_sides(
        lpi_->GetNumberOfRows().value());
    for (RowIndex row(0); row < lpi_->GetNumberOfRows(); ++row) {
      left_hand_sides[row] = lpi_->GetLeftHandSide(row);
    }
    return left_hand_sides;
  }

  absl::StrongVector<RowIndex, double> ExtractRightHandSides() {
    absl::StrongVector<RowIndex, double> right_hand_sides(
        lpi_->GetNumberOfRows().value());
    for (RowIndex row(0); row < lpi_->GetNumberOfRows(); ++row) {
      right_hand_sides[row] = lpi_->GetRightHandSide(row);
    }
    return right_hand_sides;
  }

  absl::StrongVector<ColIndex, absl::StrongVector<RowIndex, double>>
  ExtractMatrix() {
    absl::StrongVector<ColIndex, absl::StrongVector<RowIndex, double>> matrix(
        lpi_->GetNumberOfColumns().value(),
        absl::StrongVector<RowIndex, double>(lpi_->GetNumberOfRows().value()));
    for (ColIndex col(0); col < lpi_->GetNumberOfColumns(); ++col) {
      for (RowIndex row(0); row < lpi_->GetNumberOfRows(); ++row) {
        matrix[col][row] = lpi_->GetMatrixCoefficient(col, row);
      }
    }
    return matrix;
  }

  absl::Status InitEmptyProblem(int n_cols, int n_rows) {
    RETURN_IF_ERROR(lpi_->Clear());
    for (int i = 0; i < n_cols; ++i) {
      RETURN_IF_ERROR(lpi_->AddColumn(SparseCol(), -lpi_->Infinity(),
                                      lpi_->Infinity(), 0.0, ""));
    }
    for (int i = 0; i < n_rows; ++i) {
      RETURN_IF_ERROR(
          lpi_->AddRow(SparseRow(), -lpi_->Infinity(), lpi_->Infinity(), ""));
    }
    return absl::OkStatus();
  }

  // Creates a problem with a simple structure
  // obj: max sum(i*x_i) for i = 1..n_cols
  // constraints:
  // -j <= sum(((i-1)*n_rows + j)*x_i) <= j for i = 1..n_cols, j = 1..n_rows
  // i.e. the coefficients in column-major order is 1, 2, ...
  // bounds:
  // -i <= x_i <= i for i = 1..n_cols
  absl::Status InitSimpleProblem(int n_cols, int n_rows) {
    RETURN_IF_ERROR(lpi_->Clear());
    {
      absl::StrongVector<ColIndex, double> lower_bounds(n_cols);
      std::iota(lower_bounds.rbegin(), lower_bounds.rend(), -n_cols);
      absl::StrongVector<ColIndex, double> upper_bounds(n_cols);
      std::iota(upper_bounds.begin(), upper_bounds.end(), 1.0);
      absl::StrongVector<ColIndex, double> objective_coefficients(n_cols);
      std::iota(objective_coefficients.begin(), objective_coefficients.end(),
                1.0);
      absl::StrongVector<ColIndex, std::string> names(n_cols);
      RETURN_IF_ERROR(lpi_->AddColumns(
          StrongSparseMatrix(ColIndex(n_cols), RowIndex(0)), lower_bounds,
          upper_bounds, objective_coefficients, names));
    }
    {
      absl::StrongVector<RowIndex, SparseRow> matrix(n_rows);
      for (RowIndex row(0); row < n_rows; ++row) {
        for (ColIndex col(0); col < n_cols; ++col) {
          matrix[row].AddEntry(col, col.value() * n_rows + row.value() + 1);
        }
        matrix[row].CleanUpIfNeeded();
      }
      absl::StrongVector<RowIndex, double> left_hand_sides(n_rows);
      std::iota(left_hand_sides.rbegin(), left_hand_sides.rend(), -n_rows);
      absl::StrongVector<RowIndex, double> right_hand_sides(n_rows);
      std::iota(right_hand_sides.begin(), right_hand_sides.end(), 1.0);
      absl::StrongVector<RowIndex, std::string> names(n_rows, "hi");
      RETURN_IF_ERROR(
          lpi_->AddRows(matrix, left_hand_sides, right_hand_sides, names));
    }
    return absl::OkStatus();
  }

  double solver_inf_;
  std::unique_ptr<LPInterface> lpi_;
};

class ModelConstructionTest : public LPInterfaceImplementationTest {};

INSTANTIATE_TEST_SUITE_P(All, ModelConstructionTest,
                         testing::ValuesIn({LPParameters::LP_GLOP,
                                            LPParameters::LP_SOPLEX}));

TEST_P(ModelConstructionTest, AddColumnRow) {
  ASSERT_EQ(lpi_->GetNumberOfColumns(), 0);
  ASSERT_EQ(lpi_->GetNumberOfRows(), 0);

  // This should build the following model
  //
  // obj: max x1 + 2*x2 -4*x3 + 1e-3x4 + 0*x5
  // constraints:
  //                  x1 +  5*x2                          <= -1     (ct1)
  // -1       <=                              -x4         <= -3e-10 (ct2)
  // -3e-10   <=            2*x2 + 3e5*x3                 <= 0      (ct3)
  // 0        <=              x2 +  20*x3         + 10*x5 <= 1      (ct4)
  // 1        <= -1.9*x1                                  <= 3e10   (ct5)
  // 3e10     <=                          1e-2*x4                   (ct6)
  // bounds:
  // -1 <= x1 <= 10
  //       x2
  //  0 <= x3
  //       x4 <= 29.3
  //  1 <= x5 <= 1
  ASSERT_OK(lpi_->SetObjectiveSense(true));
  ASSERT_OK(lpi_->AddColumn(SparseCol(), -1, 10, 1, "x1"));
  ASSERT_OK(lpi_->AddColumn(SparseCol(), -solver_inf_, solver_inf_, 2, "x2"));
  ASSERT_OK(lpi_->AddRow(SparseRow({{ColIndex(0), 1}, {ColIndex(1), 5}}),
                         -solver_inf_, -1, "ct1"));
  ASSERT_OK(lpi_->AddRow(SparseRow(), -1, -3e-10, "ct2"));
  ASSERT_OK(lpi_->AddColumn(SparseCol(), 0, solver_inf_, -4, "x3"));
  ASSERT_OK(lpi_->AddColumn(SparseCol({{RowIndex(1), -1}}), -solver_inf_, 29.3,
                            1e-3, "x4"));
  ASSERT_OK(lpi_->AddRow(SparseRow({{ColIndex(1), 2}, {ColIndex(2), 3e5}}),
                         -3e-10, 0, "ct3"));
  ASSERT_OK(lpi_->AddRow(SparseRow({{ColIndex(1), 1}, {ColIndex(2), 20}}), 0, 1,
                         "ct4"));
  ASSERT_OK(lpi_->AddColumn(SparseCol({{RowIndex(3), 10}}), 1, 1, 0, "x5"));
  ASSERT_OK(lpi_->AddRow(SparseRow({{ColIndex(0), -1.9}}), 1, 3e10, "ct5"));
  ASSERT_OK(
      lpi_->AddRow(SparseRow({{ColIndex(3), 1e-2}}), 3e10, solver_inf_, "ct6"));

  // Check variable properties
  ASSERT_EQ(lpi_->GetNumberOfColumns(), 5);
  absl::StrongVector<ColIndex, double> lower_bounds(5);
  absl::StrongVector<ColIndex, double> upper_bounds(5);
  absl::StrongVector<ColIndex, double> obj_coeff(5);
  absl::StrongVector<ColIndex, SparseCol> columns(5);
  for (ColIndex col(0); col < 5; ++col) {
    lower_bounds[col] = lpi_->GetLowerBound(col);
    upper_bounds[col] = lpi_->GetUpperBound(col);
    obj_coeff[col] = lpi_->GetObjectiveCoefficient(col);
  }
  EXPECT_THAT(lower_bounds, ElementsAre(-1, -solver_inf_, 0, -solver_inf_, 1));
  EXPECT_THAT(upper_bounds, ElementsAre(10, solver_inf_, solver_inf_, 29.3, 1));
  EXPECT_THAT(obj_coeff, ElementsAre(1, 2, -4, 1e-3, 0));

  // Check row properties
  ASSERT_EQ(lpi_->GetNumberOfRows(), 6);
  absl::StrongVector<RowIndex, double> left_hand_sides(6);
  absl::StrongVector<RowIndex, double> right_hand_sides(6);
  for (RowIndex row(0); row < 6; ++row) {
    left_hand_sides[row] = lpi_->GetLeftHandSide(row);
    right_hand_sides[row] = lpi_->GetRightHandSide(row);
  }
  EXPECT_THAT(left_hand_sides,
              ElementsAre(-solver_inf_, -1, -3e-10, 0, 1, 3e10));
  EXPECT_THAT(right_hand_sides,
              ElementsAre(-1, -3e-10, 0, 1, 3e10, solver_inf_));

  // Check constraint matrix
  EXPECT_EQ(lpi_->GetNumberOfNonZeros(), 10);
  const absl::StrongVector<RowIndex, absl::StrongVector<ColIndex, double>>
      expected_matrix{{1, 5, 0, 0, 0},    {0, 0, 0, -1, 0},
                      {0, 2, 3e5, 0, 0},  {0, 1, 20, 0, 10},
                      {-1.9, 0, 0, 0, 0}, {0, 0, 0, 1e-2, 0}};
  for (RowIndex row(0); row < 6; ++row) {
    for (ColIndex col(0); col < 5; ++col) {
      SCOPED_TRACE(absl::StrFormat("row=%d, col=%d", row.value(), col.value()));
      EXPECT_EQ(lpi_->GetSparseRowCoefficients(row).value(col),
                expected_matrix[row][col]);
      EXPECT_EQ(lpi_->GetSparseColumnCoefficients(col).value(row),
                expected_matrix[row][col]);
      EXPECT_EQ(lpi_->GetMatrixCoefficient(col, row),
                expected_matrix[row][col]);
    }
  }
}

TEST_P(ModelConstructionTest, AddColumnsRowsBatched) {
  ASSERT_EQ(lpi_->GetNumberOfColumns(), 0);
  ASSERT_EQ(lpi_->GetNumberOfRows(), 0);

  // This should build the following model
  //
  // obj: max x1 + 2*x2 -4*x3 + 1e-3x4 + 0*x5
  // constraints:
  //                  x1 +  5*x2                          <= -1     (ct1)
  // -1       <=                              -x4         <= -3e-10 (ct2)
  // -3e-10   <=            2*x2 + 3e5*x3                 <= 0      (ct3)
  // 0        <=              x2 +  20*x3         + 10*x5 <= 1      (ct4)
  // 1        <= -1.9*x1                                  <= 3e10   (ct5)
  // 3e10     <=                          1e-2*x4                   (ct6)
  // bounds:
  // -1 <= x1 <= 10
  //       x2
  //  0 <= x3
  //       x4 <= 29.3
  //  1 <= x5 <= 1
  {
    StrongSparseMatrix matrix(ColIndex(2), RowIndex(0));
    absl::StrongVector<ColIndex, double> lower_bounds = {-1, -solver_inf_};
    absl::StrongVector<ColIndex, double> upper_bounds = {10, solver_inf_};
    absl::StrongVector<ColIndex, double> objective_coefficients = {1, 2};
    absl::StrongVector<ColIndex, std::string> names = {"x1", "x2"};
    ASSERT_OK(lpi_->AddColumns(matrix, lower_bounds, upper_bounds,
                               objective_coefficients, names));
  }
  {
    absl::StrongVector<RowIndex, SparseRow> rows = {
        SparseRow({{ColIndex(0), 1}, {ColIndex(1), 5}}), SparseRow()};
    absl::StrongVector<RowIndex, double> left_hand_sides = {-solver_inf_, -1};
    absl::StrongVector<RowIndex, double> right_hand_sides = {-1, -3e-10};
    absl::StrongVector<RowIndex, std::string> names = {"ct1", "ct2"};
    ASSERT_OK(lpi_->AddRows(rows, left_hand_sides, right_hand_sides, names));
  }
  {
    StrongSparseMatrix matrix(ColIndex(2), RowIndex(2));
    matrix.PopulateCol(ColIndex(1), SparseCol({{RowIndex(1), -1}}));
    absl::StrongVector<ColIndex, double> lower_bounds = {0, -solver_inf_};
    absl::StrongVector<ColIndex, double> upper_bounds = {solver_inf_, 29.3};
    absl::StrongVector<ColIndex, double> objective_coefficients = {-4, 1e-3};
    absl::StrongVector<ColIndex, std::string> names = {"x3", "x4"};
    ASSERT_OK(lpi_->AddColumns(matrix, lower_bounds, upper_bounds,
                               objective_coefficients, names));
  }
  {
    absl::StrongVector<RowIndex, SparseRow> rows = {
        SparseRow({{ColIndex(1), 2}, {ColIndex(2), 3e5}}),
        SparseRow({{ColIndex(1), 1}, {ColIndex(2), 20}})};
    absl::StrongVector<RowIndex, double> left_hand_sides = {-3e-10, 0};
    absl::StrongVector<RowIndex, double> right_hand_sides = {0, 1};
    absl::StrongVector<RowIndex, std::string> names = {"ct3", "ct4"};
    ASSERT_OK(lpi_->AddRows(rows, left_hand_sides, right_hand_sides, names));
  }
  {
    StrongSparseMatrix matrix(ColIndex(1), RowIndex(4));
    matrix.PopulateCol(ColIndex(0), SparseCol({{RowIndex(3), 10}}));
    absl::StrongVector<ColIndex, double> lower_bounds = {1};
    absl::StrongVector<ColIndex, double> upper_bounds = {1};
    absl::StrongVector<ColIndex, double> objective_coefficients = {0};
    absl::StrongVector<ColIndex, std::string> names = {"x5"};
    ASSERT_OK(lpi_->AddColumns(matrix, lower_bounds, upper_bounds,
                               objective_coefficients, names));
  }
  {
    absl::StrongVector<RowIndex, SparseRow> rows = {
        SparseRow({{ColIndex(0), -1.9}}), SparseRow({{ColIndex(3), 1e-2}})};
    absl::StrongVector<RowIndex, double> left_hand_sides = {1, 3e10};
    absl::StrongVector<RowIndex, double> right_hand_sides = {3e10, solver_inf_};
    absl::StrongVector<RowIndex, std::string> names = {"ct5", "ct6"};
    ASSERT_OK(lpi_->AddRows(rows, left_hand_sides, right_hand_sides, names));
  }

  // Check variable properties
  EXPECT_EQ(lpi_->GetNumberOfColumns(), 5);
  EXPECT_THAT(ExtractLowerBounds(),
              ElementsAre(-1, -solver_inf_, 0, -solver_inf_, 1));
  EXPECT_THAT(ExtractUpperBounds(),
              ElementsAre(10, solver_inf_, solver_inf_, 29.3, 1));
  EXPECT_THAT(ExtractObjectiveCoefficients(), ElementsAre(1, 2, -4, 1e-3, 0));

  // Check row properties
  EXPECT_EQ(lpi_->GetNumberOfRows(), 6);
  EXPECT_THAT(ExtractLeftHandSides(),
              ElementsAre(-solver_inf_, -1, -3e-10, 0, 1, 3e10));
  EXPECT_THAT(ExtractRightHandSides(),
              ElementsAre(-1, -3e-10, 0, 1, 3e10, solver_inf_));

  // Check constraint matrix
  EXPECT_EQ(lpi_->GetNumberOfNonZeros(), 10);
  const absl::StrongVector<ColIndex, absl::StrongVector<RowIndex, double>>
      extracted_matrix = ExtractMatrix();
  EXPECT_THAT(extracted_matrix, ElementsAre(ElementsAre(1, 0, 0, 0, -1.9, 0),
                                            ElementsAre(5, 0, 2, 1, 0, 0),
                                            ElementsAre(0, 0, 3e5, 20, 0, 0),
                                            ElementsAre(0, -1, 0, 0, 0, 1e-2),
                                            ElementsAre(0, 0, 0, 10, 0, 0)));

  // Check that the different ways to retrieve matrix coefficients are
  // consistent.
  for (ColIndex col(0); col < 5; ++col) {
    for (RowIndex row(0); row < 6; ++row) {
      SCOPED_TRACE(absl::StrFormat("row=%d, col=%d", row.value(), col.value()));
      EXPECT_EQ(lpi_->GetSparseRowCoefficients(row).value(col),
                extracted_matrix[col][row]);
      EXPECT_EQ(lpi_->GetSparseColumnCoefficients(col).value(row),
                extracted_matrix[col][row]);
      EXPECT_EQ(lpi_->GetMatrixCoefficient(col, row),
                extracted_matrix[col][row]);
    }
  }
}

TEST_P(ModelConstructionTest, Clear) {
  ASSERT_OK(InitSimpleProblem(4, 5));
  ASSERT_OK(lpi_->SolveLPWithPrimalSimplex());
  EXPECT_TRUE(lpi_->IsSolved());
  EXPECT_EQ(lpi_->GetNumberOfColumns(), 4);
  EXPECT_EQ(lpi_->GetNumberOfRows(), 5);
  EXPECT_EQ(lpi_->GetNumberOfNonZeros(), 20);

  ASSERT_OK(lpi_->Clear());
  EXPECT_FALSE(lpi_->IsSolved());
  EXPECT_EQ(lpi_->GetNumberOfColumns(), 0);
  EXPECT_EQ(lpi_->GetNumberOfRows(), 0);
  EXPECT_EQ(lpi_->GetNumberOfNonZeros(), 0);

  // Make sure the previous state is properly overwritten.
  ASSERT_OK(InitEmptyProblem(3, 2));
  EXPECT_EQ(lpi_->GetNumberOfColumns(), 3);
  EXPECT_EQ(lpi_->GetNumberOfRows(), 2);
  EXPECT_EQ(lpi_->GetNumberOfNonZeros(), 0);
  EXPECT_THAT(ExtractObjectiveCoefficients(), ElementsAre(0.0, 0.0, 0.0));
  EXPECT_THAT(ExtractLowerBounds(),
              ElementsAre(-solver_inf_, -solver_inf_, -solver_inf_));
  EXPECT_THAT(ExtractUpperBounds(),
              ElementsAre(solver_inf_, solver_inf_, solver_inf_));
  EXPECT_THAT(ExtractLeftHandSides(), ElementsAre(-solver_inf_, -solver_inf_));
  EXPECT_THAT(ExtractRightHandSides(), ElementsAre(solver_inf_, solver_inf_));
  EXPECT_THAT(ExtractMatrix(), ElementsAre(ElementsAre(0, 0), ElementsAre(0, 0),
                                           ElementsAre(0, 0)));
}

TEST_P(ModelConstructionTest, ClearStateDoesNotClearModel) {
  ASSERT_OK(InitSimpleProblem(2, 2));
  ASSERT_OK(lpi_->ClearState());
  EXPECT_EQ(lpi_->GetNumberOfColumns(), 2);
  EXPECT_EQ(lpi_->GetNumberOfRows(), 2);
  EXPECT_EQ(lpi_->GetNumberOfNonZeros(), 4);
  EXPECT_THAT(ExtractObjectiveCoefficients(), ElementsAre(1, 2));
  EXPECT_THAT(ExtractLowerBounds(), ElementsAre(-1, -2));
  EXPECT_THAT(ExtractUpperBounds(), ElementsAre(1, 2));
  EXPECT_THAT(ExtractLeftHandSides(), ElementsAre(-1, -2));
  EXPECT_THAT(ExtractRightHandSides(), ElementsAre(1, 2));
  EXPECT_THAT(ExtractMatrix(),
              ElementsAre(ElementsAre(1, 2), ElementsAre(3, 4)));
}

TEST_P(ModelConstructionTest, DeleteSingleRow) {
  ASSERT_OK(InitSimpleProblem(4, 5));
  ASSERT_OK(lpi_->SolveLPWithPrimalSimplex());
  EXPECT_TRUE(lpi_->IsSolved());
  ASSERT_OK(lpi_->DeleteRows(RowIndex(0), RowIndex(0)));
  EXPECT_FALSE(lpi_->IsSolved());
  EXPECT_EQ(lpi_->GetNumberOfColumns(), 4);
  EXPECT_EQ(lpi_->GetNumberOfRows(), 4);
  EXPECT_THAT(ExtractLeftHandSides(), ElementsAre(-2, -3, -4, -5));
  EXPECT_THAT(ExtractRightHandSides(), ElementsAre(2, 3, 4, 5));
  EXPECT_THAT(ExtractLowerBounds(), ElementsAre(-1, -2, -3, -4));
  EXPECT_THAT(ExtractUpperBounds(), ElementsAre(1, 2, 3, 4));
  EXPECT_EQ(lpi_->GetNumberOfNonZeros(), 16);
  EXPECT_THAT(
      ExtractMatrix(),
      ElementsAre(ElementsAre(2, 3, 4, 5), ElementsAre(7, 8, 9, 10),
                  ElementsAre(12, 13, 14, 15), ElementsAre(17, 18, 19, 20)));

  ASSERT_OK(InitSimpleProblem(4, 5));
  ASSERT_OK(lpi_->SolveLPWithPrimalSimplex());
  EXPECT_TRUE(lpi_->IsSolved());
  ASSERT_OK(lpi_->DeleteRows(RowIndex(2), RowIndex(2)));
  EXPECT_FALSE(lpi_->IsSolved());
  EXPECT_EQ(lpi_->GetNumberOfColumns(), 4);
  EXPECT_EQ(lpi_->GetNumberOfRows(), 4);
  EXPECT_THAT(ExtractLeftHandSides(), ElementsAre(-1, -2, -4, -5));
  EXPECT_THAT(ExtractRightHandSides(), ElementsAre(1, 2, 4, 5));
  EXPECT_THAT(ExtractLowerBounds(), ElementsAre(-1, -2, -3, -4));
  EXPECT_THAT(ExtractUpperBounds(), ElementsAre(1, 2, 3, 4));
  EXPECT_EQ(lpi_->GetNumberOfNonZeros(), 16);
  EXPECT_THAT(
      ExtractMatrix(),
      ElementsAre(ElementsAre(1, 2, 4, 5), ElementsAre(6, 7, 9, 10),
                  ElementsAre(11, 12, 14, 15), ElementsAre(16, 17, 19, 20)));

  ASSERT_OK(InitSimpleProblem(4, 5));
  ASSERT_OK(lpi_->SolveLPWithPrimalSimplex());
  EXPECT_TRUE(lpi_->IsSolved());
  ASSERT_OK(lpi_->DeleteRows(RowIndex(4), RowIndex(4)));
  EXPECT_FALSE(lpi_->IsSolved());
  EXPECT_EQ(lpi_->GetNumberOfColumns(), 4);
  EXPECT_EQ(lpi_->GetNumberOfRows(), 4);
  EXPECT_THAT(ExtractLeftHandSides(), ElementsAre(-1, -2, -3, -4));
  EXPECT_THAT(ExtractRightHandSides(), ElementsAre(1, 2, 3, 4));
  EXPECT_THAT(ExtractLowerBounds(), ElementsAre(-1, -2, -3, -4));
  EXPECT_THAT(ExtractUpperBounds(), ElementsAre(1, 2, 3, 4));
  EXPECT_EQ(lpi_->GetNumberOfNonZeros(), 16);
  EXPECT_THAT(
      ExtractMatrix(),
      ElementsAre(ElementsAre(1, 2, 3, 4), ElementsAre(6, 7, 8, 9),
                  ElementsAre(11, 12, 13, 14), ElementsAre(16, 17, 18, 19)));
}

TEST_P(ModelConstructionTest, DeleteMultipleRows) {
  ASSERT_OK(InitSimpleProblem(4, 5));
  ASSERT_OK(lpi_->SolveLPWithPrimalSimplex());
  EXPECT_TRUE(lpi_->IsSolved());
  ASSERT_OK(lpi_->DeleteRows(RowIndex(2), RowIndex(3)));
  EXPECT_FALSE(lpi_->IsSolved());
  EXPECT_EQ(lpi_->GetNumberOfColumns(), 4);
  EXPECT_EQ(lpi_->GetNumberOfRows(), 3);
  EXPECT_THAT(ExtractLeftHandSides(), ElementsAre(-1, -2, -5));
  EXPECT_THAT(ExtractRightHandSides(), ElementsAre(1, 2, 5));
  EXPECT_THAT(ExtractLowerBounds(), ElementsAre(-1, -2, -3, -4));
  EXPECT_THAT(ExtractUpperBounds(), ElementsAre(1, 2, 3, 4));
  EXPECT_EQ(lpi_->GetNumberOfNonZeros(), 12);
  EXPECT_THAT(ExtractMatrix(),
              ElementsAre(ElementsAre(1, 2, 5), ElementsAre(6, 7, 10),
                          ElementsAre(11, 12, 15), ElementsAre(16, 17, 20)));

  ASSERT_OK(InitSimpleProblem(4, 5));
  ASSERT_OK(lpi_->SolveLPWithPrimalSimplex());
  EXPECT_TRUE(lpi_->IsSolved());
  ASSERT_OK(lpi_->DeleteRows(RowIndex(0), RowIndex(4)));
  EXPECT_FALSE(lpi_->IsSolved());
  EXPECT_EQ(lpi_->GetNumberOfColumns(), 4);
  EXPECT_EQ(lpi_->GetNumberOfRows(), 0);
  EXPECT_THAT(ExtractLeftHandSides(), ElementsAre());
  EXPECT_THAT(ExtractRightHandSides(), ElementsAre());
  EXPECT_THAT(ExtractLowerBounds(), ElementsAre(-1, -2, -3, -4));
  EXPECT_THAT(ExtractUpperBounds(), ElementsAre(1, 2, 3, 4));
  EXPECT_EQ(lpi_->GetNumberOfNonZeros(), 0);
  EXPECT_THAT(ExtractMatrix(),
              ElementsAre(IsEmpty(), IsEmpty(), IsEmpty(), IsEmpty()));
}

TEST_P(ModelConstructionTest, DeleteRowSet) {
  {
    ASSERT_OK(InitSimpleProblem(4, 5));
    ASSERT_OK(lpi_->SolveLPWithPrimalSimplex());
    EXPECT_TRUE(lpi_->IsSolved());
    ASSERT_OK_AND_ASSIGN(
        (const absl::StrongVector<RowIndex, RowIndex> row_mapping),
        lpi_->DeleteRowSet({true, false, true, false, true}));
    EXPECT_THAT(row_mapping, ElementsAre(kInvalidRow, RowIndex(0), kInvalidRow,
                                         RowIndex(1), kInvalidRow));
    EXPECT_FALSE(lpi_->IsSolved());
    EXPECT_EQ(lpi_->GetNumberOfColumns(), 4);
    EXPECT_EQ(lpi_->GetNumberOfRows(), 2);
    EXPECT_THAT(ExtractLeftHandSides(), ElementsAre(-2, -4));
    EXPECT_THAT(ExtractRightHandSides(), ElementsAre(2, 4));
    EXPECT_THAT(ExtractLowerBounds(), ElementsAre(-1, -2, -3, -4));
    EXPECT_THAT(ExtractUpperBounds(), ElementsAre(1, 2, 3, 4));
    EXPECT_EQ(lpi_->GetNumberOfNonZeros(), 8);
    EXPECT_THAT(ExtractMatrix(),
                ElementsAre(ElementsAre(2, 4), ElementsAre(7, 9),
                            ElementsAre(12, 14), ElementsAre(17, 19)));
  }

  {
    ASSERT_OK(InitSimpleProblem(4, 5));
    ASSERT_OK(lpi_->SolveLPWithPrimalSimplex());
    EXPECT_TRUE(lpi_->IsSolved());
    ASSERT_OK_AND_ASSIGN(
        (const absl::StrongVector<RowIndex, RowIndex> row_mapping),
        lpi_->DeleteRowSet({false, true, true, true, false}));
    EXPECT_THAT(row_mapping, ElementsAre(RowIndex(0), kInvalidRow, kInvalidRow,
                                         kInvalidRow, RowIndex(1)));
    EXPECT_FALSE(lpi_->IsSolved());
    EXPECT_EQ(lpi_->GetNumberOfColumns(), 4);
    EXPECT_EQ(lpi_->GetNumberOfRows(), 2);
    EXPECT_THAT(ExtractLeftHandSides(), ElementsAre(-1, -5));
    EXPECT_THAT(ExtractRightHandSides(), ElementsAre(1, 5));
    EXPECT_THAT(ExtractLowerBounds(), ElementsAre(-1, -2, -3, -4));
    EXPECT_THAT(ExtractUpperBounds(), ElementsAre(1, 2, 3, 4));
    EXPECT_EQ(lpi_->GetNumberOfNonZeros(), 8);
    EXPECT_THAT(ExtractMatrix(),
                ElementsAre(ElementsAre(1, 5), ElementsAre(6, 10),
                            ElementsAre(11, 15), ElementsAre(16, 20)));
  }

  {
    ASSERT_OK(InitSimpleProblem(4, 5));
    ASSERT_OK(lpi_->SolveLPWithPrimalSimplex());
    EXPECT_TRUE(lpi_->IsSolved());
    ASSERT_OK_AND_ASSIGN(
        (const absl::StrongVector<RowIndex, RowIndex> row_mapping),
        lpi_->DeleteRowSet({true, true, true, true, true}));
    EXPECT_THAT(row_mapping, ElementsAre(kInvalidRow, kInvalidRow, kInvalidRow,
                                         kInvalidRow, kInvalidRow));
    EXPECT_FALSE(lpi_->IsSolved());
    EXPECT_EQ(lpi_->GetNumberOfColumns(), 4);
    EXPECT_EQ(lpi_->GetNumberOfRows(), 0);
    EXPECT_THAT(ExtractLeftHandSides(), IsEmpty());
    EXPECT_THAT(ExtractRightHandSides(), IsEmpty());
    EXPECT_THAT(ExtractLowerBounds(), ElementsAre(-1, -2, -3, -4));
    EXPECT_THAT(ExtractUpperBounds(), ElementsAre(1, 2, 3, 4));
    EXPECT_EQ(lpi_->GetNumberOfNonZeros(), 0);
    EXPECT_THAT(ExtractMatrix(),
                ElementsAre(IsEmpty(), IsEmpty(), IsEmpty(), IsEmpty()));
  }
}

TEST_P(ModelConstructionTest, DeleteSingleColumn) {
  {
    ASSERT_OK(InitSimpleProblem(4, 5));
    ASSERT_OK(lpi_->SolveLPWithPrimalSimplex());
    EXPECT_TRUE(lpi_->IsSolved());
    ASSERT_OK(lpi_->DeleteColumns(ColIndex(0), ColIndex(0)));
    EXPECT_FALSE(lpi_->IsSolved());
    EXPECT_EQ(lpi_->GetNumberOfColumns(), 3);
    EXPECT_EQ(lpi_->GetNumberOfRows(), 5);
    EXPECT_THAT(ExtractLeftHandSides(), ElementsAre(-1, -2, -3, -4, -5));
    EXPECT_THAT(ExtractRightHandSides(), ElementsAre(1, 2, 3, 4, 5));
    EXPECT_THAT(ExtractLowerBounds(), ElementsAre(-2, -3, -4));
    EXPECT_THAT(ExtractUpperBounds(), ElementsAre(2, 3, 4));
    EXPECT_EQ(lpi_->GetNumberOfNonZeros(), 15);
    EXPECT_THAT(ExtractMatrix(), ElementsAre(ElementsAre(6, 7, 8, 9, 10),
                                             ElementsAre(11, 12, 13, 14, 15),
                                             ElementsAre(16, 17, 18, 19, 20)));
  }

  {
    ASSERT_OK(InitSimpleProblem(4, 5));
    ASSERT_OK(lpi_->SolveLPWithPrimalSimplex());
    EXPECT_TRUE(lpi_->IsSolved());
    ASSERT_OK(lpi_->DeleteColumns(ColIndex(1), ColIndex(1)));
    EXPECT_FALSE(lpi_->IsSolved());
    EXPECT_EQ(lpi_->GetNumberOfColumns(), 3);
    EXPECT_EQ(lpi_->GetNumberOfRows(), 5);
    EXPECT_THAT(ExtractLeftHandSides(), ElementsAre(-1, -2, -3, -4, -5));
    EXPECT_THAT(ExtractRightHandSides(), ElementsAre(1, 2, 3, 4, 5));
    EXPECT_THAT(ExtractLowerBounds(), ElementsAre(-1, -3, -4));
    EXPECT_THAT(ExtractUpperBounds(), ElementsAre(1, 3, 4));
    EXPECT_EQ(lpi_->GetNumberOfNonZeros(), 15);
    EXPECT_THAT(ExtractMatrix(), ElementsAre(ElementsAre(1, 2, 3, 4, 5),
                                             ElementsAre(11, 12, 13, 14, 15),
                                             ElementsAre(16, 17, 18, 19, 20)));
  }

  {
    ASSERT_OK(InitSimpleProblem(4, 5));
    ASSERT_OK(lpi_->SolveLPWithPrimalSimplex());
    EXPECT_TRUE(lpi_->IsSolved());
    ASSERT_OK(lpi_->DeleteColumns(ColIndex(3), ColIndex(3)));
    EXPECT_FALSE(lpi_->IsSolved());
    EXPECT_EQ(lpi_->GetNumberOfColumns(), 3);
    EXPECT_EQ(lpi_->GetNumberOfRows(), 5);
    EXPECT_THAT(ExtractLeftHandSides(), ElementsAre(-1, -2, -3, -4, -5));
    EXPECT_THAT(ExtractRightHandSides(), ElementsAre(1, 2, 3, 4, 5));
    EXPECT_THAT(ExtractLowerBounds(), ElementsAre(-1, -2, -3));
    EXPECT_THAT(ExtractUpperBounds(), ElementsAre(1, 2, 3));
    EXPECT_EQ(lpi_->GetNumberOfNonZeros(), 15);
    EXPECT_THAT(ExtractMatrix(), ElementsAre(ElementsAre(1, 2, 3, 4, 5),
                                             ElementsAre(6, 7, 8, 9, 10),
                                             ElementsAre(11, 12, 13, 14, 15)));
  }
}

TEST_P(ModelConstructionTest, DeleteMultipleColumns) {
  {
    ASSERT_OK(InitSimpleProblem(4, 5));
    ASSERT_OK(lpi_->SolveLPWithPrimalSimplex());
    EXPECT_TRUE(lpi_->IsSolved());
    ASSERT_OK(lpi_->DeleteColumns(ColIndex(1), ColIndex(2)));
    EXPECT_FALSE(lpi_->IsSolved());
    EXPECT_EQ(lpi_->GetNumberOfColumns(), 2);
    EXPECT_EQ(lpi_->GetNumberOfRows(), 5);
    EXPECT_THAT(ExtractLeftHandSides(), ElementsAre(-1, -2, -3, -4, -5));
    EXPECT_THAT(ExtractRightHandSides(), ElementsAre(1, 2, 3, 4, 5));
    EXPECT_THAT(ExtractLowerBounds(), ElementsAre(-1, -4));
    EXPECT_THAT(ExtractUpperBounds(), ElementsAre(1, 4));
    EXPECT_EQ(lpi_->GetNumberOfNonZeros(), 10);
    EXPECT_THAT(ExtractMatrix(), ElementsAre(ElementsAre(1, 2, 3, 4, 5),
                                             ElementsAre(16, 17, 18, 19, 20)));
  }

  {
    ASSERT_OK(InitSimpleProblem(4, 5));
    ASSERT_OK(lpi_->SolveLPWithPrimalSimplex());
    EXPECT_TRUE(lpi_->IsSolved());
    ASSERT_OK(lpi_->DeleteColumns(ColIndex(0), ColIndex(3)));
    EXPECT_FALSE(lpi_->IsSolved());
    EXPECT_EQ(lpi_->GetNumberOfColumns(), 0);
    EXPECT_EQ(lpi_->GetNumberOfRows(), 5);
    EXPECT_THAT(ExtractLeftHandSides(), ElementsAre(-1, -2, -3, -4, -5));
    EXPECT_THAT(ExtractRightHandSides(), ElementsAre(1, 2, 3, 4, 5));
    EXPECT_THAT(ExtractLowerBounds(), IsEmpty());
    EXPECT_THAT(ExtractUpperBounds(), IsEmpty());
    EXPECT_EQ(lpi_->GetNumberOfNonZeros(), 0);
    EXPECT_THAT(ExtractMatrix(), IsEmpty());
  }
}

TEST_P(ModelConstructionTest, SetObjectiveCoefficients) {
  for (double c1 : {-1.0, 0.0, 1.0}) {
    for (double c2 : {-1.0, 0.0, 1.0}) {
      SCOPED_TRACE(absl::StrFormat("c1=%lf, c2=%lf", c1, c2));
      ASSERT_OK(InitEmptyProblem(2, 2));
      ASSERT_OK(lpi_->SolveLPWithPrimalSimplex());

      EXPECT_TRUE(lpi_->IsSolved());
      EXPECT_OK(lpi_->SetObjectiveCoefficient(ColIndex(0), c1));
      EXPECT_OK(lpi_->SetObjectiveCoefficient(ColIndex(1), c2));
      EXPECT_FALSE(lpi_->IsSolved());
      EXPECT_EQ(lpi_->GetObjectiveCoefficient(ColIndex(0)), c1);
      EXPECT_EQ(lpi_->GetObjectiveCoefficient(ColIndex(1)), c2);
    }
  }
}

TEST_P(ModelConstructionTest, SetColumnBounds) {
  const std::vector<double> lbs{-lpi_->Infinity(), -1.0, 0.0, 1.0};
  const std::vector<double> ubs{5.0, 6.0, 7.0, lpi_->Infinity()};
  for (double lb1 : lbs) {
    for (double lb2 : lbs) {
      for (double ub1 : ubs) {
        for (double ub2 : ubs) {
          SCOPED_TRACE(absl::StrFormat("lb1=%lf, lb2=%lf, ub1=%lf, ub2=%lf",
                                       lb1, lb2, ub1, ub2));
          ASSERT_OK(InitEmptyProblem(2, 2));
          ASSERT_OK(lpi_->SolveLPWithPrimalSimplex());

          EXPECT_TRUE(lpi_->IsSolved());
          ASSERT_OK(lpi_->SetColumnBounds(ColIndex(0), lb1, ub1));
          ASSERT_OK(lpi_->SetColumnBounds(ColIndex(1), lb2, ub2));
          EXPECT_EQ(lpi_->GetLowerBound(ColIndex(0)), lb1);
          EXPECT_EQ(lpi_->GetUpperBound(ColIndex(0)), ub1);
          EXPECT_EQ(lpi_->GetLowerBound(ColIndex(1)), lb2);
          EXPECT_EQ(lpi_->GetUpperBound(ColIndex(1)), ub2);
          EXPECT_FALSE(lpi_->IsSolved());
        }
      }
    }
  }
}

TEST_P(ModelConstructionTest, SetRowSides) {
  const std::vector<double> left{-lpi_->Infinity(), -1.0, 0.0, 1.0};
  const std::vector<double> right{5.0, 6.0, 7.0, lpi_->Infinity()};
  for (double l1 : left) {
    for (double l2 : left) {
      for (double r1 : right) {
        for (double r2 : right) {
          SCOPED_TRACE(absl::StrFormat("l1=%lf, l2=%lf, r1=%lf, r2=%lf", l1, l2,
                                       r1, r2));
          ASSERT_OK(InitEmptyProblem(2, 2));
          ASSERT_OK(lpi_->SolveLPWithPrimalSimplex());

          EXPECT_TRUE(lpi_->IsSolved());
          ASSERT_OK(lpi_->SetRowSides(RowIndex(0), l1, r1));
          ASSERT_OK(lpi_->SetRowSides(RowIndex(1), l2, r2));
          EXPECT_EQ(lpi_->GetLeftHandSide(RowIndex(0)), l1);
          EXPECT_EQ(lpi_->GetRightHandSide(RowIndex(0)), r1);
          EXPECT_EQ(lpi_->GetLeftHandSide(RowIndex(1)), l2);
          EXPECT_EQ(lpi_->GetRightHandSide(RowIndex(1)), r2);
          EXPECT_FALSE(lpi_->IsSolved());
        }
      }
    }
  }
}

TEST_P(ModelConstructionTest, ConstructFromMipData) {
  // This should build the following model
  //
  // obj: max x1 + 2*x2 -4*x3 + 1e-3x4 + 0*x5
  // constraints:
  //                  x1 +  5*x2                          <= -1     (ct1)
  // -1       <=                              -x4         <= -3e-10 (ct2)
  // -3e-10   <=            2*x2 + 3e5*x3                 <= 0      (ct3)
  // 0        <=              x2 +  20*x3         + 10*x5 <= 1      (ct4)
  // 1        <= -1.9*x1                                  <= 3e10   (ct5)
  // 3e10     <=                          1e-2*x4                   (ct6)
  // bounds:
  // -1 <= x1 <= 10
  //       x2
  //  0 <= x3
  //       x4 <= 29.3
  //  1 <= x5 <= 1
  MiniMipProblem problem;
  problem.is_maximization = true;
  problem.objective_offset = 0.0;
  problem.variables.push_back(MiniMipVariable{.name = "x1",
                                              .objective_coefficient = 1,
                                              .lower_bound = -1,
                                              .upper_bound = 10,
                                              .is_integer = false});
  problem.variables.push_back(MiniMipVariable{.name = "x2",
                                              .objective_coefficient = 2,
                                              .lower_bound = -solver_inf_,
                                              .upper_bound = solver_inf_,
                                              .is_integer = false});
  problem.variables.push_back(MiniMipVariable{.name = "x3",
                                              .objective_coefficient = -4,
                                              .lower_bound = 0,
                                              .upper_bound = solver_inf_,
                                              .is_integer = true});
  problem.variables.push_back(MiniMipVariable{.name = "x4",
                                              .objective_coefficient = 1e-3,
                                              .lower_bound = -solver_inf_,
                                              .upper_bound = 29.0,
                                              .is_integer = true});
  problem.variables.push_back(MiniMipVariable{.name = "x5",
                                              .objective_coefficient = 0,
                                              .lower_bound = 1,
                                              .upper_bound = 1,
                                              .is_integer = false});
  problem.constraints.push_back(
      MiniMipConstraint{.name = "ct1",
                        .var_indices = {0, 1},
                        .coefficients = {1, 5},
                        .left_hand_side = -solver_inf_,
                        .right_hand_side = -1});
  problem.constraints.push_back(MiniMipConstraint{.name = "ct2",
                                                  .var_indices = {3},
                                                  .coefficients = {-1},
                                                  .left_hand_side = -1,
                                                  .right_hand_side = -3e-10});
  problem.constraints.push_back(MiniMipConstraint{.name = "ct3",
                                                  .var_indices = {1, 2},
                                                  .coefficients = {2, 3e5},
                                                  .left_hand_side = -3e-10,
                                                  .right_hand_side = 0});
  problem.constraints.push_back(MiniMipConstraint{.name = "ct4",
                                                  .var_indices = {1, 2, 4},
                                                  .coefficients = {1, 20, 10},
                                                  .left_hand_side = 0,
                                                  .right_hand_side = 1});
  problem.constraints.push_back(MiniMipConstraint{.name = "ct5",
                                                  .var_indices = {0},
                                                  .coefficients = {-1.9},
                                                  .left_hand_side = 1,
                                                  .right_hand_side = 3e10});
  problem.constraints.push_back(
      MiniMipConstraint{.name = "ct6",
                        .var_indices = {3},
                        .coefficients = {1e-2},
                        .left_hand_side = 3e10,
                        .right_hand_side = solver_inf_});
  ASSERT_OK(lpi_->PopulateFromMipData(MipData(problem)));

  // Check variable properties
  EXPECT_EQ(lpi_->GetNumberOfColumns(), 5);
  EXPECT_THAT(ExtractLowerBounds(),
              ElementsAre(-1, -solver_inf_, 0, -solver_inf_, 1));
  EXPECT_THAT(ExtractUpperBounds(),
              ElementsAre(10, solver_inf_, solver_inf_, 29, 1));
  EXPECT_THAT(ExtractObjectiveCoefficients(), ElementsAre(1, 2, -4, 1e-3, 0));

  // Check row properties
  EXPECT_EQ(lpi_->GetNumberOfRows(), 6);
  EXPECT_THAT(ExtractLeftHandSides(),
              ElementsAre(-solver_inf_, -1, -3e-10, 0, 1, 3e10));
  EXPECT_THAT(ExtractRightHandSides(),
              ElementsAre(-1, -3e-10, 0, 1, 3e10, solver_inf_));

  // Check constraint matrix
  EXPECT_EQ(lpi_->GetNumberOfNonZeros(), 10);
  const absl::StrongVector<ColIndex, absl::StrongVector<RowIndex, double>>
      extracted_matrix = ExtractMatrix();
  EXPECT_THAT(extracted_matrix, ElementsAre(ElementsAre(1, 0, 0, 0, -1.9, 0),
                                            ElementsAre(5, 0, 2, 1, 0, 0),
                                            ElementsAre(0, 0, 3e5, 20, 0, 0),
                                            ElementsAre(0, -1, 0, 0, 0, 1e-2),
                                            ElementsAre(0, 0, 0, 10, 0, 0)));

  // Check that the different ways to retrieve matrix coefficients are
  // consistent.
  for (ColIndex col(0); col < 5; ++col) {
    for (RowIndex row(0); row < 6; ++row) {
      SCOPED_TRACE(absl::StrFormat("row=%d, col=%d", row.value(), col.value()));
      EXPECT_EQ(lpi_->GetSparseRowCoefficients(row).value(col),
                extracted_matrix[col][row]);
      EXPECT_EQ(lpi_->GetSparseColumnCoefficients(col).value(row),
                extracted_matrix[col][row]);
      EXPECT_EQ(lpi_->GetMatrixCoefficient(col, row),
                extracted_matrix[col][row]);
    }
  }
}

class FileTest : public LPInterfaceImplementationTest {
 protected:
  static absl::StatusOr<std::string> GetTemporaryFile() {
    const char* name = std::tmpnam(nullptr);
    if (name != nullptr) return name;
    return absl::InternalError("Failed to create a temporary file name");
  }

  static bool FileExists(const std::string& file_name) {
    FILE* file = std::fopen(file_name.c_str(), "r");
    if (file == nullptr) return false;
    std::fclose(file);
    return true;
  }
};

INSTANTIATE_TEST_SUITE_P(All, FileTest,
                         testing::ValuesIn({LPParameters::LP_GLOP,
                                            LPParameters::LP_SOPLEX}));

TEST_P(FileTest, WritesLPToFile) {
  ASSERT_OK(InitSimpleProblem(3, 4));
  ASSERT_OK_AND_ASSIGN(const std::string file_path, GetTemporaryFile());
  ASSERT_OK_AND_ASSIGN(const std::string actual_path,
                       lpi_->WriteLPToFile(file_path));

  EXPECT_THAT(actual_path, StartsWith(file_path));
  if (actual_path != file_path) {
    ASSERT_GT(actual_path.size(), file_path.size());
    EXPECT_EQ(actual_path[file_path.size()], '.');
    EXPECT_THAT(actual_path.substr(file_path.size()),
                testing::Not(testing::Contains(testing::AnyOf('\\', '/'))));
  }
  EXPECT_TRUE(FileExists(actual_path));
}

TEST_P(FileTest, LoadLPFromFile) {
  // max z = 2 x1 +    x2
  //           x1 +    x2 <= 1
  //           x1         <= 3/4
  //         8 x1 + 20 x2 <= 10
  //           x1,     x2 >= 0
  // with optimal solution (0.75, 0.2) and objective value 1.7
  ASSERT_OK(UploadProblem(
      /*is_maximization=*/true,
      /*lower_bounds=*/{0, 0},
      /*upper_bounds=*/{solver_inf_, solver_inf_},
      /*objective_coefficients=*/{2, 1},
      /*left_hand_sides=*/{-solver_inf_, -solver_inf_, -solver_inf_},
      /*right_hand_sides=*/{1, 0.75, 10},
      /*matrix=*/{{1, 1, 8}, {1, 0, 20}}));
  ASSERT_OK(lpi_->SolveLPWithPrimalSimplex());
  ASSERT_THAT(lpi_->GetObjectiveValue(), DoubleEq(1.7));
  {
    ASSERT_OK_AND_ASSIGN(
        (const absl::StrongVector<ColIndex, double> primal_values),
        lpi_->GetPrimalValues());
    ASSERT_THAT(primal_values, Pointwise(DoubleEq(), {0.75, 0.2}));
  }

  ASSERT_OK_AND_ASSIGN(const std::string file_path, GetTemporaryFile());
  ASSERT_OK_AND_ASSIGN(const std::string actual_path,
                       lpi_->WriteLPToFile(file_path));
  ASSERT_TRUE(FileExists(actual_path));
  ASSERT_OK(lpi_->Clear());
  ASSERT_OK(lpi_->ReadLPFromFile(actual_path));

  // The loaded model is allowed to differ from original one, as long as it
  // doesn't change the solution space. We therefore test by checking the solve
  // results, rather than checking model properties directly.
  ASSERT_OK(lpi_->SolveLPWithPrimalSimplex());
  EXPECT_THAT(lpi_->GetObjectiveValue(), DoubleEq(1.7));
  {
    ASSERT_OK_AND_ASSIGN(
        (const absl::StrongVector<ColIndex, double> primal_values),
        lpi_->GetPrimalValues());
    EXPECT_THAT(primal_values, Pointwise(DoubleEq(), {0.75, 0.2}));
  }

  ASSERT_EQ(std::remove(actual_path.c_str()), 0);
}

class SolveTest : public LPInterfaceImplementationTest {
 protected:
  // The row order of the B matrix depends on the order of the variables in the
  // basis. To make sure that it is correct, we must therefore ensure a
  // deterministic order. Specifically, we first place all rows corresponding to
  // original variables in increasing index order, then rows coresponding to
  // slack variables in increasing index order.
  absl::StrongVector<RowIndex, absl::StrongVector<ColIndex, double>>
  SortMatrixByBasis(
      absl::StrongVector<RowIndex, absl::StrongVector<ColIndex, double>>
          matrix) {
    std::vector<ColOrRowIndex> cors = lpi_->GetColumnsAndRowsInBasis();
    CHECK_EQ(matrix.size(), cors.size());
    std::vector<int> indices(cors.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&cors](int i1, int i2) {
      return std::make_tuple(cors[i1].row(), cors[i1].col()) <
             std::make_tuple(cors[i2].row(), cors[i2].col());
    });
    absl::StrongVector<RowIndex, absl::StrongVector<ColIndex, double>> res;
    res.reserve(cors.size());
    for (int i : indices) res.emplace_back(std::move(matrix[RowIndex(i)]));
    return res;
  }

  absl::StatusOr<
      absl::StrongVector<RowIndex, absl::StrongVector<ColIndex, double>>>
  ExtractBInvByColumn() {
    const int dim = lpi_->GetNumberOfRows().value();
    absl::StrongVector<RowIndex, absl::StrongVector<ColIndex, double>> binv(
        dim, absl::StrongVector<ColIndex, double>(dim, 0.0));
    for (ColIndex col_idx(0); col_idx < dim; ++col_idx) {
      ASSIGN_OR_RETURN(const SparseCol col,
                       lpi_->GetSparseColumnOfBInverted(col_idx));
      for (RowIndex row_idx(0); row_idx < dim; ++row_idx) {
        binv[row_idx][col_idx] = col.value(row_idx);
      }
    }
    return SortMatrixByBasis(binv);
  }

  absl::StatusOr<
      absl::StrongVector<RowIndex, absl::StrongVector<ColIndex, double>>>
  ExtractBInvByRow() {
    const int dim = lpi_->GetNumberOfRows().value();
    absl::StrongVector<RowIndex, absl::StrongVector<ColIndex, double>> binv(
        dim, absl::StrongVector<ColIndex, double>(dim, 0.0));
    for (RowIndex row_idx(0); row_idx < dim; ++row_idx) {
      ASSIGN_OR_RETURN(const SparseRow row,
                       lpi_->GetSparseRowOfBInverted(row_idx));
      for (ColIndex col_idx(0); col_idx < dim; ++col_idx) {
        binv[row_idx][col_idx] = row.value(col_idx);
      }
    }
    return SortMatrixByBasis(binv);
  }

  absl::StatusOr<
      absl::StrongVector<RowIndex, absl::StrongVector<ColIndex, double>>>
  ExtractBInvTimesAByColumn() {
    const int n_rows = lpi_->GetNumberOfRows().value();
    const int n_cols = lpi_->GetNumberOfColumns().value();
    absl::StrongVector<RowIndex, absl::StrongVector<ColIndex, double>> binv(
        n_rows, absl::StrongVector<ColIndex, double>(n_cols, 0.0));
    for (ColIndex col_idx(0); col_idx < n_cols; ++col_idx) {
      ASSIGN_OR_RETURN(const SparseCol col,
                       lpi_->GetSparseColumnOfBInvertedTimesA(col_idx));
      for (RowIndex row_idx(0); row_idx < n_rows; ++row_idx) {
        binv[row_idx][col_idx] = col.value(row_idx);
      }
    }
    return SortMatrixByBasis(binv);
  }

  absl::StatusOr<
      absl::StrongVector<RowIndex, absl::StrongVector<ColIndex, double>>>
  ExtractBInvTimesAByRow() {
    const int n_rows = lpi_->GetNumberOfRows().value();
    const int n_cols = lpi_->GetNumberOfColumns().value();
    absl::StrongVector<RowIndex, absl::StrongVector<ColIndex, double>> binv(
        n_rows, absl::StrongVector<ColIndex, double>(n_cols, 0.0));
    for (RowIndex row_idx(0); row_idx < n_rows; ++row_idx) {
      ASSIGN_OR_RETURN(const SparseRow row,
                       lpi_->GetSparseRowOfBInvertedTimesA(row_idx));
      for (ColIndex col_idx(0); col_idx < n_cols; ++col_idx) {
        binv[row_idx][col_idx] = row.value(col_idx);
      }
    }
    return SortMatrixByBasis(binv);
  }
};

INSTANTIATE_TEST_SUITE_P(All, SolveTest,
                         testing::ValuesIn({LPParameters::LP_GLOP,
                                            LPParameters::LP_SOPLEX}));

TEST_P(SolveTest, PrimalDualFeasible1) {
  // max 3 x1 +   x2
  //    2 x1 +   x2 <= 10
  //      x1 + 3 x2 <= 15
  //      x1,    x2 >= 0
  for (bool use_primal_simplex : {true, false}) {
    SCOPED_TRACE(absl::StrCat("Using ", use_primal_simplex ? "primal" : "dual",
                              " simplex."));
    ASSERT_OK(UploadProblem(
        /*is_maximization=*/true,
        /*lower_bounds=*/{0, 0},
        /*upper_bounds=*/{solver_inf_, solver_inf_},
        /*objective_coefficients=*/{3, 1},
        /*left_hand_sides=*/{-solver_inf_, -solver_inf_},
        /*right_hand_sides=*/{10, 15},
        /*matrix=*/{{2, 1}, {1, 3}}));
    if (use_primal_simplex) {
      ASSERT_OK(lpi_->SolveLPWithPrimalSimplex());
    } else {
      ASSERT_OK(lpi_->SolveLPWithDualSimplex());
    }

    // Check solve status.
    ASSERT_TRUE(lpi_->IsSolved());
    EXPECT_FALSE(lpi_->ObjectiveLimitIsExceeded());
    EXPECT_FALSE(lpi_->IterationLimitIsExceeded());
    EXPECT_FALSE(lpi_->TimeLimitIsExceeded());
    EXPECT_TRUE(lpi_->IsOptimal());
    EXPECT_TRUE(lpi_->IsPrimalFeasible());
    EXPECT_FALSE(lpi_->IsPrimalInfeasible());
    EXPECT_FALSE(lpi_->IsPrimalUnbounded());
    EXPECT_TRUE(lpi_->IsDualFeasible());
    EXPECT_FALSE(lpi_->IsDualInfeasible());
    EXPECT_FALSE(lpi_->IsDualUnbounded());

    // Check solution values
    EXPECT_THAT(lpi_->GetObjectiveValue(), DoubleEq(15.0));
    ASSERT_OK_AND_ASSIGN(
        (const absl::StrongVector<ColIndex, double> primal_solution),
        lpi_->GetPrimalValues());
    EXPECT_THAT(primal_solution, Pointwise(DoubleEq(), {5, 0}));
    ASSERT_OK_AND_ASSIGN(
        (const absl::StrongVector<RowIndex, double> dual_solution),
        lpi_->GetDualValues());
    EXPECT_THAT(dual_solution, Pointwise(DoubleEq(), {1.5, 0.0}));
    ASSERT_OK_AND_ASSIGN(
        (const absl::StrongVector<ColIndex, double> reduced_costs),
        lpi_->GetReducedCosts());
    EXPECT_THAT(reduced_costs, Pointwise(DoubleEq(), {0.0, -0.5}));
    ASSERT_OK_AND_ASSIGN(
        (const absl::StrongVector<RowIndex, double> activities),
        lpi_->GetRowActivities());
    EXPECT_THAT(activities, Pointwise(DoubleEq(), {10, 5}));

    // Check basis status
    ASSERT_OK_AND_ASSIGN(
        (const absl::StrongVector<ColIndex, LPBasisStatus> column_basis_status),
        lpi_->GetBasisStatusForColumns());
    EXPECT_THAT(column_basis_status, ElementsAre(LPBasisStatus::kBasic,
                                                 LPBasisStatus::kAtLowerBound));
    ASSERT_OK_AND_ASSIGN(
        (const absl::StrongVector<RowIndex, LPBasisStatus> row_basis_status),
        lpi_->GetBasisStatusForRows());
    EXPECT_THAT(row_basis_status, ElementsAre(LPBasisStatus::kAtUpperBound,
                                              LPBasisStatus::kBasic));
    ASSERT_THAT(lpi_->GetColumnsAndRowsInBasis(),
                UnorderedElementsAre(ColOrRowIndex(RowIndex(3)),
                                     ColOrRowIndex(ColIndex(0))));

    // Check B inverse
    {
      ASSERT_OK_AND_ASSIGN(
          (const absl::StrongVector<RowIndex,
                                    absl::StrongVector<ColIndex, double>>
               b_inverse),
          ExtractBInvByColumn());
      EXPECT_THAT(b_inverse, ElementsAre(Pointwise(DoubleEq(), {0.5, 0.0}),
                                         Pointwise(DoubleEq(), {-0.5, 1.0})));
    }
    {
      ASSERT_OK_AND_ASSIGN(
          (const absl::StrongVector<RowIndex,
                                    absl::StrongVector<ColIndex, double>>
               b_inverse),
          ExtractBInvByRow());
      EXPECT_THAT(b_inverse, ElementsAre(Pointwise(DoubleEq(), {0.5, 0.0}),
                                         Pointwise(DoubleEq(), {-0.5, 1.0})));
    }
    {
      ASSERT_OK_AND_ASSIGN(
          (const absl::StrongVector<RowIndex,
                                    absl::StrongVector<ColIndex, double>>
               b_inverse_times_a),
          ExtractBInvTimesAByColumn());
      EXPECT_THAT(b_inverse_times_a,
                  ElementsAre(Pointwise(DoubleEq(), {1.0, 0.5}),
                              Pointwise(DoubleEq(), {0.0, 2.5})));
    }
    {
      ASSERT_OK_AND_ASSIGN(
          (const absl::StrongVector<RowIndex,
                                    absl::StrongVector<ColIndex, double>>
               b_inverse_times_a),
          ExtractBInvTimesAByRow());
      EXPECT_THAT(b_inverse_times_a,
                  ElementsAre(Pointwise(DoubleEq(), {1.0, 0.5}),
                              Pointwise(DoubleEq(), {0.0, 2.5})));
    }

    // Check rays
    ASSERT_FALSE(lpi_->ExistsPrimalRay());
    ASSERT_FALSE(lpi_->HasPrimalRay());
    ASSERT_FALSE(lpi_->ExistsDualRay());
    ASSERT_FALSE(lpi_->HasDualRay());
  }
}

TEST_P(SolveTest, PrimalDualFeasible2) {
  // The spirit of this test is the same as the previous one, but it triggered
  // bugs in the soplex API so is kept.
  //
  // max z = 2 x1 +    x2
  //           x1 +    x2 <= 1
  //           x1         <= 3/4
  //         8 x1 + 20 x2 <= 10
  //           x1,     x2 >= 0
  for (bool use_primal_simplex : {true, false}) {
    SCOPED_TRACE(absl::StrCat("Using ", use_primal_simplex ? "primal" : "dual",
                              " simplex."));
    ASSERT_OK(UploadProblem(
        /*is_maximization=*/true,
        /*lower_bounds=*/{0, 0},
        /*upper_bounds=*/{solver_inf_, solver_inf_},
        /*objective_coefficients=*/{2, 1},
        /*left_hand_sides=*/{-solver_inf_, -solver_inf_, -solver_inf_},
        /*right_hand_sides=*/{1, 0.75, 10},
        /*matrix=*/{{1, 1, 8}, {1, 0, 20}}));
    if (use_primal_simplex) {
      ASSERT_OK(lpi_->SolveLPWithPrimalSimplex());
    } else {
      ASSERT_OK(lpi_->SolveLPWithDualSimplex());
    }

    // Check solve status.
    ASSERT_TRUE(lpi_->IsSolved());
    EXPECT_FALSE(lpi_->ObjectiveLimitIsExceeded());
    EXPECT_FALSE(lpi_->IterationLimitIsExceeded());
    EXPECT_FALSE(lpi_->TimeLimitIsExceeded());
    EXPECT_TRUE(lpi_->IsOptimal());
    EXPECT_TRUE(lpi_->IsPrimalFeasible());
    EXPECT_FALSE(lpi_->IsPrimalInfeasible());
    EXPECT_FALSE(lpi_->IsPrimalUnbounded());
    EXPECT_TRUE(lpi_->IsDualFeasible());
    EXPECT_FALSE(lpi_->IsDualInfeasible());
    EXPECT_FALSE(lpi_->IsDualUnbounded());

    // Check solution values
    EXPECT_THAT(lpi_->GetObjectiveValue(), DoubleEq(1.7));
    ASSERT_OK_AND_ASSIGN(
        (const absl::StrongVector<ColIndex, double> primal_solution),
        lpi_->GetPrimalValues());
    EXPECT_THAT(primal_solution, Pointwise(DoubleEq(), {0.75, 0.2}));
    ASSERT_OK_AND_ASSIGN(
        (const absl::StrongVector<RowIndex, double> dual_solution),
        lpi_->GetDualValues());
    EXPECT_THAT(dual_solution, Pointwise(DoubleEq(), {0.0, 1.6, 0.05}));
    ASSERT_OK_AND_ASSIGN(
        (const absl::StrongVector<ColIndex, double> reduced_costs),
        lpi_->GetReducedCosts());
    EXPECT_THAT(reduced_costs, Pointwise(DoubleEq(), {0.0, 0.0}));
    ASSERT_OK_AND_ASSIGN(
        (const absl::StrongVector<RowIndex, double> activities),
        lpi_->GetRowActivities());
    EXPECT_THAT(activities, Pointwise(DoubleEq(), {0.95, 0.75, 10.0}));

    // Check basis status
    ASSERT_OK_AND_ASSIGN(
        (const absl::StrongVector<ColIndex, LPBasisStatus> column_basis_status),
        lpi_->GetBasisStatusForColumns());
    EXPECT_THAT(column_basis_status,
                ElementsAre(LPBasisStatus::kBasic, LPBasisStatus::kBasic));
    ASSERT_OK_AND_ASSIGN(
        (const absl::StrongVector<RowIndex, LPBasisStatus> row_basis_status),
        lpi_->GetBasisStatusForRows());
    EXPECT_THAT(row_basis_status,
                ElementsAre(LPBasisStatus::kBasic, LPBasisStatus::kAtUpperBound,
                            LPBasisStatus::kAtUpperBound));
    ASSERT_THAT(lpi_->GetColumnsAndRowsInBasis(),
                UnorderedElementsAre(ColOrRowIndex(RowIndex(2)),
                                     ColOrRowIndex(ColIndex(0)),
                                     ColOrRowIndex(ColIndex(1))));

    // Check B inverse
    {
      ASSERT_OK_AND_ASSIGN(
          (const absl::StrongVector<RowIndex,
                                    absl::StrongVector<ColIndex, double>>
               b_inverse),
          ExtractBInvByColumn());
      EXPECT_THAT(b_inverse,
                  ElementsAre(Pointwise(DoubleEq(), {0.0, 1.0, 0.0}),
                              Pointwise(DoubleEq(), {0.0, -0.4, 0.05}),
                              Pointwise(DoubleEq(), {1.0, -0.6, -0.05})));
    }
    {
      ASSERT_OK_AND_ASSIGN(
          (const absl::StrongVector<RowIndex,
                                    absl::StrongVector<ColIndex, double>>
               b_inverse),
          ExtractBInvByRow());
      EXPECT_THAT(b_inverse,
                  ElementsAre(Pointwise(DoubleEq(), {0.0, 1.0, 0.0}),
                              Pointwise(DoubleEq(), {0.0, -0.4, 0.05}),
                              Pointwise(DoubleEq(), {1.0, -0.6, -0.05})));
    }
    {
      ASSERT_OK_AND_ASSIGN(
          (const absl::StrongVector<RowIndex,
                                    absl::StrongVector<ColIndex, double>>
               b_inverse_times_a),
          ExtractBInvTimesAByColumn());
      EXPECT_THAT(b_inverse_times_a,
                  ElementsAre(Pointwise(DoubleEq(), {1.0, 0.0}),
                              Pointwise(DoubleEq(), {0.0, 1.0}),
                              Pointwise(DoubleEq(), {0.0, 0.0})));
    }
    {
      ASSERT_OK_AND_ASSIGN(
          (const absl::StrongVector<RowIndex,
                                    absl::StrongVector<ColIndex, double>>
               b_inverse_times_a),
          ExtractBInvTimesAByRow());
      EXPECT_THAT(b_inverse_times_a,
                  ElementsAre(Pointwise(DoubleEq(), {1.0, 0.0}),
                              Pointwise(DoubleEq(), {0.0, 1.0}),
                              Pointwise(DoubleEq(), {0.0, 0.0})));
    }

    // Check rays
    ASSERT_FALSE(lpi_->ExistsPrimalRay());
    ASSERT_FALSE(lpi_->HasPrimalRay());
    ASSERT_FALSE(lpi_->ExistsDualRay());
    ASSERT_FALSE(lpi_->HasDualRay());
  }
}

TEST_P(SolveTest, PrimalUnboundedDualInfeasible) {
  // max 3 x1 +   x2
  //     2 x1 +   x2 <= 10
  //       x1 + 3 x2 <= 15
  //       x1, x2 free
  //
  // which is unbounded.
  for (bool use_primal_simplex : {true, false}) {
    SCOPED_TRACE(absl::StrCat("Using ", use_primal_simplex ? "primal" : "dual",
                              " simplex."));
    ASSERT_OK(UploadProblem(
        /*is_maximization=*/true,
        /*lower_bounds=*/{-solver_inf_, -solver_inf_},
        /*upper_bounds=*/{solver_inf_, solver_inf_},
        /*objective_coefficients=*/{3, 1},
        /*left_hand_sides=*/{-solver_inf_, -solver_inf_},
        /*right_hand_sides=*/{10, 15},
        /*matrix=*/{{2, 1}, {1, 3}}));
    if (use_primal_simplex) {
      ASSERT_OK(lpi_->SolveLPWithPrimalSimplex());
    } else {
      ASSERT_OK(lpi_->SolveLPWithDualSimplex());
    }

    // Check solve status.
    ASSERT_TRUE(lpi_->IsSolved());
    EXPECT_FALSE(lpi_->ObjectiveLimitIsExceeded());
    EXPECT_FALSE(lpi_->IterationLimitIsExceeded());
    EXPECT_FALSE(lpi_->TimeLimitIsExceeded());
    EXPECT_FALSE(lpi_->IsOptimal());
    EXPECT_FALSE(lpi_->IsPrimalFeasible());
    EXPECT_FALSE(lpi_->IsPrimalInfeasible());
    if (use_primal_simplex) EXPECT_TRUE(lpi_->IsPrimalUnbounded());
    EXPECT_FALSE(lpi_->IsDualFeasible());
    EXPECT_TRUE(lpi_->IsDualInfeasible());
    EXPECT_FALSE(lpi_->IsDualUnbounded());
  }
}

TEST_P(SolveTest, PrimalUnboundedRayMaximization) {
  // max 3 x1 +   x2
  //     2 x1 +   x2 <= 10
  //       x1 + 3 x2 <= 15
  //       x1, x2 free
  //
  // which is unbounded.
  ASSERT_OK(UploadProblem(
      /*is_maximization=*/true,
      /*lower_bounds=*/{-solver_inf_, -solver_inf_},
      /*upper_bounds=*/{solver_inf_, solver_inf_},
      /*objective_coefficients=*/{3, 1},
      /*left_hand_sides=*/{-solver_inf_, -solver_inf_},
      /*right_hand_sides=*/{10, 15},
      /*matrix=*/{{2, 1}, {1, 3}}));
  ASSERT_OK(lpi_->SolveLPWithPrimalSimplex());
  ASSERT_TRUE(lpi_->IsSolved());
  EXPECT_TRUE(lpi_->IsPrimalUnbounded());

  // Check rays
  EXPECT_TRUE(lpi_->ExistsPrimalRay());
  if (lpi_->HasPrimalRay()) {
    EXPECT_TRUE(lpi_->ExistsPrimalRay());
    ASSERT_OK_AND_ASSIGN(
        (const absl::StrongVector<ColIndex, double>& primal_ray),
        lpi_->GetPrimalRay());
    ASSERT_EQ(primal_ray.size(), 2);
    const SparseRow sparse_ray = SparseRow(primal_ray);

    // Ray points in direction of improving objective
    EXPECT_THAT(sparse_ray,
                Activation(CreateSparseRow({{0, 3.0}, {1, 1.0}}), Gt(0.0)));

    // Taking a large step in the ray direction still produces a valid
    // point
    const SparseRow feasible_point = CreateSparseRow({});
    const SparseRow next_point = feasible_point + 1000 * sparse_ray;
    EXPECT_THAT(next_point,
                Activation(CreateSparseRow({{0, 2.0}, {1, 1.0}}), Le(10.0)));
    EXPECT_THAT(next_point,
                Activation(CreateSparseRow({{0, 1.0}, {1, 3.0}}), Le(15.0)));
  }
  ASSERT_FALSE(lpi_->ExistsDualRay());
  ASSERT_FALSE(lpi_->HasDualRay());
}

TEST_P(SolveTest, PrimalUnboundedRayMinimization) {
  // min -3 x1 -   x2
  //      2 x1 +   x2 <= 10
  //        x1 + 3 x2 <= 15
  //        x1, x2 free
  //
  // which is unbounded.
  ASSERT_OK(UploadProblem(
      /*is_maximization=*/false,
      /*lower_bounds=*/{-solver_inf_, -solver_inf_},
      /*upper_bounds=*/{solver_inf_, solver_inf_},
      /*objective_coefficients=*/{-3, -1},
      /*left_hand_sides=*/{-solver_inf_, -solver_inf_},
      /*right_hand_sides=*/{10, 15},
      /*matrix=*/{{2, 1}, {1, 3}}));
  ASSERT_OK(lpi_->SolveLPWithPrimalSimplex());
  ASSERT_TRUE(lpi_->IsSolved());
  EXPECT_TRUE(lpi_->IsPrimalUnbounded());

  // Check rays
  EXPECT_TRUE(lpi_->ExistsPrimalRay());
  if (lpi_->HasPrimalRay()) {
    EXPECT_TRUE(lpi_->ExistsPrimalRay());
    ASSERT_OK_AND_ASSIGN(
        (const absl::StrongVector<ColIndex, double>& primal_ray),
        lpi_->GetPrimalRay());
    ASSERT_EQ(primal_ray.size(), 2);
    const SparseRow sparse_ray = SparseRow(primal_ray);

    // Ray points in direction of improving objective
    EXPECT_THAT(sparse_ray,
                Activation(CreateSparseRow({{0, -3.0}, {1, -1.0}}), Lt(0.0)));

    // Taking a large step in the ray direction still produces a valid
    // point
    const SparseRow feasible_point = CreateSparseRow({});
    const SparseRow next_point = feasible_point + 1000 * sparse_ray;
    EXPECT_THAT(next_point,
                Activation(CreateSparseRow({{0, 2.0}, {1, 1.0}}), Le(10.0)));
    EXPECT_THAT(next_point,
                Activation(CreateSparseRow({{0, 1.0}, {1, 3.0}}), Le(15.0)));
  }
  ASSERT_FALSE(lpi_->ExistsDualRay());
  ASSERT_FALSE(lpi_->HasDualRay());
}

TEST_P(SolveTest, PrimalInfeasibleDualUnbounded) {
  // min 10 y1 + 15 y2
  //     2 y1 +   y2 == 3
  //       y1 + 3 y2 == 1
  //       y1,    y2 >= 0
  //
  // This is primal infeasible, since the only solution to the constraints is
  // (y1, y2) = (8/5, -1/5), which is outside the variable bounds. It is
  // however dual unbounded. Note that it is the dual of the previous problem.
  if (GetParam() == LPParameters::LP_SOPLEX) {
    GTEST_SKIP() << "Soplex currently gets confused by infeasible/unbounded "
                    "models, so this test is temporarily disabled.";
  }
  for (bool use_primal_simplex : {true, false}) {
    SCOPED_TRACE(absl::StrCat("Using ", use_primal_simplex ? "primal" : "dual",
                              " simplex."));
    ASSERT_OK(UploadProblem(
        /*is_maximization=*/false,
        /*lower_bounds=*/{0, 0},
        /*upper_bounds=*/{solver_inf_, solver_inf_},
        /*objective_coefficients=*/{10, 15},
        /*left_hand_sides=*/{3, 1},
        /*right_hand_sides=*/{3, 1},
        /*matrix=*/{{2, 1}, {1, 3}}));
    if (use_primal_simplex) {
      ASSERT_OK(lpi_->SolveLPWithPrimalSimplex());
    } else {
      ASSERT_OK(lpi_->SolveLPWithDualSimplex());
    }

    // Check solve status.
    ASSERT_TRUE(lpi_->IsSolved());
    EXPECT_FALSE(lpi_->ObjectiveLimitIsExceeded());
    EXPECT_FALSE(lpi_->IterationLimitIsExceeded());
    EXPECT_FALSE(lpi_->TimeLimitIsExceeded());
    EXPECT_FALSE(lpi_->IsOptimal());
    EXPECT_FALSE(lpi_->IsPrimalFeasible());
    EXPECT_TRUE(lpi_->IsPrimalInfeasible());
    EXPECT_FALSE(lpi_->IsPrimalUnbounded());
    EXPECT_FALSE(lpi_->IsDualFeasible());
    EXPECT_FALSE(lpi_->IsDualInfeasible());
    if (!use_primal_simplex) EXPECT_TRUE(lpi_->IsDualUnbounded());
  }
}

TEST_P(SolveTest, DualUnboundedRayMaximization) {
  // max -10 y1 - 15 y2
  //       2 y1 +    y2 == 3
  //         y1 +  3 y2 == 1
  //         y1,     y2 >= 0
  //
  // This is primal infeasible, since the only solution to the constraints is
  // (y1, y2) = (8/5, -1/5), which is outside the variable bounds. It is
  // however dual unbounded. Note that it is the dual of the previous problem.
  if (GetParam() == LPParameters::LP_SOPLEX) {
    GTEST_SKIP() << "Soplex currently gets confused by infeasible/unbounded "
                    "models, so this test is temporarily disabled.";
  }
  ASSERT_OK(UploadProblem(
      /*is_maximization=*/true,
      /*lower_bounds=*/{0, 0},
      /*upper_bounds=*/{solver_inf_, solver_inf_},
      /*objective_coefficients=*/{-10, -15},
      /*left_hand_sides=*/{3, 1},
      /*right_hand_sides=*/{3, 1},
      /*matrix=*/{{2, 1}, {1, 3}}));
  ASSERT_OK(lpi_->SolveLPWithDualSimplex());
  ASSERT_TRUE(lpi_->IsSolved());
  EXPECT_TRUE(lpi_->IsDualUnbounded());

  // Check rays
  EXPECT_FALSE(lpi_->ExistsPrimalRay());
  EXPECT_FALSE(lpi_->HasPrimalRay());
  EXPECT_TRUE(lpi_->ExistsDualRay());
  if (lpi_->HasDualRay()) {
    EXPECT_TRUE(lpi_->ExistsDualRay());
    ASSERT_OK_AND_ASSIGN((const absl::StrongVector<RowIndex, double>& dual_ray),
                         lpi_->GetDualRay());
    ASSERT_EQ(dual_ray.size(), 2);
    const SparseCol sparse_ray = SparseCol(dual_ray);

    // Ray points in direction of improving objective
    EXPECT_THAT(sparse_ray,
                Activation(CreateSparseCol({{0, -10.0}, {1, -15.0}}), Gt(0.0)));

    // Taking a large step in the ray direction still produces a valid
    // point
    const SparseCol feasible_point = CreateSparseCol({});
    const SparseCol next_point = feasible_point + 1000 * sparse_ray;
    EXPECT_THAT(next_point,
                Activation(CreateSparseCol({{0, 2.0}, {1, 1.0}}), Le(10.0)));
    EXPECT_THAT(next_point,
                Activation(CreateSparseCol({{0, 1.0}, {1, 3.0}}), Le(15.0)));
  }
}

TEST_P(SolveTest, DualUnboundedRayMinimization) {
  // min 10 y1 + 15 y2
  //     2 y1 +   y2 == 3
  //       y1 + 3 y2 == 1
  //       y1,    y2 >= 0
  //
  // This is primal infeasible, since the only solution to the constraints is
  // (y1, y2) = (8/5, -1/5), which is outside the variable bounds. It is
  // however dual unbounded. Note that it is the dual of the previous problem.
  if (GetParam() == LPParameters::LP_SOPLEX) {
    GTEST_SKIP() << "Soplex currently gets confused by infeasible/unbounded "
                    "models, so this test is temporarily disabled.";
  }
  ASSERT_OK(UploadProblem(
      /*is_maximization=*/false,
      /*lower_bounds=*/{0, 0},
      /*upper_bounds=*/{solver_inf_, solver_inf_},
      /*objective_coefficients=*/{10, 15},
      /*left_hand_sides=*/{3, 1},
      /*right_hand_sides=*/{3, 1},
      /*matrix=*/{{2, 1}, {1, 3}}));
  ASSERT_OK(lpi_->SolveLPWithDualSimplex());
  ASSERT_TRUE(lpi_->IsSolved());
  EXPECT_TRUE(lpi_->IsDualUnbounded());

  // Check rays
  EXPECT_FALSE(lpi_->ExistsPrimalRay());
  EXPECT_FALSE(lpi_->HasPrimalRay());
  EXPECT_TRUE(lpi_->ExistsDualRay());
  if (lpi_->HasDualRay()) {
    EXPECT_TRUE(lpi_->ExistsDualRay());
    ASSERT_OK_AND_ASSIGN((const absl::StrongVector<RowIndex, double>& dual_ray),
                         lpi_->GetDualRay());
    ASSERT_EQ(dual_ray.size(), 2);
    const SparseCol sparse_ray = SparseCol(dual_ray);

    // Ray points in direction of improving objective
    EXPECT_THAT(sparse_ray,
                Activation(CreateSparseCol({{0, 10.0}, {1, 15.0}}), Lt(0.0)));

    // Taking a large step in the ray direction still produces a valid
    // point
    const SparseCol feasible_point = CreateSparseCol({});
    const SparseCol next_point = feasible_point + 1000 * sparse_ray;
    EXPECT_THAT(next_point,
                Activation(CreateSparseCol({{0, 2.0}, {1, 1.0}}), Le(10.0)));
    EXPECT_THAT(next_point,
                Activation(CreateSparseCol({{0, 1.0}, {1, 3.0}}), Le(15.0)));
  }
}

TEST_P(SolveTest, PrimalDualInfeasible) {
  // max x1 + x2
  //    x1 - x2 <= 0
  //  - x1 + x2 <= -1
  //    x1,  x2 free
  //
  // Note that the two constraints contradict each other, both in the
  // primal and in the dual.
  if (GetParam() == LPParameters::LP_SOPLEX) {
    GTEST_SKIP() << "Soplex currently gets confused by infeasible/unbounded "
                    "models, so this test is temporarily disabled.";
  }
  for (bool use_primal_simplex : {true, false}) {
    SCOPED_TRACE(absl::StrCat("Using ", use_primal_simplex ? "primal" : "dual",
                              " simplex."));
    ASSERT_OK(UploadProblem(
        /*is_maximization=*/true,
        /*lower_bounds=*/{-solver_inf_, -solver_inf_},
        /*upper_bounds=*/{solver_inf_, solver_inf_},
        /*objective_coefficients=*/{1, 1},
        /*left_hand_sides=*/{-solver_inf_, -solver_inf_},
        /*right_hand_sides=*/{0, -1},
        /*matrix=*/{{1, -1}, {-1, 1}}));
    if (use_primal_simplex) {
      ASSERT_OK(lpi_->SolveLPWithPrimalSimplex());
    } else {
      ASSERT_OK(lpi_->SolveLPWithDualSimplex());
    }

    // Check solve status.
    ASSERT_TRUE(lpi_->IsSolved());
    EXPECT_FALSE(lpi_->ObjectiveLimitIsExceeded());
    EXPECT_FALSE(lpi_->IterationLimitIsExceeded());
    EXPECT_FALSE(lpi_->TimeLimitIsExceeded());
    EXPECT_FALSE(lpi_->IsOptimal());
    EXPECT_FALSE(lpi_->IsPrimalFeasible());
    if (use_primal_simplex) EXPECT_TRUE(lpi_->IsPrimalInfeasible());
    EXPECT_FALSE(lpi_->IsPrimalUnbounded());
    EXPECT_FALSE(lpi_->IsDualFeasible());
    if (!use_primal_simplex) EXPECT_TRUE(lpi_->IsDualInfeasible());
    EXPECT_FALSE(lpi_->IsDualUnbounded());

    // Check rays.
    EXPECT_FALSE(lpi_->ExistsPrimalRay());
    EXPECT_FALSE(lpi_->HasPrimalRay());
    EXPECT_FALSE(lpi_->ExistsDualRay());
    EXPECT_FALSE(lpi_->HasDualRay());
  }
}

TEST_P(SolveTest, SolveFromGivenPrimalBasis) {
  // max x1 + x2
  // 2 x1 +   x2 <= 3
  //   x1 + 2 x2 <= 3
  // x1, x2 >= 0
  // is primal/dual feasible. The optimum is in (1, 1), where both constraints
  // are active.
  ASSERT_OK(UploadProblem(
      /*is_maximization=*/true,
      /*lower_bounds=*/{0, 0},
      /*upper_bounds=*/{solver_inf_, solver_inf_},
      /*objective_coefficients=*/{1, 1},
      /*left_hand_sides=*/{-solver_inf_, -solver_inf_},
      /*right_hand_sides=*/{3, 3},
      /*matrix=*/{{2, 1}, {1, 2}}));

  // Starting at the origin, we can take two paths to the optimum, both
  // requiring two iterations.
  ASSERT_OK(lpi_->SetBasisStatusForColumnsAndRows(
      {LPBasisStatus::kAtLowerBound, LPBasisStatus::kAtLowerBound},
      {LPBasisStatus::kBasic, LPBasisStatus::kBasic}));
  ASSERT_OK(lpi_->SolveLPWithPrimalSimplex());
  EXPECT_EQ(lpi_->GetNumIterations(), 2);

  {
    // Check that we arrived at the correct solution.
    ASSERT_TRUE(lpi_->IsSolved());
    ASSERT_TRUE(lpi_->IsOptimal());
    ASSERT_OK_AND_ASSIGN((const absl::StrongVector<ColIndex, double> solution),
                         lpi_->GetPrimalValues());
    EXPECT_THAT(solution, Pointwise(DoubleEq(), {1.0, 1.0}));
    ASSERT_OK_AND_ASSIGN(
        (const absl::StrongVector<ColIndex, LPBasisStatus> column_basis_status),
        lpi_->GetBasisStatusForColumns());
    EXPECT_THAT(column_basis_status,
                ElementsAre(LPBasisStatus::kBasic, LPBasisStatus::kBasic));
    ASSERT_OK_AND_ASSIGN(
        (const absl::StrongVector<RowIndex, LPBasisStatus> row_basis_status),
        lpi_->GetBasisStatusForRows());
    EXPECT_THAT(row_basis_status, ElementsAre(LPBasisStatus::kAtUpperBound,
                                              LPBasisStatus::kAtUpperBound));
  }

  // We now start in the point (1.5, 0), from which the optimum can be found
  // in a single iteration.
  ASSERT_OK(lpi_->SetBasisStatusForColumnsAndRows(
      {LPBasisStatus::kBasic, LPBasisStatus::kAtLowerBound},
      {LPBasisStatus::kAtUpperBound, LPBasisStatus::kBasic}));
  ASSERT_OK(lpi_->SolveLPWithPrimalSimplex());
  EXPECT_EQ(lpi_->GetNumIterations(), 1);

  {
    // Check that we arrived at the correct solution.
    ASSERT_TRUE(lpi_->IsSolved());
    ASSERT_TRUE(lpi_->IsOptimal());
    ASSERT_OK_AND_ASSIGN((const absl::StrongVector<ColIndex, double> solution),
                         lpi_->GetPrimalValues());
    EXPECT_THAT(solution, Pointwise(DoubleEq(), {1.0, 1.0}));
    ASSERT_OK_AND_ASSIGN(
        (const absl::StrongVector<ColIndex, LPBasisStatus> column_basis_status),
        lpi_->GetBasisStatusForColumns());
    EXPECT_THAT(column_basis_status,
                ElementsAre(LPBasisStatus::kBasic, LPBasisStatus::kBasic));
    ASSERT_OK_AND_ASSIGN(
        (const absl::StrongVector<RowIndex, LPBasisStatus> row_basis_status),
        lpi_->GetBasisStatusForRows());
    EXPECT_THAT(row_basis_status, ElementsAre(LPBasisStatus::kAtUpperBound,
                                              LPBasisStatus::kAtUpperBound));
  }
}

TEST_P(SolveTest, SolveFromGivenDualBasis) {
  // max x1 + x2
  // 2 x1 +   x2 <= 3
  //   x1 + 2 x2 <= 3
  // x1, x2 >= 0
  // is primal/dual feasible. The optimum is in (1, 1), where both constraints
  // are active.
  ASSERT_OK(UploadProblem(
      /*is_maximization=*/true,
      /*lower_bounds=*/{0, 0},
      /*upper_bounds=*/{solver_inf_, solver_inf_},
      /*objective_coefficients=*/{1, 1},
      /*left_hand_sides=*/{-solver_inf_, -solver_inf_},
      /*right_hand_sides=*/{3, 3},
      /*matrix=*/{{2, 1}, {1, 2}}));

  // We now start in the point (3, 0), which is dual feasible but not optimal.
  // From here the optimum can be found in a single iteration.
  ASSERT_OK(lpi_->SetBasisStatusForColumnsAndRows(
      {LPBasisStatus::kBasic, LPBasisStatus::kAtLowerBound},
      {LPBasisStatus::kBasic, LPBasisStatus::kAtUpperBound}));
  ASSERT_OK(lpi_->SolveLPWithDualSimplex());
  EXPECT_EQ(lpi_->GetNumIterations(), 1);

  {
    // Check that we arrived at the correct solution.
    ASSERT_TRUE(lpi_->IsSolved());
    ASSERT_TRUE(lpi_->IsOptimal());
    ASSERT_OK_AND_ASSIGN((const absl::StrongVector<ColIndex, double> solution),
                         lpi_->GetPrimalValues());
    EXPECT_THAT(solution, Pointwise(DoubleEq(), {1.0, 1.0}));
    ASSERT_OK_AND_ASSIGN(
        (const absl::StrongVector<ColIndex, LPBasisStatus> column_basis_status),
        lpi_->GetBasisStatusForColumns());
    EXPECT_THAT(column_basis_status,
                ElementsAre(LPBasisStatus::kBasic, LPBasisStatus::kBasic));
    ASSERT_OK_AND_ASSIGN(
        (const absl::StrongVector<RowIndex, LPBasisStatus> row_basis_status),
        lpi_->GetBasisStatusForRows());
    EXPECT_THAT(row_basis_status, ElementsAre(LPBasisStatus::kAtUpperBound,
                                              LPBasisStatus::kAtUpperBound));
  }
}

TEST_P(SolveTest, DualIncrementality) {
  //  max x1 + x2
  //  2 x1 +   x2 <= 3
  //  x1, x2 >= 0
  //  is primal/dual feasible. The optimum is in (0, 3).
  ASSERT_OK(UploadProblem(
      /*is_maximization=*/true,
      /*lower_bounds=*/{0, 0},
      /*upper_bounds=*/{solver_inf_, solver_inf_},
      /*objective_coefficients=*/{1, 1},
      /*left_hand_sides=*/{-solver_inf_},
      /*right_hand_sides=*/{3},
      /*matrix=*/{{2}, {1}}));
  ASSERT_OK(lpi_->SolveLPWithPrimalSimplex());
  ASSERT_TRUE(lpi_->IsSolved());
  ASSERT_TRUE(lpi_->IsOptimal());

  {
    // Check that we arrived at the correct solution.
    ASSERT_OK_AND_ASSIGN((const absl::StrongVector<ColIndex, double> solution),
                         lpi_->GetPrimalValues());
    EXPECT_THAT(solution, Pointwise(DoubleEq(), {0.0, 3.0}));
    ASSERT_OK_AND_ASSIGN(
        (const absl::StrongVector<ColIndex, LPBasisStatus> column_basis_status),
        lpi_->GetBasisStatusForColumns());
    EXPECT_THAT(column_basis_status, ElementsAre(LPBasisStatus::kAtLowerBound,
                                                 LPBasisStatus::kBasic));
    ASSERT_OK_AND_ASSIGN(
        (const absl::StrongVector<RowIndex, LPBasisStatus> row_basis_status),
        lpi_->GetBasisStatusForRows());
    EXPECT_THAT(row_basis_status, ElementsAre(LPBasisStatus::kAtUpperBound));
  }

  // We now add the constraint x2 <= 2, making (0.5, 2) the new optimum.
  ASSERT_OK(lpi_->AddRow(SparseRow({{ColIndex(1), 1}}), -solver_inf_, 2, ""));
  EXPECT_FALSE(lpi_->IsSolved());

  // Since the basis should be retained, it should now only take a single
  // iteration to reach the new optimum.
  ASSERT_OK(lpi_->SolveLPWithDualSimplex());
  ASSERT_TRUE(lpi_->IsSolved());
  ASSERT_TRUE(lpi_->IsOptimal());
  EXPECT_EQ(lpi_->GetNumIterations(), 1);

  {
    // Check that we arrived at the correct solution.
    ASSERT_OK_AND_ASSIGN((const absl::StrongVector<ColIndex, double> solution),
                         lpi_->GetPrimalValues());
    EXPECT_THAT(solution, Pointwise(DoubleEq(), {0.5, 2.0}));
    ASSERT_OK_AND_ASSIGN(
        (const absl::StrongVector<ColIndex, LPBasisStatus> column_basis_status),
        lpi_->GetBasisStatusForColumns());
    EXPECT_THAT(column_basis_status,
                ElementsAre(LPBasisStatus::kBasic, LPBasisStatus::kBasic));
    ASSERT_OK_AND_ASSIGN(
        (const absl::StrongVector<RowIndex, LPBasisStatus> row_basis_status),
        lpi_->GetBasisStatusForRows());
    EXPECT_THAT(row_basis_status, ElementsAre(LPBasisStatus::kAtUpperBound,
                                              LPBasisStatus::kAtUpperBound));
  }
}

}  // namespace
}  // namespace minimip
