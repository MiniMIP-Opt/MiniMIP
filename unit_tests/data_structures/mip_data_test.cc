// Copyright 2024 the MiniMIP Project
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

#include "src/data_structures/mip_data.h"

#include <limits>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "unit_tests/utils.h"

// TODO(CGraczyk): add randomized O(1000) init tests for mip data.

namespace minimip {
using ::testing::ElementsAre;
using ::testing::ElementsAreArray;
using testing::IsEmpty;
using ::testing::UnorderedElementsAreArray;

TEST(MipDataTests, CreateEmptyProblem) {
  MipData mip_data;
  EXPECT_EQ(mip_data.problem_name(), "");
  EXPECT_EQ(mip_data.is_maximization(), false);
  EXPECT_EQ(mip_data.objective_offset(), 0);

  EXPECT_TRUE(mip_data.objective().entries().empty());
  EXPECT_TRUE(mip_data.hints().empty());
  EXPECT_TRUE(mip_data.lower_bounds().empty());
  EXPECT_TRUE(mip_data.upper_bounds().empty());
  EXPECT_TRUE(mip_data.left_hand_sides().empty());
  EXPECT_TRUE(mip_data.right_hand_sides().empty());
  EXPECT_TRUE(mip_data.integer_variables().empty());
  EXPECT_TRUE(mip_data.variable_types().empty());
  EXPECT_TRUE(mip_data.variable_names().empty());
  EXPECT_TRUE(mip_data.constraint_names().empty());

  EXPECT_EQ(mip_data.matrix().num_cols(), ColIndex(0));
  EXPECT_EQ(mip_data.matrix().num_rows(), RowIndex(0));
}

TEST(MipDataTests, PopulatesProblemName) {
  const MiniMipProblem problem = {.name = "Foo"};
  const MipData mip_data(problem);
  EXPECT_EQ(mip_data.problem_name(), "Foo");
}

TEST(MipDataTests, PopulatesVariables) {
  const MiniMipVariable variable = {.name = "Bar",
                                    .objective_coefficient = 13.0,
                                    .lower_bound = 0.0,
                                    .upper_bound = 1.0,
                                    .is_integer = true};
  const MiniMipProblem problem = {.variables = {variable}};

  const MipData mip_data(problem);
  EXPECT_EQ(mip_data.lower_bounds().size(), 1);
  EXPECT_EQ(mip_data.lower_bounds()[ColIndex(0)], 0.0);
  EXPECT_EQ(mip_data.upper_bounds().size(), 1);
  EXPECT_EQ(mip_data.upper_bounds()[ColIndex(0)], 1);
  EXPECT_EQ(mip_data.objective().value(ColIndex(0)), 13.0);

  EXPECT_THAT(mip_data.integer_variables(), UnorderedElementsAreArray({0}));
  EXPECT_THAT(mip_data.variable_types(),
              ElementsAreArray({VariableType::kInteger}));
  EXPECT_THAT(mip_data.variable_names(), ElementsAreArray({"Bar"}));
}

TEST(MipDataTests, PopulatesConstraints) {
  const MiniMipVariable variable = {.name = "Bar",
                                    .objective_coefficient = 13.0,
                                    .lower_bound = 0.0,
                                    .upper_bound = 1.0,
                                    .is_integer = true};
  const MiniMipConstraint constraint = {.name = "Baz",
                                        .var_indices = {0},
                                        .coefficients = {1.5},
                                        .left_hand_side = -1.0,
                                        .right_hand_side = 0.5};
  const MiniMipProblem problem = {.variables = {variable},
                                  .constraints = {constraint}};
  const MipData mip_data = MipData(problem);

  EXPECT_THAT(mip_data.left_hand_sides(), ElementsAreArray({-1.0}));
  EXPECT_THAT(mip_data.right_hand_sides(), ElementsAreArray({0.5}));
  EXPECT_THAT(mip_data.constraint_names(), ElementsAreArray({"Baz"}));

  const StrongSparseMatrix& constraint_matrix = mip_data.matrix();
  EXPECT_EQ(constraint_matrix.GetCoefficient(ColIndex(0), RowIndex(0)), 1.5);
}

TEST(MipDataTests, IdentifiesIntegralConstraint) {
  const MiniMipVariable integer_variable = {.name = "integer",
                                            .is_integer = true};
  const MiniMipVariable continuous_variable = {.name = "continuous",
                                               .is_integer = false};

  {
    const MiniMipConstraint integral_constraint = {.name = "integral",
                                                   .var_indices = {0},
                                                   .coefficients = {1},
                                                   .left_hand_side = 0.5,
                                                   .right_hand_side = 2.5};
    const MiniMipProblem problem = {.variables = {integer_variable},
                                    .constraints = {integral_constraint}};
    const MipData mip_data = MipData(problem);
    EXPECT_THAT(mip_data.is_integral_constraint(), ElementsAre(true));
    EXPECT_THAT(mip_data.left_hand_sides(), ElementsAre(0.5));
    EXPECT_THAT(mip_data.right_hand_sides(), ElementsAre(2.5));
  }

  {
    const MiniMipConstraint integral_constraint = {.name = "integral",
                                                   .var_indices = {0, 1},
                                                   .coefficients = {1, 0.0},
                                                   .left_hand_side = 0.5,
                                                   .right_hand_side = 2.5};
    const MiniMipProblem problem = {
        .variables = {integer_variable, continuous_variable},
        .constraints = {integral_constraint}};
    const MipData mip_data = MipData(problem);
    EXPECT_THAT(mip_data.is_integral_constraint(), ElementsAre(true));
    EXPECT_THAT(mip_data.left_hand_sides(), ElementsAre(0.5));
    EXPECT_THAT(mip_data.right_hand_sides(), ElementsAre(2.5));
  }
}

TEST(MipDataTests, IdentifiesNonIntegralConstraint) {
  const MiniMipVariable integer_variable = {.name = "integer",
                                            .is_integer = true};
  const MiniMipVariable continuous_variable = {.name = "continuous",
                                               .is_integer = false};

  {
    const MiniMipConstraint constraint_with_continuous_variable = {
        .name = "nonintegral",
        .var_indices = {0, 1},
        .coefficients = {1, 1},
        .left_hand_side = 0.5,
        .right_hand_side = 2.5};
    const MiniMipProblem problem = {
        .variables = {integer_variable, continuous_variable},
        .constraints = {constraint_with_continuous_variable}};
    const MipData mip_data = MipData(problem);
    EXPECT_THAT(mip_data.is_integral_constraint(), ElementsAre(false));
    EXPECT_THAT(mip_data.left_hand_sides(), ElementsAre(0.5));
    EXPECT_THAT(mip_data.right_hand_sides(), ElementsAre(2.5));
  }

  {
    const MiniMipConstraint constraint_with_real_coefficient = {
        .name = "nonintegral",
        .var_indices = {0},
        .coefficients = {1.1},
        .left_hand_side = 0.5,
        .right_hand_side = 2.5};
    const MiniMipProblem problem = {
        .variables = {integer_variable, continuous_variable},
        .constraints = {constraint_with_real_coefficient}};
    const MipData mip_data = MipData(problem);
    EXPECT_THAT(mip_data.is_integral_constraint(), ElementsAre(false));
    EXPECT_THAT(mip_data.left_hand_sides(), ElementsAre(0.5));
    EXPECT_THAT(mip_data.right_hand_sides(), ElementsAre(2.5));
  }
}

TEST(MipDataTests, PopulatesHints) {
  const MiniMipSolutionHint hint = {.var_indices = {0}, .values = {0.5}};
  const MiniMipProblem problem = {.hints = {hint}};
  const MipData mip_data = MipData(problem);

  EXPECT_THAT(mip_data.hints()[0].values, ElementsAreArray({0.5}));
  EXPECT_THAT(mip_data.hints()[0].var_indices, ElementsAreArray({0.0}));
}

TEST(MipDataTests, PopulatesMipDataDirection) {
  const MiniMipProblem problem = {.is_maximization = true};
  const MipData mip_data = MipData(problem);

  EXPECT_EQ(mip_data.is_maximization(), true);
}

TEST(MipDataTests, PopulatesMipDataOffset) {
  const MiniMipProblem problem = {.objective_offset = 12.0};
  const MipData mip_data = MipData(problem);

  EXPECT_EQ(mip_data.objective_offset(), 12.0);
}

TEST(MipDataTests, PopulatesMipDataFromMiniMipProblemWithVarBoundConstraint) {
  // Test Problem
  // min           x1 + 2 x2 + 3 x3 +1
  // s.t.
  //     0.0 <=    x1               <= inf
  //    -inf <=  3 x1 + 5 x2        <= -1.0
  //    -2.0 <=  5 x1 + 2 x2 + 7 x3 <= 3.0
  //                 x1,x2 Integer, x3 Real

  const std::vector<MiniMipVariable> variables = {
      {.name = "x1",
       .objective_coefficient = 1.0,
       .lower_bound = 0,
       .upper_bound = std::numeric_limits<double>::infinity(),
       .is_integer = true},
      {.name = "x2",
       .objective_coefficient = 2.0,
       .lower_bound = -std::numeric_limits<double>::infinity(),
       .upper_bound = std::numeric_limits<double>::infinity(),
       .is_integer = true},
      {.name = "x3",
       .objective_coefficient = 3.0,
       .lower_bound = -std::numeric_limits<double>::infinity(),
       .upper_bound = std::numeric_limits<double>::infinity(),
       .is_integer = false}};

  const std::vector<MiniMipConstraint> constraints = {
      {.name = "first",
       .var_indices = {0, 1},
       .coefficients = {3.0, 5.0},
       .left_hand_side = -std::numeric_limits<double>::infinity(),
       .right_hand_side = -1},
      {.name = "second",
       .var_indices = {0, 1, 2},
       .coefficients = {5.0, 2.0, 7.0},
       .left_hand_side = -2,
       .right_hand_side = 3}};

  const MiniMipProblem problem = {.name = "Test",
                                  .variables = variables,
                                  .constraints = constraints,
                                  .hints = {},
                                  .is_maximization = false,
                                  .objective_offset = 1.0};

  const MipData mip_data = MipData(problem);

  EXPECT_EQ(mip_data.problem_name(), "Test");
  EXPECT_EQ(mip_data.is_maximization(), false);
  EXPECT_EQ(mip_data.objective_offset(), 1.0);

  const SparseRow& objective = mip_data.objective();
  EXPECT_EQ(objective.entries().size(), 3);

  EXPECT_THAT(
      objective.entries(),
      ElementsAreArray({RowEntry(ColIndex(0), 1.0), RowEntry(ColIndex(1), 2.0),
                        RowEntry(ColIndex(2), 3.0)}));

  ASSERT_TRUE(mip_data.hints().empty());
  EXPECT_THAT(mip_data.lower_bounds(),
              ElementsAreArray({0.0, -std::numeric_limits<double>::infinity(),
                                -std::numeric_limits<double>::infinity()}));
  EXPECT_THAT(mip_data.upper_bounds(),
              ElementsAreArray({std::numeric_limits<double>::infinity(),
                                std::numeric_limits<double>::infinity(),
                                std::numeric_limits<double>::infinity()}));
  EXPECT_THAT(
      mip_data.left_hand_sides(),
      ElementsAreArray({-std::numeric_limits<double>::infinity(), -2.0}));
  EXPECT_THAT(mip_data.right_hand_sides(), ElementsAreArray({-1.0, 3.0}));
  EXPECT_THAT(mip_data.integer_variables(), UnorderedElementsAreArray({0, 1}));
  EXPECT_THAT(mip_data.variable_types(),
              ElementsAreArray({VariableType::kInteger, VariableType::kInteger,
                                VariableType::kFractional}));
  EXPECT_THAT(mip_data.variable_names(), ElementsAreArray({"x1", "x2", "x3"}));
  EXPECT_THAT(mip_data.constraint_names(),
              ElementsAreArray({"first", "second"}));

  const StrongSparseMatrix& constraint_matrix = mip_data.matrix();
  EXPECT_EQ(constraint_matrix.num_rows(), RowIndex(2));
  EXPECT_EQ(constraint_matrix.num_cols(), ColIndex(3));

  EXPECT_THAT(
      constraint_matrix.row(RowIndex(0)).entries(),
      ElementsAre(RowEntry(ColIndex(0), 3.0), RowEntry(ColIndex(1), 5.0)));
  EXPECT_THAT(
      constraint_matrix.row(RowIndex(1)).entries(),
      ElementsAre(RowEntry(ColIndex(0), 5.0), RowEntry(ColIndex(1), 2.0),
                  RowEntry(ColIndex(2), 7.0)));

  EXPECT_THAT(
      constraint_matrix.col(ColIndex(0)).entries(),
      ElementsAre(ColEntry(RowIndex(0), 3.0), ColEntry(RowIndex(1), 5.0)));
  EXPECT_THAT(
      constraint_matrix.col(ColIndex(1)).entries(),
      ElementsAre(ColEntry(RowIndex(0), 5.0), ColEntry(RowIndex(1), 2.0)));
  EXPECT_THAT(constraint_matrix.col(ColIndex(2)).entries(),
              ElementsAre(ColEntry(RowIndex(1), 7.0)));

  for (int row_idx = 0; row_idx < constraints.size(); ++row_idx) {
    for (int col_idx = 0; col_idx < constraints[row_idx].coefficients.size();
         ++col_idx) {
      EXPECT_EQ(constraint_matrix.GetCoefficient(ColIndex(col_idx),
                                                 RowIndex(row_idx)),
                constraints[row_idx].coefficients[col_idx]);
    }
  }
}

TEST(MipDataTests, PopulatesMipDataFrom3x3MiniMipProblem) {
  // Test Problem
  // min           x1 + 2 x2 + 3 x3 +1
  // s.t.
  //     0.0 <=  1 x1 +  .   + 1 x3 <= inf
  //    -inf <=  3 x1 + 5 x2    .   <= -1.0
  //    -2.0 <=  5 x1 + 2 x2 + 7 x3 <= 3.0
  //                 x1,x2 Integer, x3 Real

  std::vector<MiniMipVariable> variables = {
      {.name = "x1",
       .objective_coefficient = 1.0,
       .lower_bound = -std::numeric_limits<double>::infinity(),
       .upper_bound = std::numeric_limits<double>::infinity(),
       .is_integer = true},
      {.name = "x2",
       .objective_coefficient = 2.0,
       .lower_bound = -std::numeric_limits<double>::infinity(),
       .upper_bound = std::numeric_limits<double>::infinity(),
       .is_integer = true},
      {.name = "x3",
       .objective_coefficient = 3.0,
       .lower_bound = -std::numeric_limits<double>::infinity(),
       .upper_bound = std::numeric_limits<double>::infinity(),
       .is_integer = false}};

  std::vector<MiniMipConstraint> constraints = {
      {.name = "first",
       .var_indices = {0, 2},
       .coefficients = {1.0, 1.0},
       .left_hand_side = 0,
       .right_hand_side = std::numeric_limits<double>::infinity()},
      {.name = "second",
       .var_indices = {0, 1},
       .coefficients = {3.0, 5.0},
       .left_hand_side = -std::numeric_limits<double>::infinity(),
       .right_hand_side = -1},
      {.name = "third",
       .var_indices = {0, 1, 2},
       .coefficients = {5.0, 2.0, 7.0},
       .left_hand_side = -2,
       .right_hand_side = 3}};

  MiniMipProblem problem = {.name = "Test",
                            .variables = variables,
                            .constraints = constraints,
                            .hints = {},
                            .is_maximization = false,
                            .objective_offset = 1.0};

  MipData mip_data = MipData(problem);

  EXPECT_EQ(mip_data.problem_name(), "Test");
  EXPECT_EQ(mip_data.is_maximization(), false);
  EXPECT_EQ(mip_data.objective_offset(), 1.0);

  const SparseRow& objective = mip_data.objective();
  EXPECT_EQ(objective.entries().size(), 3);
  EXPECT_THAT(
      objective.entries(),
      ElementsAreArray({RowEntry(ColIndex(0), 1.0), RowEntry(ColIndex(1), 2.0),
                        RowEntry(ColIndex(2), 3.0)}));

  ASSERT_TRUE(mip_data.hints().empty());
  EXPECT_THAT(mip_data.lower_bounds(),
              ElementsAreArray({-std::numeric_limits<double>::infinity(),
                                -std::numeric_limits<double>::infinity(),
                                -std::numeric_limits<double>::infinity()}));
  EXPECT_THAT(mip_data.upper_bounds(),
              ElementsAreArray({std::numeric_limits<double>::infinity(),
                                std::numeric_limits<double>::infinity(),
                                std::numeric_limits<double>::infinity()}));
  EXPECT_THAT(
      mip_data.left_hand_sides(),
      ElementsAreArray({0.0, -std::numeric_limits<double>::infinity(), -2.0}));
  EXPECT_THAT(
      mip_data.right_hand_sides(),
      ElementsAreArray({std::numeric_limits<double>::infinity(), -1.0, 3.0}));
  EXPECT_THAT(mip_data.integer_variables(), UnorderedElementsAreArray({0, 1}));
  EXPECT_THAT(mip_data.variable_types(),
              ElementsAreArray({VariableType::kInteger, VariableType::kInteger,
                                VariableType::kFractional}));
  EXPECT_THAT(mip_data.variable_names(), ElementsAreArray({"x1", "x2", "x3"}));
  EXPECT_THAT(mip_data.constraint_names(),
              ElementsAreArray({"first", "second", "third"}));

  const StrongSparseMatrix& constraint_matrix = mip_data.matrix();
  EXPECT_EQ(constraint_matrix.num_rows(), RowIndex(3));
  EXPECT_EQ(constraint_matrix.num_cols(), ColIndex(3));

  EXPECT_THAT(
      constraint_matrix.row(RowIndex(0)).entries(),
      ElementsAre(RowEntry(ColIndex(0), 1.0), RowEntry(ColIndex(2), 1.0)));
  EXPECT_THAT(
      constraint_matrix.row(RowIndex(1)).entries(),
      ElementsAre(RowEntry(ColIndex(0), 3.0), RowEntry(ColIndex(1), 5.0)));
  EXPECT_THAT(
      constraint_matrix.row(RowIndex(2)).entries(),
      ElementsAre(RowEntry(ColIndex(0), 5.0), RowEntry(ColIndex(1), 2.0),
                  RowEntry(ColIndex(2), 7.0)));

  EXPECT_THAT(
      constraint_matrix.col(ColIndex(0)).entries(),
      ElementsAre(ColEntry(RowIndex(0), 1.0), ColEntry(RowIndex(1), 3.0),
                  ColEntry(RowIndex(2), 5.0)));
  EXPECT_THAT(
      constraint_matrix.col(ColIndex(1)).entries(),
      ElementsAre(ColEntry(RowIndex(1), 5.0), ColEntry(RowIndex(2), 2.0)));
  EXPECT_THAT(
      constraint_matrix.col(ColIndex(2)).entries(),
      ElementsAre(ColEntry(RowIndex(0), 1.0), ColEntry(RowIndex(2), 7.0)));

  for (int row_idx = 0; row_idx < constraints.size(); ++row_idx) {
    for (int col_idx = 0; col_idx < constraints[row_idx].coefficients.size();
         ++col_idx) {
      EXPECT_EQ(constraint_matrix.GetCoefficient(
                    ColIndex(constraints[row_idx].var_indices[col_idx]),
                    RowIndex(row_idx)),
                constraints[row_idx].coefficients[col_idx])
          << absl::StrFormat("row_idx: %i, col_idx: %i", row_idx, col_idx);
    }
  }
}

}  // namespace minimip
