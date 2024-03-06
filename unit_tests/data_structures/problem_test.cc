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

#include "src/data_structures/problem.h"

#include <algorithm>
#include <fstream>
#include <ios>
#include <limits>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "unit_tests/utils.h"

// TODO(CGraczyk): add randomized O(1000) init tests for mip data.

namespace minimip {
namespace {

absl::StatusOr<std::string> ReadFileToString(const std::string& file_name) {
  std::fstream fs(file_name, std::ios_base::in);
  if (!fs.good()) {
    return util::InvalidArgumentErrorBuilder()
           << "File not found: " << file_name;
  }
  std::string data;
  std::string line;
  while (std::getline(fs, line)) absl::StrAppend(&data, line, "\n");
  return data;
}

}  // namespace

using testing::Each;
using testing::Eq;
using testing::Field;
using testing::Ge;
using testing::Gt;
using testing::IsEmpty;
using testing::Le;
using testing::SizeIs;

TEST(MiniMipProblemTest, LoadProblemFromMPSFile) {
  constexpr double kInf = std::numeric_limits<double>::infinity();

  // We rely on ortools to provide the core reader, so we assume it is well
  // tested. We therefore only test our own reader by loading 50v-10 from
  // MIPLIB2017, which is a medium-size problem with varying types of variables.
  // https://miplib.zib.de/instance_details_50v-10.html
  ASSERT_OK_AND_ASSIGN(const MiniMipProblem problem,
                       ReadProblemFromMPSFile(
                           "unit_tests/data_structures/test_data/50v-10.mps"));
  EXPECT_EQ(problem.name, "50v-10");
  EXPECT_THAT(problem.hints, IsEmpty());
  EXPECT_EQ(problem.objective_offset, 0.0);

  // Strangely, the MPS format doesn't specify optimization direction.
  // Therefore, we assume minimization as default.
  EXPECT_EQ(problem.is_maximization, false);

  EXPECT_THAT(problem.variables, SizeIs(2013));
  EXPECT_EQ(
      std::count_if(problem.variables.begin(), problem.variables.end(),
                    [](const MiniMipVariable& var) { return var.is_integer; }),
      1647);

  // In this problem, all variables are non-negative. 366 of them are
  // unbounded and don't participate in the objective, the rest are bounded and
  // participate with positive coefficient in the objective.
  EXPECT_THAT(problem.variables,
              Each(Field(&MiniMipVariable::lower_bound, Eq(0.0))));
  EXPECT_THAT(problem.variables,
              Each(Field(&MiniMipVariable::upper_bound, Ge(0.0))));
  EXPECT_THAT(problem.variables,
              Each(Field(&MiniMipVariable::objective_coefficient, Ge(0.0))));
  EXPECT_EQ(std::count_if(problem.variables.begin(), problem.variables.end(),
                          [](const MiniMipVariable& variable) {
                            return variable.upper_bound == kInf;
                          }),
            366);
  EXPECT_EQ(std::count_if(problem.variables.begin(), problem.variables.end(),
                          [](const MiniMipVariable& variable) {
                            return variable.objective_coefficient == 0.0;
                          }),
            366);

  // There are 233 rows, of which 50 are ==, the rest being <=.
  EXPECT_THAT(problem.constraints, SizeIs(233));
  EXPECT_EQ(
      std::count_if(problem.constraints.begin(), problem.constraints.end(),
                    [](const MiniMipConstraint& constraint) {
                      return constraint.left_hand_side > -kInf &&
                             constraint.right_hand_side < kInf &&
                             constraint.left_hand_side ==
                                 constraint.right_hand_side;
                    }),
      50);
  EXPECT_EQ(
      std::count_if(problem.constraints.begin(), problem.constraints.end(),
                    [](const MiniMipConstraint& constraint) {
                      return constraint.left_hand_side == -kInf &&
                             constraint.right_hand_side < kInf;
                    }),
      183);
}

// This is the exact same test as above, but where the file contents is given in
// a data string.
TEST(MiniMipProblemTest, LoadProblemFromMPSData) {
  constexpr double kInf = std::numeric_limits<double>::infinity();

  // We rely on ortools to provide the core reader, so we assume it is well
  // tested. We therefore only test our own reader by loading 50v-10 from
  // MIPLIB2017, which is a medium-size problem with varying types of variables.
  // https://miplib.zib.de/instance_details_50v-10.html
  ASSERT_OK_AND_ASSIGN(
      const std::string data,
      ReadFileToString("unit_tests/data_structures/test_data/50v-10.mps"));
  ASSERT_OK_AND_ASSIGN(const MiniMipProblem problem,
                       ReadProblemFromMPSData(data));

  EXPECT_EQ(problem.name, "50v-10");
  EXPECT_THAT(problem.hints, IsEmpty());
  EXPECT_EQ(problem.objective_offset, 0.0);

  // Strangely, the MPS format doesn't specify optimization direction.
  // Therefore, we assume minimization as default.
  EXPECT_EQ(problem.is_maximization, false);

  EXPECT_THAT(problem.variables, SizeIs(2013));
  EXPECT_EQ(
      std::count_if(problem.variables.begin(), problem.variables.end(),
                    [](const MiniMipVariable& var) { return var.is_integer; }),
      1647);

  // In this problem, all variables are non-negative and bounded.
  std::vector<double> lower_bounds(problem.variables.size());
  std::transform(problem.variables.begin(), problem.variables.end(),
                 lower_bounds.begin(),
                 [](const MiniMipVariable& var) { return var.lower_bound; });
  EXPECT_THAT(lower_bounds, Each(0.0));
  std::vector<double> upper_bounds(problem.variables.size());
  std::transform(problem.variables.begin(), problem.variables.end(),
                 upper_bounds.begin(),
                 [](const MiniMipVariable& var) { return var.lower_bound; });
  EXPECT_THAT(upper_bounds, Each(Ge(0.0)));
  EXPECT_THAT(upper_bounds, Each(Le(kInf)));

  // There are 233 rows, of which 50 are ==, the rest being <=.
  EXPECT_THAT(problem.constraints, SizeIs(233));
  EXPECT_EQ(
      std::count_if(problem.constraints.begin(), problem.constraints.end(),
                    [](const MiniMipConstraint& constraint) {
                      return constraint.left_hand_side > -kInf &&
                             constraint.right_hand_side < kInf &&
                             constraint.left_hand_side ==
                                 constraint.right_hand_side;
                    }),
      50);
  EXPECT_EQ(
      std::count_if(problem.constraints.begin(), problem.constraints.end(),
                    [](const MiniMipConstraint& constraint) {
                      return constraint.left_hand_side == -kInf &&
                             constraint.right_hand_side < kInf;
                    }),
      183);
}

// TODO(CGraczyk): add tests for FindErrorInMiniMipProblem.

}  // namespace minimip
