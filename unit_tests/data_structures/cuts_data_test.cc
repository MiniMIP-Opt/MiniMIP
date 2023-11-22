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

#include "src/data_structures/cuts_data.h"

#include "gmock/gmock.h"
#include "google/protobuf/text_format.h"
#include "gtest/gtest.h"
#include "src/data_structures/mip_data.h"
#include "unit_tests/utils.h"

namespace minimip {

TEST(CutDataTests, CreateCut) {
  // Create a dummy problem where all variables are non-negative integers.
  // The problem consists of only the variables and their objective
  // coefficients. The mip_data object is created to access the relevant fields
  // for the cut creation.

  // Note: The objective function is not explicitly declared in the problem.
  //       Instead, we declare the objective coefficients of the variables.
  //       This is a valid way to encode the objective function in a MIP
  //       problem.

  CutRegistry cut_registry;
  MiniMipProblem problem;

  problem.variables.push_back(MiniMipVariable{.name = "x1",
                                              .objective_coefficient = 1.0,
                                              .lower_bound = 0,
                                              .upper_bound = kInf,
                                              .is_integer = true});
  problem.variables.push_back(MiniMipVariable{.name = "x2",
                                              .objective_coefficient = 0.0,
                                              .lower_bound = 0,
                                              .upper_bound = kInf,
                                              .is_integer = true});

  MipData mip_data(problem);

  // Create a cutting plane a^Tx <= b, a = (0,1) and b = 1, x^*_lp = (0,2).
  CutData cut =
      cut_registry.CreateCut(mip_data, CreateSparseRow({{0, 0.0}, {1, 2.0}}),
                             CreateSparseRow({{1, 1.0}}), 1.0, "cut");

  EXPECT_EQ(cut.name(), "cut");
  EXPECT_EQ(cut.number_of_non_zeros(), 1.0);
  EXPECT_EQ(cut.number_of_integer_variables(), 1.0);
  EXPECT_EQ(cut.right_hand_side(), 1.0);
  EXPECT_EQ(cut.efficacy(), 1.0);
  EXPECT_EQ(cut.objective_parallelism(), 0.0);
  EXPECT_EQ(cut.is_forced(), false);
  EXPECT_EQ(cut.is_active(), false);
}

TEST(CutDataTests, CreateDifferentCut) {
  // Create a dummy problem where all variables are non-negative integers.
  // The problem consists of only the variables and their objective
  // coefficients. The mip_data object is created to access the relevant fields
  // for the cut creation.

  // Note: The objective function is not explicitly declared in the problem.
  //       Instead, we declare the objective coefficients of the variables.
  //       This is a valid way to encode the objective function in a MIP
  //       problem.

  CutRegistry cut_registry;

  MiniMipProblem problem;

  problem.variables.push_back(MiniMipVariable{.name = "x1",
                                              .objective_coefficient = 1.0,
                                              .lower_bound = 0,
                                              .upper_bound = kInf,
                                              .is_integer = true});
  problem.variables.push_back(MiniMipVariable{.name = "x2",
                                              .objective_coefficient = 1.0,
                                              .lower_bound = 0,
                                              .upper_bound = kInf,
                                              .is_integer = false});

  MipData mip_data(problem);

  // Create a cutting plane a^Tx <= b, a = (1,1) and b = 2, x^*_lp = (2,2).
  CutData cut = cut_registry.CreateCut(
      mip_data, CreateSparseRow({{0, 2.0}, {1, 2.0}}),
      CreateSparseRow({{0, 1.0}, {1, 1.0}}), 2.0, "cut", true);

  EXPECT_EQ(cut.name(), "cut");
  EXPECT_EQ(cut.number_of_non_zeros(), 2.0);
  EXPECT_EQ(cut.number_of_integer_variables(), 1.0);
  EXPECT_EQ(cut.right_hand_side(), 2.0);
  EXPECT_NEAR(cut.efficacy(), sqrt(2), 1e-9);
  EXPECT_EQ(cut.objective_parallelism(), 1.0);
  EXPECT_EQ(cut.is_forced(), true);
  EXPECT_EQ(cut.is_active(), false);
}

TEST(CutRegistryTest, InitTest) {
  CutRegistry cut_registry;

  EXPECT_EQ(cut_registry.cuts().size(), 0);
  EXPECT_EQ(cut_registry.active_cuts().size(), 0);
}

TEST(CutRegistryTest, AddCuts) {
  CutRegistry cut_registry;
  const std::vector<CutData>& cuts = cut_registry.cuts();

  MiniMipProblem problem;

  problem.variables.push_back(MiniMipVariable{.name = "x1",
                                              .objective_coefficient = 1.0,
                                              .lower_bound = 0,
                                              .upper_bound = kInf,
                                              .is_integer = true});
  problem.variables.push_back(MiniMipVariable{.name = "x2",
                                              .objective_coefficient = 0.0,
                                              .lower_bound = 0,
                                              .upper_bound = kInf,
                                              .is_integer = true});
  MipData mip_data(problem);

  EXPECT_EQ(cuts.size(), 0);

  // Create a cutting plane a^Tx <= b, a = (0,1) and b = 1, x^*_lp = (0,2).
  cut_registry.AddCut(
      cut_registry.CreateCut(mip_data, CreateSparseRow({{0, 0.0}, {1, 2.0}}),
                             CreateSparseRow({{1, 1.0}}), 1.0, "cut1"));

  EXPECT_EQ(cuts.size(), 1);
  EXPECT_EQ(cuts[0].index(), 0);  // Check the index of the first cut.

  // Create a cutting plane a^Tx <= b, a = (0,1) and b = 1.5, x^*_lp = (0,2).
  cut_registry.AddCut(
      cut_registry.CreateCut(mip_data, CreateSparseRow({{0, 0.0}, {1, 2.0}}),
                             CreateSparseRow({{1, 1.0}}), 1.5, "cut2"));

  EXPECT_EQ(cuts.size(), 2);
  EXPECT_EQ(cuts[0].index(), 0);
  EXPECT_EQ(cuts[0].name(), "cut1");
  EXPECT_EQ(cuts[1].index(), 1);
  EXPECT_EQ(cuts[1].name(), "cut2");
}

TEST(CutRegistryTest, ActivateAndRemoveCuts) {
  CutRegistry cut_registry;
  const std::vector<CutData>& cuts = cut_registry.cuts();
  const std::vector<int>& active_cuts = cut_registry.active_cuts();

  MiniMipProblem problem;

  problem.variables.push_back(MiniMipVariable{.name = "x1",
                                              .objective_coefficient = 1.0,
                                              .lower_bound = 0,
                                              .upper_bound = kInf,
                                              .is_integer = true});
  problem.variables.push_back(MiniMipVariable{.name = "x2",
                                              .objective_coefficient = 0.0,
                                              .lower_bound = 0,
                                              .upper_bound = kInf,
                                              .is_integer = true});

  MipData mip_data(problem);

  // Create a cutting plane a^Tx <= b, a = (0,1) and b = 1, x^*_lp = (0,2).
  cut_registry.AddCut(
      cut_registry.CreateCut(mip_data, CreateSparseRow({{0, 0.0}, {1, 2.0}}),
                             CreateSparseRow({{1, 1.0}}), 1.0, "cut1"));

  // Create a cutting plane a^Tx <= b, a = (0,1) and b = 1.5, x^*_lp = (0,2).
  cut_registry.AddCut(
      cut_registry.CreateCut(mip_data, CreateSparseRow({{0, 0.0}, {1, 2.0}}),
                             CreateSparseRow({{1, 1.0}}), 1.5, "cut2"));

  EXPECT_EQ(active_cuts.size(), 0);
  cut_registry.ActivateCut(1);
  EXPECT_EQ(active_cuts.size(), 1);
  EXPECT_EQ(active_cuts[0], 1);

  cut_registry.RemoveCut(0);
  EXPECT_EQ(active_cuts[0], 0);
  EXPECT_EQ(cuts.size(), 1);
  EXPECT_EQ(cuts[0].name(), "cut2");
  EXPECT_EQ(cuts[0].index(), 0);
}

TEST(CutRegistryTest, UpdateCutEfficacy) {
  CutRegistry cut_registry;
  const std::vector<CutData>& cuts = cut_registry.cuts();

  MiniMipProblem problem;

  problem.variables.push_back(MiniMipVariable{.name = "x1",
                                              .objective_coefficient = 1.0,
                                              .lower_bound = 0,
                                              .upper_bound = kInf,
                                              .is_integer = true});
  problem.variables.push_back(MiniMipVariable{.name = "x2",
                                              .objective_coefficient = 0.0,
                                              .lower_bound = 0,
                                              .upper_bound = kInf,
                                              .is_integer = true});

  MipData mip_data(problem);

  // Create a cutting plane a^Tx <= b, a = (0,1) and b = 1.5, x^*_lp = (0,2).
  cut_registry.AddCut(
      cut_registry.CreateCut(mip_data, CreateSparseRow({{0, 0.0}, {1, 2.0}}),
                             CreateSparseRow({{1, 1.0}}), 1.5, "cut"));

  cut_registry.UpdateCutEfficacy(CreateSparseRow({{0, 1.0}, {1, 1.0}}));
  CutData cut = cut_registry.GetCut(0);

  EXPECT_EQ(cut.efficacy(), -0.5);
  cut_registry.RemoveCut(0);
  EXPECT_EQ(cuts.size(), 0);
}

}  // namespace minimip
