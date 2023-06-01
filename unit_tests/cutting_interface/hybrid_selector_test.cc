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

#include "src/cutting_interface/hybrid_selector.h"

#include <algorithm>

#include "gmock/gmock.h"
#include "google/protobuf/text_format.h"
#include "gtest/gtest.h"
#include "src/parameters.pb.h"
#include "src/solver.h"
#include "unit_tests/utils.h"

namespace minimip {
namespace {

using testing::Gt;
using testing::IsEmpty;
using testing::Le;
using testing::Not;

class CutSelectorTest : public ::testing::Test {
 protected:
  void SetUp() final {
    // Here we initialize the cut selector.
    SelectorParameters params;
    cut_selector_ = std::make_unique<HybridSelector>(params);

    std::vector<CutData> cuts;
  }
  static std::pair<MiniMipProblem, SparseRow> CreateCuttingplanes(
      int switch_variable_index) {
    switch (switch_variable_index) {
      case 0:
        // Create simple Cutting planes.
        for (ColIndex col(0); col < mip_data.matrix().num_cols(); ++col) {
          RETURN_IF_ERROR(AddColumn(
              mip_data.matrix().col(col), mip_data.lower_bounds()[col],
              mip_data.upper_bounds()[col], mip_data.objective().value(col),
              mip_data.variable_names()[col]));
        }
        SparseRow objective = solver.mip_data().objective();

        double objective_parallelism =
            row.DotProduct(objective) /
            sqrt(row.DotProduct(row) * objective.DotProduct(objective));

        bool is_forced = false;
        if (row.entries().size() == 1) {
          is_forced = true;
        }

        int number_of_non_zeros = row.entries().size();
        int number_of_integer_variables = row.entries().size();
        double efficacy =
            (row.DotProduct(lp_optimum) - rounded_row->right_hand_side) /
            row.DotProduct(row);

        CutData cut(std::move(rounded_row->variable_coefficients),
                    rounded_row->right_hand_side, number_of_non_zeros,
                    number_of_integer_variables, objective_parallelism,
                    efficacy, cutname, is_forced);
        cutting_planes.push_back(std::move(cut));
    }

    std::unique_ptr<HybridSelector> cut_selector_;
  };

}  // namespace
}  // namespace minimip
