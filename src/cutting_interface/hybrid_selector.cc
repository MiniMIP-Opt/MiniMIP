// Copyright 2022 the MiniMIP Project
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

#include "hybrid_selector.h"

#include <limits>

#include "cuts_selector.h"
#include "cuts_separator.h"
#include "src/solver.h"

namespace minimip {

namespace {

void scoring_function(const MiniMipSolver& solver,
                      HybridSelectorParameters params, CutData& cut) {
  // TODO: implement incumbent solution for directed cutoff distance.

  double efficacy_weight =
      params.efficacy_weight() + params.directed_cutoff_distance_weight();
  double integer_support_weight = params.integer_support_weight();
  double objective_parallelism_weight = params.objective_parallelism_weight();
  double directed_cutoff_distance_weight =
      params.directed_cutoff_distance_weight();

  double efficacy =
      (efficacy_weight > 0) ? efficacy_weight * cut.efficacy : 0.0;

  double integer_support =
      (integer_support_weight > 0)
          ? integer_support_weight *
                (cut.number_of_integer_variables / cut.number_of_non_zeros)
          : 0.0;
  double objective_parallelism =
      (objective_parallelism_weight > 0)
          ? objective_parallelism_weight * cut.objective_parallelism
          : 0.0;

  cut.current_score = efficacy + integer_support + objective_parallelism;
}

bool compute_row_parallelism(CutData& cut_reference, CutData& cut, double maximum_parallelism,
                             bool signed_orthogonality = false) {

  std::vector<double> reference_row_values = cut_reference.row.values();
  std::vector<double> cut_row_values = cut.row.values();

  double squared_norm_reference = std::inner_product(reference_row_values.begin(), reference_row_values.end(), reference_row_values.begin(), 0.0);
  double squared_norm_cut = std::inner_product(cut_row_values.begin(), cut_row_values.end(), cut_row_values.begin(), 0.0);

  double dot_product = cut_reference.row.DotProduct(cut.row);
  double cos_angle = dot_product / (std::sqrt(squared_norm_reference) * std::sqrt(squared_norm_cut));

  // TODO: check unequal sign for signed_orthogonality
  return signed_orthogonality ? std::abs(cos_angle) > maximum_parallelism : cos_angle > maximum_parallelism;
}

int select_best_cut(const MiniMipSolver& solver, std::vector<CutData>& cuts) {
  double max_score = std::numeric_limits<double>::lowest();
  int best_cut_index = -1;
  for (CutData cut : cuts) {
    if (cut.current_score > max_score) {
      max_score = cut.current_score;
      best_cut_index = cut.cut_index;
    }
  }
  return best_cut_index;
}

std::vector<CutData> filter_cuts(const MiniMipSolver& solver,
                 HybridSelectorParameters params, CutData& cut_reference,
                 std::vector<CutData>& cuts) {
  const double parallel_cutoff = 1.0 - params.minimum_orthogonality();

  std::vector<CutData> filtered_cuts;

  filtered_cuts.push_back(cut_reference);

  for (CutData cut : cuts) {
    bool is_parallel = compute_row_parallelism(cut_reference, cut, parallel_cutoff,
                                               params.signed_orthogonality());
        if(!is_parallel) {
          filtered_cuts.push_back(cut);
        }
  }

  return filtered_cuts;
}


}  // namespace

absl::StatusOr<std::vector<CutData>> HybridSelector::SelectCuttingPlanes(
    const MiniMipSolver& solver, std::vector<CutData>& cuts) {
  const int max_cuts = params_.max_num_cuts();

  // todo: implement the cutselector for hybrid selection

}

}  // namespace minimip
