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

#include "hybrid_selector.h"

#include <limits>

#include "cuts_selector.h"
#include "cuts_separator.h"
#include "src/solver.h"

namespace minimip {

namespace {

void scoring_function(const MiniMipSolver& solver,
                      const HybridSelectorParameters params, Cut& cut) {
  // TODO: implement incumbent solution for directed cutoff distance.

  double efficacy_weight =
      params.efficacy_weight() + params.directed_cutoff_distance_weight();
  double integer_support_weight = params.integer_support_weight();
  double objective_parallelism_weight = params.objective_parallelism_weight();

  double efficacy =
      (efficacy_weight > 0) ? efficacy_weight * cut.efficacy() : 0.0;

  double integer_support =
      (integer_support_weight > 0)
          ? integer_support_weight *
                (cut.getConstData().number_of_integer_variables /
                 cut.getConstData().number_of_non_zeros)
          : 0.0;
  double objective_parallelism =
      (objective_parallelism_weight > 0)
          ? objective_parallelism_weight *
                cut.getConstData().objective_parallelism
          : 0.0;

  cut.setScore(efficacy + integer_support + objective_parallelism);
}

int select_best_cut(const MiniMipSolver& solver, std::vector<Cut>& cuts) {
  double max_score = std::numeric_limits<double>::lowest();
  int best_cut_index = -1;
  for (Cut cut : cuts) {
    if (cut.score() > max_score) {
      max_score = cut.score();
      best_cut_index = cut.index();
    }
  }
  return best_cut_index;
}

bool compute_row_parallelism(Cut& reference_cut, Cut& cut,
                             double maximum_parallelism,
                             bool signed_orthogonality = false) {
  const SparseRow& reference_cut_row = reference_cut.getConstData().row;
  const SparseRow& cut_row = cut.getConstData().row;
  std::vector<double> reference_row_values = reference_cut_row.values();
  std::vector<double> cut_row_values = cut_row.values();

  double squared_norm_reference = std::inner_product(
      reference_row_values.begin(), reference_row_values.end(),
      reference_row_values.begin(), 0.0);
  double squared_norm_cut =
      std::inner_product(cut_row_values.begin(), cut_row_values.end(),
                         cut_row_values.begin(), 0.0);

  double dot_product = reference_cut_row.DotProduct(cut_row);
  double cos_angle = dot_product / (std::sqrt(squared_norm_reference) *
                                    std::sqrt(squared_norm_cut));

  // TODO: check unequal sign for signed_orthogonality
  return signed_orthogonality ? std::abs(cos_angle) > maximum_parallelism
                              : cos_angle > maximum_parallelism;
}

std::vector<Cut> filter_cuts(const MiniMipSolver& solver,
                             HybridSelectorParameters params,
                             Cut& cut_reference, std::vector<Cut>& cuts) {
  const double parallel_cutoff = 1.0 - params.minimum_orthogonality();

  // Todo: fix memory leak
  std::vector<Cut> filtered_cuts;

  filtered_cuts.push_back(cut_reference);

  for (Cut cut : cuts) {
    bool is_parallel = compute_row_parallelism(
        cut_reference, cut, parallel_cutoff, params.signed_orthogonality());
    if (!is_parallel) {
      filtered_cuts.push_back(cut);
    }
  }

  return filtered_cuts;
}

}  // namespace

absl::StatusOr<std::vector<Cut>> HybridSelector::SelectCuttingPlanes(
    const MiniMipSolver& solver, std::vector<Cut>& cuts) {
  const int max_cuts = params_.max_num_cuts();

  // 1. compute the score for each cut
  for (Cut& cut : cuts) {
    scoring_function(solver, params_.hybrid_selector_parameters(), cut);
  }

  int selected_cuts = 0;

  while (!cuts.empty()) {
    // 2. select the best cut
    int best_cut_index = select_best_cut(solver, cuts);

    selected_cuts++;

    if (cuts[best_cut_index].score() <
        params_.hybrid_selector_parameters().score_threshold()) {
      break;
    }

    // Todo: fix memory leak
    //  3. filter the cuts
    std::vector<Cut> filtered_cuts =
        filter_cuts(solver, params_.hybrid_selector_parameters(),
                    cuts[best_cut_index], cuts);
    cuts = std::move(filtered_cuts);

    if (selected_cuts == max_cuts) {
      break;
    }
  }

  return cuts;
}

}  // namespace minimip

// TODO: refactor the whole cutting interface, remove the runner and instead
// implement
//  the generator and incorporate the filtering into the cut_storage.
//  The generator should be able to generate cuts and store them in the
//  cut_storage. The cut_storage should be able to filter the cuts and return
//  the best cuts. The runner should be able to run the generator and apply the
//  cut_storage to the LP.
