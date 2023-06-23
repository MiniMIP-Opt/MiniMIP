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

double ComputeCutScore(const HybridSelectorParameters& params, CutData& cut) {
  double efficacy = (params.efficacy_weight() > 0)
                        ? params.efficacy_weight() * cut.efficacy()
                        : 0.0;

  double integer_support =
      (params.integer_support_weight() > 0)
          ? params.integer_support_weight() *
                (cut.number_of_integer_variables() / cut.number_of_non_zeros())
          : 0.0;

  double objective_parallelism =
      (params.objective_parallelism_weight() > 0)
          ? params.objective_parallelism_weight() * cut.objective_parallelism()
          : 0.0;

  return efficacy + integer_support + objective_parallelism;
}

bool TwoCutsAreOrthogonal(const CutData& cut_reference, const CutData& cut,
                          const double minimum_orthogonality,
                          const bool signed_orthogonality = false) {
  if (cut.is_forced()) {
    return false;
  }
  std::vector<double> reference_row_values = cut_reference.row().values();
  std::vector<double> cut_row_values = cut.row().values();

  const double squared_norm_reference = std::inner_product(
      reference_row_values.begin(), reference_row_values.end(),
      reference_row_values.begin(), 0.0);
  const double squared_norm_cut =
      std::inner_product(cut_row_values.begin(), cut_row_values.end(),
                         cut_row_values.begin(), 0.0);

  const double cos_angle =
      cut_reference.row().DotProduct(cut.row()) /
      (std::sqrt(squared_norm_reference) * std::sqrt(squared_norm_cut));

  // Return false if rows are insufficiently orthogonal.
  return signed_orthogonality ? cos_angle > 1 - minimum_orthogonality
                              : std::abs(cos_angle) > 1 - minimum_orthogonality;
}

}  // namespace

absl::StatusOr<std::vector<CutData>> HybridSelector::SelectCuttingPlanes(
    const Solver& solver, std::vector<CutData>& cuts) {
  // 1. Compute the score for each cut
  for (CutData& cut : cuts) {
    cut.SetScore(ComputeCutScore(params_, cut));
  }

  // 2. Select the best cut and filter the remaining cuts.
  std::vector<CutData> selected_cuts;

  while (!cuts.empty()) {
    // 2.1. Sort the cuts by score.
    std::sort(cuts.begin(), cuts.end(),
              [](const CutData& cut1, const CutData& cut2) {
                return cut1.score() > cut2.score();
              });

    // The remaining cut with the highest score.
    CutData& cut_reference = cuts[0];

    if (cut_reference.score() < params_.score_threshold()) {
      break;
    }

    // 2.2 Add the current best cut to the selected cuts.
    selected_cuts.push_back(cut_reference);

    if (selected_cuts.size() == max_num_cuts_) {
      break;
    }

    // 3. Filter the cuts
    auto predicate = [cut_reference, this](const CutData& cut) {
      return TwoCutsAreOrthogonal(cut_reference, cut,
                                  params_.minimum_orthogonality(),
                                  params_.signed_orthogonality());
    };

    // Remove all elements that match the predicate.
    cuts.erase(std::remove_if(cuts.begin() + 1, cuts.end(), predicate),
               cuts.end());

    // Remove the first element (since we've already processed it).
    cuts.erase(cuts.begin());
  }

  return selected_cuts;
}

}  // namespace minimip
