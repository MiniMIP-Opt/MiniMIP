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

#ifndef SRC_CUTTING_INTERFACE_CUT_SELECTORS_DYNAMIC_ORTHOGONALITY_H_
#define SRC_CUTTING_INTERFACE_CUT_SELECTORS_DYNAMIC_ORTHOGONALITY_H_

#include "src/cutting_interface/selector.h"

namespace minimip {

class DynamicOrthogonality : public Selector {
  // TODO: Implement Dynamic Orthogonality Filtering based cut selector.

 public:
  DynamicOrthogonality();

  ~DynamicOrthogonality() override;
 protected:
  // Extract the relevant data needed for the following cut selection loop.
  // This must include some sort of pre-scoring notion, refreshing the current score of a cut.
  // TODO: add "isCutFresh()" like function corresponding to its current_score.
  absl::Status PrepareSelection(const LPInterface *lpi,
                                const MipData &mip_data,
                                const std::vector<CuttingPlane> &cuts) const final;

  // Prepare the next iteration of filtering cuts and any preprocessing needed
  absl::Status PrepareIteration(std::vector<CuttingPlane>& cuts) final;

  // Computing the cutting plane from the current iterative and any additional data needed.
  absl::Status Filtering(std::vector<CuttingPlane>& cuts) final;

  // Compute new scores after filtering if necessary.
  absl::Status Scoring(std::vector<CuttingPlane>& cuts) final;
};

}  // namespace minimip
#endif  // SRC_CUTTING_INTERFACE_CUT_SELECTORS_DYNAMIC_ORTHOGONALITY_H_
