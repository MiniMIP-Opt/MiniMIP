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

#ifndef MINIMIP_CUTTING_INTERFACE_HYBRID_SELECTOR_H_
#define MINIMIP_CUTTING_INTERFACE_HYBRID_SELECTOR_H_

#include <iostream>
#include <optional>

#include "minimip/cutting_interface/cuts_selector.h"
#include "minimip/solver_context_interface.h"
#include "ortools/base/status_macros.h"

// TODO(Cgraczyk): Add a check that duplicates are removed.
// The cutselector should not be responsible for removing duplicates.
// Maybe this can be avoided by using std::set instead of std::vector in the
// future.

// The assumption here is that all cuts are initialized correctly and are valid.
// Note:
// - make sure that cuts are updated if returned from the cut_registry.

namespace minimip {

class HybridSelector : public CutSelectorInterface {
 public:
  explicit HybridSelector(CutSelectorParameters params)
      : params_(std::move(params.hybrid_selector_parameters())) {}

  absl::StatusOr<std::vector<CutData>> SelectCuttingPlanes(
      const SolverContextInterface& context, std::vector<CutData>& cuts) final;

 private:
  const HybridSelectorParameters params_;
};

}  // namespace minimip

#endif  // MINIMIP_CUTTING_INTERFACE_HYBRID_SELECTOR_H_
