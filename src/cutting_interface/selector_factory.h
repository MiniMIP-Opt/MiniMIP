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

#ifndef SRC_CUTTING_INTERFACE_SELECTOR_FACTORY_H_
#define SRC_CUTTING_INTERFACE_SELECTOR_FACTORY_H_

#include "src/cutting_interface/cuts_selector.h"
#include "src/cutting_interface/hybrid_selector.h"

namespace minimip {

inline absl::StatusOr<std::unique_ptr<CutSelectorInterface>>
ConfigureSelectorFromProto(const CutSelectorParameters& selector_parameters) {
  if (selector_parameters.has_hybrid_selector_parameters()) {
    return std::make_unique<HybridSelector>(selector_parameters);
  }
  return absl::InvalidArgumentError("No cut-selector-specific parameters set.");
}

}  // namespace minimip

#endif  // SRC_CUTTING_INTERFACE_SELECTOR_FACTORY_H_
