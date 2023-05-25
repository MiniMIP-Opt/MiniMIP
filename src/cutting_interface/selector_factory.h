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

#include "src/cutting_interface/cuts_selector.h"
#include "src/parameters.pb.h"

namespace minimip {

inline absl::StatusOr<std::unique_ptr<Selector>> ConfigureSelectorFromProto(
    const SelectorParameters& selector_parameters) {
  if (selector_parameters.has_some_cut_selector_parameters()) {
    // Here we would create and return the cut runner.
  }
  return absl::InvalidArgumentError("No separator-specific parameters set.");
}

}  // namespace minimip
