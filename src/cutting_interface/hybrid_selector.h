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

#ifndef MINIMIP_SRC_CUTTING_INTERFACE_HYBRID_SELECTOR_H_
#define MINIMIP_SRC_CUTTING_INTERFACE_HYBRID_SELECTOR_H_

#include <iostream>
#include <optional>

#include "ortools/base/status_macros.h"
#include "src/cutting_interface/cuts_selector.h"
#include "src/solver.h"

namespace minimip {

class HybridSelector : public Selector {
 public:
  explicit HybridSelector(SelectorParameters params)
      : params_(std::move(params)) {}

  absl::StatusOr<std::vector<CutData>> SelectCuttingPlanes(
      const MiniMipSolver& solver, std::vector<CutData>& cuts) final;

 private:
  const SelectorParameters params_;
};

}  // namespace minimip

#endif  // MINIMIP_SRC_CUTTING_INTERFACE_HYBRID_SELECTOR_H_
