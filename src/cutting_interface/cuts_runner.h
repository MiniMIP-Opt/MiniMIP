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

#ifndef SRC_CUTTING_INTERFACE_CUTS_RUNNER_H_
#define SRC_CUTTING_INTERFACE_CUTS_RUNNER_H_

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "ortools/base/logging.h"
#include "ortools/base/status_macros.h"
#include "src/cutting_interface/cuts_selector.h"
#include "src/cutting_interface/cuts_separator.h"

namespace minimip {

// Forward declaration of Solver. This is required to break the circular
// dependency between the solver and the cutting interface.
class Solver;

class CuttingInterface {
 public:
  virtual ~CuttingInterface() = default;

  // TODO: add searchtree etc.
  virtual absl::Status SeparateCurrentLPSolution(
      const Solver& solver, CutStorage& mutable_cut_storage) = 0;

  void AddSeparator(std::unique_ptr<Separator> separator) {
    separators_.push_back(std::move(separator));
  }

  void AddSelector(std::unique_ptr<Selector> selector) {
    selector_ = std::move(selector);
  }

 protected:
  std::vector<std::unique_ptr<Separator>> separators_;
  std::unique_ptr<Selector> selector_;
};

class CutRunner : CuttingInterface {
 public:
  virtual ~CutRunner() = default;

  virtual absl::Status SeparateCurrentLPSolution(
      const Solver& solver, CutStorage& mutable_cut_storage) final;
};

}  // namespace minimip

#endif  // SRC_CUTTING_INTERFACE_CUTS_RUNNER_H_
