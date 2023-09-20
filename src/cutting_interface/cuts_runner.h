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

#ifndef SRC_CUTTING_INTERFACE_CUTS_RUNNER_H_
#define SRC_CUTTING_INTERFACE_CUTS_RUNNER_H_

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "ortools/base/logging.h"
#include "ortools/base/status_macros.h"
#include "src/cutting_interface/cuts_generator.h"
#include "src/cutting_interface/cuts_selector.h"
#include "src/parameters.pb.h"

namespace minimip {

// Forward declaration of Solver. This is required to break the circular
// dependency between the solver and the cutting interface.
class Solver;

class CuttingInterface {
 public:
  virtual ~CuttingInterface() = default;

  // TODO(CG): add searchtree etc.
  virtual absl::Status SeparateCurrentLPSolution(
      const Solver& solver, LPInterface* mutable_lpi,
      CutRegistry& mutable_cut_registry) = 0;

  virtual bool CutCondition(const Solver& solver) = 0;

  void AddGenerator(std::unique_ptr<CutGenerator> generator) {
    generators_.push_back(std::move(generator));
  }

  void AddSelector(std::unique_ptr<CutSelector> selector) {
    selector_ = std::move(selector);
  }

 protected:
  std::vector<std::unique_ptr<CutGenerator>> generators_;
  std::unique_ptr<CutSelector> selector_;
};

class CutRunner : public CuttingInterface {
 public:
  explicit CutRunner(const RunnerParameters& params)
      : params_(params.default_runner_parameters()) {}

  bool CutCondition(const Solver& solver) final;

  absl::Status SeparateCurrentLPSolution(
      const Solver& solver, LPInterface* mutable_lpi,
      CutRegistry& mutable_cut_registry) final;

 private:
  const DefaultRunnerParameters& params_;
};

}  // namespace minimip

#endif  // SRC_CUTTING_INTERFACE_CUTS_RUNNER_H_
