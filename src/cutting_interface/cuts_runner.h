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

// Forward declaration of MiniMipSolver. This is required to break the circular
// dependency between the solver and the cutting interface.
class MiniMipSolver;

class CuttingInterface {
 public:
  virtual ~CuttingInterface() = default;

  // TODO: add searchtree etc.
  virtual absl::Status SeparateCurrentLPSolution(
      const MiniMipSolver& solver) = 0;

  void AddSeparator(std::unique_ptr<Separator> separator) {
    separators_.push_back(std::move(separator));
  }

 protected:
  std::vector<std::unique_ptr<Separator>> separators_;
};

class CutRunner : CuttingInterface {
  // TODO: Implement abstract cut runner class to allow the implementation of
  //       different cut runners to control the hierarchy and application of cut
  //       generators.
 public:
  virtual ~CutRunner() = default;

  absl::Status SeparateCurrentLPSolution(const MiniMipSolver& solver) final {
    ASSIGN_OR_RETURN(int max_num_cuts, PrepareSeparation(solver));

    ASSIGN_OR_RETURN(std::vector<CutData> cuts,
                     separator_->GenerateCuttingPlanes(solver));

    RETURN_IF_ERROR(PrepareSelection(solver, cuts));

    ASSIGN_OR_RETURN(int number_of_cuts_selected,
                     selector_->SelectCuttingPlanes(solver, cuts));

    RETURN_IF_ERROR(StoreCutsInArchive(solver));
  };

 private:
  Separator* separator_;
  Selector* selector_;

 protected:
  virtual absl::StatusOr<int> PrepareSeparation(
      const MiniMipSolver& solver) = 0;

  virtual absl::Status PrepareSelection(const MiniMipSolver& solver,
                                        std::vector<CutData> cuts) = 0;

  virtual absl::Status StoreCutsInArchive(const MiniMipSolver& solver) = 0;
};

}  // namespace minimip

#endif  // SRC_CUTTING_INTERFACE_CUTS_RUNNER_H_
