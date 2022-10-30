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

#ifndef SRC_CUTTING_INTERFACE_CUT_RUNNER_H_
#define SRC_CUTTING_INTERFACE_CUT_RUNNER_H_

#include "absl/status/status.h"
#include "absl/status/statusor.h"

#include "ortools/base/status_macros.h"
#include "ortools/base/logging.h"

#include "src/data_structures/cuts_data.h"
#include "src/data_structures/mip_data.h"
#include "src/lp_interface/lpi.h"

namespace minimip {

class CuttingInterface {
 public:
  // TODO: add searchtree etc.
  virtual absl::Status SeparateCurrentLPSolution(const LPInterface *lpi,
                                                 const MipData &mip_data,
                                                 const CutArchive& cut_archive) = 0;
};


class CutRunner : CuttingInterface {
  // TODO: Implement abstract cut runner class to allow the implementation of
  //       different cut runners to control the hierarchy and application of cut
  //       generators.
 public:
  virtual ~CutRunner() = default;

  absl::Status SeparateCurrentLPSolution(const LPInterface *lpi,
                                         const MipData &mip_data,
                                         const CutArchive& cut_archive) final {
    RETURN_IF_ERROR(PrepareSeparation(lpi, mip_data, cut_archive));

    ASSIGN_OR_RETURN(std::vector<CuttingPlane> cuts,
                     separator_->GenerateCuttingPlanes(lpi,
                                                       mip_data,
                                                       cut_archive));

    ASSIGN_OR_RETURN(int number_of_cuts_selected,
                     selector_->SelectCuttingPlanes(lpi,
                                                    mip_data,
                                                    cuts));

    RETURN_IF_ERROR(StoreCutsInArchive(lpi, mip_data, cut_archive));
    RETURN_IF_ERROR(SolveLP());
  };
 private:
  SeparatorInterface separator_;
  SelectorInterface selector_;

 protected:
  virtual absl::Status PrepareSeparation(
      const LPInterface* lpi, const MipData& mip_data, const CutArchive& cut_archive) = 0;

  virtual absl::Status StoreCutsInArchive(
      const LPInterface* lpi, const MipData& mip_data, const CutArchive& cut_archive) = 0;


};

}  // namespace minimip

#endif  // SRC_CUTTING_INTERFACE_CUT_RUNNER_H_
