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

#ifndef SRC_CUTTING_INTERFACE_SEPARATOR_H_
#define SRC_CUTTING_INTERFACE_SEPARATOR_H_

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "ortools/base/status_macros.h"
#include "ortools/base/logging.h"

#include "src/data_structures/cuts_data.h"
#include "src/data_structures/mip_data.h"
#include "src/lp_interface/lpi.h"

namespace minimip {

class SeparatorInterface{
 public:
  // Generate up to `max_num_cuts` cutting planes. Returns an empty vector if no
  // cuts could be generated given the current state. If called multiple times,
  // the generator is expected to continue where it left of, if appropriate.
  virtual absl::StatusOr<std::vector<CuttingPlane>> GenerateCuttingPlanes(
      const LPInterface* lpi, const MipData& mip_data, int max_num_cuts) = 0;
};

class Separator : SeparatorInterface {
 public:
  virtual ~Separator() = default;

  // Generate up to `max_num_cuts` cutting planes. Returns an empty vector if no
  // cuts could be generated given the current state. If called multiple times,
  // the generator is expected to continue where it left of, if appropriate.
  absl::StatusOr<std::vector<CuttingPlane>> GenerateCuttingPlanes(const LPInterface* lpi, const MipData& mip_data, int max_num_cuts) final {
    RETURN_IF_ERROR( PrepareSeparationRound(lpi, mip_data, max_num_cuts) );

    std::vector<CuttingPlane> cuts;

    while (MayBeRun(max_num_cuts)) {
      RETURN_IF_ERROR(PrepareIteration());

      ASSIGN_OR_RETURN(std::vector<CuttingPlane> current_cuts, ComputeCuts());

      RETURN_IF_ERROR(PostProcessing(current_cuts));

      cuts.insert(cuts.end(), current_cuts.begin(), current_cuts.end());
    }

    RETURN_IF_ERROR(LocallyFilter(cuts, max_num_cuts));
    return cuts;
  };

 protected:
  // Extract the relevant data needed for the cut generation loop iteratives.
  // Including the MipData address and the LPInterface into the Separator Data
  // is fine if behaviour-dependent information is required for the generator.
  // If this function is not implemented, it simply passes the MipData address
  // as well as the LPI Pointer to the SeparatorData.
  virtual absl::Status PrepareSeparationRound(
      const LPInterface* lpi, const MipData& mip_data, int max_num_cuts) = 0;

  // Indicate if this generator can be run in the current state. This is not a
  // guarantee that `GenerateNextCuttingPlanes` will be able to find any cuts.
  // TODO: Make sure the while loop is easy to exit.
  virtual bool MayBeRun(int max_num_cuts) const {
    return true;
  }

  // Modify the iterative used to generate the next cut, e.g. tightening constraints,
  // modifying coefficients etc. If this function is not implemented, it simply
  // returns the iterative.
  virtual absl::Status PrepareIteration() = 0;

  // Computing the cutting plane from the current iterative and any additional data needed.
  virtual absl::StatusOr<std::vector<CuttingPlane>> ComputeCuts() = 0;

  // Modifying the newly extracted cutting plane, e.g. strengthening or sparsification.
  // If this function is not implemented, it simply returns the cutting plane.
  virtual absl::Status PostProcessing(std::vector<CuttingPlane>& current_cuts) {
    return absl::OkStatus();
  }

  // Apply a local selection heuristic to the generator specific cuts of the separation round.
  // If this function is not implemented, it simply returns the vector of cuts.
  virtual absl::Status LocallyFilter(std::vector<CuttingPlane>& cuts, int max_num_cuts) const {
    return absl::OkStatus();
    }

};

}  // namespace minimip

#endif  // SRC_CUTTING_INTERFACE_SEPARATOR_H_
