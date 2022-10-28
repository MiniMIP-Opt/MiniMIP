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

#ifndef SRC_CUTTING_INTERFACE_CUT_GENERATOR_H_
#define SRC_CUTTING_INTERFACE_CUT_GENERATOR_H_

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "ortools/base/logging.h"

#include "src/data_structures/cuts_data.h"
#include "src/data_structures/mip_data.h"
#include "src/lp_interface/lpi.h"

namespace minimip {

// This datastructure should carry the data that is used by the generator to
// iterate over or computes an iterative from. Can store the MipData adress and
// LPInterface pointer by default.
class SeparatorData{
 public:
  virtual ~SeparatorData() = default;
  MipData& mip_data = null_ptr;
  LPInterface* lpi = null_ptr;
};

// The object that is used by the generator to compute the next cutting plane.
class Iterative{
 public:
  virtual ~Iterative() = default;
};

class CutGenerator {
 public:
  virtual ~CutGenerator() = default;

  // Indicate if this generator can be run in the current state. This is not a
  // guarantee that `GenerateNextCuttingPlanes` will be able to find any cuts.
  virtual bool MayBeRun(const LPInterface* lpi,
                        const MipData& mip_data) const = 0;

  // Extract the relevant data needed for the cut generation loop iteratives.
  // Including the MipData adress and the LPInterface into the SepaData is fine
  // if behaviour-dependent information is required for the generator.
  // If this function is not implemented, it simply passes the MipData address
  // as well as the LPI Pointer to the SeparatorData.
  virtual absl::StatusOr<SeparatorData> InitializeSeparatorData(
    const LPInterface* lpi, const MipData& mip_data, int max_num_cuts) = 0;

  // Acts both as the while condition and the modification of the SepaData
  // to allow refinement of information used to compute iteratives in the
  // generation loop.
  // If this function is not implemented, it simply returns true.
  virtual bool PrepareNextIteration(SeparatorData SepaData) = 0;

  // Construct the Iterative object from the SepaData that the generator acts on.
  // Usually this would be a constraint or a column, depending on the generator.
  virtual absl::StatusOr<Iterative> GetItFrom(const SeparatorData& SepaData) = 0;

  // Modify the iterative used to generate the next cut, e.g. tightening constraints,
  // modifying coefficients etc. If this function is not implemented, it simply
  // returns the iterative.
  virtual absl::StatusOr<Iterative> PreProcessing(Iterative It, const SeparatorData& SepaData) = 0;

  // Computing the cutting plane from the current iterative and any additinal data needed.
  virtual absl::StatusOr<CuttingPlane> ComputeCutFrom(Iterative It, const SeparatorData& SepaData) = 0;

  // Modifying the newly extracted cutting plane, e.g. strengthening or sparsification.
  // If this function is not implemented, it simply returns the cutting plane.
  virtual absl::StatusOr<CuttingPlane> PostProcessing(CuttingPlane cut, const SeparatorData& SepaData) = 0;

  // Apply a local selection heuristic to the generator specific cuts of the separation round.
  // If this function is not implemented, it simply returns the vector of cuts.
  virtual absl::StatusOr<std::vector<CuttingPlane>> LocalSelect(std::vector<CuttingPlane> cuts, int max_num_cuts) = 0;

  // Generate up to `max_num_cuts` cutting planes. Returns an empty vector if no
  // cuts could be generated given the current state. If called multiple times,
  // the generator is expected to continue where it left of, if appropriate.
  absl::StatusOr<std::vector<CuttingPlane>> GenerateCuttingPlanes(
      const LPInterface* lpi, const MipData& mip_data, int max_num_cuts)
      {
        std::vector<CuttingPlane> cuts;
        cuts.reserve(max_num_cuts);

        absl::StatusOr<SeparatorData> SepaData = InitializeSeparatorData(lpi, mip_data, max_num_cuts);
        if(SepaData.ok()) {
          while( PrepareNextIteration(*SepaData) ) {
            absl::StatusOr<Iterative> It = GetItFrom(*SepaData);

            if(It.ok())
              It = PreProcessing(*It, *SepaData);

            if(!It.ok())
              return It.status();

            absl::StatusOr<CuttingPlane> cut = ComputeCutFrom(*It, *SepaData);

            if(cut.ok())
              cut = PostProcessing(*cut, *SepaData);

            if(!cut.ok())
              return cut.status();

            cuts.push_back(*cut);
          }
          return LocalSelect(cuts,max_num_cuts);
        }
        else
          return SepaData.status();
      };
};

}  // namespace minimip

#endif  // SRC_CUTTING_INTERFACE_CUT_GENERATOR_H_
