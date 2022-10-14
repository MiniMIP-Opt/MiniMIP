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
#include "src/data_structures/cuts_data.h"
#include "src/data_structures/mip_data.h"
#include "src/lp_interface/lpi.h"

namespace minimip {

class CutGenerator {
 public:
  virtual ~CutGenerator() = default;

  // Indicate if this generator can be run in the current state. This is not a
  // guarantee that `GenerateNextCuttingPlanes` will be able to find any cuts.
  virtual bool MayBeRun(const LPInterface* lpi,
                        const MipData& mip_data) const = 0;

  // Generate up to `max_num_cuts` cutting planes. Returns an empty vector if no
  // cuts could be generated given the current state. If called multiple times,
  // the generator is expected to continue where it left of, if appropriate.
  virtual absl::StatusOr<std::vector<CuttingPlane>> GenerateNextCuttingPlanes(
      const LPInterface* lpi, const MipData& mip_data, int max_num_cuts) = 0;
};

}  // namespace minimip

#endif  // SRC_CUTTING_INTERFACE_CUT_GENERATOR_H_
