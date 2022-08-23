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

#ifndef SRC_CUTTING_INTERFACE_CUT_GENERATORS_GOMORY_MIXED_INTEGER_H_
#define SRC_CUTTING_INTERFACE_CUT_GENERATORS_GOMORY_MIXED_INTEGER_H_

#include "src/cutting_interface/cut_generator.h"

namespace minimip {

class GomoryMixedInteger : public CutGenerator {
  // TODO: Implement Gomory Mixed-Integer cutting plane generator.

 public:
  GomoryMixedInteger();

  ~GomoryMixedInteger() override;

  absl::Status MyCutGeneratorFunction() final;
};

}  // namespace minimip

#endif  // SRC_CUTTING_INTERFACE_CUT_GENERATORS_GOMORY_MIXED_INTEGER_H_
