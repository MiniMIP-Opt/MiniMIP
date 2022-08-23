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

namespace minimip {

class CutGenerator {
  // TODO: Implement abstract cut generator class to allow the implementation
  //       of different cut generators for the solver.

 public:
  virtual ~CutGenerator() = default;

  virtual absl::Status MyCutGeneratorFunction() = 0;
};

}  // namespace minimip

#endif  // SRC_CUTTING_INTERFACE_CUT_GENERATOR_H_
