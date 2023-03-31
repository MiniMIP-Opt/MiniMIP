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

#include "cuts_runner.h"

namespace minimip {
// Remove duplicate cuts and other preliminary operations to be done before
// calling the cut selector itself.
absl::Status CutRunner::PrepareSelection(const MiniMipSolver& solver,
                                         std::vector<CutData> cuts) {
  // remove duplicate cuts
  for (CutData& cut : cuts) {

  }
}

}  // namespace minimip
