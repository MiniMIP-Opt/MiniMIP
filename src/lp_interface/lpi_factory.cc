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

#include "src/lp_interface/lpi_factory.h"

#include <memory>

#include "ortools/base/logging.h"
#include "src/lp_interface/lpi_glop.h"
#include "src/lp_interface/lpi_soplex.h"

namespace minimip {

std::unique_ptr<LPInterface> CreateLPInterface(LPInterfaceCode interface_code) {
  switch (interface_code) {
    case LPInterfaceCode::kGlop:
      return std::make_unique<LPGlopInterface>();
    case LPInterfaceCode::kSoplex:
      return std::make_unique<LPSoplexInterface>();
    default:
      LOG(FATAL) << "Unsupported interface code " << interface_code;
  }
}

std::ostream& operator<<(std::ostream& os, const LPInterfaceCode& code) {
  os << static_cast<int>(code);
  return os;
}

}  // namespace minimip
