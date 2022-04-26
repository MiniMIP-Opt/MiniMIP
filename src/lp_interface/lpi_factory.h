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

#ifndef SRC_LP_INTERFACE_LPI_FACTORY_H_
#define SRC_LP_INTERFACE_LPI_FACTORY_H_

#include "src/lp_interface/lpi_glop.h"
#include "src/lp_interface/lpi_soplex.h"

namespace minimip {

enum class InterfaceCode {
  kGlop   = 0,
  kSoplex = 1,
};

class LPInterfaceFactory {
 public:
  LPInterface* CreateLPInterface(InterfaceCode interface_code) {
    switch (interface_code) {
      case InterfaceCode::kSoplex:
        return new LPSoplexInterface();
      default:
        return new LPGlopInterface();
    }
  }
};

}  // namespace minimip
#endif  // SRC_LP_INTERFACE_LPI_FACTORY_H_
