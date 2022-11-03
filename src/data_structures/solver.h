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

#ifndef SRC_DATA_STRUCTURES_SOLVER_H_
#define SRC_DATA_STRUCTURES_SOLVER_H_

#include "mip_data.h"
#include "cuts_data.h"
#include "src/lp_interface/lpi.h"

namespace minimip {

class Solver {
 public:
  explicit Solver(MipData &mip_data,
                  CutStorage& cut_storage,
                  LPInterface *lp_interface) :
      mip_data_(mip_data),
      cut_storage_(cut_storage),
      lp_interface_(lp_interface) {}

  const MipData &mip_data() const { return mip_data_; };
  const CutStorage &cut_storage() const { return cut_storage_; };
  const LPInterface *lp_interface() const { return lp_interface_; };

 private:
  MipData &mip_data_;
  CutStorage &cut_storage_;
  LPInterface *lp_interface_;

};

}  // namespace minimip
#endif //SRC_DATA_STRUCTURES_SOLVER_H_
