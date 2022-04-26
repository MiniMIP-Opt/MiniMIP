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
// @TODO (cgraczyk): Remove unnecessary files like this / messagehandler before
// push
#ifndef SRC_MINIMIP_MINIMIP_DEF_H_
#define SRC_MINIMIP_MINIMIP_DEF_H_

#include "absl/status/status.h"
#include "src/lp_interface/lp_types.h"
#include "src/messagehandler/message_handler.h"
#include "src/messagehandler/message_macros.h"
#include "src/minimip/sparse_types.h"

// #define MINIMIP_VERSION 000   // MiniMIP version number (multiplied by 100 to
// get integer number)
//  #define MINIMIP_SUBVERSION               0 // MiniMIP sub version number
//  #define MINIMIP_APIVERSION              00 // MiniMIP API version number
//  #define MINIMIP_COPYRIGHT   // PLACEHOLDER:"~~Copyright (C) 2021
//  Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)~~"

#define EPSCEIL(x, eps) (ceil((x) - (eps)))
#define EPSFLOOR(x, eps) (floor((x) + (eps)))
#define REALABS(x) (fabs(x))

#define MINIMIP_CALL(x)                                       \
  {                                                           \
    absl::Status _restat_ = (x);                              \
    if (!_restat_.ok()) {                                     \
      MiniMIPerrorMessage("Error <%d> in function call\n",    \
                          static_cast<int>(_restat_.code())); \
      return _restat_;                                        \
    }                                                         \
  }

#endif  // SRC_MINIMIP_MINIMIP_DEF_H_
