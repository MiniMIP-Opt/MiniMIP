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

#ifndef SRC_PARAMETERS_FACTORY_H_
#define SRC_PARAMETERS_FACTORY_H_

#include "src/parameters.pb.h"

namespace minimip {
// Unified factory method to create default parameters
inline MiniMipParameters DefaultMiniMipParameters() {
  MiniMipParameters default_params;

  // initialize LP parameters.
  default_params.mutable_lp_parameters();
  
  return default_params;
}

// TODO(cgraczy): Add intelligent merging of user parameters and default
// parameters. Factory method to create fully initialized user modified
// parameters.
inline MiniMipParameters UserCustomizedParameters(
    const MiniMipParameters& user_params) {
  MiniMipParameters params = DefaultMiniMipParameters();

  // The user's settings will overwrite the defaults where they're provided.
  params.MergeFrom(user_params);

  return params;
}

}  // namespace minimip
#endif  // SRC_PARAMETERS_FACTORY_H_
