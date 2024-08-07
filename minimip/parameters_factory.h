// Copyright 2024 the MiniMIP Project
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

#ifndef MINIMIP_PARAMETERS_FACTORY_H_
#define MINIMIP_PARAMETERS_FACTORY_H_

#include "minimip/parameters.pb.h"

namespace minimip {

// Create the default CutGenerator parameters.
// As of 26.10.23: The TableauRoundingGenerator as the default.
inline CutGeneratorParameters DefaultCutGeneratorParameters() {
  VLOG(10) << "calling DefaultCutGeneratorParameters().";
  CutGeneratorParameters params;
  params.mutable_tableau_rounding_generator_parameters();
  return params;
}

// Create the default CutSelector parameters.
// As of 26.10.23: The HybridSelector as the default.
inline CutSelectorParameters DefaultCutSelectorParameters() {
  VLOG(10) << "calling DefaultCutSelectorParameters().";
  CutSelectorParameters params;
  params.mutable_hybrid_selector_parameters();
  return params;
}

// Create the default CutRunnerParameters
// As of 26.10.23: The DefaultRunner is the default.
inline CutRunnerParameters DefaultCutRunnerParameters() {
  VLOG(10) << "calling DefaultCutRunnerParameters().";
  CutRunnerParameters runner_params;
  // Set the default runner parameters.
  runner_params.mutable_default_runner_parameters();

  // Set the default generator and selector parameters.
  *runner_params.add_generator_parameters() = DefaultCutGeneratorParameters();
  *runner_params.add_selector_parameters() = DefaultCutSelectorParameters();

  return runner_params;
}

// Create the default BranchingParameters
// As of 30.04.24: The MaxFractionalBranching is the default.
inline BranchingParameters DefaultBranchingParameters() {
  BranchingParameters params = BranchingParameters();
  params.mutable_max_fractional_branching_parameters();
  return params;
}

// Unified factory method to create default parameters
inline MiniMipParameters DefaultMiniMipParameters() {
  VLOG(10) << "calling DefaultMiniMipParameters().";
  MiniMipParameters default_params;

  // initialize LP parameters.
  default_params.mutable_lp_parameters();

  // initialize branching parameters.
  *default_params.mutable_branching_parameters() = DefaultBranchingParameters();

  // Initialize the default CutRunner parameters.
  *default_params.mutable_cut_runner() = DefaultCutRunnerParameters();

  return default_params;
}

// TODO(cgraczy): Add intelligent merging of user parameters and default
// parameters. Factory method to create fully initialized user modified
// parameters.
inline MiniMipParameters UserCustomizedParameters(
    const MiniMipParameters& user_params) {
  VLOG(10) << "calling UserCustomizedParameters().";
  MiniMipParameters params = DefaultMiniMipParameters();

  // The user's settings will overwrite the defaults where they're provided.
  params.MergeFrom(user_params);

  return params;
}

}  // namespace minimip

#endif  // MINIMIP_PARAMETERS_FACTORY_H_
