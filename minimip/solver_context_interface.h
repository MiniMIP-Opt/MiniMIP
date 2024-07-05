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

#ifndef MINIMIP_SOLVER_CONTEXT_INTERFACE_H_
#define MINIMIP_SOLVER_CONTEXT_INTERFACE_H_

#include "minimip/cutting_interface/cuts_runner.h"
#include "minimip/branching_interface/branching_interface.h"
#include "minimip/data_structures/cuts_data.h"
#include "minimip/data_structures/mip_data.h"
#include "minimip/data_structures/mip_tree.h"
#include "minimip/data_structures/problem.h"
#include "minimip/lp_interface/lpi.h"

namespace minimip {

class SolverContextInterface {
 public:
  virtual ~SolverContextInterface() = default;

  virtual const MipData& mip_data() const = 0;
  virtual MipData& mutable_mip_data() = 0;

  virtual const MiniMipResult& result() const = 0;
  virtual MiniMipResult& mutable_result() = 0;

  virtual const MipTree& mip_tree() const = 0;
  virtual MipTree& mutable_mip_tree() = 0;

  virtual const CutRegistry& cut_registry() const = 0;
  virtual CutRegistry& mutable_cut_registry() = 0;

  virtual CutRunnerInterface* mutable_cut_runner() const = 0;

  virtual const BranchingRuleInterface* branching_interface() const = 0;

  virtual const LpInterface* lpi() const = 0;
  virtual LpInterface* mutable_lpi() = 0;

  virtual bool IsIntegerWithinTolerance(double d) const = 0;

  virtual bool IsEqualToWithinTolerance(double d, double b) const = 0;

  virtual double FloorWithTolerance(double d) const = 0;
};

}  // namespace minimip
#endif  // MINIMIP_SOLVER_CONTEXT_INTERFACE_H_
