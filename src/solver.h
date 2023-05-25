// Copyright 2023 the MiniMIP Project
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

#ifndef SRC_SOLVER_H
#define SRC_SOLVER_H

#include <memory>
#include <vector>

#include "absl/status/status.h"
#include "ortools/base/status_builder.h"
#include "ortools/base/status_macros.h"
#include "src/cutting_interface/cuts_runner.h"
#include "src/cutting_interface/cuts_separator.h"
#include "src/data_structures/cuts_data.h"
#include "src/data_structures/mip_data.h"
#include "src/data_structures/mip_tree.h"
#include "src/data_structures/problem.h"
#include "src/lp_interface/lpi.h"
#include "src/lp_interface/lpi_factory.h"
#include "src/parameters.pb.h"

namespace minimip {

// This is the main solver class. It serves as both the main point of contact
// between MiniMIP and client code. It also owns all global data structures and
// modules.
class Solver {
 public:
  // Factory method to create a Solver object from a set of parameters.
  static absl::StatusOr<std::unique_ptr<Solver>> Create(
      const MiniMipParameters& params, const MiniMipProblem& problem) {
    const std::string problem_error = FindErrorInMiniMipProblem(problem);
    if (!problem_error.empty()) {
      return util::InvalidArgumentErrorBuilder()
             << "Error found in problem: " << problem_error;
    }
    MipData mip_data(problem);
    MipTree mip_tree;
    CutStorage cut_storage;
    ASSIGN_OR_RETURN(std::unique_ptr<LPInterface> lpi,
                     ConfigureLPSolverFromProto(params.lp_parameters()));
    auto solver = std::unique_ptr<Solver>(
        new Solver(params, std::move(mip_data), std::move(mip_tree),
                   std::move(cut_storage), std::move(lpi)));
    return solver;
  }

  absl::StatusOr<MiniMipResult> Solve() {
    return absl::UnimplementedError("Solve isn't implemented yet.");
  }

  const MipData& mip_data() const { return mip_data_; }
  MipData& mutable_mip_data() { return mip_data_; }

  const MipTree& mip_tree() const { return mip_tree_; }
  MipTree& mutable_mip_tree() { return mip_tree_; }

  const CutStorage& cut_storage() const { return cut_storage_; }
  CutStorage& mutable_cut_storage() { return cut_storage_; }

  const LPInterface* lpi() const { return lpi_.get(); }
  LPInterface* mutable_lpi() { return lpi_.get(); }

  bool IsIntegerWithinTolerance(double d) const {
    return std::abs(d - std::round(d)) <= params_.integrality_tolerance();
  }

  double FloorWithTolerance(double d) const {
    return IsIntegerWithinTolerance(d) ? std::round(d) : std::floor(d);
  }

 private:
  // Parameters specifying the MiniMip configuration.
  const MiniMipParameters params_;

  // The canonical copy of the MIP data. This contains all information from the
  // original problem, but not cuts etc. For the currently used LP information
  // (including cuts), see the LPI.
  MipData mip_data_;

  // Contains all open nodes and bound changes.
  MipTree mip_tree_;

  // Contains all currently active and stored cuts.
  CutStorage cut_storage_;
  // Entry-point for the cut API, used to create and activate cuts in the main
  // cut-and-price loop (not yet implemented).
  std::unique_ptr<CuttingInterface> cut_runner_;
  // Handle to an LP solver.
  std::unique_ptr<LPInterface> lpi_;

  // Protected constructor, use Create() instead.
  Solver(const MiniMipParameters& params, MipData mip_data, MipTree mip_tree,
         CutStorage cut_storage,
         std::unique_ptr<CuttingInterface> cut_runner, std::unique_ptr<LPInterface> lpi)
      : params_{std::move(params)},
        mip_data_{std::move(mip_data)},
        mip_tree_{std::move(mip_tree)},
        cut_storage_(std::move(cut_storage)),
        cut_runner_{std::move(cut_runner)},
        lpi_{std::move(lpi)} {}
};

// Convenience function to create a solver and solve the given problem.
inline absl::StatusOr<MiniMipResult> Solve(const MiniMipParameters& parameters,
                                           const MiniMipProblem& problem) {
  ASSIGN_OR_RETURN(std::unique_ptr<Solver> solver,
                   Solver::Create(parameters, problem));
  return solver->Solve();
}

}  // namespace minimip

#endif  // SRC_SOLVER_H
