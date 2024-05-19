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

#ifndef SRC_SOLVER_H
#define SRC_SOLVER_H

#include <memory>
#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "solver_context_interface.h"
#include "src/cutting_interface/runner_factory.h"
#include "src/branching_interface/branching_factory.h"
#include "src/lp_interface/lpi_factory.h"
#include "src/parameters.pb.h"
#include "src/parameters_factory.h"

namespace minimip {

// This is the main solver class. It serves as both the main point of contact
// between MiniMIP and client code. It also owns all global data structures and
// modules.
class Solver : public SolverContextInterface {
 public:
  // Factory method to create a Solver object from a set of parameters.
  static absl::StatusOr<std::unique_ptr<Solver>> Create(
      const MiniMipProblem& problem = MiniMipProblem(),
      const MiniMipParameters& user_params = MiniMipParameters()) {
    VLOG(10) << "calling Create().";
    const std::string problem_error = FindErrorInMiniMipProblem(problem);
    if (!problem_error.empty()) {
      return util::InvalidArgumentErrorBuilder()
             << "Error found in problem: " << problem_error;
    }
    // The user's settings will overwrite the defaults where they're provided.
    MiniMipParameters params = UserCustomizedParameters(user_params);

    MipData mip_data(problem);
    MipTree mip_tree;
    CutRegistry cut_registry;
    MiniMipResult result;

    ASSIGN_OR_RETURN(std::unique_ptr<LpInterface> lpi,
                     CreateLpSolver(params.lp_parameters()));

    ASSIGN_OR_RETURN(std::unique_ptr<CutRunnerInterface> cut_runner,
                     CreateCutRunner(params.cut_runner()));

    ASSIGN_OR_RETURN(std::unique_ptr<BranchingInterface> branching_interface,
                     CreateBranching(params.branching_parameters()));

    auto solver = std::unique_ptr<Solver>(new Solver(
        params, std::move(mip_data), std::move(result), std::move(mip_tree),
        std::move(cut_registry), std::move(cut_runner), std::move(branching_interface),
        std::move(lpi)));

    // Populate LP with initial problem data
    CHECK_OK(solver->mutable_lpi()->PopulateFromMipData(solver->mip_data_));

    return solver;
  }

  absl::Status Solve();

  const MipData& mip_data() const override { return mip_data_; }
  MipData& mutable_mip_data() override { return mip_data_; }

  const MiniMipResult& result() const override { return result_; }
  MiniMipResult& mutable_result() override { return result_; }

  const MipTree& mip_tree() const override { return mip_tree_; }
  MipTree& mutable_mip_tree() override { return mip_tree_; }

  const CutRegistry& cut_registry() const override { return cut_registry_; }
  CutRegistry& mutable_cut_registry() override { return cut_registry_; }

  CutRunnerInterface* mutable_cut_runner() const override {
    return cut_runner_.get();
  }

  const BranchingInterface* branching_interface() const override {
    return branching_interface_.get();
  }

  const LpInterface* lpi() const override { return lpi_.get(); }
  LpInterface* mutable_lpi() override { return lpi_.get(); }

  bool IsIntegerWithinTolerance(double d) const override {
    return std::abs(d - std::round(d)) <= params_.integrality_tolerance();
  }

  bool IsEqualToWithinTolerance(double d, double b) const override {
    return std::abs(d - b) <= params_.numerical_tolerance();
  }

  double FloorWithTolerance(double d) const override {
    return IsIntegerWithinTolerance(d) ? std::round(d) : std::floor(d);
  }

 private:
  // Parameters specifying the MiniMip configuration.
  const MiniMipParameters params_;

  // The canonical copy of the MIP data. This contains all information from the
  // original problem, but not cuts etc. For the currently used LP information
  // (including cuts), see the LPI.
  MipData mip_data_;

  // The result struct of the optimization process. This is updated as the
  // solver progresses and stores incumbents, bounds, etc.
  MiniMipResult result_;

  // Contains all open nodes and bound changes.
  MipTree mip_tree_;

  // Contains all currently active and stored cuts.
  CutRegistry cut_registry_;

  // Entry-point for the cut API, used to create and activate cuts in the main
  // cut-and-price loop (not yet implemented).
  std::unique_ptr<CutRunnerInterface> cut_runner_;

  // Handle Branching Candidate selection.
  std::unique_ptr<BranchingInterface> branching_interface_;

  // Handle to an LP solver.
  std::unique_ptr<LpInterface> lpi_;

  // Protected constructor, use Create() instead.
  Solver(MiniMipParameters params, MipData mip_data, MiniMipResult result,
         MipTree mip_tree, CutRegistry cut_registry,
         std::unique_ptr<CutRunnerInterface> cut_runner,
         std::unique_ptr<BranchingInterface> branching_interface,
         std::unique_ptr<LpInterface> lpi)
      : params_{std::move(params)},
        mip_data_{std::move(mip_data)},
        result_{std::move(result)},
        mip_tree_{std::move(mip_tree)},
        cut_registry_{std::move(cut_registry)},
        cut_runner_{std::move(cut_runner)},
        branching_interface_{std::move(branching_interface)},
        lpi_{std::move(lpi)} {}
};

// Convenience function to create a solver and solve the given problem.
inline absl::StatusOr<MiniMipResult> Solve(
    const MiniMipProblem& problem, const MiniMipParameters& user_params) {
  ASSIGN_OR_RETURN(std::unique_ptr<Solver> solver,
                   Solver::Create(problem, user_params));
  return solver->Solve();
}

}  // namespace minimip

#endif  // SRC_SOLVER_H
