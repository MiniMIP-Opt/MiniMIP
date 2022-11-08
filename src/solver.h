#ifndef SRC_SOLVER_H
#define SRC_SOLVER_H

#include <memory>
#include <vector>

#include "absl/status/status.h"
#include "ortools/base/status_builder.h"
#include "ortools/base/status_macros.h"
#include "src/cutting_interface/cut_runner.h"
#include "src/cutting_interface/separator.h"
#include "src/data_structures/cuts_data.h"
#include "src/data_structures/mip_data.h"
#include "src/data_structures/problem.h"
#include "src/lp_interface/lpi.h"
#include "src/lp_interface/lpi_factory.h"
#include "src/parameters.pb.h"

namespace minimip {

// This is the main solver class. It serves as both the main point of contact
// between MiniMIP and client code. It also owns all global data structures and
// modules.
class MiniMipSolver {
 public:
  // Factory method to create a Solver object from a set of parameters.
  static absl::StatusOr<std::unique_ptr<MiniMipSolver>> Create(
      const MiniMipParameters& params, const MiniMipProblem& problem);

  absl::StatusOr<MiniMipResult> Solve() {
    return absl::UnimplementedError("Solve isn't implemented yet.");
  }

  const MipData& mip_data() const { return mip_data_; }
  MipData& mutable_mip_data() { return mip_data_; }

  const CutStorage& cut_storage() const { return cut_storage_; }
  CutStorage& mutable_cut_storage() { return cut_storage_; }

  const LPInterface* lpi() const { return lpi_.get(); }
  LPInterface* mutable_lpi() { return lpi_.get(); }

  bool IsIntegerWithinTolerance(double d) {
    return std::abs(d - std::round(d)) <= params_.integrality_tolerance();
  }

 private:
  const MiniMipParameters params_;
  MipData mip_data_;
  CutStorage cut_storage_;
  std::unique_ptr<LPInterface> lpi_;
  std::unique_ptr<CuttingInterface> cut_runner_;

  // Protected constructor, use Create() instead.
  MiniMipSolver(const MiniMipParameters& params, MipData mip_data,
                CutStorage cut_storage, std::unique_ptr<LPInterface> lpi,
                std::unique_ptr<CuttingInterface> cut_runner)
      : params_{std::move(params)},
        mip_data_{std::move(mip_data)},
        cut_storage_(std::move(cut_storage)),
        lpi_{std::move(lpi)},
        cut_runner_{std::move(cut_runner)} {}
};

// Convenience function to create a solver and solve the given problem.
inline absl::StatusOr<MiniMipResult> Solve(const MiniMipParameters& parameters,
                                           const MiniMipProblem& problem) {
  ASSIGN_OR_RETURN(std::unique_ptr<MiniMipSolver> solver,
                   MiniMipSolver::Create(parameters, problem));
  return solver->Solve();
}

}  // namespace minimip

#endif  // SRC_SOLVER_H
