#ifndef SRC_SOLVER_H
#define SRC_SOLVER_H

#include <memory>
#include <vector>

#include "absl/status/status.h"
#include "ortools/base/status_builder.h"
#include "ortools/base/status_macros.h"
#include "src/data_structures/mip_data.h"
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
    ASSIGN_OR_RETURN(std::unique_ptr<LPInterface> lpi,
                     ConfigureLPSolverFromProto(params.lp_parameters()));
    auto solver = std::unique_ptr<Solver>(
        new Solver(params, std::move(mip_data), std::move(lpi)));
    return solver;
  }

  absl::StatusOr<MiniMipResult> Solve() {
    return absl::UnimplementedError("Solve isn't implemented yet.");
  }

  const MipData& mip_data() const { return mip_data_; }
  MipData& mutable_mip_data() { return mip_data_; }

  const LPInterface* lpi() const { return lpi_.get(); }
  LPInterface* mutable_lpi() { return lpi_.get(); }

  bool IsIntegerWithinTolerance(double d) {
    return std::abs(d - std::round(d)) <= params_.integrality_tolerance();
  }

 private:
  const MiniMipParameters params_;
  MipData mip_data_;
  std::unique_ptr<LPInterface> lpi_;

  // Protected constructor, use Create() instead.
  Solver(MiniMipParameters params, MipData mip_data,
         std::unique_ptr<LPInterface> lpi)
      : params_{std::move(params)},
        mip_data_{std::move(mip_data)},
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