#include "absl/status/statusor.h"
#include "src/data_structures/problem.h"
#include "src/solver.h"

namespace minimip {

// Factory function that configures
absl::StatusOr<std::unique_ptr<MiniMipSolver>> ConfigureMiniMipSolverFromProto(
    const MiniMipParameters& params, const MiniMipProblem& problem);

}  // namespace minimip