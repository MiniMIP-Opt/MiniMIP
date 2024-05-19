#include "src/lp_interface/lpi.h"

namespace minimip {

absl::Status LpParametersAreValid(const LpParameters& params) {
  VLOG(10) << "calling LpParametersAreValid().";
  if (params.feasibility_tolerance() < 0.0 &&
      params.feasibility_tolerance() != -1.0) {
    return absl::InvalidArgumentError(
        "Feasibility tolerance must be non-negative or -1 for default.");
  }
  if (params.optimality_tolerance() < 0.0 &&
      params.optimality_tolerance() != -1.0) {
    return absl::InvalidArgumentError(
        "Optimality tolerance must be non-negative or -1 for default.");
  }
  if (params.min_markowitz_threshold() < 0.0 &&
      params.min_markowitz_threshold() != -1.0) {
    return absl::InvalidArgumentError(
        "Markowitz threshold must be non-negative or -1 for default.");
  }
  if (params.refactorization_interval() < 0) {
    return absl::InvalidArgumentError(
        "Refactorization interval must be positive or 0 for default.");
  }
  if (params.timing_mode() != LpParameters::TIMING_OFF &&
      params.time_limit() <= 0) {
    return absl::InvalidArgumentError("Time limit must be positive.");
  }
  if (params.iteration_limit() < 0) {
    return absl::InvalidArgumentError(
        "Iteration limit must be positive or 0 for default.");
  }
  if (params.num_threads() < 0) {
    return absl::InvalidArgumentError(
        "Number of threads must be positive or 0 for default.");
  }
  if (params.random_seed() < 0) {
    return absl::InvalidArgumentError(
        "Random seed must be positive or 0 for default.");
  }

  return absl::OkStatus();
}

}  // namespace minimip
