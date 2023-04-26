#include "src/cutting_interface/cut_runner.h"
#include "src/parameters.pb.h"

namespace minimip {

inline absl::StatusOr<std::unique_ptr<CuttingInterface>>
ConfigureCutInterfaceFromProto(const CutRunnerParameters& parameters) {
  if (parameters.has_some_cut_runner_parameters()) {
    // Here we would create and return the cut runner.
  }
  return nullptr;
  /* TODO(cgraczy): Return an error instead once we have at least one runner
     implementation. We currently return nullptr to allow unit tests to create a
     solver object.
  return absl::InvalidArgumentError("No cut runner implementation was chosen.");
  */
}

}  // namespace minimip