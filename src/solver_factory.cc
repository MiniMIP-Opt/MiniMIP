#include "src/solver_factory.h"

#include "src/cutting_interface/cut_runner_factory.h"
#include "src/cutting_interface/separator_factory.h"
#include "src/lp_interface/lpi_factory.h"

namespace minimip {

absl::StatusOr<std::unique_ptr<MiniMipSolver>> ConfigureMiniMipSolverFromProto(
    const MiniMipParameters& params, const MiniMipProblem& problem) {
  const std::string problem_error = FindErrorInMiniMipProblem(problem);
  if (!problem_error.empty()) {
    return util::InvalidArgumentErrorBuilder()
           << "Error found in problem: " << problem_error;
  }
  auto mip_data = std::make_unique<MipData>(problem);
  auto cut_storage = std::make_unique<CutStorage>();
  ASSIGN_OR_RETURN(std::unique_ptr<LPInterface> lpi,
                   ConfigureLPSolverFromProto(params.lp_parameters()));
  ASSIGN_OR_RETURN(std::unique_ptr<CuttingInterface> cut_runner,
                   ConfigureCutInterfaceFromProto(params.cut_runner()));
  for (const SeparatorParameters& separator_params : params.separators()) {
    ASSIGN_OR_RETURN(std::unique_ptr<Separator> separator,
                     ConfigureSeparatorFromProto(separator_params));
    cut_runner->AddSeparator(std::move(separator));
  }
  // TODO(cgraczy): Configure the selector(s) in a similar way.
  auto solver = std::make_unique<MiniMipSolver>(
      params, std::move(mip_data), std::move(cut_storage), std::move(lpi),
      std::move(cut_runner));
  return solver;
}

}  // namespace minimip