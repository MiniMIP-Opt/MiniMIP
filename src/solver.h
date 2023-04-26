#ifndef SRC_SOLVER_H
#define SRC_SOLVER_H

#include <memory>
#include <vector>

#include "absl/status/status.h"
#include "cutting_interface/cut_runner.h"
#include "lp_interface/lpi.h"
#include "ortools/base/status_builder.h"
#include "ortools/base/status_macros.h"
#include "src/cutting_interface/cuts_runner.h"
#include "src/cutting_interface/cuts_separator.h"
#include "src/data_structures/cuts_data.h"
#include "src/data_structures/mip_data.h"
#include "src/lp_interface/lpi.h"
#include "src/parameters.pb.h"

namespace minimip {

// This is the main solver class. It serves as both the main point of contact
// between MiniMIP and client code. It also owns all global data structures and
// modules.
class MiniMipSolver {
public:
  absl::StatusOr<MiniMipResult> Solve() {
    return absl::UnimplementedError("Solve isn't implemented yet.");
  }

  const MipData &mip_data() const { return *mip_data_; }
  MipData &mutable_mip_data() { return *mip_data_; }

  const CutStorage &cut_storage() const { return *cut_storage_; }
  CutStorage &mutable_cut_storage() { return *cut_storage_; }

  const LPInterface &lpi() const { return *lpi_; }
  LPInterface &mutable_lpi() { return *lpi_; }

  bool IsIntegerWithinTolerance(double d) const {
    return std::abs(d - std::round(d)) <= params_.integrality_tolerance();
  }

  double FloorWithTolerance(double d) const {
    return IsIntegerWithinTolerance(d) ? std::round(d) : std::floor(d);
  }

  // Directly configure the solver. For most use cases, the factory function in
  // `solver_factory.h` is preferred.
  MiniMipSolver(MiniMipParameters params, std::unique_ptr<MipData> mip_data,
                std::unique_ptr<CutStorage> cut_storage,
                std::unique_ptr<LPInterface> lpi,
                std::unique_ptr<CuttingInterface> cut_runner)
      : params_{std::move(params)}, mip_data_{std::move(mip_data)},
        cut_storage_(std::move(cut_storage)), lpi_{std::move(lpi)},
        cut_runner_{std::move(cut_runner)} {}

private:
  // Parameters specifying the MiniMip configuration.
  const MiniMipParameters params_;

  // The canonical copy of the MIP data. This contains all information from the
  // original problem, but not cuts etc. For the currently used LP information
  // (including cuts), see the LPI.
  std::unique_ptr<MipData> mip_data_;

  // Contains all currently active and stored cuts.
  std::unique_ptr<CutStorage> cut_storage_;

  // Handle to an LP solver.
  std::unique_ptr<LPInterface> lpi_;

  // Entry-point for the cut API, used to create and activate cuts in the main
  // cut-and-price loop (not yet implemented).
  std::unique_ptr<CuttingInterface> cut_runner_;
};

} // namespace minimip

#endif // SRC_SOLVER_H
