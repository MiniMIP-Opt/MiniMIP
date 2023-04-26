#include "src/cutting_interface/aggregating_separator.h"
#include "src/cutting_interface/separator.h"

namespace minimip {

inline absl::StatusOr<std::unique_ptr<Separator>> ConfigureSeparatorFromProto(
    const SeparatorParameters& separator_parameters) {
  if (separator_parameters.has_tableau_rounding_separator_parameters()) {
    return std::make_unique<TableauRoundingSeparator>(separator_parameters);
  }
  return absl::InvalidArgumentError("No separator-specific parameters set.");
}

}  // namespace minimip