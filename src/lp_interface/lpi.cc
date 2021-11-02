#ifndef SRC_LP_INTERFACE_LPI_C_
#define SRC_LP_INTERFACE_LPI_C_

#include "src/lp_interface/lpi.h"

namespace minimip {

absl::Status LPInterface::SetObjectiveCoefficients(
    const std::vector<int>&
        indices,  // column indices to change objective value for
    const std::vector<double>&
        objective_coefficients  // new objective values for columns
) {
  for (size_t idx = 0; idx < indices.size(); ++idx) {
    MINIMIP_CALL(this->SetObjectiveCoefficient(indices[idx],
                                           objective_coefficients[idx]));
  }
  return absl::OkStatus();
}

}  // namespace minimip
#endif  // SRC_LP_INTERFACE_LPI_H_
