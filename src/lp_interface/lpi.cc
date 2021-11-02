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
    absl::Status status = absl::OkStatus();
    for (size_t idx = 0; idx < indices.size(); ++idx) {
        status = this->SetObjectiveCoefficient(indices[idx], objective_coefficients[idx]);
        if (status != absl::OkStatus()) {
            return status;
        }
    }
    return status;
}

} // namespace minimip
#endif  // SRC_LP_INTERFACE_LPI_H_
