#ifndef SRC_CUTTING_INTERFACE_CUT_GENERATORS_AGGREGATING_CUT_GENERATORS_H_
#define SRC_CUTTING_INTERFACE_CUT_GENERATORS_AGGREGATING_CUT_GENERATORS_H_

#include <iostream>
#include <optional>

#include "ortools/base/status_macros.h"
#include "separator.h"
#include "src/cutting_interface/separator.h"
#include "src/solver.h"

// The cut generators in this file creates cutting planes by rounding row
// aggregations. Assume that row i in the LP is given by
//
// lhs[i] <= sum_j{ a[i, j]*x[j] } <= rhs[i]
//
// By introducing a slack variable s[i] >= 0, this inequality can be rewritten
// as the equality
//
// sum_j{ a[i, j]*x[j] } + s[i] = rhs[i]
// OR
// sum_j{ a[i, j]*x[j] } - s[i] = lhs[i]
//
// or, for unified notation
//
// sum_j{ a[i, j]*x[j] } + slack_sign[i]*s[i] = side_value[i]             (*)
//
// These equalities are valid for the entire LP, and consequently any linear
// combination is also valid for the entire LP. Hence we can derive new valid
// equalities
//
// sum{ variable_coefficients[j]*x[j] } + sum{ slack_coefficients[i]*s[i] }
//                                                          = right_hand_side
//
// where
//
// variable_coefficients[j] = sum_i{ weight[i]*a[i, j] }
// slack_coefficients[i] = slack_sign[i] * weight[i]
// right_hand_side = sum_i{ weight[i]*side_value[i] }
//
// The aggregated row can now be rounded in various ways to produce a cutting
// plane that is valid in the MIP but (ideally) not valid in the LP. Equation
// (*) is then used to remove the slack variables before the resulting cut is
// returned.
//
// The parameters in these cut generators can be boiled down to
// 1. Which side to use for each row.
// 2. How to weight each row in the aggregation.
// 3. How to round the obtained aggregation.
// Steps 1-2 are handled by the ConstraintAggregator classes, while step 3 is
// implemented in RoundingCutGenerator and controlled via parameters.

namespace minimip {

// Represents an aggregated row. See the top-level comment for the semantics
// of each field.
struct AggregatedRow {
  SparseRow variable_coefficients;
  SparseCol slack_coefficients;
  SparseCol slack_signs;
  SparseCol side_values;
  double right_hand_side = 0.0;

  std::string DebugString() const {
    std::stringstream ss;
    ss << "variables: " << variable_coefficients << "\n";
    ss << "slacks: " << slack_coefficients << "\n";
    ss << "signs: " << slack_signs << "\n";
    ss << "side values: " << side_values << "\n";
    ss << "rhs: " << right_hand_side;
    return ss.str();
  }
};

class Rounder {
 public:
  virtual ~Rounder() = default;

  virtual std::optional<AggregatedRow> RoundAggregatedRow(
      const MiniMipSolver& solver, AggregatedRow aggregated_row) const = 0;
};

// Create a cut from the aggregated row using Mixed Integer Rounding.
class MIRRounder : public Rounder {
 public:
  std::optional<AggregatedRow> RoundAggregatedRow(
      const MiniMipSolver& solver, AggregatedRow aggregated_row) const final;
};

// Create a cut from the aggregated row using strong Chvatal-Gomory rounding.
// Fails if the row contains any continuous variables with negative weight, in
// which case std::nullopt is returned.
class StrongCGRounder : public Rounder {
 public:
  std::optional<AggregatedRow> RoundAggregatedRow(
      const MiniMipSolver& solver, AggregatedRow aggregated_row) const final;
};

class TableauRoundingSeparator : public Separator {
 public:
  explicit TableauRoundingSeparator(SeparatorParameters params)
      : params_(std::move(params)) {
    DCHECK(params_.has_tableau_rounding_separator_parameters());
    if (params_.tableau_rounding_separator_parameters()
            .use_mixed_integer_rounding()) {
      rounders_.push_back(std::make_unique<MIRRounder>());
    }
    if (params_.tableau_rounding_separator_parameters()
            .use_strong_cg_rounding()) {
      rounders_.push_back(std::make_unique<StrongCGRounder>());
    }
  }

  absl::StatusOr<std::vector<CutData>> GenerateCuttingPlanes(
      const MiniMipSolver& solver) final;

 private:
  std::vector<std::unique_ptr<Rounder>> rounders_;
  const SeparatorParameters params_;
};

}  // namespace minimip

#endif  // SRC_CUTTING_INTERFACE_CUT_GENERATORS_AGGREGATING_CUT_GENERATORS_H_