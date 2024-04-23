// Copyright 2024 the MiniMIP Project
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef SRC_DATA_STRUCTURES_MIP_DATA_H_
#define SRC_DATA_STRUCTURES_MIP_DATA_H_

#include <limits>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "src/data_structures/problem.h"
#include "src/data_structures/strong_sparse_matrix.h"

namespace minimip {

enum class VariableType {
  kFractional = 0,
  kBinary = 1,
  kInteger = 2,
  kImpliedInteger = 3,
};

// ============================================================================
// Class storing the internal problem representation
// ============================================================================

class MipData {
 public:
  // ==========================================================================
  // Constructors
  // ==========================================================================

  MipData() : MipData(MiniMipProblem{}) {};

  explicit MipData(const MiniMipProblem& problem);

  bool SolutionIsIntegral(const absl::StrongVector<ColIndex, double>& primal_values, double tolerance) const;
  // ==========================================================================
  // Problem getters.
  // ==========================================================================

  std::string problem_name() const { return problem_name_; }

  bool is_maximization() const { return is_maximization_; }

  double objective_offset() const { return objective_offset_; }

  const SparseRow& objective() const { return objective_; }

  const std::vector<MiniMipSolutionHint>& hints() const {
    return solution_hints_;
  }

  const absl::StrongVector<ColIndex, double>& lower_bounds() const {
    return lower_bounds_;
  }

  const absl::StrongVector<ColIndex, double>& upper_bounds() const {
    return upper_bounds_;
  }

  const absl::StrongVector<RowIndex, double>& left_hand_sides() const {
    return left_hand_sides_;
  }

  const absl::StrongVector<RowIndex, double>& right_hand_sides() const {
    return right_hand_sides_;
  }

  const absl::StrongVector<RowIndex, bool>& is_integral_constraint() const {
    return is_integral_constraint_;
  }

  const StrongSparseMatrix& matrix() const { return constraint_matrix_; }

  const absl::StrongVector<ColIndex, VariableType>& variable_types() const {
    return variable_types_;
  }

  const absl::StrongVector<ColIndex, std::string>& variable_names() const {
    return variable_names_;
  }

  const absl::StrongVector<RowIndex, std::string>& constraint_names() const {
    return constraint_names_;
  }

  const absl::flat_hash_set<ColIndex>& integer_variables() const {
    return integer_variables_;
  }

 private:
  // TODO(lpawel): Clean up and comment all fields.
  std::string problem_name_;
  bool is_maximization_ = false;
  double objective_offset_ = 0.0;
  SparseRow objective_;
  std::vector<MiniMipSolutionHint> solution_hints_;

  absl::StrongVector<ColIndex, std::string> variable_names_;
  absl::StrongVector<ColIndex, double> lower_bounds_;
  absl::StrongVector<ColIndex, double> upper_bounds_;
  absl::flat_hash_set<ColIndex> integer_variables_;
  absl::StrongVector<ColIndex, VariableType> variable_types_;

  absl::StrongVector<RowIndex, std::string> constraint_names_;
  absl::StrongVector<RowIndex, double> left_hand_sides_;
  absl::StrongVector<RowIndex, double> right_hand_sides_;
  absl::StrongVector<RowIndex, bool> is_integral_constraint_;
  StrongSparseMatrix constraint_matrix_;
};

}  // namespace minimip
#endif  // SRC_DATA_STRUCTURES_MIP_DATA_H_
