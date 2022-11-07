// Copyright 2022 the MiniMIP Project
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

#include "problem.h"

#include <fstream>

#include "ortools/linear_solver/linear_solver.pb.h"
#include "ortools/lp_data/mps_reader.h"

namespace minimip {
namespace {

absl::StatusOr<MiniMipProblem> MPModelProtoToMiniMipProblem(
    const operations_research::MPModelProto& model_proto) {
  using operations_research::MPConstraintProto;
  using operations_research::MPVariableProto;

  if (model_proto.general_constraint_size() != 0) {
    return absl::InvalidArgumentError(
        "MiniMIP only supports linear constraints.");
  }
  if (model_proto.has_quadratic_objective()) {
    return absl::InvalidArgumentError(
        "MiniMIP only supports linear objectives.");
  }
  if (model_proto.annotation_size()) {
    return absl::InvalidArgumentError("MiniMIP doesn't support annotations.");
  }

  MiniMipProblem problem;
  problem.name = model_proto.name();
  problem.is_maximization = model_proto.maximize();
  problem.objective_offset = model_proto.objective_offset();

  problem.variables.reserve(model_proto.variable_size());
  for (const MPVariableProto& variable_proto : model_proto.variable()) {
    problem.variables.push_back(
        {.name = variable_proto.name(),
         .objective_coefficient = variable_proto.objective_coefficient(),
         .lower_bound = variable_proto.lower_bound(),
         .upper_bound = variable_proto.upper_bound(),
         .is_integer = variable_proto.is_integer()});
  }
  problem.constraints.reserve(model_proto.constraint_size());
  for (const MPConstraintProto& constraint_proto : model_proto.constraint()) {
    problem.constraints.push_back(
        {.name = constraint_proto.name(),
         .var_indices = {constraint_proto.var_index().begin(),
                         constraint_proto.var_index().end()},
         .coefficients = {constraint_proto.coefficient().begin(),
                          constraint_proto.coefficient().end()},
         .left_hand_side = constraint_proto.lower_bound(),
         .right_hand_side = constraint_proto.upper_bound()});
  }
  if (model_proto.has_solution_hint()) {
    problem.hints.push_back(
        {.var_indices = {model_proto.solution_hint().var_index().begin(),
                         model_proto.solution_hint().var_index().end()},
         .values = {model_proto.solution_hint().var_value().begin(),
                    model_proto.solution_hint().var_value().end()}});
  }

  const std::string error_string = FindErrorInMiniMipProblem(problem);
  if (!error_string.empty()) {
    return util::InvalidArgumentErrorBuilder()
           << "MiniMipProblem converted from MPModelProto was invalid: "
           << error_string;
  }
  return problem;
}

}  // namespace

// Returns the first error encountered while validating MiniMipProblem.
// Returns empty string if `problem` is valid.
std::string FindErrorInMiniMipProblem(const MiniMipProblem& problem) {
  for (const MiniMipVariable& variable : problem.variables) {
    if (variable.name.empty()) {
      return "No variable name given.";
    }
    if (variable.upper_bound < variable.lower_bound) {
      return absl::StrFormat(
          "Variable %s's lower variable bound is less than its upper variable "
          "bound.",
          variable.name);
    }
  }

  for (const MiniMipConstraint& constraint : problem.constraints) {
    if (constraint.name.empty()) {
      return "No constraint name given.";
    }
    absl::flat_hash_set<int> var_indices_set;
    for (int var_index : constraint.var_indices) {
      if (var_index < 0) {
        return absl::StrFormat("Negative variable index in constraint %s.",
                               constraint.name);
      }
      if (var_index > problem.variables.size()) {
        return absl::StrFormat("Variable index %d does not exist.", var_index);
      }
      if (!var_indices_set.insert(var_index).second) {
        return absl::StrFormat(
            "Variable indices in constraint %s are not unique.",
            constraint.name);
      }
    }
    if (constraint.var_indices.size() != constraint.coefficients.size()) {
      return absl::StrFormat(
          "Constraint %s has unequal number of indices and coefficients.",
          constraint.name);
    }
    if (constraint.right_hand_side < constraint.left_hand_side) {
      return absl::StrFormat(
          "Constraint %s has lesser left- than right-hand side.",
          constraint.name);
    }
  }
  return "";
}

absl::StatusOr<MiniMipProblem> ReadProblemFromMPSFile(
    const std::string& file_name) {
  // If the file isn't found, ortools only logs a warning and returns an empty
  // problem. This may lead to silent failures, so we issue a proper error.
  std::fstream fs(file_name, std::ios_base::in);
  if (!fs.good()) {
    return util::InvalidArgumentErrorBuilder()
           << "File not found: " << file_name;
  }

  ASSIGN_OR_RETURN(const operations_research::MPModelProto model_proto,
                   operations_research::glop::MpsFileToMPModelProto(file_name));
  return MPModelProtoToMiniMipProblem(model_proto);
}

absl::StatusOr<MiniMipProblem> ReadProblemFromMPSData(
    const std::string& mps_data) {
  ASSIGN_OR_RETURN(const operations_research::MPModelProto model_proto,
                   operations_research::glop::MpsDataToMPModelProto(mps_data));
  return MPModelProtoToMiniMipProblem(model_proto);
}

}  // namespace minimip
