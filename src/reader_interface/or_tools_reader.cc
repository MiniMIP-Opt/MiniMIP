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

#include "src/reader_interface/or_tools_reader.h"

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
  if (model_proto.annotation_size() != 0) {
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

absl::StatusOr<MiniMipProblem> OrToolsReader::ReadProblemDataFromFile(
    const std::string& file_path) const {
  // If the file isn't found, ortools only logs a warning and returns an empty
  // problem. This may lead to silent failures, so we issue a proper error.
  std::fstream fs(file_path, std::ios_base::in);
  if (!fs.good()) {
    return util::InvalidArgumentErrorBuilder()
           << "File not found: " << file_path;
  } else if (file_path.rfind(".mps") != file_path.length() - 4) {
    return absl::InvalidArgumentError("File format not supported");
  }

  ASSIGN_OR_RETURN(const operations_research::MPModelProto model_proto,
                   operations_research::glop::MpsFileToMPModelProto(file_path));

  return MPModelProtoToMiniMipProblem(model_proto);
}

absl::StatusOr<MiniMipProblem> OrToolsReader::ReadProblemDataFromString(
    const std::string& file_data) const {
  // Make sure that the data is in MPS format. Otherwise, undefined behaviour is
  // expected.
  ASSIGN_OR_RETURN(const operations_research::MPModelProto model_proto,
                   operations_research::glop::MpsDataToMPModelProto(file_data));

  return MPModelProtoToMiniMipProblem(model_proto);
}

// absl::StatusOr<operations_research::MPModelProto> ReadLpFile(
//     const std::string& file_path) const {
//   absl::StatusOr<operations_research::MPModelProto> problem_proto =
//       operations_research::glop::ParseLP(file_path);
//
//   if (!problem_proto.ok()) {
//     return problem_proto.status();
//   }
//
//   return problem_proto.value();
// }

}  // namespace minimip
