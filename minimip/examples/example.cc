#include "minimip/examples/example.h"

#include <utility>

#include "ortools/base/logging.h"
#include "ortools/linear_solver/linear_solver.h"
#include "ortools/linear_solver/linear_solver.pb.h"

namespace minimip {

using ::operations_research::MPModelProto;
using ::operations_research::MPModelRequest;
using ::operations_research::MPSolutionResponse;
using ::operations_research::MPSolver;

MPSolutionResponse SolveMPModelProtoWithGlop(MPModelProto model) {
  LOG(INFO) << "Will solve a model with Glop! Yay.";
  MPModelRequest request;
  *request.mutable_model() = std::move(model);
  request.set_solver_type(MPModelRequest::GLOP_LINEAR_PROGRAMMING);
  MPSolutionResponse response;
  MPSolver::SolveWithProto(request, &response);
  return response;
}

}  // namespace minimip
