#ifndef MINIMIP_EXAMPLES_EXAMPLE_H_
#define MINIMIP_EXAMPLES_EXAMPLE_H_

#include "ortools/linear_solver/linear_solver.pb.h"

namespace minimip {

operations_research::MPSolutionResponse SolveMPModelProtoWithGlop(
    operations_research::MPModelProto model);

}  // namespace minimip

#endif  // MINIMIP_EXAMPLES_EXAMPLE_H_
