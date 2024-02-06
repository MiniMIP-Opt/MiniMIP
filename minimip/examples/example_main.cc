#include <google/protobuf/text_format.h>

#include "minimip/examples/example.h"
#include "ortools/base/logging.h"
#include "ortools/linear_solver/linear_solver.pb.h"

namespace minimip {
namespace {

void RealMain() {
  operations_research::MPModelProto model;
  CHECK(google::protobuf::TextFormat::ParseFromString(
      R"pb(
        maximize: true
        variable { upper_bound: 13.0 objective_coefficient: 2.0 }
        variable { lower_bound: 42.0 objective_coefficient: -0.5 }
      )pb",
      &model));
  LOG(INFO) << "Solving:\n" << model.DebugString();
  const operations_research::MPSolutionResponse response =
      SolveMPModelProtoWithGlop(std::move(model));
  LOG(INFO) << "Response:\n" << response.DebugString();
}

}  // namespace
}  // namespace minimip

int main(int argc, char* argv[]) {
  absl::ParseCommandLine(argc, argv);
  google::InitGoogleLogging(argv[0]);
  minimip::RealMain();
  return EXIT_SUCCESS;
}
