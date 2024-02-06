#include "minimip/examples/example.h"

#include "google/protobuf/text_format.h"
#include "google/protobuf/util/message_differencer.h"
#include "gtest/gtest.h"
#include "ortools/base/logging.h"
#include "ortools/linear_solver/linear_solver.pb.h"

namespace minimip {
namespace {

using ::google::protobuf::TextFormat;
using ::google::protobuf::util::MessageDifferencer;
using ::operations_research::MPModelProto;
using ::operations_research::MPSolutionResponse;

TEST(ExampleTest, SolvesSimpleLP) {
  // We optimize a simple LP:
  //   maximize 2.0 x - 0.5 y
  //   x <= 13
  //   y >= 42
  MPModelProto model;
  ASSERT_TRUE(TextFormat::ParseFromString(
      R"pb(
        maximize: true
        variable { upper_bound: 13.0 objective_coefficient: 2.0 }
        variable { lower_bound: 42.0 objective_coefficient: -0.5 }
      )pb",
      &model));

  const MPSolutionResponse observed =
      SolveMPModelProtoWithGlop(std::move(model));

  // NOTE(pawell): Modify the expected response (e.g., change a variable_value)
  // to see how the tests' failures (and proto differences) are reported.
  MPSolutionResponse expected;
  ASSERT_TRUE(TextFormat::ParseFromString(R"pb(
                                            status: MPSOLVER_OPTIMAL
                                            objective_value: 5
                                            variable_value: [ 13, 42 ]
                                            reduced_cost: [ 2, -0.5 ]
                                          )pb",
                                          &expected));

  std::string differences;
  MessageDifferencer differencer;
  differencer.ReportDifferencesToString(&differences);
  // TODO(ondrasej): Include `EqualsProto` from EXEgesis open-source code.
  EXPECT_TRUE(differencer.Compare(observed, expected)) << differences;
}

}  // namespace
}  // namespace minimip
