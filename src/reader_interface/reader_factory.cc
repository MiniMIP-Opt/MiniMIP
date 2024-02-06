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

#include "src/reader_interface/reader_factory.h"

#include "ortools/base/status_macros.h"
#include "src/reader_interface/or_tools_reader.h"

namespace minimip {

absl::StatusOr<std::unique_ptr<ReaderInterface>> CreateReader(
    const ReaderParameters& params) {
  std::unique_ptr<ReaderInterface> reader;
  switch (params.reader_type()) {
    case ReaderParameters::OR_TOOLS:
      reader = std::make_unique<OrToolsReader>();
      break;
    default:
      return util::InvalidArgumentErrorBuilder()
             << "Invalid lp solver type: " << params.reader_type();
  }

  return reader;
}

}  // namespace minimip
