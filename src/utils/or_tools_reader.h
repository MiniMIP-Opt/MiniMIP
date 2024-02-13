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

#ifndef SRC_OR_TOOLS_READER_H_
#define SRC_OR_TOOLS_READER_H_

#include "src/data_structures/problem.h"

namespace minimip {

class OrToolsReader {
 public:
  OrToolsReader() = default;

  static absl::StatusOr<MiniMipProblem> ReadMipProblemFromMPSFile(
      const std::string& file_path);

  static absl::StatusOr<MiniMipProblem> ReadMipProblemFromString(
      const std::string& file_data);
};

}  // namespace minimip
#endif  // SRC_OR_TOOLS_READER_H_
