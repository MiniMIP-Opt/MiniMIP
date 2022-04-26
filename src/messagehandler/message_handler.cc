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

#include "src/messagehandler/message_handler.h"

#include <cstdio>

namespace minimip {

void messagehandler::PrintDebugHeader(
    const char*
        source_file,  // name of the source file that called the function
    int source_line   // line in the source file where the function was called
) {
  (void)printf("[%s:%d] DEBUG: ", source_file, source_line);
}
void messagehandler::PrintDebugMessage() {
  // to be implemented
}

void messagehandler::PrintErrorHeader(
    const char*
        source_file,  // name of the source file that called the function
    int source_line   // line in the source file where the function was called
) {
  (void)printf("[%s:%d] ERROR: ", source_file, source_line);
}
void messagehandler::PrintError() {
  // to be implemented
}

}  // namespace minimip
