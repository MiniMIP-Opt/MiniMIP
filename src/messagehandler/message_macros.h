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

#ifndef SRC_MESSAGEHANDLER_MESSAGE_MACROS_H_
#define SRC_MESSAGEHANDLER_MESSAGE_MACROS_H_

#include "src/messagehandler/message_handler.h"

#ifdef MINIMIP_DEBUG
// prints a debugging message if MINIMIP_DEBUG flag is set
#define MiniMIPdebugMessage                                 \
  messagehandler::PrintDebugHeader(__FILENAME__, __LINE__); \
  messagehandler::PrintDebugMessage
#else
// prints a debugging message if MINIMIP_DEBUG flag is set
#define MiniMIPdebugMessage \
  while (false) printf
#endif

#ifdef MINIMIP_DEBUG
// prints an error message
#define MiniMIPerrorMessage                                 \
  messagehandler::PrintErrorHeader(__FILENAME__, __LINE__); \
  messagehandler::PrintError
#else
// prints a debugging message if MINIMIP_DEBUG flag is set - also consider using
// MiniMIPdebugMsg/MiniMIPsetDebugMsg
#define MiniMIPerrorMessage \
  while (false) printf
#endif

#endif  // SRC_MESSAGEHANDLER_MESSAGE_MACROS_H_
