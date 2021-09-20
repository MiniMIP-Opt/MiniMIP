
#ifndef MINIMIP_SRC_MESSAGEHANDLER_MESSAGE_MACROS_H_
#define MINIMIP_SRC_MESSAGEHANDLER_MESSAGE_MACROS_H_

#include "src/messagehandler/message_handler.h"

#ifdef MINIMIP_DEBUG
/** prints a debugging message if MINIMIP_DEBUG flag is set */
#define MiniMIPdebugMessage                                 \
  messagehandler::PrintDebugHeader(__FILENAME__, __LINE__); \
  messagehandler::PrintDebugMessage
#else
/** prints a debugging message if MINIMIP_DEBUG flag is set */
#define MiniMIPdebugMessage      \
  while (false) /*lint -e{530}*/ \
  printf
#endif

#ifdef MINIMIP_DEBUG
/** prints an error message */
#define MiniMIPerrorMessage                                 \
  messagehandler::PrintErrorHeader(__FILENAME__, __LINE__); \
  messagehandler::PrintError
#else
/** prints a debugging message if MINIMIP_DEBUG flag is set - also consider using MiniMIPdebugMsg/MiniMIPsetDebugMsg */
#define MiniMIPerrorMessage      \
  while (false) /*lint -e{530}*/ \
  printf
#endif

#endif /* MINIMIP_SRC_MESSAGEHANDLER_MESSAGE_MACROS_H_ */