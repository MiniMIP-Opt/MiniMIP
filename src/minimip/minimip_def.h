
#ifndef SRC_MINIMIP_MINIMIP_DEF_H_
#define SRC_MINIMIP_MINIMIP_DEF_

#include "absl/status/status.h"
#include "src/lp_interface/lp_types.h"
#include "src/messagehandler/message_handler.h"
#include "src/messagehandler/message_macros.h"

// #define MINIMIP_VERSION 000   // MiniMIP version number (multiplied by 100 to
// get integer number)
//  #define MINIMIP_SUBVERSION               0 // MiniMIP sub version number
//  #define MINIMIP_APIVERSION              00 // MiniMIP API version number
//  #define MINIMIP_COPYRIGHT   // PLACEHOLDER:"~~Copyright (C) 2021
//  Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)~~"

#define EPSCEIL(x, eps) (ceil((x) - (eps)))
#define EPSFLOOR(x, eps) (floor((x) + (eps)))
#define REALABS(x) (fabs(x))

#define MINIMIP_CALL(x)                                    \
  {                                                        \
    absl::Status _restat_ = (x);                           \
    if (!_restat_.ok()) {                                  \
      MiniMIPerrorMessage("Error <%d> in function call\n", \
                          (int)_restat_.code());           \
      return _restat_;                                     \
    }                                                      \
  }

#endif  // SRC_MINIMIP_MINIMIP_DEF_H_
