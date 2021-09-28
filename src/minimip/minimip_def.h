
#ifndef MINIMIP_SRC_MINIMIP_DEF_H_
#define MINIMIP_SRC_MINIMIP_DEF_H_

#include "src/messagehandler/message_handler.h"
#include "src/messagehandler/message_macros.h"
#include "src/lp_interface/lp_types.h"

/* #define MINIMIP_VERSION 000 / **< MiniMIP version number (multiplied by 100 to get integer number) * /
  #define MINIMIP_SUBVERSION               0 / **< MiniMIP sub version number * /
  #define MINIMIP_APIVERSION              00 / **< MiniMIP API version number * /
  #define MINIMIP_COPYRIGHT   PLACEHOLDER:"~~Copyright (C) 2002-2020 Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)~~" */

#define EPSCEIL(x, eps) (ceil((x) - (eps)))
#define EPSFLOOR(x, eps) (floor((x) + (eps)))
#define REALABS(x)        (fabs(x))

#define MINIMIP_CALL(x)                                               \
  do {                                                                \
    minimip::RetCode _restat_;                                         \
    if ((_restat_ = (x)) != minimip::RetCode::kOkay) {                 \
    MiniMIPerrorMessage("Error <%d> in function call\n", static_cast<int>(_restat_)); \
      return _restat_;                                                \
    }                                                                 \
  } while (false)

#endif /* MINIMIP_SRC_MINIMIP_DEF_H_ */
