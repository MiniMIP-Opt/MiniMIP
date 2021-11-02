#ifndef SRC_LP_INTERFACE_LPI_FACTORY_H_
#define SRC_LP_INTERFACE_LPI_FACTORY_H_

#include "src/lp_interface/lpi_glop.h"
#include "src/lp_interface/lpi_soplex.h"

namespace minimip {

enum class InterfaceCode {
  kGlop   = 0,
  kSoplex = 1,
};

class LPInterfaceFactory {
 public:
  LPInterface* CreateLPInterface(InterfaceCode interface_code) {
    switch (interface_code) {
      case InterfaceCode::kSoplex:
        return new LPSoplexInterface();
      default:
        return new LPGlopInterface();
    }
  }
};

}  // namespace minimip
#endif  // SRC_LP_INTERFACE_LPI_FACTORY_H_
