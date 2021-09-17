
#include "src/messagehandler/message_handler.h"
#include <cstdio> /** TODO: printf should become cout and  include -> iostream (c++ conform)
                    *
                    * */
namespace minimip {

  void messagehandler::PrintDebugHeader(
    const char* source_file, /**< name of the source file that called the function */
    int source_line          /**< line in the source file where the function was called */
    ) {
    (void) printf("[%s:%d] DEBUG: ", source_file, source_line);
  }
  void messagehandler::PrintDebugMessage(){
    //to be implemented
  }

  void messagehandler::PrintErrorHeader(
    const char* source_file, /**< name of the source file that called the function */
    int source_line          /**< line in the source file where the function was called */
    ) {
    (void) printf("[%s:%d] ERROR: ", source_file, source_line);
  }
  void messagehandler::PrintError(){
    // to be implemented
  }

} // namespace minimip
