
#ifndef MINIMIP_SRC_MESSAGEHANDLER_MESSAGE_HANDLER_H_
#define MINIMIP_SRC_MESSAGEHANDLER_MESSAGE_HANDLER_H_

namespace minimip {

class messagehandler {
 public:
  void PrintDebugHeader(
    const char* source_file, /**< name of the source file that called the function */
    int source_line          /**< line in the source file where the function was called */
  );

  void PrintDebugMessage();

  void PrintErrorHeader(
    const char* source_file, /**< name of the source file that called the function */
    int source_line          /**< line in the source file where the function was called */
  );

  void PrintError();
};

} /* namespace minimip */

#endif /* MINIMIP_SRC_MESSAGEHANDLER_MESSAGE_HANDLER_H_ */
