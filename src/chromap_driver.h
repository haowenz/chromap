#ifndef CHROMAP_DRIVER_H_
#define CHROMAP_DRIVER_H_

namespace chromap {

class ChromapDriver {
 public:
  ChromapDriver() = default;
  ~ChromapDriver() = default;
  void ParseArgsAndRun(int argc, char *argv[]);
};

}  // namespace chromap

#endif  // CHROMAP_DRIVER_H_
