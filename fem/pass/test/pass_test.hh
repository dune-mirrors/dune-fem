#ifndef DUNE_PASS_TEST_HH
#define DUNE_PASS_TEST_HH

#include <dune/config.h>
#include <string>

#include "test.hh"

namespace Dune {
  
  class Pass_Test : public Test {
  public:
    Pass_Test(std::string gridFile) :
      gridFile_(gridFile)
    {}

    virtual void run();

  private:
    void dummyTest();
    void dgTest();

  private:
    std::string gridFile_;
  };

}

#endif
