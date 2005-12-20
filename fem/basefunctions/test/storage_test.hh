#ifndef DUNE_STORAGE_TEST_HH
#define DUNE_STORAGE_TEST_HH

#include <string>

#include "../../misc/test.hh"

namespace Dune {

  class Storage_Test : public Test {
  public:
    Storage_Test(std::string gridFile) :
      gridFile_(gridFile)
    {}

    virtual void run();

    void storageComparisonTest();

  private:
    std::string gridFile_;
  };

}

#endif
