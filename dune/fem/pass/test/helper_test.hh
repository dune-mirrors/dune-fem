#ifndef DUNE_PASSHELPER_TEST_HH
#define DUNE_PASSHELPER_TEST_HH

#include <dune/fem/misc/test.hh>

#include "../selection.hh"
#include "../callerutility.hh"

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/common/utility.hh>

#include "lagrange_fixture.hh"

namespace Dune {

  class PassHelper_Test : public Test 
  {
  public:
    PassHelper_Test(std::string name) : gridFile_(name) {}

    virtual void run();

  private:
    void selectorTest();
    void filterTest();
    void tupleConverterTest();
  
  private:
    std::string gridFile_;
  };

}

#endif
