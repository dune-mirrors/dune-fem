#ifndef DUNE_PASSHELPER_TEST_HH
#define DUNE_PASSHELPER_TEST_HH

#include "test.hh"

#include <dune/fem/pass/selection.hh>
#include <dune/fem/pass/caller.hh>

#include <dune/fem/dfadapt.hh>
#include <dune/common/utility.hh>

#include "lagrange_fixture.hh"

namespace Dune {

  class PassHelper_Test : public Test 
  {
  public:
    PassHelper_Test(std::string name) : gridFile_(name) {}

    virtual void run();

  private:
    void filterTest();
    void tupleConverterTest();
    void functorTest();
    void callerTest();

  private:
    std::string gridFile_;
  };

}

#endif
