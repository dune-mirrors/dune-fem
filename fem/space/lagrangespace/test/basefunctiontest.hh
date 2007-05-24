#ifndef DUNE_LAGRANGESPACE_BASEFUNCTIONTEST_HH
#define DUNE_LAGRANGESPACE_BASEFUNCTIONTEST_HH

#include <string>

#define TEST_SECOND_ORDER

#include <dune/fem/misc/test.hh>

namespace Dune {

  class LagrangeBase_Test : public Test {
  public:
    LagrangeBase_Test(std::string gridFile) :
      gridFile_(gridFile)
    {}

    virtual void run();

    void testBaseFunctions();

  private:
    template <class SpaceType>
    void checkLagrangeBase(const SpaceType&);
    
    std::string gridFile_;
  };

}

#endif
