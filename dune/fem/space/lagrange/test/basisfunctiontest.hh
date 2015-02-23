#ifndef DUNE_FEM_LAGRANGESPACE_BASEFUNCTIONTEST_HH
#define DUNE_FEM_LAGRANGESPACE_BASEFUNCTIONTEST_HH

#include <string>

#define TEST_SECOND_ORDER

namespace Dune
{

  namespace Fem
  {

  class LagrangeBasis_Test {
  public:
    LagrangeBasis_Test(std::string gridFile) :
      gridFile_(gridFile)
    {}

    virtual void run();

    void testBasisFunctions();

  private:
    template <class SpaceType>
    void checkLagrangeBasis(const SpaceType&);

    std::string gridFile_;
  };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_LAGRANGESPACE_BASEFUNCTIONTEST_HH
