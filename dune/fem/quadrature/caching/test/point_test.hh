#ifndef DUNE_POINT_TEST_HH
#define DUNE_POINT_TEST_HH

#include "../pointprovider.hh"

namespace Dune {
  namespace Fem {

  class PointProvider_Test
  {
  public:
    virtual void run();

    void codim0Test();
    void sameOutputTest();
    void transformationTest();

  private:
    template <class PointType>
    bool findPoint(const PointType& p, const std::vector<PointType>& vec);
  };

  } // end namespace Fem
} // end namespace Dune

#endif
