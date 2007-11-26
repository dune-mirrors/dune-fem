#ifndef DUNE_POINT_TEST_HH
#define DUNE_POINT_TEST_HH

#include <dune/fem/misc/test.hh>
#include "../pointprovider.hh"

namespace Dune {

  class PointProvider_Test : public Test 
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

} // end namespace Dune

#endif
