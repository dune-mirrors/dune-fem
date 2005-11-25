#include "test.hh"

#include "../quadrature.hh"
#include <dune/quadrature/fixedorder.hh>

namespace Dune {
  
  class Quad_Test : public Test 
  {
  public:
    virtual void run();

    template <class Quad, class Fixed>
    void fixedOrderComparisonExec(Quad& quad, Fixed& fixed);
    template <class Quad>
    void weightSummationExec(Quad& quad);
    template <class Quad>
    void integrationExec(Quad& quad);
    template <class Quad>
    void orderExec(Quad& quad, int order);
    template <class Quad>
    void sameGeometryExec(Quad& quad1, Quad& quad2);
    void indicesTest();

  private:
  };

}
