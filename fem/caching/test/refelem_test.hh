#ifndef DUNE_REFELEM_TEST_HH
#define DUNE_REFELEM_TEST_HH

#include <dune/grid/common/referenceelements.hh>
#include <dune/common/fvector.hh>

#include "../../misc/test.hh"


namespace Dune {

  class ReferenceElement_Test : public Test 
  {
  public:
    ReferenceElement_Test() :
      refCube3_(ReferenceElements<double, 3>::general(cube)),
      refSimplex3_(ReferenceElements<double, 3>::general(simplex)),
      refCube2_(ReferenceElements<double, 2>::general(cube)),
      refSimplex2_(ReferenceElements<double, 2>::general(simplex)),
      refLine_(ReferenceElements<double, 1>::general(simplex))
    {}

    virtual void run();

    void globalTest();

  private:
    template <int dim, int codim>
    void checkCorners(const ReferenceElement<double, dim>& refElem);

  private:
    const ReferenceElement<double, 3>& refCube3_;
    const ReferenceElement<double, 3>& refSimplex3_;
    const ReferenceElement<double, 2>& refCube2_;
    const ReferenceElement<double, 2>& refSimplex2_;
    const ReferenceElement<double, 1>& refLine_;
  };

}

#endif
