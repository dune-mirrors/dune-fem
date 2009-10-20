#ifndef DUNE_REFELEM_TEST_HH
#define DUNE_REFELEM_TEST_HH

#include <dune/grid/common/referenceelements.hh>
#include <dune/common/fvector.hh>

#include <dune/fem/misc/test.hh>


namespace Dune {

  class ReferenceElement_Test : public Test 
  {
    GeometryType cube3_;
    GeometryType cube2_;
    GeometryType simplex3_;
    GeometryType simplex2_;
    GeometryType simplex1_;
  public:
    ReferenceElement_Test() :
      cube3_(GeometryType(GeometryType::cube,3)),
      cube2_(GeometryType(GeometryType::cube,2)),
      simplex3_(GeometryType(GeometryType::simplex,3)),
      simplex2_(GeometryType(GeometryType::simplex,2)),
      simplex1_(GeometryType(GeometryType::simplex,1)),
      refCube3_(ReferenceElements<double, 3>::general(cube3_)),
      refSimplex3_(ReferenceElements<double, 3>::general(simplex3_)),
      refCube2_(ReferenceElements<double, 2>::general(cube2_)),
      refSimplex2_(ReferenceElements<double, 2>::general(simplex2_)),
      refLine_(ReferenceElements<double, 1>::general(simplex1_))
    {}

    virtual void run();

    void globalTest();
    void allGeometriesTest();

  private:
    template <int dim, int codim>
    void checkCorners(const ReferenceElement<double, dim>& refElem);

    template <typename RefElemType, int codim>
    void checkSingle(const RefElemType& ref);

  private:
    const ReferenceElement<double, 3>& refCube3_;
    const ReferenceElement<double, 3>& refSimplex3_;
    const ReferenceElement<double, 2>& refCube2_;
    const ReferenceElement<double, 2>& refSimplex2_;
    const ReferenceElement<double, 1>& refLine_;
  };

} // end namespace Dune

#endif
