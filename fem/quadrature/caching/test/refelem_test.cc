#include <config.h>

#include "refelem_test.hh"

namespace Dune {
  
  void ReferenceElement_Test::run() 
  {
    globalTest();
    allGeometriesTest();
  }

  void ReferenceElement_Test::globalTest() 
  {
    //- Tests for Codim 1
    //- map all corners to reference element and check identity
    checkCorners<3, 1>(refCube3_);
    checkCorners<3, 1>(refSimplex3_);
    checkCorners<3, 2>(refCube3_);
    checkCorners<3, 2>(refSimplex3_);
    checkCorners<2, 1>(refCube2_);
    checkCorners<2, 1>(refSimplex2_);

    //- map point (0.1, 0.2) resp (0.1) into referece element and check correctness
    FieldVector<double, 1> vec1(0.1);
    FieldVector<double, 2> vec2(0.1); vec2[1] = 0.2;

    FieldVector<double, 3> result3;
    FieldVector<double, 2> result2;

    result3 = refCube3_.global<1>(vec2, 1, 1);
    _floatTest(result3[0], 1.0);
    _floatTest(result3[1], 0.1);
    _floatTest(result3[2], 0.2);
    
    result3 = refSimplex3_.global<1>(vec2, 0, 1);
    _floatTest(result3[0], 0.7);
    _floatTest(result3[1], 0.1);
    _floatTest(result3[2], 0.2);
 
    result3 = refCube3_.global<2>(vec1, 6, 2);
    _floatTest(result3[0], 0.0);
    _floatTest(result3[1], 0.1);
    _floatTest(result3[2], 1.0);
  
    result3 = refSimplex3_.global<2>(vec1, 4, 2);
    _floatTest(result3[0], 0.9);
    _floatTest(result3[1], 0.0);
    _floatTest(result3[2], 0.1);
    
    result2 = refCube2_.global<1>(vec1, 1, 1);
    _floatTest(result2[0], 1.0);
    _floatTest(result2[1], 0.1);
 
    result2 = refSimplex2_.global<1>(vec1, 2, 1);
    _floatTest(result2[0], 0.1);
    _floatTest(result2[1], 0.0);
 
  }

  void ReferenceElement_Test::allGeometriesTest() {
    const ReferenceCube<double, 3>& rc3 = 
      ReferenceElements<double, 3>::cube(cube3_);
    checkSingle<ReferenceCube<double, 3>, 1>(rc3);
    checkSingle<ReferenceCube<double, 3>, 2>(rc3);

    const ReferenceCube<double, 2>& rc2 = 
      ReferenceElements<double, 2>::cube(cube2_);
    checkSingle<ReferenceCube<double, 2>, 1>(rc2);

    const ReferenceSimplex<double, 3>& rs3 =
      ReferenceElements<double, 3>::simplices(simplex3_);
    checkSingle<ReferenceSimplex<double, 3>, 1>(rs3);
    checkSingle<ReferenceSimplex<double, 3>, 2>(rs3);

    const ReferenceSimplex<double, 2>& rs2 =
      ReferenceElements<double, 2>::simplices(simplex2_);
    checkSingle<ReferenceSimplex<double, 2>, 1>(rs2);

    GeometryType prism(GeometryType::prism,3);
    const ReferencePrism<double, 3>& rpr = 
      ReferenceElements<double, 3>::prism(prism);
    checkSingle<ReferencePrism<double, 3>, 1>(rpr);
    checkSingle<ReferencePrism<double, 3>, 2>(rpr);

    //const ReferencePyramid<double, 3>& rpy =
    //  ReferenceElements<double, 3>::pyram(pyramid);
    //checkSingle<ReferencePyramid<double, 3>, 1>(rpy);
    //checkSingle<ReferencePyramid<double, 3>, 2>(rpy);
      
  }

  template <int dim, int codim>
  void ReferenceElement_Test::
  checkCorners(const ReferenceElement<double, dim>& refElem) 
  {
    const ReferenceElement<double, dim-codim>& refCd = 
      ReferenceElements<double, dim-codim>::general(refElem.type(0, codim));

    for (int subEn = 0; subEn < refElem.size(codim); ++subEn) {
      for (int corner = 0; corner < refCd.size(dim-codim); ++corner) {
        int globalCorner = refElem.subEntity(subEn, codim, corner, dim);
        FieldVector<double, dim> cornerVec = 
          refElem.position(globalCorner, dim);
        FieldVector<double, dim> global = 
          refElem.template global<codim>(refCd.position(corner, dim-codim), subEn, codim);
        for (int d = 0; d < dim; ++d) {
          _floatTest(cornerVec[d], global[d]);
        }
        //std::cout << cornerVec << " == " << global << std::endl;
      }
    }
    
  }

  template <typename RefElemType, int codim>
  void ReferenceElement_Test::checkSingle(const RefElemType& ref)
  {
    const int dim = RefElemType::d;
    
    const ReferenceElement<double, dim-codim>& refCd = 
      ReferenceElements<double, dim-codim>::general(ref.type(0, codim));

    int globalCorner = ref.subEntity(0, codim, 0, dim);
    FieldVector<double, dim> cornerVec = ref.position(globalCorner, dim);
    FieldVector<double, dim> global = 
      ref.template global<codim>(refCd.position(0, dim-codim), 0, codim);

    for (int d = 0; d < dim; ++d) {
      _floatTest(cornerVec[d], global[d]);
    }
  }
  
} // end namespace Dune
