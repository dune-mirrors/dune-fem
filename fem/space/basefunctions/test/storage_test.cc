#include "storage_test.hh"

#include <dune/grid/alu3dgrid.hh>
#include <dune/grid/utility/twistutility.hh>
#include <dune/common/functionspace.hh>

#include "basefunctionstub.hh"
#include "../basefunctionstorage.hh"
#include "../../quadrature/cachequad.hh"
#include "alugrid_fixture.hh"

namespace Dune {

  void Storage_Test::run() 
  {
    storageComparisonTest();
  }

  void Storage_Test::storageComparisonTest() 
  {
    const int dim = 3;
    const int codim = 1;

    typedef FunctionSpace<double, double, dim, 1> FunctionSpaceType;
    typedef FunctionSpaceType::RangeType RangeType;
    typedef FunctionSpaceType::JacobianRangeType JacobianRangeType;
    typedef FunctionSpaceType::DomainType DomainType;
    typedef CachingStorage<FunctionSpaceType> CachingStorageType;
    typedef SimpleStorage<FunctionSpaceType> SimpleStorageType;
    typedef BaseFunctionStubFactory<FunctionSpaceType> FactoryType;
 
    typedef ALUGridFixture<hexa> FixtureType;
    typedef FixtureType::GridType GridType;
    typedef GridType::Traits::IntersectionIterator IntersectionIteratorType;
    typedef CachingQuadrature<GridType, codim> QuadratureType;
    typedef CachingQuadrature<GridType, 0> VolumeQuadratureType;

    RangeType phiSimple;
    RangeType phiCaching;

    JacobianRangeType gradPhiSimple;
    JacobianRangeType gradPhiCaching;

    FieldVector<int, 0> diffVar0;
    FieldVector<int, 1> diffVar1(dim-1);

    FactoryType factory(simplex);
    SimpleStorageType simple(factory);
    CachingStorageType caching(factory);

    FixtureType fix(gridFile_);
    GridType& grid = fix.grid();

    const int numBaseFunctions = factory.numBaseFunctions();

    VolumeQuadratureType volQuad(*grid.leafbegin<0>(), 2);
    for (int q = 0; q < volQuad.nop(); ++q) {
      for (int i = 0; i < numBaseFunctions; ++i) {
        simple.evaluate(i, diffVar0, volQuad, q, phiSimple);
        caching.evaluate(i, diffVar0, volQuad, q, phiCaching);
        
        _floatTest(phiSimple[0], phiCaching[0]);
        std::cout << phiSimple[0] << " =v= " << phiCaching[0] << std::endl;
        
        simple.evaluate(i, diffVar1, volQuad, q, phiSimple);
        caching.evaluate(i, diffVar1, volQuad, q, phiCaching);
        
        _floatTest(phiSimple[0], phiCaching[0]);
        std::cout << phiSimple[0] << " =d= " << phiCaching[0] << std::endl;
        
        simple.jacobian(i, volQuad, q, gradPhiSimple);
        caching.jacobian(i, volQuad, q, gradPhiCaching);
        
        for (int i = 0; i < dim; ++i) {
          _floatTest(gradPhiSimple[i][0], gradPhiCaching[i][0]);
          std::cout << gradPhiSimple[i][0] << " =j= " 
                    << gradPhiCaching[i][0] << std::endl;
        }
      }
    }

    IntersectionIteratorType endit = grid.leafbegin<0>()->iend();
    IntersectionIteratorType it = grid.leafbegin<0>()->ibegin();
    for (; it != endit; ++it) {
      TwistUtility<GridType> twistUtil(grid);
      int twist = twistUtil.twistInSelf(it);
      //std::cout << "Twist == " << twist << std::endl;
      QuadratureType quad(it, 4, twist, ElementQuadrature<GridType,1>::INSIDE);
  
      for (int q = 0; q < quad.nop(); ++q) {
        std::cout << "x == " << quad.point(q) << std::endl;

        for (int i = 0; i < numBaseFunctions; ++i) {
          simple.evaluate(i, diffVar0, quad, q, phiSimple);
          caching.evaluate(i, diffVar0, quad, q, phiCaching);

          _floatTest(phiSimple[0], phiCaching[0]);
          std::cout << phiSimple[0] << " =v= " << phiCaching[0] << std::endl;

          simple.evaluate(i, diffVar1, quad, q, phiSimple);
          caching.evaluate(i, diffVar1, quad, q, phiCaching);

          _floatTest(phiSimple[0], phiCaching[0]);
          std::cout << phiSimple[0] << " =d= " << phiCaching[0] << std::endl;

          simple.jacobian(i, quad, q, gradPhiSimple);
          caching.jacobian(i, quad, q, gradPhiCaching);
          
          for (int i = 0; i < dim; ++i) {
            _floatTest(gradPhiSimple[i][0], gradPhiCaching[i][0]);
            std::cout << gradPhiSimple[i][0] << " =j= " 
                      << gradPhiCaching[i][0] << std::endl;
          }
        }
      
      }
    


    }

  }

} // end namespace Dune
