#include "storage_test.hh"

#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/quadrature/cachequad.hh>

#include "basefunctionstub.hh"
#include "../basefunctionstorage.hh"
#include "prismgrid_fixture.hh"
#include <dune/fem/quadrature/caching/twistutility.hh>
#include <dune/fem/gridpart/gridpart.hh>

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
 
    typedef PrismGridFixture<dim> FixtureType;
    typedef FixtureType::GridType GridType;
    typedef LeafGridPart< GridType > GridPartType;
    typedef GridPartType :: IntersectionIteratorType IntersectionIteratorType;
    typedef GridPartType :: Codim<0> :: IteratorType IteratorType; 
    typedef GridType :: Codim<0> :: Entity EntityType;
    
    typedef CachingQuadrature<GridPartType, codim> QuadratureType;
    typedef CachingQuadrature<GridPartType, 0> VolumeQuadratureType;

    RangeType phiSimple;
    RangeType phiCaching;

    JacobianRangeType gradPhiSimple;
    JacobianRangeType gradPhiCaching;

    FieldVector<int, 0> diffVar0;
    FieldVector<int, 1> diffVar1(dim-1);

    GeometryType geoType (GeometryType::prism, dim );
    
    FactoryType factory(geoType);
    SimpleStorageType simple(factory);
    CachingStorageType caching(factory);

    FixtureType fix(gridFile_);
    GridType& grid = fix.grid();

    GridPartType gridPart( grid );

    const int numBaseFunctions = factory.numBaseFunctions();

#if 0
    VolumeQuadratureType volQuad( *gridPart.begin<0>(), 2);
    for (int q = 0; q < volQuad.nop(); ++q) 
    {
      for (int i = 0; i < numBaseFunctions; ++i) 
      {
        simple.evaluate(i, diffVar0, volQuad[ q ], phiSimple);
        caching.evaluate(i, diffVar0, volQuad[ q ], phiCaching);
        
        _floatTest(phiSimple[0], phiCaching[0]);
        std::cout << phiSimple[0] << " =v= " << phiCaching[0] << std::endl;
        
        simple.jacobian(i, volQuad[ q ], gradPhiSimple);
        caching.jacobian(i, volQuad[ q ], gradPhiCaching);
        
        for (int i = 0; i < dim; ++i) {
          _floatTest(gradPhiSimple[i][0], gradPhiCaching[i][0]);
          std::cout << gradPhiSimple[i][0] << " =j= " 
                    << gradPhiCaching[i][0] << std::endl;
        }
      }
    }
#endif

    IteratorType grit = gridPart.begin<0>();
    EntityType & en = *grit ;
    
#if 1
    IntersectionIteratorType endit = gridPart.iend( en );
    for (IntersectionIteratorType it = gridPart.ibegin( en ); 
         it != endit; ++it) 
    {
      //std::cout << "Twist == " << twist << std::endl;
      QuadratureType quad(gridPart, it, 4, QuadratureType :: INSIDE);
  
      /*
      for (int q = 0; q < quad.nop(); ++q) 
      {
        std::cout << "x == " << quad.point(q) << std::endl;

        for (int i = 0; i < numBaseFunctions; ++i) 
        {
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
          
          for (int i = 0; i < dim; ++i) 
          {
            _floatTest(gradPhiSimple[i][0], gradPhiCaching[i][0]);
            std::cout << gradPhiSimple[i][0] << " =j= " 
                      << gradPhiCaching[i][0] << std::endl;
          }
        }
      
      }
      */
    }
#endif 
  }

} // end namespace Dune
