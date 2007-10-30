#include <config.h>

#include <dune/grid/alugrid.hh>
#include <dune/grid/common/gridpart.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>

#include "cachequad_test.hh"
#include "albertagrid_fixture.hh"
#include "alugrid_fixture.hh"
#include "sgrid_fixture.hh"

namespace Dune {

  void CachingQuadrature_Test::run() 
  {
    codim0Test();
    codim1AlbertaTest();
    codim1ALUHexaTest();
    codim1ALUTetraTest();
    codim1SGridTest();
  }

  void CachingQuadrature_Test::codim0Test() 
  {
    const int dim = GRIDDIM;
    const int codim = 0;
    GeometryType triangle ( GeometryType::simplex, 2);

    typedef AlbertaGridFixture<dim, dim> GridFixtureType;
    typedef GridFixtureType::GridType GridType;
    typedef LeafGridPart< GridType > GridPartType;
    typedef CachingQuadrature<GridPartType, codim> QuadratureType;
    typedef PointProvider<double, dim, codim> PointProviderType;
    typedef PointProviderType::GlobalPointVectorType PointVectorType;

    GridFixtureType fix(albertaGridFile_);
    GridType& grid = fix.grid();
    GridPartType gridPart( grid );
    
    QuadratureType quad(*grid.leafbegin<0>(), 3);
    const PointVectorType& points = 
      PointProviderType::getPoints(quad.id(), triangle);

    _test(points.size() == (size_t)quad.nop());
    
    for (int i = 0; i < quad.nop(); ++i) {
      for (int d = 0; d < dim; ++d) {
        _floatTest(points[i][d], quad.point(quad.cachingPoint(i))[d]);
      }
    }
  } 

  void CachingQuadrature_Test::codim1AlbertaTest() 
  {
    const int dim = GRIDDIM ;
    const int codim = 1;
    const GeometryType tetrahedron( GeometryType::simplex, dim);

    typedef AlbertaGridFixture<dim, dim> GridFixtureType;
    typedef GridFixtureType::GridType GridType;
    typedef LeafGridPart< GridType > GridPartType;
    typedef GridPartType::IntersectionIteratorType IntersectionIterator;
    typedef GridPartType :: Codim<0> :: IteratorType IteratorType;
    typedef CachingQuadrature<GridPartType, codim> QuadratureType;
    typedef PointProvider<GridType::ctype, dim, codim> PointProviderType;
    typedef PointProviderType::GlobalPointVectorType PointVectorType;
    typedef IntersectionIterator::LocalGeometry LocalGeometryType;
    typedef IntersectionIterator::Geometry GlobalGeometryType;
    typedef FieldVector<double,dim> DomainType;

    GridFixtureType fix(albertaGridFile_);
    GridType& grid = fix.grid();
    GridPartType gridPart( grid );

    const int quadOrd = 4;

    for(int l=0; l<3; ++l) 
    {
      IteratorType enditer = gridPart.end<0> ();
      for(IteratorType eiter = gridPart.begin<0> (); 
          eiter != enditer; ++eiter)
      {
        IntersectionIterator endit = gridPart.iend( *eiter );
        for (IntersectionIterator it = gridPart.ibegin(*eiter );
             it != endit; ++it) 
        {
          const LocalGeometryType& geo = it.intersectionSelfLocal();
          QuadratureType quad(gridPart, it, quadOrd , QuadratureType::INSIDE);

          const PointVectorType& points = 
            PointProviderType::getPoints(quad.id(), tetrahedron);

          _test(points.size() == (size_t)((dim +1) * quad.nop()));

          typedef TwistUtility<GridType> TwistUtilityType; 
          //std::cout << "New Intersection: Twists: ";
          //std::cout << TwistUtilityType :: twistInSelf( grid, it ) << " ";
          //std::cout << TwistUtilityType :: twistInNeighbor( grid, it ) << "\n";

          for (int i = 0; i < quad.nop(); ++i) 
          {
            for (int d = 0; d < dim; ++d) 
            {
              _floatTest(points[quad.cachingPoint(i)][d],
                         geo.global(quad.localPoint(i))[d]);
            }
            //std::cout << "nis: " << it.numberInSelf() << " nin: " << it.numberInNeighbor();
            //std::cout << " pt " << i << ": " << points[quad.cachingPoint(i)]
            //          << " == " << geo.global(quad.localPoint(i)) << std::endl;
          }

          if( it.neighbor ())
          {
            if( TwistUtilityType::conforming( grid, it ))
            {
              const GlobalGeometryType& gGeo = it.intersectionGlobal();
              QuadratureType outerQuad(gridPart, it, quadOrd , QuadratureType::OUTSIDE);
              
              for (int i = 0; i < quad.nop(); ++i) 
              {
                DomainType p = gGeo.global( quad.localPoint(i) );
                DomainType q = gGeo.global( outerQuad.localPoint(i) );

                for (int d = 0; d < dim; ++d) 
                {
                  _floatTest(p[d],q[d]);
                }
              }
            }
          }

        } // end iterator loop
      }
      grid.globalRefine(1);
    }
  }

  
  void CachingQuadrature_Test::codim1ALUHexaTest() 
  {
    const int dim = 3;
    const int codim = 1;
    const GeometryType hexahedron( GeometryType::cube, dim);

    typedef ALUCubeGridFixture GridFixtureType;
    typedef GridFixtureType::GridType GridType;
    typedef LeafGridPart< GridType > GridPartType;
    typedef GridPartType::IntersectionIteratorType IntersectionIterator;
    typedef GridPartType:: Codim<0> :: IteratorType IteratorType;
    typedef CachingQuadrature<GridPartType, codim> QuadratureType;
    typedef PointProvider<GridType::ctype, dim, codim> PointProviderType;
    typedef PointProviderType::GlobalPointVectorType PointVectorType;
    typedef IntersectionIterator::LocalGeometry LocalGeometryType;
    typedef IntersectionIterator::Geometry GlobalGeometryType;
    typedef FieldVector<double,dim> DomainType;

    GridFixtureType fix(aluGridHexaFile_);
    GridType& grid = fix.grid();
    GridPartType gridPart( grid );

    const int quadOrd = 4;

    for(int l=0; l<3; ++l) 
    {
      IteratorType enditer = gridPart.end<0> ();
      for(IteratorType eiter = gridPart.begin<0> (); 
          eiter != enditer; ++eiter)
      {
        IntersectionIterator endit = gridPart.iend( *eiter );
        for (IntersectionIterator it = gridPart.ibegin(*eiter);
            it != endit; ++it) 
        {
          const LocalGeometryType& geo = it.intersectionSelfLocal();
          //std::cout << "Twist == " << twist << std::endl;
          QuadratureType quad(gridPart, it, quadOrd,  QuadratureType::INSIDE);
          const PointVectorType& points = 
            PointProviderType::getPoints(quad.id(), hexahedron);

          _test(points.size() == (size_t)(dim*2*quad.nop()));
          
          for (int i = 0; i < quad.nop(); ++i) {
            for (int d = 0; d < dim; ++d) {
              // geo has twist
              _floatTest(points[quad.cachingPoint(i)][d],
                         geo.global(quad.localPoint(i))[d]);
            }
            //std::cout << "pt " << i << ": " << points[quad.cachingPoint(i)]
            //          << " == " << geo.global(quad.localPoint(i)) << std::endl;
          }

          if( it.neighbor())
          {
            typedef TwistUtility<GridType> TwistUtilityType; 
            if( TwistUtilityType::conforming( grid, it ))
            {
              const GlobalGeometryType& gGeo = it.intersectionGlobal();
              QuadratureType outerQuad(gridPart, it, quadOrd , QuadratureType::OUTSIDE);
              
              for (int i = 0; i < quad.nop(); ++i) 
              {
                DomainType p = gGeo.global( quad.localPoint(i) );
                DomainType q = gGeo.global( outerQuad.localPoint(i) );

                for (int d = 0; d < dim; ++d) 
                {
                  _floatTest(p[d],q[d]);
                }
              }
            }
          }
        }
      } // end iterator loop
      grid.globalRefine(1);
    }
  }

  void CachingQuadrature_Test::codim1ALUTetraTest() 
  {
    const int dim = 3;
    const int codim = 1;
    const GeometryType tetrahedron( GeometryType::simplex, dim);

    typedef ALUSimplexGridFixture GridFixtureType;
    typedef GridFixtureType::GridType GridType;
    typedef LeafGridPart< GridType > GridPartType;
    typedef GridPartType::IntersectionIteratorType IntersectionIterator;
    typedef GridPartType :: Codim<0> :: IteratorType IteratorType;
    typedef CachingQuadrature<GridPartType, codim> QuadratureType;
    typedef PointProvider<GridType::ctype, dim, codim> PointProviderType;
    typedef PointProviderType::GlobalPointVectorType PointVectorType;
    typedef IntersectionIterator::LocalGeometry LocalGeometryType;
    typedef IntersectionIterator::Geometry GlobalGeometryType;
    typedef FieldVector<double,dim> DomainType;

    GridFixtureType fix(aluGridTetraFile_);
    GridType& grid = fix.grid();
    GridPartType gridPart( grid );

    const int quadOrd = 4;

    for(int l=0; l<3; ++l) 
    {
      IteratorType enditer = gridPart.end<0> ();
      for(IteratorType eiter = gridPart.begin<0> (); 
          eiter != enditer; ++eiter)
      {
        IntersectionIterator endit = gridPart.iend( *eiter );
        for (IntersectionIterator it = gridPart.ibegin(*eiter );
             it != endit; ++it) 
        {
          const LocalGeometryType& geo = it.intersectionSelfLocal();
          QuadratureType quad(gridPart, it, quadOrd , QuadratureType::INSIDE);

          const PointVectorType& points = 
            PointProviderType::getPoints(quad.id(), tetrahedron);

          _test(points.size() == (size_t)((dim +1) * quad.nop()));

          typedef TwistUtility<GridType> TwistUtilityType; 
          //std::cout << "New Intersection: Twists: ";
          //std::cout << TwistUtilityType :: twistInSelf( grid, it ) << " ";
          //std::cout << TwistUtilityType :: twistInNeighbor( grid, it ) << "\n";

          for (int i = 0; i < quad.nop(); ++i) 
          {
            for (int d = 0; d < dim; ++d) 
            {
              _floatTest(points[quad.cachingPoint(i)][d],
                         geo.global(quad.localPoint(i))[d]);
            }
            //std::cout << "nis: " << it.numberInSelf() << " nin: " << it.numberInNeighbor();
            //std::cout << " pt " << i << ": " << points[quad.cachingPoint(i)]
            //          << " == " << geo.global(quad.localPoint(i)) << std::endl;
          }

          if( it.neighbor())
          {
            if( TwistUtilityType::conforming( grid, it ))
            {
              const GlobalGeometryType& gGeo = it.intersectionGlobal();
              QuadratureType outerQuad(gridPart, it, quadOrd , QuadratureType::OUTSIDE);
              
              for (int i = 0; i < quad.nop(); ++i) 
              {
                DomainType p = gGeo.global( quad.localPoint(i) );
                DomainType q = gGeo.global( outerQuad.localPoint(i) );

                for (int d = 0; d < dim; ++d) 
                {
                  _floatTest(p[d],q[d]);
                }
              }
            }
          }

        } // end iterator loop
      }
      grid.globalRefine(1);
    }
  }

  void CachingQuadrature_Test::codim1SGridTest() 
  {
    const int dim = 2;
    const int codim = 1;
    GeometryType cube ( GeometryType::cube, dim);

    typedef SGridFixture<dim, dim> GridFixtureType;
    typedef GridFixtureType::GridType GridType;
    typedef LeafGridPart< GridType > GridPartType;
    typedef GridPartType::IntersectionIteratorType IntersectionIterator;
    typedef CachingQuadrature<GridPartType, codim> QuadratureType;
    typedef PointProvider<GridType::ctype, dim, codim> PointProviderType;
    typedef PointProviderType::GlobalPointVectorType PointVectorType;
    typedef IntersectionIterator::LocalGeometry LocalGeometryType;

    const int N = 2;

    GridFixtureType fix(N);
    GridType& grid = fix.grid();
    GridPartType gridPart( grid );

    //_fail("if you set the order to 4 here, the cache provider will find the wrong quadrature, since there is no distinction between line quadrature for triangles and quadrilaterals right now");


    IntersectionIterator endit = gridPart.iend( *grid.leafbegin<0> ());
    for (IntersectionIterator it = gridPart.ibegin(*grid.leafbegin<0>());
        it != endit; ++it) {
      const LocalGeometryType& geo = it.intersectionSelfLocal();
    
      QuadratureType quad(gridPart, it, 5, QuadratureType::INSIDE);
      const PointVectorType& points = 
        PointProviderType::getPoints(quad.id(), cube);

      _test(points.size() == (size_t)2*dim*quad.nop());
      
      for (int i = 0; i < quad.nop(); ++i) {
        for (int d = 0; d < dim; ++d) {
          _floatTest(points[quad.cachingPoint(i)][d],
                     geo.global(quad.localPoint(i))[d]);
        }
        //std::cout << "pt " << i << ": " << points[quad.cachingPoint(i)]
        //          << " == " << geo.global(quad.localPoint(i)) << std::endl;
      }
      
    } // end iterator loop
  }


} // end namespace Dune
 
