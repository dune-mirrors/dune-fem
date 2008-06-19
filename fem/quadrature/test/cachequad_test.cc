#include <config.h>

#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#ifdef ENABLE_UG
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#endif

#include "cachequad_test.hh"
#include "albertagrid_fixture.hh"
#include "alugrid_fixture.hh"
#include "sgrid_fixture.hh"

namespace Dune {

  void CachingQuadrature_Test::run() 
  {
    codim0Test();
    codim1AlbertaTest();
    codim1SGridTest();
    codim1UGTest();
    codim1ALUHexaTest();
    codim1ALUTetraTest();
    codim1YaspGridTest();
  }

  void CachingQuadrature_Test::codim0Test() 
  {
#if ENABLE_ALBERTA
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
#endif
  } 

  template <class GridPartType> 
  void CachingQuadrature_Test::checkLeafsCodim1(GridPartType& gridPart,
                  const int quadOrd)
  {
    typedef typename GridPartType :: GridType GridType; 
    enum { dim = GridType :: dimension };
    enum { codim = 1 };
    typedef typename GridType::ctype ctype;

    typedef typename GridPartType::IntersectionIteratorType IntersectionIterator;
    typedef typename IntersectionIterator::Intersection Intersection;
    typedef typename GridPartType :: template Codim<0> :: IteratorType IteratorType;
    typedef CachingQuadrature<GridPartType, codim> QuadratureType;
    typedef PointProvider<ctype, dim, codim> PointProviderType;
    typedef typename PointProviderType::GlobalPointVectorType PointVectorType;
    typedef typename Intersection::LocalGeometry LocalGeometryType;
    typedef typename Intersection::Geometry GlobalGeometryType;
    typedef FieldVector<ctype,dim> DomainType;

    IteratorType enditer = gridPart.template end<0> ();
    for(IteratorType eiter = gridPart.template begin<0> (); 
        eiter != enditer; ++eiter)
    {
      const GeometryType geomType = eiter->geometry().type();
      const ReferenceElement< ctype, dim > & refElem =
                    ReferenceElements< ctype, dim >::general(geomType);
      const int numFaces = refElem.size(codim);
      //std::cout << "For type " << geomType << " got " << numFaces << " numFaces\n";

      //int twist = -4;
      const IntersectionIterator endit = gridPart.iend( *eiter );
      for (IntersectionIterator it = gridPart.ibegin( *eiter );
           it != endit; ++it) 
      {
        const Intersection& inter=*it;
        if( dim > 2 )
        {
          checkLocalIntersectionConsistency( *inter.inside(),
              inter.intersectionSelfLocal(), inter.numberInSelf() );
          if( inter.neighbor() ) 
          {
            checkLocalIntersectionConsistency( *inter.outside(),
                inter.intersectionNeighborLocal(), inter.numberInNeighbor() );
          }
        }

        const LocalGeometryType& geo = inter.intersectionSelfLocal();
        typedef TwistUtility<GridType> TwistUtilityType; 

        QuadratureType quad(gridPart, inter, quadOrd , QuadratureType :: INSIDE);

        const PointVectorType& points = 
          PointProviderType::getPoints(quad.id(), geomType);

        _test((int) points.size() == numFaces * quad.nop());
        //std::cout << points.size() << " ps | qnop " << numFaces * quad.nop() << "\n";

        //std::cout << "New Intersection: Twists: ";
        //std::cout << TwistUtilityType :: twistInSelf( grid, it ) << " ";
        //std::cout << TwistUtilityType :: twistInNeighbor( grid, it ) << "\n";

        for (int i = 0; i < quad.nop(); ++i) 
        {
          for (int d = 0; d < dim; ++d) 
          {
            //std::cout << quad.cachingPoint(i) << " cp | size " << points.size() << "\n";
            assert( quad.cachingPoint(i) < points.size() );
            _floatTest(points[quad.cachingPoint(i)][d],
                       geo.global(quad.localPoint(i))[d]);
          }
          //std::cout << "nis: " << it.numberInSelf();
          //std::cout << " pt " << i << ": " << points[quad.cachingPoint(i)]
          //          << " == " << geo.global(quad.localPoint(i)) << std::endl;
        }

        if( inter.neighbor ())
        {
          if( TwistUtilityType::conforming( gridPart.grid(), inter ))
          {
            const LocalGeometryType& nGeo = inter.intersectionNeighborLocal();
            QuadratureType outerQuad(gridPart, inter, quadOrd , QuadratureType::OUTSIDE);
            
            for (int i = 0; i < outerQuad.nop(); ++i) 
            {
              for (int d = 0; d < dim; ++d) 
              {
                assert( outerQuad.cachingPoint(i) < points.size() );
                _floatTest(points[outerQuad.cachingPoint(i)][d],
                           nGeo.global(outerQuad.localPoint(i))[d]);
              }
              //std::cout << "nin: " << it.numberInNeighbor();
              //std::cout << " nis: " << it.numberInSelf();
              //std::cout << " pt " << i << ": " << points[outerQuad.cachingPoint(i)]
              //          << " == " << nGeo.global(outerQuad.localPoint(i)) << std::endl;
            }
          }
        }
      } // end iterator loop
    }
  }

  void CachingQuadrature_Test::codim1AlbertaTest() 
  {
#ifdef ENABLE_ALBERTA
    const int dim = GRIDDIM ;

    typedef AlbertaGridFixture<dim, dim> GridFixtureType;
    typedef GridFixtureType::GridType GridType;
    typedef LeafGridPart< GridType > GridPartType;

    GridFixtureType fix(albertaGridFile_);
    GridType& grid = fix.grid();
    GridPartType gridPart( grid );

    const int quadOrd = 4;

    for(int l=0; l<5; ++l) 
    {
      checkLeafsCodim1(gridPart,quadOrd);
      grid.globalRefine(1);
    }
#endif
  }

  void CachingQuadrature_Test::codim1UGTest() 
  {
#ifdef ENABLE_UG 
    std::cout << "CachingQuadrature_Test checking UGGrid \n";
    // 2d test 
    {
      const int dim = 2;

      typedef UGGrid<dim> GridType;
      typedef LeafGridPart< GridType > GridPartType;

      GridPtr<GridType> gridPtr(dgf2DGridFile_);
      GridType& grid = *gridPtr;
      GridPartType gridPart( grid );

      const int quadOrd = 4;

      for(int l=0; l<3; ++l) 
      {
        checkLeafsCodim1(gridPart,quadOrd);
        grid.globalRefine(1);
      }
    }
    /*
    // 3d test 
    {
      const int dim = 3;

      typedef UGGrid<dim> GridType;
      typedef LeafGridPart< GridType > GridPartType;

      GridPtr<GridType> gridPtr(dgf3DGridFile_);
      GridType& grid = *gridPtr;
      GridPartType gridPart( grid );

      const int quadOrd = 4;

     // for(int l=0; l<3; ++l) 
      {
        checkLeafsCodim1(gridPart, quadOrd);
        grid.globalRefine(1);
      }
    }
    */
#endif
  }

  void CachingQuadrature_Test::codim1ALUHexaTest() 
  {
#ifdef ENABLE_ALUGRID 
    typedef ALUCubeGridFixture GridFixtureType;
    typedef GridFixtureType::GridType GridType;
    typedef LeafGridPart< GridType > GridPartType;

    GridFixtureType fix(aluGridHexaFile_);
    GridType& grid = fix.grid();
    GridPartType gridPart( grid );

    const int quadOrd = 4;

    for(int l=0; l<3; ++l) 
    {
      checkLeafsCodim1(gridPart,quadOrd);
      grid.globalRefine(1);
    }
#endif
  }

  void CachingQuadrature_Test::codim1ALUTetraTest() 
  {
#ifdef ENABLE_ALUGRID 
    typedef ALUSimplexGridFixture GridFixtureType;
    typedef GridFixtureType::GridType GridType;
    typedef LeafGridPart< GridType > GridPartType;

    GridFixtureType fix(aluGridTetraFile_);
    GridType& grid = fix.grid();
    GridPartType gridPart( grid );

    const int quadOrd = 4;

    for(int l=0; l<3; ++l) 
    {
      checkLeafsCodim1(gridPart,quadOrd);
      grid.globalRefine(1);
    }
#endif
  }

  void CachingQuadrature_Test::codim1SGridTest() 
  {
    std::cout << "CachingQuadrature_Test checking SGrid \n";
    const int dim = 2;

    typedef SGridFixture<dim, dim> GridFixtureType;
    typedef GridFixtureType::GridType GridType;
    typedef LeafGridPart< GridType > GridPartType;

    const int N = 2;

    GridFixtureType fix(N);
    GridType& grid = fix.grid();
    GridPartType gridPart( grid );

    const int quadOrd = 5;

    for(int l=0; l<2; ++l) 
    {
      checkLeafsCodim1(gridPart,quadOrd);
      grid.globalRefine( 1 );
    }
  }

  void CachingQuadrature_Test::codim1YaspGridTest() 
  {
    std::cout << "CachingQuadrature_Test checking YaspGrid \n";
    const int dim = 3;

    typedef YaspGrid<dim, dim> GridType;
    typedef LeafGridPart< GridType > GridPartType;

    GridPtr<GridType> gridPtr ( dgf3DGridFile_ );
    GridType& grid = *gridPtr;
    GridPartType gridPart( grid );

    const int quadOrd = 5;

    for(int l=0; l<2; ++l) 
    {
      checkLeafsCodim1(gridPart,quadOrd);
      grid.globalRefine( 1 );
    }
  }

  template <class EntityType, class LocalGeometryType>
  void CachingQuadrature_Test::checkLocalIntersectionConsistency(
      const EntityType& en, const LocalGeometryType& localGeom, 
      const int face, const bool output) const
  {
    enum { dim = EntityType :: dimension };
    typedef typename EntityType :: ctype ctype;

    // get reference element 
    const ReferenceElement< ctype , dim > & refElem = 
      ReferenceElements< ctype , dim >::general(en.geometry().type());

    const int vxSize = refElem.size( face, 1, dim );
    std::vector<int> vx( vxSize ,-1);
    for(int i=0; i<vxSize; ++i) 
    {
      // get face vertices of number in self face 
      vx[i] = refElem.subEntity( face, 1 , i, dim);
    }

    // debugging output 
    if( output )
    {
      std::cout << "Found face["<< face << "] vx = {";
      for(size_t i=0; i<vx.size(); ++i) 
      {
        std::cout << vx[i] << ",";
      }
      std::cout << "} \n";
    }

    bool allRight = true;
    std::vector< int > faceMap ( vxSize , -1 );

    typedef  FieldVector<ctype,dim> CoordinateVectorType;

    for(int i=0; i<vxSize; ++i) 
    {
      // standard face map is identity 
      faceMap[i] = i;

      // get position in reference element of vertex i
      CoordinateVectorType refPos = refElem.position( vx[i], dim );

      // get position as we get it from intersectionSelfLocal 
      // in the best case this should be the same 
      // at least the orientatation should be the same 
      CoordinateVectorType localPos = localGeom[i];

      if( (refPos - localPos).infinity_norm() > 1e-8 )
      {
        allRight = false;
        if( output )
          std::cout << "RefPos (" << refPos << ") != (" << localPos << ") localPos !\n";
      }
    }

    if( !allRight ) 
    {
      for(int i=0; i<vxSize; ++i) 
      {
        // get position in reference element of vertex i
        CoordinateVectorType refPos = refElem.position( vx[i], dim );

        for( int j=1; j<vxSize; ++j) 
        {
          int newVx = (i+j)% vxSize;
          CoordinateVectorType localPos = localGeom[newVx];
          if( (refPos - localPos).infinity_norm() < 1e-8 )
          {
            faceMap[i] = newVx;
          }
        }
      }

      // calculate twist 
      const int twist = (faceMap[1] == (faceMap[0]+1)%vxSize) ? faceMap[0] : faceMap[1] - vxSize;
      std::cout << "Got twist = "<< twist << "\n";
      // off set 
      const int offset = (2 * vxSize) + 1;

      // now check mapping with twist 
      for(int i=0; i<vxSize; ++i) 
      {
        // get position in reference element of vertex i
        CoordinateVectorType refPos = refElem.position( vx[i], dim );

        int newVx = (twist < 0) ? (offset - i + twist)% vxSize : (twist + i)%vxSize ;
        CoordinateVectorType localPos = localGeom[newVx];
        if( (refPos - localPos).infinity_norm() > 1e-8 )
        {
          std::cout << "RefPos (" << refPos << ") != (" << localPos << ") localPos !\n";
          DUNE_THROW(GridError,"LocalGeometry has wrong mapping !");
        }
      }
    }
  }




} // end namespace Dune
 
