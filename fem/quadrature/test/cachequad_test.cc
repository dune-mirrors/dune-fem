#include <config.h>

#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#ifdef ENABLE_UG
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#endif

#ifdef ENABLE_PRISMGRID
#include <dune/prismgrid/grid.hh>
#endif

#include <dune/grid/alugrid/3d/topology.hh>

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
    codim1PrismTest();
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
              inter.intersectionSelfLocal(), inter.numberInSelf() , false );
          if( inter.neighbor() ) 
          {
            checkLocalIntersectionConsistency( *inter.outside(),
                inter.intersectionNeighborLocal(), inter.numberInNeighbor(), true );
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

  void CachingQuadrature_Test::codim1PrismTest()
  {
#ifdef ENABLE_PRISMGRID 
    std::cout << "\n**********************************************\n";
    std::cout << "CachingQuadrature_Test checking PrismGrid \n";
    {
      const int dim = 3;

      typedef PrismGrid<dim,dim> GridType;
      typedef LeafGridPart< GridType > GridPartType;

      GridType grid(dgf2DGridFile_, 0, 3);
      GridPartType gridPart( grid );

      const int quadOrd = 4;

      // check grid 
      checkLeafsCodim1(gridPart,quadOrd);
    }
#endif
  }

  void CachingQuadrature_Test::codim1UGTest() 
  {
#ifdef ENABLE_UG 
    std::cout << "\n**********************************************\n";
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

    // 3d test 
    {
      const int dim = 3;

      typedef UGGrid<dim> GridType;
      typedef LeafGridPart< GridType > GridPartType;

      GridPtr<GridType> gridPtr(dgf3DGridFile_);
      GridType& grid = *gridPtr;
      GridPartType gridPart( grid );

      const int quadOrd = 4;

      for(int l=0; l<3; ++l) 
      {
        checkLeafsCodim1(gridPart, quadOrd);
        grid.globalRefine(1);
      }
    }
#endif
  }

  void CachingQuadrature_Test::codim1ALUHexaTest() 
  {
#ifdef ENABLE_ALUGRID 
    std::cout << "\n**********************************************\n";
    std::cout << "CachingQuadrature_Test checking ALUCubeGrid\n";

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
    std::cout << "\n**********************************************\n";
    std::cout << "CachingQuadrature_Test checking ALUSimplexGrid\n";

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
    std::cout << "\n**********************************************\n";
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
      const int face, const bool neighbor, const bool output ) const
  {
    enum { dim = EntityType :: dimension };
    typedef typename EntityType :: ctype ctype;

    typedef FaceTopologyMapping<tetra> SimplexFaceMapping;
    typedef FaceTopologyMapping<hexa>  CubeFaceMapping;

    // get reference element 
    const ReferenceElement< ctype , dim > & refElem = 
      ReferenceElements< ctype , dim >::general(en.geometry().type());

    const int vxSize = refElem.size( face, 1, dim );
    std::vector<int> vx( vxSize ,-1);
    for(int i=0; i<vxSize; ++i) 
    {
      //const int idx = i;
      const int idx = ( localGeom.type().isCube() ) ?
            CubeFaceMapping::dune2aluVertex( i ) :
            SimplexFaceMapping::dune2aluVertex( i );

      // get face vertices of number in self face 
      vx[i] = refElem.subEntity( face, 1 , idx , dim);
    }

    // debugging output 
    if( output )
    {
      std::string neighout ((neighbor)?"outside":"inside");
      std::cout << "\n******************************************\n";
      std::cout << "Found ("<<neighout<<") face["<< face << "] vx = {";
      for(size_t i=0; i<vx.size(); ++i) 
      {
        std::cout << vx[i] << ",";
      }
      std::cout << "} \n";
    }

    bool faceTwisted = false;
    std::vector< int > faceMap ( vxSize , -1 );

    typedef  FieldVector<ctype,dim> CoordinateVectorType;

    for(int i=0; i<vxSize; ++i) 
    {
      //const int idx = i;
      const int idx = ( localGeom.type().isCube() ) ?
            CubeFaceMapping::dune2aluVertex( i ) :
            SimplexFaceMapping::dune2aluVertex( i );

      // standard face map is identity 
      faceMap[i] = idx;

      // get position in reference element of vertex i
      CoordinateVectorType refPos = refElem.position( vx[i], dim );

      // get position as we get it from intersectionSelfLocal 
      // in the best case this should be the same 
      // at least the orientatation should be the same 
      CoordinateVectorType localPos = localGeom[idx];

      if( (refPos - localPos).infinity_norm() > 1e-8 )
      {
        faceTwisted = true;
        if( output )
          std::cout << "RefPos (" << refPos << ") != (" << localPos << ") localPos !\n";
      }
    }

    if( faceTwisted ) 
    {
      if( output ) 
      {
        std::string neighout ((neighbor)?"outside":"inside");
        std::cout <<"Face "<< face << " ("<<neighout<< ") is twisted! \n";
      }
      
      // generate new twisted mapping 
      for(int i=0; i<vxSize; ++i) 
      {
        // get position in reference element of vertex i
        CoordinateVectorType refPos = refElem.position( vx[i], dim );

        for( int j=1; j<vxSize; ++j) 
        {
          int newVx = (i+j)% vxSize;
          newVx = ( localGeom.type().isCube() ) ?
                         CubeFaceMapping::dune2aluVertex( newVx ) :
                         SimplexFaceMapping::dune2aluVertex( newVx );
          
          CoordinateVectorType localPos = localGeom[newVx];
          if( (refPos - localPos).infinity_norm() < 1e-8 )
          {
            faceMap[i] = newVx;
          }
        }
      }

      // consistency check, face mapping should map all points 
      for(int i=0; i<vxSize; ++i) 
      {
        // get position in reference element of vertex i
        CoordinateVectorType refPos   = refElem.position( vx[i], dim );
        // get corner from local geometry 
        CoordinateVectorType localPos = localGeom[ faceMap[i] ];

        if( (refPos - localPos).infinity_norm() > 1e-8 )
        {
          DUNE_THROW(GridError,"Inconsistent face mapping");
        }
      }

      // now calculate twist by trial and error for all possible twists 
      // the calculated twist is with respect to the ALUGrid 
      // reference face, see twistprovider.cc  
      int twistFound = -66;
      for(int twist = -vxSize; twist<vxSize; ++twist) 
      {
        bool twistOk = true; 
        // now check mapping with twist 
        for(int i=0; i<vxSize; ++i) 
        {
          // get face vertices of number in self face 
          int vxIdx = refElem.subEntity( face, 1 , i , dim);
          
          // get position in reference element of vertex i
          CoordinateVectorType refPos = refElem.position( vxIdx, dim );

          int twistedDuneIndex = -1;
          if( localGeom.type().isCube() ) 
          {
            const int aluIndex = CubeFaceMapping::dune2aluVertex( i );
            twistedDuneIndex = CubeFaceMapping::alu2duneVertex(aluIndex, twist);
          }
          else 
          {
            const int aluIndex = SimplexFaceMapping::dune2aluVertex( i );
            twistedDuneIndex = SimplexFaceMapping::alu2duneVertex(aluIndex, twist);
          }
          
          // check coordinates again 
          CoordinateVectorType localPos = localGeom[ twistedDuneIndex ];
          if( (refPos - localPos).infinity_norm() > 1e-8 )
          {
            twistOk = false;
            break; 
          }
        }

        if( twistOk ) 
        {
          twistFound = twist;
          break ;
        }
      }

      // if no twist found, then something is wrong 
      if( twistFound == -66 ) 
      {
        DUNE_THROW(GridError,"Not matching twist found");
      }
      
      if( output )
      {
        std::string twistIn( (neighbor) ? "twistInNeighbor()" : "twistInSelf" );
        std::string numberIn( (neighbor) ? "numberInNeighbor()" : "numberInSelf" );
        std::cout << "Face "<< face << " : twist = "<< twistFound << std::endl;
        std::cout << "\nPut twist = "<< twistFound << " In TwistUtility::"<< twistIn << " for " << numberIn << " = " << face << " ! \n";
        std::cout << "******************************************\n";
      }
    }
  }




} // end namespace Dune
 
