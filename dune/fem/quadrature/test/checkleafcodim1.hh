#ifndef DUNE_CHECKLEAFCODIM1_HH
#define DUNE_CHECKLEAFCODIM1_HH

#include <dune/fem/io/parameter.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

namespace Dune {
  namespace Fem {

struct CachingQuadratureTest
{
  template <class GridPartType>
  static void checkLeafsCodimOne(GridPartType& gridPart,
                                 const int quadOrd)
  {
    CachingQuadratureTest tester;
    tester.checkLeafsCodim1( gridPart, quadOrd );
  }

  static void doTest( const double& a, const double& b )
  {
    if( std::abs( a - b ) > 1e-12 )
    {
      assert( false );
      std::abort();
    }
  }

  static void doTest( const bool value )
  {
    if( ! value )
    {
      assert( false );
      std::abort();
    }
  }

  void run (){};

protected:
  template <class EntityType, class LocalGeometryType>
  int aluTwistCheck(const EntityType& en, const LocalGeometryType& localGeom,
                    const int face, const bool neighbor, const bool output ) const
  {
    enum { dim = EntityType :: dimension };
    typedef typename EntityType :: Geometry :: ctype ctype;

    typedef FaceTopologyMapping<tetra> SimplexFaceMapping;
    typedef FaceTopologyMapping<hexa>  CubeFaceMapping;

    // get reference element
    const auto refElem = Dune::referenceElement< ctype, dim >( en.type() );

    const int vxSize = refElem.size( face, 1, dim );
    typedef  FieldVector<ctype,dim> CoordinateVectorType;

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
        int twistedDuneIndex = -1;
        if( localGeom.type().isCube() )
        {
          twistedDuneIndex = CubeFaceMapping::twistedDuneIndex( i, twist );
        }
        else
        {
          twistedDuneIndex = SimplexFaceMapping::twistedDuneIndex( i, twist );
        }

        // get face vertices of number in self face
        int vxIdx = refElem.subEntity( face, 1 , twistedDuneIndex , dim);

        // get position in reference element of vertex i
        CoordinateVectorType refPos = refElem.position( vxIdx, dim );

        // check coordinates again
        CoordinateVectorType localPos = localGeom.corner( i );
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
      assert(false);
      DUNE_THROW(GridError,"Not matching twist found");
    }

    if( output )
    {
      std::string twistIn( (neighbor) ? "twistInNeighbor()" : "twistInSelf" );
      std::string numberIn( (neighbor) ? "indexInOutside()" : "indexInInside" );
      std::cout << "ERROR: Face "<< face << " : twist = "<< twistFound << std::endl;
      std::cout << "\nPut twist = "<< twistFound << " In TwistUtility::"<< twistIn << " for " << numberIn << " = " << face << " ! \n";
      std::cout << "******************************************\n";
    }

    return twistFound;
  }

  template <class GridPartType>
  void checkLeafsCodim1(GridPartType& gridPart,
                        const int quadOrd)
  {
    enum { dim = GridPartType :: dimension };
    enum { codim = 1 };
    typedef typename GridPartType ::ctype ctype;

    typedef typename GridPartType::IntersectionIteratorType IntersectionIterator;
    typedef typename IntersectionIterator::Intersection Intersection;
    typedef typename GridPartType :: template Codim<0> :: IteratorType IteratorType;
    typedef CachingQuadrature<GridPartType, codim> QuadratureType;
    typedef PointProvider<ctype, dim, codim> PointProviderType;
    typedef typename PointProviderType::GlobalPointVectorType PointVectorType;
    typedef typename Intersection::LocalGeometry LocalGeometryType;
    typedef typename  GridPartType :: template Codim<0> :: EntityType EntityType;

    IteratorType enditer = gridPart.template end<0> ();
    for(IteratorType eiter = gridPart.template begin<0> ();
        eiter != enditer; ++eiter)
    {
      const EntityType& entity = (*eiter);
      const GeometryType geomType = entity.geometry().type();
      const auto refElem = Dune::referenceElement< ctype, dim >( geomType );
      const int numFaces = refElem.size(codim);
      //std::cout << "For type " << geomType << " got " << numFaces << " numFaces\n";

      //int twist = -4;
      const IntersectionIterator endit = gridPart.iend( entity );
      for (IntersectionIterator it = gridPart.ibegin( entity );
           it != endit; ++it)
      {
        const Intersection& inter=*it;
        typedef typename GridPartType::TwistUtilityType TwistUtilityType;

        // set this flag to true for output of twists that have been calculated
        const bool output = Parameter :: verbose() ;

        if( dim > 2 )
        {
          //const int twistFound = checkLocalIntersectionConsistency( entity,
          //                      inter.geometryInInside(), inter.indexInInside() , false, false);
          const int twistFound = aluTwistCheck( entity,
                                inter.geometryInInside(), inter.indexInInside() , false, false);
          const int twistInside = TwistUtilityType::twistInSelf( gridPart.grid(), inter);
          if( output && twistFound != twistInside )
          {
            std::cout << "Twist inconsistent: calculated twist " << twistFound << "  not equal to inside " << twistInside << "\n";
          }

          if( inter.neighbor() )
          {
            EntityType neighbor = inter.outside();
            //const int twstF = checkLocalIntersectionConsistency( neighbor,
            //              inter.geometryInOutside(), inter.indexInOutside(), true, false);
            const int twstF = aluTwistCheck( neighbor,
                          inter.geometryInOutside(), inter.indexInOutside(), true, false);

            const int twistOutside = TwistUtilityType::twistInNeighbor( gridPart.grid(), inter);
            if( output && twstF != twistOutside )
            {
              std::cout << "Twist inconsistent: calculated twist " << twstF << "  not equal to outside " << twistOutside << "\n";
            }
          }
        }

        const LocalGeometryType& geo = inter.geometryInInside();

        QuadratureType quad(gridPart, inter, quadOrd , QuadratureType :: INSIDE);

        const PointVectorType& points =
          PointProviderType::getPoints(quad.id(), geomType);

        doTest( points.size() == numFaces * quad.nop());
        if ( output )
        {
          std::cout << points.size() << " ps | qnop " << numFaces * quad.nop() << "\n";
          std::cout << "New Intersection: Twists: ";
          std::cout << TwistUtilityType :: twistInSelf( gridPart.grid(), inter);
          if( inter.neighbor() )
            std::cout << " " << TwistUtilityType :: twistInNeighbor( gridPart.grid(), inter);
          std::cout << std::endl;
        }

        // this check is only ok if we are on the inside of
        // the smaller element of the non-conforming intersection
        // or if the intersection is conforming
        if( inter.conforming() ||
            (inter.neighbor() && entity.level() > inter.outside().level()) )
        {
          for (size_t i = 0; i < quad.nop(); ++i)
          {
            /*
            {
              std::cout << "nis: " << inter.indexInInside();
              std::cout << " pt " << i << ": " << points[quad.cachingPoint(i)]
                        << " == " << geo.global(quad.localPoint(i)) << std::endl;
            }
            */
            for (int d = 0; d < dim; ++d)
            {
              assert( quad.cachingPoint(i) < points.size() );
              doTest(points[quad.cachingPoint(i)][d],
                     geo.global(quad.localPoint(i))[d]);
            }

          }
        }

        if( inter.neighbor ())
        {
          if( inter.conforming() )
          {
            const LocalGeometryType& nGeo = inter.geometryInOutside();
            QuadratureType outerQuad(gridPart, inter, quadOrd , QuadratureType::OUTSIDE);

            for (size_t i = 0; i < outerQuad.nop(); ++i)
            {
              for (int d = 0; d < dim; ++d)
              {
                assert( outerQuad.cachingPoint(i) < points.size() );
                doTest(points[outerQuad.cachingPoint(i)][d],
                           nGeo.global(outerQuad.localPoint(i))[d]);
              }

              /*
              {
                std::cout << "nin: " << inter.indexInOutside();
                std::cout << " nis: " << inter.indexInInside();
                std::cout << " pt " << i << ": " << points[outerQuad.cachingPoint(i)]
                          << " == " << nGeo.global(outerQuad.localPoint(i)) << std::endl;
              }
              */
            }
          }
        }

      } // end iterator loop
    }
  }
#if 0
  template <class EntityType, class LocalGeometryType>
  int checkLocalIntersectionConsistency(
      const EntityType& en, const LocalGeometryType& localGeom,
      const int face, const bool neighbor, const bool output ) const
  {
    enum { dim = EntityType :: dimension };
    typedef typename EntityType :: ctype ctype;

    typedef FaceTopologyMapping<tetra> SimplexFaceMapping;
    typedef FaceTopologyMapping<hexa>  CubeFaceMapping;

    // get reference element
    const  auto refElem = Dune::referenceElement< ctype, dim >( en.type() );

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
      CoordinateVectorType localPos = localGeom.corner( idx );

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
          int twistedDuneIndex = -1;
          if( localGeom.type().isCube() )
          {
            twistedDuneIndex = CubeFaceMapping::twistedDuneIndex( i, twist );
          }
          else
          {
            twistedDuneIndex = SimplexFaceMapping::twistedDuneIndex( i, twist );
          }

          // get face vertices of number in self face
          int vxIdx = refElem.subEntity( face, 1 , twistedDuneIndex , dim);

          // get position in reference element of vertex i
          CoordinateVectorType refPos = refElem.position( vxIdx, dim );

          // check coordinates again
          CoordinateVectorType localPos = localGeom.corner( i );
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
        assert(false);
        DUNE_THROW(GridError,"Not matching twist found");
      }

      if( output )
      {
        std::string twistIn( (neighbor) ? "twistInNeighbor()" : "twistInSelf" );
        std::string numberIn( (neighbor) ? "indexInOutside()" : "indexInInside" );
        std::cout << "Face "<< face << " : twist = "<< twistFound << std::endl;
        std::cout << "\nPut twist = "<< twistFound << " In TwistUtility::"<< twistIn << " for " << numberIn << " = " << face << " ! \n";
        std::cout << "******************************************\n";
      }

      return twistFound;
    }

    assert( ! faceTwisted );
    return 0;
  }
#endif
};

  } // end namespace Fem
} // namespace Dune

#endif
