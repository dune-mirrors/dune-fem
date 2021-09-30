// C++ includes
#include <cassert>

// dune-geometry includes
#include <dune/geometry/referenceelements.hh>

// dune-fem includes
#include <dune/fem/quadrature/caching/registry.hh>
#include <dune/fem/misc/mpimanager.hh>

#include <dune/fem/space/localfiniteelement/quadratureinterpolation.hh>

namespace Dune
{

  namespace Fem
  {

    template <class ct, int dim>
    void PointProvider<ct, dim, 0>::
    registerQuadrature(const QuadratureType& quad)
    {
      QuadratureKeyType key( quad.geometryType(), quad.id() );

      PointContainerType& points_ = points();

      if (points_.find( key ) == points_.end() )
      {
        // only register when in single thread mode
        if( ! Fem :: MPIManager :: singleThreadMode() )
        {
          DUNE_THROW(SingleThreadModeError, "PointProvider::registerQuadrature: only call in single thread mode!");
        }

        PointIteratorType it =
          points_.insert(std::make_pair
                         (key,
                          GlobalPointVectorType(quad.nop()))
                         ).first;
        for (size_t i = 0; i < quad.nop(); ++i)
          it->second[i] = quad.point(i);

        // register quadrature to existing storages
        QuadratureStorageRegistry::registerQuadrature( quad );
      }
    }

    template <class ct, int dim>
    const typename PointProvider<ct, dim, 0>::GlobalPointVectorType&
    PointProvider<ct, dim, 0>::getPoints(const size_t id, const GeometryType& elementGeo)
    {
      QuadratureKeyType key( elementGeo, id );

      PointContainerType& points_ = points();

      typename PointContainerType::const_iterator pos = points_.find( key );
#ifndef NDEBUG
      if( pos == points_.end() )
      {
        std::cerr << "Unable to find quadrature points in list (key = " << key << ")." << std::endl;
        for( pos = points_.begin(); pos != points_.end(); ++pos )
          std::cerr << "found key: " << pos->first << std::endl;
        std::cerr << "Aborting..." << std::endl;
        abort();
      }
#endif
      return pos->second;
    }

    template <class ct, int dim>
    const typename PointProvider<ct, dim, 1>::MapperVectorPairType&
    PointProvider<ct, dim, 1>::getMappers(const QuadratureType& quad,
                                          const GeometryType& elementGeo)
    {
      QuadratureKeyType key( elementGeo , quad.id() );

      MapperContainerType& mappers_ = mappers();

      MapperIteratorType it = mappers_.find( key );
      if (it == mappers_.end())
      {
        std::vector<LocalPointType> pts(quad.nop());
        for (size_t i = 0; i < quad.nop(); ++i)
        {
          pts[i] = quad.point(i);
        }
        it = addEntry(quad, pts, elementGeo);
      }

      return it->second;
    }

    template <class ct, int dim>
    const typename PointProvider<ct, dim, 1>::MapperVectorPairType&
    PointProvider<ct, dim, 1>::getMappers(const QuadratureType& quad,
                                          const LocalPointVectorType& pts,
                                          const GeometryType& elementGeo)
    {
      QuadratureKeyType key( elementGeo, quad.id() );

      MapperContainerType& mappers_ = mappers();

      MapperIteratorType it = mappers_.find( key );
      if (it == mappers_.end()) {
        it = addEntry(quad, pts, elementGeo);
      }

      return it->second;
    }

    template <class ct, int dim>
    const typename PointProvider<ct, dim, 1>::GlobalPointVectorType&
    PointProvider<ct, dim, 1>::getPoints(const size_t id, const GeometryType& elementGeo)
    {
      QuadratureKeyType key( elementGeo, id );

      PointContainerType& points_ = points();
      assert(points_.find(key) != points_.end());
      return points_.find(key)->second;
    }

    template <class ct, int dim>
    typename PointProvider<ct, dim, 1>::MapperIteratorType
    PointProvider<ct, dim, 1>::addEntry(const QuadratureType& quad,
                                        const LocalPointVectorType& pts,
                                        GeometryType elementGeo)
    {
      // only addEntry when in single thread mode
      if( ! Fem :: MPIManager :: singleThreadMode() )
      {
        DUNE_THROW(SingleThreadModeError, "PointProvider::addEntry: only call in single thread mode!");
      }

      // generate key
      QuadratureKeyType key ( elementGeo, quad.id() );

      const auto &refElem = Dune::ReferenceElements<ct, dim>::general(elementGeo);

      const int numLocalPoints = pts.size();
      const int numFaces = refElem.size(codim);
      const int numGlobalPoints = numFaces*numLocalPoints;


      PointContainerType& points_ = points();
      PointIteratorType pit =
        points_.insert(std::make_pair(key,
                                      GlobalPointVectorType(numGlobalPoints))).first;
      MapperContainerType& mappers_ = mappers();
      MapperIteratorType mit =
        mappers_.insert(std::make_pair(key,
                                       std::make_pair(MapperVectorType(numFaces), MapperVectorType(numFaces) ))).first;

      MapperIteratorType iit;
      std::vector< GlobalPointType > itps = quad.interpolationPoints( dim );

      int globalNum = 0;
      const size_t nItp = itps.size();
      for (int face = 0; face < numFaces; ++face)
      {
        // Assumption: all faces have the same type
        // (not true for pyramids and prisms)
        MapperType pMap(numLocalPoints);
        MapperType itpMap(numLocalPoints);

        MapperVectorPairType& map = mit->second;

        for (int pt = 0; pt < numLocalPoints; ++pt, ++globalNum)
        {
          // Store point on reference element
          pit->second[globalNum] =
            refElem.template geometry<codim>(face).global( pts[pt] );

          if( nItp > 0 )
          {
            // compare interpolation points with created point to get mapping
            for( size_t i=0; i<nItp; ++i )
            {
              if( (itps[ i ] - pit->second[globalNum]).two_norm() < 1e-12 )
              {
                //std::cout << "found match for point " << itps[ i ] << " = ("<<  i << "," << pt << ")" << std::endl;
                itpMap[ pt ] = i;
              }
            }
          }

          // Store point mapping
          pMap[pt] = globalNum;
        }

        map.first[face] = pMap;   // = face*numLocalPoints+pt
        if( nItp > 0 )
        {
          map.second[face] = itpMap; // = face*numLocalPoints+pt
        }

      } // end for all faces

      // register quadrature to existing storages
      QuadratureStorageRegistry::registerQuadrature( quad, elementGeo, 1 );

      return mit;
    }

  } // namespace Fem

} // namespace Dune
