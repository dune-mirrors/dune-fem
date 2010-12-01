#include <cassert>
#include <dune/grid/common/genericreferenceelements.hh>

#include <dune/fem/space/basefunctions/storageinterface.hh>
#include <dune/fem/misc/ompmanager.hh>

namespace Dune {

  template <class ct, int dim>
  typename PointProvider<ct, dim, 0>::PointContainerType
  PointProvider<ct, dim, 0>::points_;

  template <class ct, int dim>
  void PointProvider<ct, dim, 0>::
  registerQuadrature(const QuadratureType& quad)
  {
    // only register on thread 0 
    if( Fem :: OMPManager :: thread() == 0 ) 
    {
      QuadratureKeyType key( quad.geometry(), quad.id() );
      
      if (points_.find( key ) == points_.end() ) 
      {
        PointIteratorType it =
          points_.insert(std::make_pair
                         (key,
                          GlobalPointVectorType(quad.nop()))
                         ).first;
        for (size_t i = 0; i < quad.nop(); ++i) {
          it->second[i] = quad.point(i);
        }
        Fem::StorageInterface<dim>::registerQuadratureToStorages(quad);
      }
    }

    Fem :: OMPManager :: barrier();
  }

  template <class ct, int dim>
  const typename PointProvider<ct, dim, 0>::GlobalPointVectorType&
  PointProvider<ct, dim, 0>::getPoints(const size_t id, const GeometryType& elementGeo) 
  {
    QuadratureKeyType key( elementGeo, id );
    
    assert(points_.find(key) != points_.end());
    return points_.find(key)->second;
  }

  template <class ct, int dim>
  typename PointProvider<ct, dim, 1>::PointContainerType
  PointProvider<ct, dim, 1>::points_;

  template <class ct, int dim>
  typename PointProvider<ct, dim, 1>::MapperContainerType
  PointProvider<ct, dim, 1>::mappers_;

  template <class ct, int dim>
  const typename PointProvider<ct, dim, 1>::MapperVectorType&
  PointProvider<ct, dim, 1>::getMappers(const QuadratureType& quad,
                                        const GeometryType& elementGeo) 
  {
    QuadratureKeyType key( elementGeo , quad.id() );
    
    MapperIteratorType it = mappers_.find( key );
    if (it == mappers_.end()) {
      std::vector<LocalPointType> pts(quad.nop());
      for (size_t i = 0; i < quad.nop(); ++i) {
        pts[i] = quad.point(i);
      }
      it = addEntry(quad, pts, elementGeo);
    }

    return it->second;
  }

  template <class ct, int dim>
  const typename PointProvider<ct, dim, 1>::MapperVectorType&
  PointProvider<ct, dim, 1>::getMappers(const QuadratureType& quad,
                                        const LocalPointVectorType& pts,
                                        const GeometryType& elementGeo)
  {
    QuadratureKeyType key( elementGeo, quad.id() );
    
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
    
    assert(points_.find(key) != points_.end());
    return points_.find(key)->second;
  }

  template <class ct, int dim>
  typename PointProvider<ct, dim, 1>::MapperIteratorType
  PointProvider<ct, dim, 1>::addEntry(const QuadratureType& quad,
                                      const LocalPointVectorType& points,
                                      GeometryType elementGeo)
  {
    // generate key 
    QuadratureKeyType key ( elementGeo, quad.id() );
      
    MapperIteratorType mit ;
    if( Fem::OMPManager::thread() == 0 )
    {
      const GenericReferenceElement<ct, dim>& refElem = 
        GenericReferenceElements<ct, dim>::general(elementGeo);

      const int numLocalPoints = points.size();
      const int numFaces = refElem.size(codim);
      const int numGlobalPoints = numFaces*numLocalPoints;
      
      PointIteratorType pit = 
        points_.insert(std::make_pair(key, 
                                      GlobalPointVectorType(numGlobalPoints))).first;
      mit =
        mappers_.insert(std::make_pair(key,
                                       MapperVectorType(numFaces))).first;
      int globalNum = 0;
      for (int face = 0; face < numFaces; ++face) 
      {
        // Assumption: all faces have the same type
        // (not true for pyramids and prisms)
        //assert(sameGeometry(quad.geometry(), refElem.type(face, codim)));
        MapperType pMap(numLocalPoints);
          
        for (int pt = 0; pt < numLocalPoints; ++pt, ++globalNum) {
          // Store point on reference element
          pit->second[globalNum] = 
            refElem.template global<codim>(points[pt], face, codim);
          
          // Store point mapping
          pMap[pt] = globalNum;
        }
        mit->second[face] = pMap;  // = face*numLocalPoints+pt
      } // end for all faces

      Fem::StorageInterface<dim>::registerQuadratureToStorages(quad,1);
    }
    
    Fem::OMPManager::barrier();

    // assign iterator for higher threads 
    if( Fem::OMPManager::thread() > 0 )
      mit = mappers_.find( key );

    return mit;
  }

  template <class ct, int dim>
  bool PointProvider<ct, dim, 1>::sameGeometry(GeometryType geo1, 
                                               GeometryType geo2)
  {
    return geo1 == geo2;
  }  
} // end namespace Dune
