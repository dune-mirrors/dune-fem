#include <cassert>
#include <dune/grid/common/referenceelements.hh>

#include <dune/fem/space/basefunctions/storageinterface.hh>

namespace Dune {

  template <class ct, int dim>
  typename PointProvider<ct, dim, 0>::PointContainerType
  PointProvider<ct, dim, 0>::points_;

  template <class ct, int dim>
  void PointProvider<ct, dim, 0>::
  registerQuadrature(const QuadratureType& quad)
  {
    if (points_.find(quad.id()) == points_.end()) 
    {
      PointIteratorType it =
        points_.insert(std::make_pair
                       (quad.id(),
                        GlobalPointVectorType(quad.nop()))
                       ).first;
      for (size_t i = 0; i < quad.nop(); ++i) {
        it->second[i] = quad.point(i);
      }
      StorageInterface<dim>::registerQuadratureToStorages(quad);
    }
  }

  template <class ct, int dim>
  const typename PointProvider<ct, dim, 0>::GlobalPointVectorType&
  PointProvider<ct, dim, 0>::getPoints(size_t id, GeometryType elementGeo) 
  {
    assert(points_.find(id) != points_.end());
    return points_.find(id)->second;
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
                                        GeometryType elementGeo) 
  {
    MapperIteratorType it = mappers_.find(quad.id());
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
                                        GeometryType elementGeo)
  {
    MapperIteratorType it = mappers_.find(quad.id());
    if (it == mappers_.end()) {
      it = addEntry(quad, pts, elementGeo);
    }

    return it->second;
  }

  template <class ct, int dim>
  const typename PointProvider<ct, dim, 1>::GlobalPointVectorType&
  PointProvider<ct, dim, 1>::getPoints(size_t id, GeometryType elementGeo)
  {
    assert(points_.find(id) != points_.end());
    return points_.find(id)->second;
  }

  template <class ct, int dim>
  typename PointProvider<ct, dim, 1>::MapperIteratorType
  PointProvider<ct, dim, 1>::addEntry(const QuadratureType& quad,
                                      const LocalPointVectorType& points,
                                      GeometryType elementGeo)
  {
    const ReferenceElementContainer<ct, dim> refContainer;
    const ReferenceElement<ct, dim>& refElem = refContainer(elementGeo);

    const int numLocalPoints = points.size();
    const int numFaces = refElem.size(codim);
    const int numGlobalPoints = numFaces*numLocalPoints;
    
    PointIteratorType pit = 
      points_.insert(std::make_pair(quad.id(), 
                                    GlobalPointVectorType(numGlobalPoints))).first;
    MapperIteratorType mit =
      mappers_.insert(std::make_pair(quad.id(),
                                     MapperVectorType(numFaces))).first;

    int globalNum = 0;
    for (int face = 0; face < numFaces; ++face) {
      // Assumption: all faces have the same type
      // (not true for pyramids and prisms)
      assert(sameGeometry(quad.geometry(), refElem.type(face, codim)));
      MapperType pMap(numLocalPoints);
        
      for (int pt = 0; pt < numLocalPoints; ++pt, ++globalNum) {
        // Store point on reference element
        pit->second[globalNum] = 
          refElem.template global<codim>(points[pt], face, codim);
        
        // Store point mapping
        pMap[pt] = globalNum;
      }
      mit->second[face] = pMap;
    } // end for all faces

    StorageInterface<dim>::registerQuadratureToStorages(quad,1);
    return mit;
  }

  template <class ct, int dim>
  bool PointProvider<ct, dim, 1>::sameGeometry(GeometryType geo1, 
                                               GeometryType geo2)
  {
    return geo1 == geo2;
  }  
} // end namespace Dune
