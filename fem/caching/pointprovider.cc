#include <cassert>

namespace Dune {

  template <class ct, int dim>
  typename PointProvider<ct, dim, 0>::PointContainerType
  PointProvider<ct, dim, 0>::points_;

  template <class ct, int dim>
  void PointProvider<ct, dim, 0>::
  registerQuadrature(const QuadratureType& quad)
  {
    if (points_.find(quad.id()) != points_.end()) {
      PointIteratorType it =
        points_.insert(std::make_pair
                       (quad.id(),
                        GlobalPointVectorType(quad.nop()))
                       ).first;
      for (int i = 0; i < quad.nop(); ++i) {
        it->second[i] = quad.point(i);
      }
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
      for (int i = 0; i < quad.nop(); ++i) {
        pts[i] = quad.point(i);
      }
      it = addEntry(quad, pts, elementGeo);
    }
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
                                      const LocalPointVectorType& pts,
                                      GeometryType elementGeo)
  {
    const ReferenceElement<ct, dim>& refElem =
      ReferenceElements<ct, dim>::general(elementGeo);

    const int numLocalPoints = pts.size();
    const int numFaces = refElem.size(codim);
    const int numGlobalPoints = numFaces*numLocalPoints;
    
    PointIteratorType pit = 
      points_.insert(std::make_pair(quad.id(), 
                                    GlobalPointVectorType(numGlobalPoints)));
    MapperIteratorType mit =
      mappers_.insert(std::make_pair(quad.id(),
                                     MapperStorageType(numFaces)));


    int globalNum = 0;
    for (int face = 0; face < numFaces; ++face) {
      MapperType pMap(numLocalPoints);
        
      for (int pt = 0; pt < numLocalPoints; ++pts, ++globalNum) {
        // Store point on reference element
        pit->second[globalNum] = refElem.global(pts[pt], face, codim);
        
        // Store point mapping
        pMap[pt] = globalNum;
      }
      mit->second.addMapper(pMap, face);
    } // end for all faces

    return mit->second;
  }
    
} // end namespace Dune
