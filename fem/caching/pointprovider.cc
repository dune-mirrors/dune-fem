#include <cassert>
#include <dune/grid/common/referenceelements.hh>

namespace Dune {

  template <class ct, int dim>
  typename PointProvider<ct, dim, 0>::PointContainerType
  PointProvider<ct, dim, 0>::points_;

  template <class ct, int dim>
  void PointProvider<ct, dim, 0>::
  registerQuadrature(const QuadratureType& quad)
  {
    if (points_.find(quad.id()) == points_.end()) {
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
    // std::cout << "Add entry for id " << quad.id() << " called\n";
    const ReferenceElement<ct, dim>& refElem =
      ReferenceElements<ct, dim>::general(elementGeo);
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

    return mit;
  }

  template <class ct, int dim>
  bool PointProvider<ct, dim, 1>::sameGeometry(GeometryType geo1, 
                                               GeometryType geo2)
  {
    return geo1 == geo2;
    /*
    // Assume here that the geometries belong to the same dimension

    switch (geo1) {
    case vertex:
      return (geo2 == vertex || geo2 == simplex || geo2 == cube);
    case line:
      return (geo2 == line || geo2 == simplex || geo2 == cube);
    case triangle:
      return (geo2 == triangle || geo2 == simplex);
    case quadrilateral:
      return (geo2 == quadrilateral || geo2 == cube);
    case tetrahedron:
      return (geo2 == tetrahedron || geo2 == simplex);
    case pyramid:
      return geo2 == pyramid;
    case prism:
      return geo2 == prism;
    case hexahedron:
      return (geo2 == hexahedron || geo2 == cube);
    case simplex:
      return (geo2 == simplex || geo2 == triangle || geo2 == tetrahedron
              || geo2 == line || (dim-1 == 1 && geo2 == cube));
    case cube:
      return (geo2 == cube || geo2 == quadrilateral || geo2 == hexahedron
              || geo2 == line || (dim-1 == 1 && geo2 == simplex));
    default:
      assert(false);
    }
    return false;
    */
  }  
} // end namespace Dune
