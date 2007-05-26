#ifndef DUNE_POINTMAPPER_HH
#define DUNE_POINTMAPPER_HH

//- system includes 
#include <vector>

//- Dune includes 
#include <dune/fem/quadrature/quadrature.hh>

namespace Dune {

  template <class ct, int dim>
  struct CachingTraits {
    //! type of integration point list implementation, fix type here 
    typedef IntegrationPointListImp<ct, dim> QuadratureType;
    //! extracted types from integration point list 
    typedef typename QuadratureType::CoordinateType PointType;
    typedef std::vector<PointType>    PointVectorType;
    typedef std::vector<size_t>       MapperType;
    typedef std::vector<MapperType>   MapperVectorType;

    static const int twistOffset_ = 5;
  };

} // end namespace Dune

#endif
