#ifndef DUNE_POINTMAPPER_HH
#define DUNE_POINTMAPPER_HH

//- system includes 
#include <vector>

//- Dune includes 
#include <dune/common/version.hh>
#include <dune/fem/quadrature/quadrature.hh>

namespace Dune
{

  struct QuadratureKey  
  {
    QuadratureKey ( const GeometryType &geoType, const size_t id )
    : id_( ((topologyId( geoType ) >> 1) << 16) + id )
    {
      assert( id < (1 << 16) );
    }
    
    bool operator< ( const QuadratureKey &other ) const
    {
      return (id_ < other.id_);
    }
    
    bool operator== ( const QuadratureKey &other ) const
    {
      return (id_ == other.id_);
    }

    friend std::ostream &operator<< ( std::ostream &out, const QuadratureKey &key )
    {
      return out << "(topologyId " << ((key.id_ >> 16) << 1) << ", quadId " << (key.id_ & ((1u << 16)-1)) << ")";
    }

  protected:
    static unsigned int topologyId ( const GeometryType &type )
    {
      return type.id();
    }

    const size_t id_;
  };



  template <class ct, int dim>
  struct CachingTraits {
    //! type of integration point list implementation, fix type here 
    typedef IntegrationPointListImp<ct, dim> QuadratureType;
    //! extracted types from integration point list 
    typedef typename QuadratureType::CoordinateType PointType;
    typedef std::vector<PointType>    PointVectorType;
    typedef std::vector<size_t>       MapperType;
    typedef std::vector<MapperType>   MapperVectorType;

    typedef QuadratureKey QuadratureKeyType;

    // minimal twist is -4 for hexahedrons 
    // so we add 4 to start from zero 
    enum { twistOffset_ = 4 };

  };
  

} // end namespace Dune

#endif
