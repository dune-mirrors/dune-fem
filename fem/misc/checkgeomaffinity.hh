#ifndef DUNE_CHECKGEOMETRYAFFINITY_HH
#define DUNE_CHECKGEOMETRYAFFINITY_HH

namespace Dune {

/*! @addtogroup HelperClasses
 ** @{
 */
  
//! Helper class to check affinity of the grid's geometries  
template <class QuadratureType>
struct GeometryAffinityCheck 
{
  //! check whether all geometry mappings are affine 
  template <class IteratorType>
  static inline bool checkAffinity(const IteratorType& begin,
                                   const IteratorType& endit, 
                                   const int quadOrd)  
  {
    bool affinity = true ;
    typedef typename IteratorType :: Entity :: Geometry Geometry;
    for(IteratorType it = begin; it != endit; ++it)
    {
      // get quadrature of desired order 
      QuadratureType volQuad( *it, quadOrd );
      const int nop = volQuad.nop();
      const Geometry& geo = it->geometry();

      // check all integration elements against the first 
      const double oldIntel = geo.integrationElement( volQuad.point(0) );
      for(int l=0; l<nop; ++l)
      {
        const double intel = geo.integrationElement( volQuad.point(l) );
        if( std::abs( oldIntel - intel ) > 1e-12 ) affinity = false;
      }
    }
    return affinity;
  }
};

//! @}  
} // end namespace Dune 
#endif
