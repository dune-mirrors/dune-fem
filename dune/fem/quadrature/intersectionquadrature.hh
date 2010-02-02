#ifndef DUNE_INTERSECTIONQUADRATURE_HH
#define DUNE_INTERSECTIONQUADRATURE_HH

//- Dune includes
#include <dune/common/misc.hh>

//- Local includes
#include "elementquadrature.hh"
#include "caching/twistutility.hh"
#include "caching/pointmapper.hh"
#include "caching/cacheprovider.hh"

#include "elementquadrature.hh"
#include "cachingquadrature.hh"

namespace Dune
{
  
  /** \brief IntersectionQuadrature is a helper class for creating the appropriate face quadratures 
             for integrating over intersections. */
  template< typename FaceQuadrature, bool conforming  >
  class IntersectionQuadrature 
  {
    template < typename FaceQuadratureImp, bool isConforming > 
    struct QuadSelector
    {
      // use given quadrature 
      typedef FaceQuadratureImp  FaceQuadratureType;
    };

    template < typename FaceQuadratureImp > 
    struct QuadSelector<FaceQuadratureImp, false>
    {
      // in this case non conforming type is used 
      typedef typename  FaceQuadratureImp :: 
        NonConformingQuadratureType  FaceQuadratureType;
    };

  public:
    //! type of grid partition
    typedef typename FaceQuadrature :: GridPartType GridPartType;
      
    //! type of the grid
    typedef typename GridPartType :: GridType GridType;

    //! Type of the intersection iterator
    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename IntersectionIteratorType::Intersection IntersectionType;

    //! codimension of the element quadrature
    enum { codimension = FaceQuadrature :: codimension };
 
    //! type of intersection quadrature implementation 
    typedef typename QuadSelector<FaceQuadrature, conforming> :: FaceQuadratureType FaceQuadratureType; 
    
    //! Dimension of the world.
    enum { dimension = FaceQuadratureType ::dimension };

    //! Just another name for double...
    typedef typename FaceQuadratureType :: RealType RealType;
    //! The type of the coordinates in the codim-0 reference element.
    typedef typename FaceQuadratureType :: CoordinateType CoordinateType;

    //! for compatibility
    typedef typename GridType::template Codim< 0 >::Entity EntityType;
    
    /** \brief Constructor creating an inside and an outside face quadrature for
               integrating over an intersection. 
            
        \param[in]  gridPart      grid partition
        \param[in]  intersection  intersection
        \param[in]  order         desired order of the quadrature
 
     */ 
    IntersectionQuadrature( const GridPartType &gridPart, 
                            const IntersectionType &intersection,
                            const int order) 
    : inside_ ( gridPart, intersection, order, FaceQuadratureType :: INSIDE ),
      outside_( gridPart, intersection, order, FaceQuadratureType :: OUTSIDE ) 
    {}

    //! \brief return reference to inside face quadrature 
    const FaceQuadratureType& inside()  const { return inside_;  }

    //! \brief return reference to outside face quadrature 
    const FaceQuadratureType& outside() const { return outside_; }

  private:
    // prohibit copying 
    IntersectionQuadrature( const IntersectionQuadrature & );

  protected:  
    const FaceQuadratureType inside_;
    const FaceQuadratureType outside_;
  };

} // end namespace Dune 

#endif // #ifndef DUNE_INTERSECTIONQUADRATURE_HH
