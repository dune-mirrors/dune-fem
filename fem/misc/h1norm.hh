#ifndef DUNE_FEM_H1NORM_HH
#define DUNE_FEM_H1NORM_HH

#include <dune/fem/quadrature/cachequad.hh>

namespace Dune
{

  template< class GridPartImp >
  class H1Norm
  {
  public:
    typedef GridPartImp GridPartType;

  private:
    typedef H1Norm< GridPartType > ThisType;

  protected:
    typedef typename GridPartType :: template Codim< 0 > :: IteratorType
      GridIteratorType;
    typedef typename GridIteratorType :: Entity EntityType;
    typedef typename EntityType :: Geometry GeometryType;
    typedef CachingQuadrature< GridPartType, 0 > QuadratureType;

  protected:
    const GridPartType &gridPart_;

  public:
    inline explicit H1Norm ( const GridPartType &gridPart );
    inline H1Norm ( const ThisType &other );

  private:
    // prohibit assignment
    ThisType operator= ( const ThisType &other );

  public:
    template< class DiscreteFunctionType >
    inline typename DiscreteFunctionType :: RangeFieldType
    norm ( const DiscreteFunctionType &u ) const;
    
    template< class UDiscreteFunctionType, class VDiscreteFunctionType >
    inline typename UDiscreteFunctionType :: RangeFieldType
    distance ( const UDiscreteFunctionType &u,
               const VDiscreteFunctionType &v ) const;
  };

}

#include "h1norm_inline.hh"

#endif
