#ifndef DUNE_FEM_L2NORM_HH
#define DUNE_FEM_L2NORM_HH

#include <dune/fem/quadrature/cachequad.hh>

namespace Dune
{

  template< class GridPartImp >
  class L2Norm
  {
  public:
    typedef GridPartImp GridPartType;

  private:
    typedef L2Norm< GridPartType > ThisType;

  protected:
    typedef typename GridPartType :: template Codim< 0 > :: IteratorType
      GridIteratorType;
    typedef typename GridIteratorType :: Entity EntityType;
    typedef typename EntityType :: Geometry GeometryType;
    typedef CachingQuadrature< GridPartType, 0 > QuadratureType;

  protected:
    const GridPartType &gridPart_;

  public:
    inline explicit L2Norm ( const GridPartType &gridPart );
    inline L2Norm ( const ThisType &other );

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


  
  template< class WeightFunctionImp >
  class WeightedL2Norm
  {
  public:
    typedef WeightFunctionImp WeightFunctionType;
    
  private:
    typedef WeightedL2Norm< WeightFunctionType > ThisType;

  public:
    typedef typename WeightFunctionType :: DiscreteFunctionSpaceType
      WeightFunctionSpaceType;
    typedef typename WeightFunctionSpaceType :: GridPartType GridPartType;

  protected:
    typedef typename WeightFunctionType :: LocalFunctionType
      LocalWeightFunctionType;
    typedef typename WeightFunctionType :: RangeType WeightType;
    typedef typename GridPartType :: template Codim< 0 > :: IteratorType
      GridIteratorType;
    typedef typename GridIteratorType :: Entity EntityType;
    typedef typename EntityType :: Geometry GeometryType;
    typedef CachingQuadrature< GridPartType, 0 > QuadratureType;

  protected:
    const WeightFunctionType &weightFunction_;
    const GridPartType &gridPart_;

  public:
    inline explicit WeightedL2Norm ( const WeightFunctionType &weightFunction );
    inline WeightedL2Norm ( const ThisType &other );

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

#include "l2norm_inline.hh"

#endif
