#ifndef DUNE_FEM_L2NORM_HH
#define DUNE_FEM_L2NORM_HH

#include <dune/fem/quadrature/cachequad.hh>
#include <dune/fem/quadrature/integrator.hh>

namespace Dune
{

  template< class GridPart >
  class L2Norm
  {
  public:
    typedef GridPart GridPartType;

  private:
    typedef L2Norm< GridPartType > ThisType;

  protected:
    template< class Function >
    class FunctionSquare;

    template< class UFunction, class VFunction >
    class FunctionDistance;

  protected:
    typedef typename GridPartType :: template Codim< 0 > :: IteratorType
      GridIteratorType;
    typedef typename GridIteratorType :: Entity EntityType;
    typedef CachingQuadrature< GridPartType, 0 > QuadratureType;
    typedef Integrator< QuadratureType > IntegratorType;

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

  protected:
    inline const GridPartType &gridPart () const
    {
      return gridPart_;
    }

    inline typename GridPartType :: GridType
      :: template Codim< 0 > :: CollectiveCommunication 
      comm () const
    {
      return gridPart().grid().comm();
    }
  };


  
  template< class WeightFunction >
  class WeightedL2Norm
  : public L2Norm
    < typename WeightFunction :: DiscreteFunctionSpaceType :: GridPartType >
  {
  public:
    typedef WeightFunction WeightFunctionType;

   public:
    typedef typename WeightFunctionType :: DiscreteFunctionSpaceType
      WeightFunctionSpaceType;
    typedef typename WeightFunctionSpaceType :: GridPartType GridPartType;
   
  private:
    typedef WeightedL2Norm< WeightFunctionType > ThisType;
    typedef L2Norm< GridPartType > BaseType;

  protected:
    template< class Function >
    class WeightedFunctionSquare;
    
  protected:
    typedef typename WeightFunctionType :: LocalFunctionType
      LocalWeightFunctionType;
    typedef typename WeightFunctionType :: RangeType WeightType;
    
    typedef typename BaseType :: GridIteratorType GridIteratorType;
    typedef typename BaseType :: IntegratorType IntegratorType;

    typedef typename GridIteratorType :: Entity EntityType;

  protected:
    const WeightFunctionType &weightFunction_;

  protected:
    using BaseType :: gridPart;
    using BaseType :: comm;

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
