#ifndef DUNE_FEM_REDUCEDBASISSPACE_SINESPACE_HH
#define DUNE_FEM_REDUCEDBASISSPACE_SINESPACE_HH

#include <dune/fem/space/reducedbasisspace.hh>

namespace Dune
{

  template< class FunctionSpaceImp >
  class SineBaseFunction
  : public Function< FunctionSpaceImp, SineBaseFunction< FunctionSpaceImp > >
  {
  public:
    typedef FunctionSpaceImp FunctionSpaceType;

  private:
    typedef SineBaseFunction< FunctionSpaceType > ThisType;
    typedef Function< FunctionSpaceType, ThisType > BaseType;

  public:
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;

    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

    enum { dimDomain = FunctionSpaceType :: dimDomain };
    enum { dimRange = FunctionSpaceType :: dimRange };

    typedef FieldVector< int, dimDomain > CoefficientType;
    
  protected:
    const CoefficientType coefficient_;
    
  public:
    inline SineBaseFunction ( const FunctionSpaceType &functionSpace,
                              const CoefficientType coefficient )
    : BaseType( functionSpace ),
      coefficient_( coefficient )
    {
    }

    inline void evaluate ( const DomainType &x, RangeType &y ) const
    {
      y = 1;
      for( unsigned int i = 0; i < dimDomain; ++i )
      {
        /*
        if( coefficient_[ i ] < 0 )
          y *= sqrt( 2 ) * cos( 2 * M_PI * coefficient_[ i ] * x[ i ] );
        else if( coefficient_[ i ] > 0 )
          y *= sqrt( 2 ) * sin( 2 * M_PI * coefficient_[ i ] * x[ i ] );
        */
        y *= sqrt( 2 ) * sin( M_PI * coefficient_[ i ] * x[ i ] );
      }
    }

    inline void evaluate ( const DomainType &x, const RangeFieldType t, RangeType &y ) const
    {
      evaluate( x, y );
    }
  };



  template< class BaseFunctionSpaceImp, unsigned int maxCoefficient >
  class SineReducedBasisSpace
  : public ReducedBasisSpace< AdaptiveDiscreteFunction< BaseFunctionSpaceImp > >
  {
  public:
    typedef BaseFunctionSpaceImp BaseFunctionSpaceType;

    typedef AdaptiveDiscreteFunction< BaseFunctionSpaceType > BaseFunctionType;

  private:
    typedef SineReducedBasisSpace< BaseFunctionSpaceType, maxCoefficient > ThisType;
    typedef ReducedBasisSpace< BaseFunctionType > BaseType;

    using BaseType :: baseFunctionSpace_;

  public:
    typedef SineBaseFunction< typename BaseFunctionSpaceType :: FunctionSpaceType >
      ContinuousBaseFunctionType;

  private:
    typedef typename ContinuousBaseFunctionType :: CoefficientType CoefficientType;

  public:
    inline explicit SineReducedBasisSpace ( BaseFunctionSpaceType &baseFunctionSpace )
    : BaseType( baseFunctionSpace )
    {
      // CoefficientType coefficient( -maxCoefficient );
      CoefficientType coefficient( 1 );
      while( true )
      {
        addBaseFunction( coefficient );

        ++coefficient[ 0 ];
        for( unsigned int i = 0; coefficient[ i ] > (int)maxCoefficient; ++i )
        {
          // coefficient[ i ] = -maxCoefficient;
          coefficient[ i ] = 1;
          if( i+1 < CoefficientType :: dimension )
            ++coefficient[ i+1 ];
          else
            return;
        }
      }
    }

  private:
    static inline int abs( const CoefficientType &coefficient )
    {
      int value = 0;
      for( unsigned int i = 0; i < CoefficientType :: dimension; ++i )
        value += (coefficient[ i ] < 0 ? -coefficient[ i ] : coefficient[ i ]);
      return value;
    }
    
    inline void addBaseFunction( const CoefficientType &coefficient )
    {
      BaseFunctionType discreteBaseFunction( "base function", baseFunctionSpace_ );
      ContinuousBaseFunctionType continuousBaseFunction( baseFunctionSpace_, coefficient );
      LagrangeInterpolation< BaseFunctionType >
        :: interpolateFunction( continuousBaseFunction, discreteBaseFunction );
      BaseType :: addBaseFunction( discreteBaseFunction );
    }
  };

}

#endif
