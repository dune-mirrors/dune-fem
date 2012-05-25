#ifndef DUNE_FEM_REDUCEDBASISSPACE_SINESPACE_HH
#define DUNE_FEM_REDUCEDBASISSPACE_SINESPACE_HH

#include <dune/fem/space/reducedbasisspace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/lagrangeinterpolation.hh>

namespace Dune
{

  template< class FunctionSpace >
  class SineBaseFunction
  : public Fem::Function< FunctionSpace, SineBaseFunction< FunctionSpace > >
  {
    typedef SineBaseFunction< FunctionSpace > ThisType;
    typedef Fem::Function< FunctionSpace, ThisType > BaseType;

  public:
    typedef FunctionSpace FunctionSpaceType;

    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;

    typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

    static const int dimDomain = FunctionSpaceType::dimDomain;
    static const int dimRange = FunctionSpaceType::dimRange;

    typedef FieldVector< int, dimDomain > CoefficientType;
    
    explicit SineBaseFunction ( const CoefficientType coefficient )
    : coefficient_( coefficient )
    {}

    void evaluate ( const DomainType &x, RangeType &y ) const
    {
      y = 1;
      for( int i = 0; i < dimDomain; ++i )
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

    void evaluate ( const DomainType &x, const RangeFieldType t, RangeType &y ) const
    {
      evaluate( x, y );
    }

  protected:
    const CoefficientType coefficient_;
  };



  template< class BaseFunctionSpaceImp, unsigned int maxCoefficient >
  class SineReducedBasisSpace
  : public ReducedBasisSpace< AdaptiveDiscreteFunction< BaseFunctionSpaceImp > >
  {
    typedef SineReducedBasisSpace< BaseFunctionSpaceImp, maxCoefficient > ThisType;
    typedef ReducedBasisSpace< AdaptiveDiscreteFunction< BaseFunctionSpaceImp > > BaseType;

  public:
    typedef BaseFunctionSpaceImp BaseFunctionSpaceType;

    typedef AdaptiveDiscreteFunction< BaseFunctionSpaceType > BaseFunctionType;

    typedef SineBaseFunction< typename BaseFunctionSpaceType::FunctionSpaceType >
      ContinuousBaseFunctionType;

  private:
    typedef typename ContinuousBaseFunctionType::CoefficientType CoefficientType;

  public:
    explicit SineReducedBasisSpace ( BaseFunctionSpaceType &baseFunctionSpace )
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
    static int abs( const CoefficientType &coefficient )
    {
      int value = 0;
      for( unsigned int i = 0; i < CoefficientType :: dimension; ++i )
        value += (coefficient[ i ] < 0 ? -coefficient[ i ] : coefficient[ i ]);
      return value;
    }
    
    void addBaseFunction( const CoefficientType &coefficient )
    {
      BaseFunctionType discreteBaseFunction( "base function", baseFunctionSpace_ );
      ContinuousBaseFunctionType continuousBaseFunction( coefficient );
      LagrangeInterpolation< BaseFunctionType >
        ::interpolateFunction( continuousBaseFunction, discreteBaseFunction );
      BaseType::addBaseFunction( discreteBaseFunction );
    }

    using BaseType::baseFunctionSpace_;
  };

}

#endif // #ifndef DUNE_FEM_REDUCEDBASISSPACE_SINESPACE_HH
