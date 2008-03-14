#ifndef DUNE_CORNERPROBLEM_HH
#define DUNE_CORNERPROBLEM_HH

namespace Dune
{

  // right hand side of governing problem 
  template< class FunctionSpace >
  class RHSFunction
  : public Function< FunctionSpace, RHSFunction< FunctionSpace > >
  {
    typedef RHSFunction< FunctionSpace > ThisType;
    typedef Function< FunctionSpace, ThisType > BaseType;

  public:
    typedef FunctionSpace FunctionSpaceType;

    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;

    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

  public:
    inline RHSFunction ( const FunctionSpaceType &functionSpace )
    : BaseType( functionSpace )
    {}
     
    // f( x, y, z ) = 2 \sum_{i=1}^{dimworld} \prod_{j \neq i} (x_j-x_j^2)
    inline void evaluate( const DomainType &x , RangeType &phi ) const
    {
      phi = 0;
    }
  }; // end class RHSFunction



  template< class FunctionSpace >
  class ExactSolution
    : public Function< FunctionSpace, ExactSolution< FunctionSpace > >
  {
    typedef ExactSolution< FunctionSpace > ThisType;
    typedef Function< FunctionSpace, ThisType > BaseType;

  public:
    typedef FunctionSpace FunctionSpaceType;

    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;
    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

  public:
    inline ExactSolution ( const FunctionSpaceType &functionSpace )
    : BaseType( functionSpace )
    {}

    inline void evaluate ( const DomainType &x, RangeType &ret ) const
    {
      const double r = x.two_norm();
      const double fac = 2.0 / 3.0;
 
      if( r > 1e-8 )
      {
        const double phi
          = (x[ 1 ] >= 0 ? acos( x[ 0 ] / r ) : acos( -x[ 0 ] / r ) + M_PI);
      
        ret = sin( fac * phi );
      
        ret *= pow( r, fac );
      }
      else
        ret = 0;
    }

    inline void jacobian ( const DomainType &x, JacobianRangeType &ret ) const
    {
      const double r = x.two_norm();
      const double fac = 2.0 / 3.0;

      const double phi
        = (x[ 1 ] >= 0 ? acos( x[ 0 ] / r ) : acos( -x[ 0 ] / r ) + M_PI);

      ret[ 0 ][ 0 ] = sin( fac * phi ) * x[ 0 ] - cos( fac * phi ) * x[ 1 ];
      ret[ 0 ][ 1 ] = sin( fac * phi ) * x[ 1 ] + cos( fac * phi ) * x[ 0 ];
      
      ret *= fac * pow( r, -2 * fac );
    }

    inline void evaluate( const DomainType &x, RangeFieldType t, RangeType &ret ) const
    {
      evaluate( x, ret );
    }
  };



  // diffusion coefficient for this problem
  template< class FunctionSpace >
  class Tensor
  : public Function< FunctionSpace, Tensor< FunctionSpace > >
  {
    typedef Tensor< FunctionSpace > ThisType;
    typedef Function< FunctionSpace, ThisType > BaseType;

  public:
    typedef FunctionSpace FunctionSpaceType;

    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;

    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

  public:
    inline Tensor ( const FunctionSpaceType &functionSpace )
    : BaseType( functionSpace )
    {}
    
    inline void evaluate ( int i, int j, const DomainType &x, RangeType &phi ) const
    {
      evaluate( x, phi );
    }
    
    inline void evaluate ( const DomainType &x, RangeType &phi ) const
    {
      phi = 1;
    }

    inline void evaluate ( const DomainType &x, RangeFieldType time, RangeType &phi ) const
    {
      evaluate( x, phi );
    }
  }; //end class Tensor

}

#endif
