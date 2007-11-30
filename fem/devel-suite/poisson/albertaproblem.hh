#ifndef DUNE_ALBERTAPROBLEM_HH
#define DUNE_ALBERTAPROBLEM_HH

/* This problem is taken from the Alberta demo */

namespace Dune
{

  // right hand side of governing problem 
  template< class FunctionSpaceImp >
  class RHSFunction
  : public Function< FunctionSpaceImp, RHSFunction< FunctionSpaceImp > >
  {
  public:
    typedef FunctionSpaceImp FunctionSpaceType;

  private:
    typedef RHSFunction< FunctionSpaceType > ThisType;
    typedef Function< FunctionSpaceType, ThisType > BaseType;

  public:
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;

    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

  public:
    inline RHSFunction ( const FunctionSpaceType &functionSpace )
    : BaseType( functionSpace )
    {
    }
     
    // f( x, y, z ) = 2 \sum_{i=1}^{dimworld} \prod_{j \neq i} (x_j-x_j^2)
    inline void evaluate( const DomainType &x , RangeType &phi ) const
    {
      enum { dimension = DomainType :: dimension };
      
      const RangeFieldType xsqr = x*x;

      phi = -(400.0 * xsqr - 20.0 * dimension) * exp( -10.0 * xsqr );
    }
  }; // end class RHSFunction



  //! the exact solution to the problem for EOC calculation 
  template< class FunctionSpaceImp >
  class ExactSolution
  : public Function < FunctionSpaceImp, ExactSolution< FunctionSpaceImp > >
  {
  public:
    typedef FunctionSpaceImp FunctionSpaceType;

  private:
    typedef ExactSolution< FunctionSpaceType > ThisType;
    typedef Function< FunctionSpaceType, ThisType > BaseType;

  public:
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;
    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

  public:
    inline ExactSolution ( const FunctionSpaceType &functionSpace )
    : BaseType( functionSpace )
    {
    }
   
    // u( x, y, z ) = \prod_{i=1}^{dimworld} (x_i - x_i^2).
    inline void evaluate ( const DomainType &x, RangeType &phi ) const
    {
      phi = exp( -10.0 * (x*x) );
    }

    inline void jacobian ( const DomainType &x, JacobianRangeType &ret ) const
    {
      ret[ 0 ] = x;
      ret[ 0 ] *= -20.0 * exp( -10.0 * (x*x) );

    }
   
    inline void evaluate ( const DomainType &x , RangeFieldType time, RangeType &phi ) const
    {
      evaluate( x, phi );
    }
  }; // end class ExactSolution



  // diffusion coefficient for this problem
  template< class FunctionSpaceImp >
  class Tensor
  : public Function< FunctionSpaceImp, Tensor< FunctionSpaceImp > >
  {
  public:
    typedef FunctionSpaceImp FunctionSpaceType;

  private:
    typedef Tensor< FunctionSpaceType > ThisType;
    typedef Function< FunctionSpaceType, ThisType > BaseType;

  public:
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;

    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

  public:
    inline Tensor ( const FunctionSpaceType &functionSpace )
    : BaseType( functionSpace )
    {
    }
    
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
