#ifndef DUNE_EXAMPLES_LAPLACE_QUADRATICPROBLEM_HH
#define DUNE_EXAMPLES_LAPLACE_QUADRATICPROBLEM_HH

namespace Dune
{

  template< class FunctionSpace >
  class QuadraticProblem
  {
  public:
    typedef FunctionSpace FunctionSpaceType;

  private:
    class RightHandSide;
    class ExactSolution;

  public:
    typedef RightHandSide RightHandSideType;
    typedef ExactSolution ExactSolutionType;
  };



  template< class FunctionSpace >
  class QuadraticProblem< FunctionSpace > :: RightHandSide
  : public Function< FunctionSpace, RightHandSide >
  {
  public:
    typedef FunctionSpace FunctionSpaceType;
  
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;

    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;
    
  private:
    typedef Function< FunctionSpaceType, RightHandSide > BaseType;

  public:
    inline RightHandSide ( const FunctionSpaceType &functionSpace )
    : BaseType( functionSpace )
    {
    }

    inline void evaluate ( const DomainType &x, RangeType &y ) const
    {
      enum { dimension = DomainType :: dimension };

      y = 0;
      for( int i = 0; i < dimension; ++i )
      {
        RangeType z = 2;
        for( int j = 0; j < dimension; ++j )
        {
          if( i == j )
            continue;
          
          const DomainFieldType &xj = x[ j ];
          z *= xj - xj * xj;
        }
        y += z;
      }
    }

    inline void evaluate ( const DomainType &x, RangeFieldType t, RangeType &y ) const
    {
      evaluate( x, y );
    }
  };



  template< class FunctionSpace >
  class QuadraticProblem< FunctionSpace > :: ExactSolution
  : public Function< FunctionSpace, ExactSolution >
  {
  public:
    typedef FunctionSpace FunctionSpaceType;
    
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;
    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

  private:
    typedef Function< FunctionSpaceType, ExactSolution > BaseType;

  public:
    inline ExactSolution ( FunctionSpaceType &functionSpace )
    : BaseType( functionSpace )
    {
    }

    inline void evaluate ( const DomainType &x, RangeType &y ) const
    {
      enum { dimension = DomainType :: dimension };
      
      y = 1;
      for( int i = 0; i < dimension; ++i )
      {
        const DomainFieldType &xi = x[ i ];
        y *= xi - xi * xi;
      }
    }

    inline void jacobian ( const DomainType &x, JacobianRangeType &ret ) const
    {
      enum { dimension = DomainType :: dimension };
     
      for( int i = 0; i < dimension; ++i )
      {
        ret[ 0 ][ i ] = 1;
        for( int j = 0; j < dimension; ++j )
        {
          const DomainFieldType &xj = x[ j ];
          ret[ 0 ][ i ] *= ((i == j) ? (1 - 2 * xj) : xj - xj * xj);
        }
      }
    }

    inline void evaluate ( const DomainType &x, RangeFieldType t, RangeType &y ) const
    {
      evaluate( x, y );
    }
  };

}

#endif
