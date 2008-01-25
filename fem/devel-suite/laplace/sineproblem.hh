#ifndef DUNE_EXAMPLES_LAPLACE_SINEPROBLEM_HH
#define DUNE_EXAMPLES_LAPLACE_SINEPROBLEM_HH

namespace Dune
{

  template< class FunctionSpace >
  class SineProblem
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
  class SineProblem< FunctionSpace > :: RightHandSide
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
        RangeType z = M_PI * M_PI;
        for( int j = 0; j < dimension; ++j )
        {
          const DomainFieldType &xj = x[ j ];
          z *= sin( M_PI * xj );
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
  class SineProblem< FunctionSpace > :: ExactSolution
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
        y *= sin( M_PI * xi );
      }
    }

    inline void jacobian ( const DomainType &x, JacobianRangeType &ret ) const
    {
      enum { dimension = DomainType :: dimension };
     
      for( int i = 0; i < dimension; ++i )
      {
        ret[ 0 ][ i ] = M_PI;
        for( int j = 0; j < dimension; ++j )
        {
          const DomainFieldType &xj = x[ j ];
          ret[ 0 ][ i ] *= ((i == j) ? cos( M_PI * xj ) : sin( M_PI * xj ));
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
