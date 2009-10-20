#ifndef DUNE_FEM_TEST_EXACTSOLUTION_HH
#define DUNE_FEM_TEST_EXACTSOLUTION_HH

#include <dune/fem/function/common/function.hh>

namespace Dune
{

  template< class FunctionSpaceImp >
  class ExactSolution
  : public Function< FunctionSpaceImp, ExactSolution< FunctionSpaceImp > >
  {
  public:
    typedef FunctionSpaceImp FunctionSpaceType;

  private:
    typedef ExactSolution< FunctionSpaceType > ThisType;
    typedef Function< FunctionSpaceType, ThisType > BaseType;
    
  public:
    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;
    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

  public:
    ExactSolution( FunctionSpaceType &functionSpace )
    : BaseType( functionSpace )
    {
    }

    void evaluate( const DomainType &x, RangeType &phi ) const
    {
      phi = 1;
      for(int r = 0; r < RangeType :: dimension; ++r )
        for( int i = 0; i < DomainType :: dimension; ++i )
          phi[ r ] += pow(sin( M_PI * x[ i ] ),double(r+1)); 
    }

    void evaluate( const DomainType &x, RangeFieldType t, RangeType &phi ) const
    {
      evaluate( x, phi );
    }

    void jacobian( const DomainType &x, JacobianRangeType &Dphi ) const
    {
      Dphi = 1;
      for (int r = 0; r < RangeType :: dimension; ++r) 
        for( int i = 0; i < DomainType :: dimension; ++i )
          for( int j = 0; j < DomainType :: dimension; ++j )
            Dphi[ r ][ j ] += double(r+1)*pow(sin( M_PI * x[ i ] ),double(r))*
              ((i != j) ? 0 : M_PI * cos( M_PI * x[ i ] ));
    }

    void jacobian( const DomainType &x, RangeFieldType t, JacobianRangeType &Dphi ) const
    {
      jacobian( x, Dphi );
    }
  };

}
  
#endif
