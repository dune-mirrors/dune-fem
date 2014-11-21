#ifndef DUNE_FEM_TEST_WEIGHTFUNCTION_HH
#define DUNE_FEM_TEST_WEIGHTFUNCTION_HH

#include <dune/fem/function/common/function.hh>

namespace Dune
{

  template< class FunctionSpaceImp >
  class WeightFunction
  : public Fem::Function< FunctionSpaceImp, WeightFunction< FunctionSpaceImp > >
  {
    typedef WeightFunction< FunctionSpaceImp > ThisType;
    typedef Fem::Function< FunctionSpaceImp, ThisType > BaseType;

  public:
    typedef FunctionSpaceImp FunctionSpaceType;

    typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  public:
    void evaluate( const DomainType &x, RangeType &phi ) const
    {
      phi = 1;
      if( x[0] > 0.5 )
        phi *= 0.5;
    }

    void evaluate( const DomainType &x, RangeFieldType t, RangeType &phi ) const
    {
      evaluate( x, phi );
    }
  };

}

#endif // #ifndef DUNE_FEM_TEST_WEIGHTFUNCTION_HH
