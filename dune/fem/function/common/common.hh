#ifndef DUNE_FEM_FUNCTION_COMMON_COMMON_HH
#define DUNE_FEM_FUNCTION_COMMON_COMMON_HH

#include <dune/fem/function/common/function.hh>
#include <dune/fem/space/common/functionspace.hh>

namespace Dune
{

  namespace Fem
  {

    // Identity
    // --------

    template< class FunctionSpace >
    class Identity;

    template< class ct, int dimworld >
    class Identity< FunctionSpace< ct, ct, dimworld, dimworld > >
      : public Function< FunctionSpace< ct, ct, dimworld, dimworld >, Identity< FunctionSpace< ct, ct, dimworld, dimworld > > >
    {
      typedef Function< FunctionSpace< ct, ct, dimworld, dimworld >, Identity< FunctionSpace< ct, ct, dimworld, dimworld > > > BaseType;

    public:
      typedef typename BaseType::DomainType DomainType;
      typedef typename BaseType::RangeType RangeType;
      typedef typename BaseType::JacobianRangeType JacobianRangeType;
      typedef typename BaseType::HessianRangeType HessianRangeType;

      void evaluate ( const DomainType &x, RangeType &value ) const
      {
        value = x;
      }

      void jacobian ( const DomainType &x, JacobianRangeType &jacobian ) const
      {
        jacobian = JacobianRangeType( ct( 0 ) );
        for( int i = 0; i < dimworld; ++i )
          jacobian[ i ][ i ] = 1;
      }

      void hessian ( const DomainType &x, HessianRangeType &hessian ) const
      {
        hessian = HessianRangeType( ct( 0 ) );
      }
    };


  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_COMMON_COMMON_HH
