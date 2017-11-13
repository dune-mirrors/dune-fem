#ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_DIFFERENCE_HH
#define DUNE_FEM_FUNCTION_LOCALFUNCTION_DIFFERENCE_HH

#include <algorithm>
#include <functional>

#include "continuous.hh"

namespace Dune
{

  namespace Fem
  {

    // LocalFunctionDifference
    // -----------------------

    template< class LocalFunction, class LocalDiscreteFunction >
    class LocalFunctionDifference
    : public Dune::Fem::DefaultContinuousLocalFunction< typename LocalFunction::EntityType, typename LocalFunction::RangeType, LocalFunctionDifference< LocalFunction, LocalDiscreteFunction > >
    {
      using BaseType = Dune::Fem::DefaultContinuousLocalFunction< typename LocalFunction::EntityType, typename LocalFunction::RangeType, LocalFunctionDifference< LocalFunction, LocalDiscreteFunction > >;

    public:
      using RangeType = typename BaseType::RangeType;
      using JacobianRangeType = typename BaseType::JacobianRangeType;
      using HessianRangeType = typename BaseType::HessianRangeType;

      using EntityType = typename BaseType::EntityType;

      /** \name Construction
       *  \{
       */

      LocalFunctionDifference ( const LocalFunction &u, const LocalDiscreteFunction &uh )
        : u_( u ),
          uh_( uh )
      {}

      /** \} */

      /** \brief Public member methods
       *  \{
       */

      int order () const { return std::max( u().order(), uh().order() ); }

      template< class Point >
      void evaluate ( const Point &x, RangeType &value ) const
      {
        u().evaluate( x, value );
        RangeType tmp;
        uh().evaluate( x, tmp );
        value -= tmp;
      }

      template< class Point >
      void jacobian ( const Point &x, JacobianRangeType &jacobian ) const
      {
        u().jacobian( x, jacobian );
        JacobianRangeType tmp;
        uh().jacobian( x, tmp );
        jacobian -= tmp;
      }

      template< class Point >
      void hessian ( const Point &x, HessianRangeType &hessian ) const
      {
        u().hessian( x, hessian );
        JacobianRangeType tmp;
        uh().hessian( x, tmp );
        hessian -= tmp;
      }

      const EntityType &entity () const { return u().entity(); }

      /** \} */

    private:
      const LocalFunction &u () const { return u_.get(); }

      const LocalDiscreteFunction &uh () const { return uh_.get(); }

      std::reference_wrapper< const LocalFunction > u_;
      std::reference_wrapper< const LocalDiscreteFunction > uh_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_DIFFERENCE_HH
