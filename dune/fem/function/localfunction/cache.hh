#ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_CACHE_HH
#define DUNE_FEM_FUNCTION_LOCALFUNCTION_CACHE_HH

#include <cassert>

#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/fem/quadrature/quadrature.hh>
#include <dune/fem/space/common/derivatives.hh>

namespace Dune
{

  namespace Fem
  {

    // LocalFunctionCache
    // ------------------

    template< class FunctionSpace >
    class LocalFunctionCache
    {
      typedef LocalFunctionCache< FunctionSpace > ThisType;

    public:
      //! type of functionspace
      typedef FunctionSpace FunctionSpaceType;

      //! field type of the domain
      typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
      //! field type of the range
      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
      //! type of domain vectors, i.e., type of coordinates
      typedef typename FunctionSpaceType::DomainType DomainType;
      //! type of range vectors, i.e., type of function values
      typedef typename FunctionSpaceType::RangeType RangeType;
      //! type of the Jacobian, i.e., type of evaluated Jacobian matrix
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      //! type of the Hessian
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      //! dimension of the domain
      static const int dimDomain = FunctionSpaceType::dimDomain;
      //! dimension of the range
      static const int dimRange = FunctionSpaceType::dimRange;

      //! default constructor, calls default ctor of BasisFunctionSetType and LocalDofVectorType
      LocalFunctionCache () {}

      //! \todo please doc me
      template< class LocalFunction, class Quadrature, int... order >
      void init ( const LocalFunction &localFunction, const Quadrature &quadrature, Derivatives::Type< order... > )
      {
        quadratureId_ = quadrature.id();
        values_.clear();
        jacobians_.clear();
        hessians_.clear();
        std::ignore = std::make_tuple( (init( localFunction, quadrature, std::integral_constant< order >() ), order)... );
      }

      /**
       * \brief evaluate the local function
       *
       * \param[in]   x      evaluation point in local coordinates
       * \param[out]  value  value of the function in the given point
       **/
      template< class Quadrature >
      void evaluate ( const QuadraturePointWrapper< Quadrature > &x, RangeType &value ) const
      {
        assert( x.quadrature().id() == quadratureId_ );
        assert( x.index() < values_.size() );
        value = values_[ x.index() ];
      }

      /**
       * \brief evaluate Jacobian of the local function
       *
       * \note Though the Jacobian is evaluated on the reference element, the
       *       return value is the Jacobian with respect to the actual entity.
       *
       * \param[in]   x         evaluation point in local coordinates
       * \param[out]  jacobian  Jacobian of the function in the evaluation point
       */
      template< class Quadrature >
      void jacobian ( const QuadraturePointWrapper< Quadrarture > &x, JacobianRangeType &jacobian ) const
      {
        assert( x.quadrature().id() == quadratureId_ );
        assert( x.index() < jacobians_.size() );
        jacobian = jacobians_[ x.index() ];
      }

      /**
       * \brief evaluate Hessian of the local function
       *
       * \note Though the Hessian is evaluated on the reference element, the
       *       return value is the Hessian with respect to the actual entity.
       *
       * \param[in]   x        evaluation point in local coordinates
       * \param[out]  ret  Hessian of the function in the evaluation point
       */
      template< class Quadrature >
      void hessian ( const QuadraturePointWrapper< Quadrature > &x, HessianRangeType &hessian ) const
      {
        assert( x.quadrature().id() == quadratureId_ );
        assert( x.index() < hessians_.size() );
        hessian = hessians_[ x.index() ];
      }

    protected:
      template< class LocalFunction, class Quadrature >
      void init ( const LocalFunction &localFunction, const Quadrature &quadrature, std::integral_constant< 0 > )
      {
        values_.resize( quadrature.size() );
        localFunction.evaluateQuadrature( quadrature, values_ );
      }

      template< class LocalFunction, class Quadrature >
      void init ( const LocalFunction &localFunction, const Quadrature &quadrature, std::integral_constant< 1 > )
      {
        jacobians_.resize( quadrature.size() );
        localFunction.evaluateQuadrature( quadrature, jacobians_ );
      }

      template< class LocalFunction, class Quadrature >
      void init ( const LocalFunction &localFunction, const Quadrature &quadrature, std::integral_constant< 2 > )
      {
        hessians_.resize( quadrature.size() );
        for( const auto &qp : quadrature )
          localFunction.hessian( qp, hessians_[ qp.index() ] );
      }

      std::vector< RangeType > values_;
      std::vector< JacobianRangeType > jacobians_;
      std::vector< HessianRangeType > hessians_;
    };

    /** \} */

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_CACHE_HH
