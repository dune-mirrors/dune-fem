#ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_FUNCTOR_HH
#define DUNE_FEM_FUNCTION_LOCALFUNCTION_FUNCTOR_HH

// dune-common includes
#include <dune/common/fvector.hh>

/**
  @file
  @author Christoph Gersbacher
  @brief  Provides functors for evaluating local functions
*/


namespace Dune
{

  namespace Fem
  {

    // LocalFunctionEvaluateFunctor
    // ----------------------------

    template< class LocalFunction, int diffOrder = 0 >
    struct LocalFunctionEvaluateFunctor;

    template< class LocalFunction >
    struct LocalFunctionEvaluateFunctor< LocalFunction, 0 >
    {
      typedef LocalFunction LocalFunctionType;

      typedef typename LocalFunctionType::LocalCoordinateType LocalCoordinateType;
      typedef typename LocalFunctionType::RangeType RangeType;

      LocalFunctionEvaluateFunctor ( RangeType &value )
      : value_( value )
      {}

      void operator() ( const LocalCoordinateType &x, const LocalFunction &localFunction )
      {
        localFunction.evaluate( x, value_ );
      }

    private:
      RangeType &value_;
    };

    template< class LocalFunction, int diffOrder >
    struct LocalFunctionEvaluateFunctor
    {
      typedef LocalFunction LocalFunctionType;

      typedef typename LocalFunctionType::LocalCoordinateType LocalCoordinateType;
      typedef typename LocalFunctionType::RangeType RangeType;

      LocalFunctionEvaluateFunctor ( const Dune::FieldVector< int, diffOrder > &diffVariable,
                                     RangeType &value )
      : value_( value ),
        diffVariable_( &diffVariable )
      {}

      void operator() ( const LocalCoordinateType &x, const LocalFunction &localFunction )
      {
        localFunction.evaluate( x, value_ );
      }

    private:
      RangeType &value_;
      const Dune::FieldVector< int, diffOrder > *diffVariable_;
    };



    // LocalFunctionJacobianFunctor
    // ----------------------------

    template< class LocalFunction >
    struct LocalFunctionJacobianFunctor
    {
      typedef LocalFunction LocalFunctionType;

      typedef typename LocalFunctionType::LocalCoordinateType LocalCoordinateType;
      typedef typename LocalFunctionType::JacobianRangeType JacobianRangeType;

      LocalFunctionJacobianFunctor ( JacobianRangeType &jacobian )
      : jacobian_( jacobian )
      {}

      void operator() ( const LocalCoordinateType &x, const LocalFunction &localFunction )
      {
        localFunction.jacobian( x, jacobian_);
      }

    private:
      JacobianRangeType &jacobian_;
    };



    // LocalFunctionHessianFunctor
    // ---------------------------

    template< class LocalFunction >
    struct LocalFunctionHessianFunctor
    {
      typedef LocalFunction LocalFunctionType;

      typedef typename LocalFunctionType::LocalCoordinateType LocalCoordinateType;
      typedef typename LocalFunctionType::HessianRangeType HessianRangeType;

      LocalFunctionHessianFunctor ( HessianRangeType &hessian )
      : hessian_( hessian )
      {}

      void operator() ( const LocalCoordinateType &x, const LocalFunction &localFunction )
      {
        localFunction.hessian( x, hessian_ );
      }

    private:
      HessianRangeType &hessian_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_FUNCTOR_HH
