#ifndef DUNE_FEM_SHAPEFUNCTIONSET_SHAPEFUNCTIONSET_HH
#define DUNE_FEM_SHAPEFUNCTIONSET_SHAPEFUNCTIONSET_HH

//- C++ includes
#include <cstddef>

/**
  @file
  @brief Interface for shape function sets
*/


namespace Dune
{

  namespace Fem
  {

    // ShapeFunctionSet
    // ----------------

    /**
      \brief Interface class for shape function sets

       This class cannot be used itself, it is for documentation purposes
       only.

       \note Constructor signatures are explicitly not specified by this
             interface.
     */
    template< class FunctionSpace >
    class ShapeFunctionSet
    {
    public:
      //! \brief function space type
      typedef FunctionSpace FunctionSpaceType;

      //! \brief domain type
      typedef typename FunctionSpaceType::DomainType DomainType;
      //! \brief range type
      typedef typename FunctionSpaceType::RangeType RangeType;
      //! \brief jacobian range type
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      //! \brief hessian range type
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      //! \brief return order of shape functions
      int order () const;

      //! \brief return number of shape functions
      std::size_t size () const;

      /**
         \brief evalute each shape function

          \param[in]  x        coordinate or quadrature point
          \param[in]  functor  functor call for evaluating each shape function

          The functor has to be a copyable object satisfying the following
          interface:
          \code
          struct Functor
          {
            template< class Value >
            void operator() ( const int shapeFunction, const Value &value );
          };
          \endcode
       */
      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor functor ) const;

      /**
         \brief evalute jacobian of each shape function

          \param[in]  x        coordinate or quadrature point
          \param[in]  functor  functor call for evaluating the jacobian of each shape function

          The functor has to be a copyable object satisfying the following
          interface:
          \code
          struct Functor
          {
            template< class Jacobian >
            void operator() ( const int shapeFunction, const Jacobian &jacobian );
          };
          \endcode
       */
      template< class Point, class Functor >
      void jacobianEach ( const Point &x, Functor functor ) const;

      /**
         \brief evalute hessian of each shape function

          \param[in]  x        coordinate or quadrature point
          \param[in]  functor  functor call for evaluating the hessian of each shape function

          The functor has to be a copyable object satisfying the following
          interface:
          \code
          struct Functor
          {
            template< class Hessian >
            void operator() ( const int shapeFunction, const Hessian &hessian );
          };
          \endcode
       */
      template< class Point, class Functor >
      void hessianEach ( const Point &x, Functor functor ) const;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SHAPEFUNCTIONSET_SHAPEFUNCTIONSET_HH
