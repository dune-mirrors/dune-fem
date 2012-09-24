#ifndef DUNE_FEM_SHAPEFUNCTIONSET_SHAPEFUNCTIONSET_HH
#define DUNE_FEM_SHAPEFUNCTIONSET_SHAPEFUNCTIONSET_HH

//- C++ includes
#include <cstddef>

//- dune-common includes
#include <dune/common/documentation.hh>

/**
  \file
  \brief Interface for shape function sets
*/


namespace Dune
{

  namespace Fem
  {

    // ShapeFunctionSet
    // ----------------

    /**
     * \brief Interface class for shape function sets
     *
     * This class cannot be used itself, it is for documentation purposes
     * only.
     *
     * \note Constructor signatures are explicitly not specified by this
     *       interface.
     */
    class ShapeFunctionSet
    {
    public:
      //! \brief function space type
      typedef ImplementationDefined FunctionSpaceType;
      
      //! \brief domain type
      typedef ImplementationDefined DomainType;
      //! \brief range type
      typedef ImplementationDefined RangeType;
      //! \brief jacobian range type
      typedef ImplementationDefined JacobianRangeType;
      //! \brief hessian range type
      typedef ImplementationDefined HessianRangeType;

      //! \brief return number of shape functions
      std::size_t size () const;

      /**
       * \brief evalute each shape function
       *
       *  \param[in]  x        coordinate or quadrature point
       *  \param[in]  functor  functor call for evaluating each shape function
       *
       *  The functor has to be a copyable object satisfying the following
       *  interface:
       *  \code
       *  struct Functor
       *  {
       *    template< class Value >
       *    void operator() ( const int shapeFunction, const Value &value );
       *  };
       *  \endcode
       */
      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor functor ) const;

      /**
       * \brief evalute jacobian of each shape function
       *
       *  \param[in]  x        coordinate or quadrature point
       *  \param[in]  functor  functor call for evaluating the jacobian of each shape function
       *
       *  The functor has to be a copyable object satisfying the following
       *  interface:
       *  \code
       *  struct Functor
       *  {
       *    template< class Jacobian >
       *    void operator() ( const int shapeFunction, const Jacobian &jacobian );
       *  };
       *  \endcode
       */
      template< class Point, class Functor > 
      void jacobianEach ( const Point &x, Functor functor ) const;

      /**
       * \brief evalute hessian of each shape function
       *
       *  \param[in]  x        coordinate or quadrature point
       *  \param[in]  functor  functor call for evaluating the hessian of each shape function
       *
       *  The functor has to be a copyable object satisfying the following
       *  interface:
       *  \code
       *  struct Functor
       *  {
       *    template< class Hessian >
       *    void operator() ( const int shapeFunction, const Hessian &hessian );
       *  };
       *  \endcode
       */
      template< class Point, class Functor > 
      void hessianEach ( const Point &x, Functor functor ) const;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SHAPEFUNCTIONSET_SHAPEFUNCTIONSET_HH
