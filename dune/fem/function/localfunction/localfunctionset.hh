#ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_LOCALFUNCTIONSET_HH
#define DUNE_FEM_FUNCTION_LOCALFUNCTION_LOCALFUNCTIONSET_HH

#include <cstddef>

#include <dune/fem/space/common/functionspace.hh>

namespace Dune
{

  namespace Fem
  {

    // LocalFunctionSet
    // ----------------

    /** \class LocalFunctionSet
     *
     *  \brief Local basis functions
     *
     *  This class documents the local basis functions interface.
     *
     *  \tparam  FunctionSpace  function space
     */
    template< class Entity, class Range >
    struct LocalFunctionSet 
    {
      //! \brief entity type
      typedef Entity EntityType;

      //! \brief function space type
      typedef FunctionSpace< typename Entity::ctype, typename Range::value_type,
                             Entity::dimensionworld, Range::dimension > FunctionSpaceType;
      
      //! \brief domain type
      typedef typename FunctionSpaceType::DomainType DomainType;
      //! \brief range type
      typedef typename FunctionSpaceType::RangeType RangeType;
      //! \brief jacobian range type
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      //! \brief hessian range type
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      //! \brief return order of basis functions
      int order () const;

      //! \brief return entity
      const EntityType &entity () const;

      //! \brief return number of basis functions
      std::size_t size () const;

      /**
       * \brief evalute each basis function
       *
       *  \param[in]  x        local coordinate or quadrature point
       *  \param[in]  functor  functor call for evaluating each basis function
       *
       *  The functor has to be a copyable object satisfying the following
       *  interface:
\code
  struct Functor
  {
    template< class Value >
    void operator() ( const int basisFunction, const Value &value );
  };
\endcode
       */
      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor functor ) const;

      /**
       * \brief evalute jacobian of each basis function
       *
       *  \param[in]  x        local coordinate or quadrature point
       *  \param[in]  functor  functor call for evaluating the jacobian of each basis function
       *
       *  The functor has to be a copyable object satisfying the following
       *  interface:
\code
  struct Functor
  {
    template< class Jacobian >
    void operator() ( const int basisFunction, const Jacobian &jacobian );
  };
\endcode
       */
      template< class Point, class Functor >
      void jacobianEach ( const Point &x, Functor functor ) const;

      /**
       * \brief evalute hessian of each basis function
       *
       *  \param[in]  x        local coordinate or quadrature point
       *  \param[in]  functor  functor call for evaluating the hessian of each basis function
       *
       *  The functor has to be a copyable object satisfying the following
       *  interface:
\code
  struct Functor
  {
    template< class Hessian >
    void operator() ( const int basisFunction, const Hessian &hessian );
  };
\endcode
       */
      template< class Point, class Functor >
      void hessianEach ( const Point &x, Functor functor ) const;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_LOCALFUNCTIONSET_HH
