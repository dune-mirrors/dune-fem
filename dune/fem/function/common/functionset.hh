#ifndef DUNE_FEM_FUNCTION_COMMON_FUNCTIONSET_HH
#define DUNE_FEM_FUNCTION_COMMON_FUNCTIONSET_HH

#include <cstddef>

namespace Dune
{

  namespace Fem
  {

    // FunctionSet 
    // -----------

    /** \class FunctionSet 
     *
     *  \brief Global basis functions.
     *
     *  This class documents the function set interface.
     *
     *  \tparam  FunctionSpace  function space
     */
    template< class FunctionSpace >
    struct FunctionSet 
    {
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

      //! \brief return order of basis functions
      int order () const;

      //! \brief return number of basis functions
      std::size_t size () const;

      /**
       * \brief evalute each basis function
       *
       *  \param[in]  x        global coordinate
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
      template< class Functor >
      void evaluateEach ( const DomainType &x, Functor functor ) const;

      /**
       * \brief evalute jacobian of each basis function
       *
       *  \param[in]  x        global coordinate
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
      template< class Functor >
      void jacobianEach ( const DomainType &x, Functor functor ) const;

      /**
       * \brief evalute hessian of each basis function
       *
       *  \param[in]  x        global coordinate
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
      template< class Functor >
      void hessianEach ( const DomainType &x, Functor functor ) const;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_COMMON_FUNCTIONSET_HH
