#ifndef DUNE_FEM_SHAPEFUNCTIONSET_SHAPEFUNCTIONSET_HH
#define DUNE_FEM_SHAPEFUNCTIONSET_SHAPEFUNCTIONSET_HH

//- C++ includes
#include <cstddef>

//- dune-common includes
#include <dune/common/bartonnackmanifcheck.hh>

//- dune-geometry includes
#include <dune/geometry/type.hh>


namespace Dune
{

  namespace Fem
  {

    // ShapeFunctionSet
    // ----------------

    /**
     * \brief interface class for shape function sets
     *
     * \tparam  FunctionSpace   Function space
     * \tparam  Implementation  Implementation of this interface
     *
     * \note FunctionSpace::dimDomain has to be equal to type().dim().
     */
    template< class FunctionSpace, class Implementation >
    class ShapeFunctionSet
    {
      // this type
      typedef ShapeFunctionSet< FunctionSpace, Implementation > ThisType;

    public:
      //! \brief function space type
      typedef FunctionSpace FunctionSpaceType;
      //! \brief type of real implementation
      typedef Implementation ImplementationType;

      //! \brief return reference to implementation
      ImplementationType &impl () { return static_cast< ImplementationType & >( *this ); }
      //! \brief return const reference to implementation
      const ImplementationType &impl () const { return static_cast< const ImplementationType & >( *this ); }
      
      //! \brief domain type
      typedef typename FunctionSpaceType::DomainType DomainType;
      //! \brief range type
      typedef typename FunctionSpaceType::RangeType RangeType;
      //! \brief jacobian range type
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      //! \brief hessian range type
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      //! \brief return number of shape functions
      std::size_t size () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( impl().size() );
        return impl().size();
      }

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
      void evaluateEach ( const Point &x, Functor functor ) const
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( impl().evaluateEach( x, functor ) );
      }

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
      void jacobianEach ( const Point &x, Functor functor ) const
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( impl().jacobianEach( x, functor ) );
      }

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
      void hessianEach ( const Point &x, Functor functor ) const
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( impl().hessianEach( x, functor ) );
      }

    private:
      // forbid copy constructor
      ShapeFunctionSet ( const ThisType & );
      // forbid assignment operator 
      ThisType &operator= ( const ThisType & );
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SHAPEFUNCTIONSET_SHAPEFUNCTIONSET_HH
