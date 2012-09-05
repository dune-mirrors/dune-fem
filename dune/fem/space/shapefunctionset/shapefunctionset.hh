#ifndef DUNE_FEM_SHAPEFUNCTIONSET_SHAPEFUNCTIONSET_HH
#define DUNE_FEM_SHAPEFUNCTIONSET_SHAPEFUNCTIONSET_HH

//- dune-geometry includes
#include <dune/geometry/type.hh>

//- dune-fem includes
#include <dune/fem/misc/bartonnackmaninterface.hh>


namespace Dune
{

  namespace Fem
  {

    // SimpleShapeFunction
    // -------------------

    template< class Traits >
    class ShapeFunctionSet
    {
    public:
      //! \brief type of real implementation
      typedef typename Traits::ShapeFunctionSetType ImplementationType;
      //! \brief return reference to implementation
      ImplementationType &impl () { return static_cast< ImplementationType & >( *this ); }
      //! \brief return const reference to implementation
      const ImplementationType &impl () const { return static_cast< const ImplementationType & >( *this ); }
      
      //! \brief function space type
      typedef typename Traits::FunctionSpaceType FunctionSpaceType;
      //! \brief range type
      typedef typename FunctionSpaceType::RangeType RangeType;
      //! \brief jacobian range type
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      //! \brief hessian range type
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      //! \brief return geometry type
      GeometryType type () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( impl().type() );
        return impl().type();
      }

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
       *    void operator() ( const int shapeFunction, RangeType &value );
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
       *    void operator() ( const int shapeFunction, JacobianRangeType &jacobian );
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
       *    void operator() ( const int shapeFunction, HessianRangeType &hessian );
       *  };
       *  \endcode
       */
      template< class Point, class Functor > 
      void hessianEach ( const Point &x, Functor functor ) const
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( impl().hessianEach( x, functor ) );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SHAPEFUNCTIONSET_SHAPEFUNCTIONSET_HH
