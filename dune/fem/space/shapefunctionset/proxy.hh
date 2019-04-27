#ifndef DUNE_FEM_SHAPEFUNCTIONSET_PROXY_HH
#define DUNE_FEM_SHAPEFUNCTIONSET_PROXY_HH

// C++ includes
#include <cassert>
#include <cstddef>

/**
  @file
  @author Christoph Gersbacher
  @brief Provides a proxy class for pointers to a shape function set
*/


namespace Dune
{

  namespace Fem
  {

    // ShapeFunctionSetProxy
    // ---------------------

    /*
     * \brief A proxy object converting a pointer to a shape function set to a object
     *
     * \tparam  ShapeFunctionSet  An implementation of Dune::Fem::ShapeFunctionSet
     *
     * \note This class has an implicit constructor from a pointer to a shape function set.
     */
    template< class ShapeFunctionSet >
    class ShapeFunctionSetProxy
    {
      typedef ShapeFunctionSetProxy< ShapeFunctionSet > ThisType;

    public:
      typedef ShapeFunctionSet ImplementationType;

      typedef typename ImplementationType::FunctionSpaceType FunctionSpaceType;

      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      const ImplementationType &impl () const
      {
        assert( shapeFunctionSet_ );
        return *shapeFunctionSet_;
      }

      ShapeFunctionSetProxy ()
      : shapeFunctionSet_( nullptr )
      {}

      ShapeFunctionSetProxy ( const ShapeFunctionSet *shapeFunctionSet )
      : shapeFunctionSet_( shapeFunctionSet )
      {}

      int order () const { return impl().order(); }

      std::size_t size () const { return impl().size(); }

      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor functor ) const
      {
        impl().evaluateEach( x, functor );
      }

      template< class Point, class Functor >
      void jacobianEach ( const Point &x, Functor functor ) const
      {
        impl().jacobianEach( x, functor );
      }

      template< class Point, class Functor >
      void hessianEach ( const Point &x, Functor functor ) const
      {
        impl().hessianEach( x, functor );
      }

    private:
      const ShapeFunctionSet *shapeFunctionSet_;
    };

    // Same as above but instead of a raw pointer use a shared_ptr so that
    // deletion of the pointer is not an issue
    template< class ShapeFunctionSet >
    class SmartShapeFunctionSetProxy
    {
      typedef SmartShapeFunctionSetProxy< ShapeFunctionSet > ThisType;

    public:
      typedef ShapeFunctionSet ImplementationType;

      typedef typename ImplementationType::FunctionSpaceType FunctionSpaceType;

      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      const ImplementationType &impl () const
      {
        assert( shapeFunctionSet_ );
        return *shapeFunctionSet_;
      }

      SmartShapeFunctionSetProxy ()
      : shapeFunctionSet_( nullptr )
      {}

      template <class... Args>
      SmartShapeFunctionSetProxy (const Args&... args )
      : shapeFunctionSet_( std::make_shared<ImplementationType>(args...) )
      {}
      SmartShapeFunctionSetProxy (const std::shared_ptr<ImplementationType> &uptr )
      : shapeFunctionSet_( uptr )
      {}

      int order () const { return impl().order(); }

      std::size_t size () const { return impl().size(); }

      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor functor ) const
      {
        impl().evaluateEach( x, functor );
      }

      template< class Point, class Functor >
      void jacobianEach ( const Point &x, Functor functor ) const
      {
        impl().jacobianEach( x, functor );
      }

      template< class Point, class Functor >
      void hessianEach ( const Point &x, Functor functor ) const
      {
        impl().hessianEach( x, functor );
      }

    private:
      std::shared_ptr<ShapeFunctionSet> shapeFunctionSet_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SHAPEFUNCTIONSET_PROXY_HH
