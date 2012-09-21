#ifndef DUNE_FEM_SHAPEFUNCTIONSET_PROXY_HH
#define DUNE_FEM_SHAPEFUNCTIONSET_PROXY_HH

// C++ includes
#include <cassert>

// dune-common includes
#include <dune/common/nullptr.hh>

// dune-geometry includes
#include <dune/geometry/type.hh>

// dune-fem includes
#include <dune/fem/space/shapefunctionset/shapefunctionset.hh>

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
     * \tparam  ShapeFunctionSetImpl  An implementation of Dune::Fem::ShapeFunctionSet
     *
     * \note This class has an implicit constructor from a pointer to a shape function set.
     */
    template< class ShapeFunctionSetImpl >
    class ShapeFunctionSetProxy
    : public ShapeFunctionSet< typename ShapeFunctionSetImpl::FunctionSpace, ShapeFunctionSetProxy< ShapeFunctionSetImpl > >
    {
      typedef ShapeFunctionSetProxy< ShapeFunctionSetImpl > ThisType;
      typedef ShapeFunctionSet< typename ShapeFunctionSetImpl::FunctionSpace, ThisType > BaseType;

    protected:
      typedef ShapeFunctionSetImpl ImplementationType;

    public:
      ShapeFunctionSetProxy ()
      : implementation_( nullptr )
      {}

      ShapeFunctionSetProxy ( const ImplementationType *implementation )
      : implementation_( implementation )
      {}

      GeometryType type () const { return implementation().type(); }

      std::size_t size () const { return implementation().size(); }

      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor functor ) const
      {
        implementation().evaluateEach( x, functor );
      }

      template< class Point, class Functor > 
      void jacobianEach ( const Point &x, Functor functor ) const
      {
        implementation().jacobianEach( x, functor );
      }

      template< class Point, class Functor > 
      void hessianEach ( const Point &x, Functor functor ) const
      {
        implementation().hessianEach( x, functor );
      }

    protected:
      const ImplementationType &implementation () const
      {
        assert( implementation_ );
        return *implementation_;
      }

    private:
      const ImplementationType *implementation_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SHAPEFUNCTIONSET_PROXY_HH
