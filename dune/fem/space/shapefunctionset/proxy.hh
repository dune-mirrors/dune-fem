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
     * \tparam  ShapeFunctionSet  An implementation of Dune::Fem::ShapeFunctionSet
     *
     * \note This class has an implicit constructor from a pointer to a shape function set.
     */
    template< class ShapeFunctionSet >
    class ShapeFunctionSetProxy
    : public Dune::Fem::ShapeFunctionSet< typename ShapeFunctionSet::FunctionSpaceType, ShapeFunctionSetProxy< ShapeFunctionSet > >
    {
      typedef ShapeFunctionSetProxy< ShapeFunctionSet > ThisType;
      typedef Dune::Fem::ShapeFunctionSet< typename ShapeFunctionSet::FunctionSpaceType, ThisType > BaseType;

    public:
      ShapeFunctionSetProxy ()
      : shapeFunctionSet_( nullptr )
      {}

      ShapeFunctionSetProxy ( const ShapeFunctionSet *shapeFunctionSet )
      : shapeFunctionSet_( shapeFunctionSet )
      {}

      GeometryType type () const { return shapeFunctionSet().type(); }

      std::size_t size () const { return shapeFunctionSet().size(); }

      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor functor ) const
      {
        shapeFunctionSet().evaluateEach( x, functor );
      }

      template< class Point, class Functor > 
      void jacobianEach ( const Point &x, Functor functor ) const
      {
        shapeFunctionSet().jacobianEach( x, functor );
      }

      template< class Point, class Functor > 
      void hessianEach ( const Point &x, Functor functor ) const
      {
        shapeFunctionSet().hessianEach( x, functor );
      }

    private:
      const ShapeFunctionSet &shapeFunctionSet () const
      {
        assert( shapeFunctionSet_ );
        return *shapeFunctionSet_;
      }

      const ShapeFunctionSet *shapeFunctionSet_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SHAPEFUNCTIONSET_PROXY_HH
