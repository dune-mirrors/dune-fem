#ifndef DUNE_FEM_SPACE_BREZZIDOUGLASMARINI_LOCALFINITEELEMEMT_HH
#define DUNE_FEM_SPACE_BREZZIDOUGLASMARINI_LOCALFINITEELEMEMT_HH

#include <dune/common/fvector.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/type.hh>

#if HAVE_DUNE_LOCALFUNCTIONS

#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini1cube2d.hh>
#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini1cube3d.hh>
#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini1simplex2d.hh>
#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini2cube2d.hh>
#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini2simplex2d.hh>

namespace Dune
{

  namespace Fem
  {

    // BDMLocalFiniteElement
    // ---------------------

    template< unsigned int id, class DomainField, class RangeField, int dimension, int order >
    struct BDMLocalFiniteElement
    {
      static_assert( AlwaysFalse< DomainField >::value, "BDMLocalFiniteElement not implemented for your choice." );
    };

    // The following local finite elements are implemented

    // 2d, Cube, first order
    template< class D >
    struct BDMLocalFiniteElement< Dune::Impl::CubeTopology< 2 >::type::id, D, D, 2, 1 >
      : public BDM1Cube2DLocalFiniteElement< D, D >
    {
      static const int numOrientations = 16;
      template< class ... Args >
      BDMLocalFiniteElement ( Args && ... args )
        : BDM1Cube2DLocalFiniteElement< D, D >( std::forward< Args >( args ) ... ) {}
    };

    // 3d, Cube, first order
    template< class D >
    struct BDMLocalFiniteElement< Dune::Impl::CubeTopology< 3 >::type::id, D, D, 3, 1 >
      : public BDM1Cube3DLocalFiniteElement< D, D >
    {
      static const int numOrientations = 64;
      template< class ... Args >
      BDMLocalFiniteElement ( Args && ... args )
        : BDM1Cube3DLocalFiniteElement< D, D >( std::forward< Args >( args ) ... ) {}
    };

    // 2d, Cube, second order
    template< class D >
    struct BDMLocalFiniteElement< Dune::Impl::CubeTopology< 2 >::type::id, D, D, 2, 2 >
      : public BDM2Cube2DLocalFiniteElement< D, D >
    {
      static const int numOrientations = 16;
      template< class ... Args >
      BDMLocalFiniteElement ( Args && ... args )
        : BDM2Cube2DLocalFiniteElement< D, D >( std::forward< Args >( args ) ... ) {}
    };


    // 2d, simplex, first order
    template< class D >
    struct BDMLocalFiniteElement< Dune::Impl::SimplexTopology< 2 >::type::id, D, D, 2, 1 >
      : public BDM1Simplex2DLocalFiniteElement< D, D >
    {
      static const int numOrientations = 8;
      template< class ... Args >
      BDMLocalFiniteElement ( Args && ... args )
        : BDM1Simplex2DLocalFiniteElement< D, D >( std::forward< Args >( args ) ... ) {}
    };

    // 2d, simplex, second order
    template< class D >
    struct BDMLocalFiniteElement< Dune::Impl::SimplexTopology< 2 >::type::id, D, D, 2, 2 >
      : public BDM2Simplex2DLocalFiniteElement< D, D >
    {
      static const int numOrientations = 8;
      template< class ... Args >
      BDMLocalFiniteElement ( Args && ... args )
        : BDM2Simplex2DLocalFiniteElement< D, D >( std::forward< Args >( args ) ... ) {}
    };

  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_DUNE_LOCALFUNCTIONS

#endif // #ifndef DUNE_FEM_SPACE_BREZZIDOUGLASMARINI_LOCALFINITEELEMEMT_HH
