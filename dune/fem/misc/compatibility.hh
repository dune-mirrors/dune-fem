#ifndef DUNE_FEM_MISC_COMPATIBILITY_HH
#define DUNE_FEM_MISC_COMPATIBILITY_HH

#warning DO NOT INCLUDE 'COMPATIBILITY': the function make_entity is not needed anymore

#include <utility>

#include <dune/grid/common/entity.hh>
#include <dune/fem/version.hh>

namespace Dune
{

  namespace Fem
  {

    // make_entity
    // -----------

    template< int codim, int dim, class Grid, template< int, int, class > class Implementation >
    [[deprecate("make_entity")]]
    typename Dune::Entity< codim, dim, Grid, Implementation >
    make_entity ( Dune::Entity< codim, dim, Grid, Implementation > entity )
    {
      return std::move( entity );
    }

  } // namespace Fem

} // end namespace Dune

#endif // #ifndef DUNE_FEM_MISC_COMPATIBILITY_HH
