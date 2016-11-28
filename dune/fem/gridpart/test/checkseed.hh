#ifndef DUNE_FEM_GRIDPART_TEST_CHECKSEED_HH
#define DUNE_FEM_GRIDPART_TEST_CHECKSEED_HH

//- C++ includes
#include <limits>
#include <type_traits>

//- dune-common includes
#include <dune/fem/common/forloop.hh>
#include <dune/common/exceptions.hh>

//- dune-fem includes
#include <dune/fem/gridpart/common/capabilities.hh>


namespace Dune
{

  namespace Fem
  {

    /** \brief Iterate over all entities of all codimensions available in the grid part,
     *         and perform a number of tests on each entity.
     */

    template< class GridPartType >
    class CheckEntitySeed
    {
      template< class GridPart, bool >
      struct If
      {
        template< int codim >
        static void check ( const GridPart & )
        { }
      };

      template< class GridPart >
      struct If< GridPart, true >
      {
        template< const int codim >
        static void check ( const GridPart &gridPart )
        {
          typedef typename GridPart::template Codim< codim >::EntityType EntityType;
          typedef typename GridPart::template Codim< codim >::EntitySeedType EntitySeedType;
          static_assert( std::is_same< EntitySeedType, typename EntityType::EntitySeed >::value, "Types of EntitySeed do not coincide" );

          for( auto it = gridPart.template begin< codim >(); it != gridPart.template end< codim >(); ++it )
          {
            auto entity = *it;
            auto seed = entity.seed();
            if( gridPart.entity( seed ) != entity )
              DUNE_THROW( Dune::Exception, "Could not recover entity from seed" );
          }
        }
      };

      template< int codim >
      struct CheckCodim
      {
        template< class GridPart >
        static void apply ( const GridPart &gridPart )
        {
          // check whether grid part has entity of given codimension
          constexpr bool hasEntity = Dune::Fem::GridPartCapabilities::hasEntity< GridPart, codim >::v;
          If< GridPart, hasEntity >::template check< codim >( gridPart );
        }
      };

    public:
      static const int dimension = GridPartType::dimension;

      static void check ( const GridPartType &gridPart )
      {
        // generic for loop over all codimensions
        Dune::Fem::ForLoop< CheckCodim, 0, dimension >::apply( gridPart );
      }
    };

  } // end namespace Fem

} // end namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_TEST_CHECKSEED_HH
