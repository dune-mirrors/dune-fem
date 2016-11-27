#ifndef DUNE_FEM_GRIDPART_TEST_CHECKGEOMETRY_HH
#define DUNE_FEM_GRIDPART_TEST_CHECKGEOMETRY_HH

//- dune-common includes
#include <dune/fem/common/forloop.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/stdstreams.hh>

//- dune-grid includes
#include <dune/grid/test/checkgeometry.hh>

//- dune-fem includes
#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/gridpart/test/failure.hh>


/** \brief Check geometries for all entities in the
 *         grid part. The tests are implemented in
 *         \code
<dune/grid/test/checkgeometry.cc>
 *         \endcode
 *         Thus, failure handlers passed to this test
 *         will be ignored.
 *
 */

namespace Dune
{
  namespace Fem
  {

    template< class GridPartType, class FailureHandler >
    class CheckGeometry
    {
      template< class GridPart, bool >
      struct If
      {
        template< int codim, class Entity  >
        static void check ( const Entity&, FailureHandler & )
        { }
      };

      template< class GridPart >
      struct If< GridPart, true >
      {
        template< int codim, class Entity  >
        static void check ( const Entity &entity, FailureHandler & )
        {
          const int subEntities = entity.subEntities( codim );
          for ( int i = 0; i < subEntities; ++i )
          {
            auto subEntity = entity.template subEntity< codim >( i );

            auto geometry = subEntity.geometry();
            if( subEntity.type() != geometry.type() )
              std::cerr << "Error: Entity and geometry report different geometry types on codimension " << codim << "." << std::endl;
            Dune::checkGeometry( geometry );
          }
        }
      };

      template< int codim >
      struct CheckSubEntityGeometry
      {
        template< class Entity >
        static void apply ( const Entity &entity, FailureHandler &failureHandler )
        {
          // check whether grid part has entity of given codimension
          constexpr bool hasEntity = Dune::Fem::GridPartCapabilities::hasEntity< GridPartType, codim >::v;
          If< GridPartType, hasEntity >::template check< codim, Entity >( entity, failureHandler );
        }
      };

    public:
      static const int dimension = GridPartType::dimension;

      static void check ( const GridPartType &gridPart, FailureHandler &failureHandler )
      {
        // tests will be forwarded to dune-grid, where failures are not supported
        std::cout << "Note: FailureHandler will be ignored" << std::endl;

        for( const auto& entity : elements( gridPart ) )
          Fem::ForLoop< CheckSubEntityGeometry, 0, dimension >::apply( entity, failureHandler );
      }
    };

  } // end namespace Fem

} // end namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_TEST_CHECKGEOMETRY_HH
