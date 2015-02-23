#ifndef DUNE_FEM_GRIDPART_TEST_CHECKSEED_HH
#define DUNE_FEM_GRIDPART_TEST_CHECKSEED_HH

//- C++ includes
#include <limits>

//- dune-common includes
#include <dune/common/forloop.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/stdstreams.hh>

//- dune-fem includes
#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/gridpart/test/failure.hh>


namespace Dune
{
  namespace Fem
  {

    /** \brief check conversion from entity to entity pointer
     *         and entity seed, and from entity seed to entity pointer
     */
    template< class EntityType, class GridPartType, class FailureHandler >
    int conversion ( const EntityType &entity, const GridPartType &gridPart, FailureHandler &failureHandler )
    {
      const int codimension = EntityType::codimension;

      typedef typename GridPartType::template Codim< codimension >::EntityPointerType EntityPointerType;
      typedef typename GridPartType::template Codim< codimension >::EntitySeedType EntitySeedType;

      EntityPointerType ep1 ( entity );

      EntitySeedType seed = entity.seed();
      EntityPointerType ep2 = gridPart.entityPointer( seed );


      if( ep1 != ep2 )
      {
        struct ConversionFailure
        : public Failure
        {
          virtual void writeTo ( std::ostream &out ) const
          {
            out <<  __FILE__
                << ":" << __LINE__ << ": Failure :"
                << "Conversion failed";
          }
        };
        Failure *failure = new ConversionFailure();
        failureHandler( *failure );
      }
      return 0;
    }



    /** \brief Iterate over all entities of all codimensions available in the grid part,
     *         and perform a number of tests on each entity.
     */

    template< class GridPartType, class FailureHandler >
    class CheckEntitySeed
    {
      template< class GridPart, bool >
      struct If
      {
        template< int codim >
        static void check ( const GridPart &, FailureHandler & )
        { }
      };

      template< class GridPart >
      struct If< GridPart, true >
      {
        template< const int codim >
        static void check ( const GridPart &gridPart, FailureHandler &failureHandler )
        {
          typedef typename GridPart::template Codim< codim >::IteratorType IteratorType;
          typedef typename GridPart::template Codim< codim >::EntityType EntityType;

          const IteratorType end = gridPart.template end< codim >();
          for( IteratorType it = gridPart.template begin< codim >(); it != end; ++it )
          {
            const EntityType &entity = *it;
            conversion( entity, gridPart, failureHandler );
          }
        }
      };

      template< int codim >
      struct CheckCodim
      {
        template< class GridPart >
        static void apply ( const GridPart &gridPart, FailureHandler &failureHandler )
        {
          // check whether grid part has entity of given codimension
          const bool hasEntity = Dune::Fem::GridPartCapabilities::hasEntity< GridPart, codim >::v;
          If< GridPart, hasEntity >::template check< codim >( gridPart, failureHandler );
        }
      };

    public:
      static const int dimension = GridPartType::dimension;

      static void check ( const GridPartType &gridPart, FailureHandler &failureHandler )
      {
        // generic for loop over all codimensions
        Dune::ForLoop< CheckCodim, 0, dimension >::apply( gridPart, failureHandler );
      }
    };

  } // end namespace Fem

} // end namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_TEST_CHECKSEED_HH
