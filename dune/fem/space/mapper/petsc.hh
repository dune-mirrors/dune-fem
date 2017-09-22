#ifndef DUNE_FEM_SPACE_MAPPER_PETSC_HH
#define DUNE_FEM_SPACE_MAPPER_PETSC_HH

#include <cassert>

#include <memory>
#include <tuple>
#include <utility>

#include <dune/fem/space/mapper/ghost.hh>
#include <dune/fem/space/mapper/parallel.hh>
#include <dune/fem/storage/singletonlist.hh>

namespace Dune
{

  namespace Fem
  {

    // PetscMappers
    // ------------

#if HAVE_PETSC
    template< class DiscreteFunctionSpace >
    class PetscMappers
    {
      typedef PetscMappers< DiscreteFunctionSpace > PetscMappersType;

    public:
      typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

      typedef typename DiscreteFunctionSpaceType::BlockMapperType BlockMapperType;
      typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

      typedef GhostDofMapper< GridPartType, BlockMapperType, PetscInt > GhostMapperType;
      typedef ParallelDofMapper< GridPartType, GhostMapperType, PetscInt > ParallelMapperType;

    private:
      template< class Mapper >
      struct MapperFactory
      {
        static std::pair< Mapper, int > *createObject ( std::pair< GridPartType *, typename Mapper::BaseMapperType * > key, int sequence )
        {
          return new std::pair< Mapper, int >( std::piecewise_construct, std::tie( *key.first, *key.second ), std::make_tuple( sequence ) );
        }

        static void deleteObject ( std::pair< Mapper, int > *object ) { delete object; }
      };

      template< class Mapper >
      using MapperProvider = SingletonList< std::pair< GridPartType *, typename Mapper::BaseMapperType * >, std::pair< Mapper, int >, MapperFactory< Mapper > >;

      typedef MapperProvider< GhostMapperType > GhostMapperProviderType;
      typedef MapperProvider< ParallelMapperType > ParallelMapperProviderType;

    public:
      //! constructor creating petsc mapper
      PetscMappers ( const DiscreteFunctionSpaceType &space )
        : space_( space )
      {
        const int sequence = space_.sequence();

        ghostMapper_.reset( &GhostMapperProviderType::getObject( std::make_pair( &space_.gridPart(), &space_.blockMapper() ), sequence ) );
        update( *ghostMapper_, sequence );

        parallelMapper_.reset( &ParallelMapperProviderType::getObject( std::make_pair( &space_.gridPart(), &ghostMapper_->first ), sequence ) );
        update( *parallelMapper_, sequence );
      }

      //! copy constructor obtaining pointers for mapper objects
      PetscMappers ( const PetscMappers &other )
        : space_( other.space_ )
      {
        const int sequence = space_.sequence();

        ghostMapper_.reset( &GhostMapperProviderType::getObject( std::make_pair( &space_.gridPart(), &space_.blockMapper() ), sequence ) );
        parallelMapper_.reset( &ParallelMapperProviderType::getObject( std::make_pair( &space_.gridPart(), &ghostMapper_->first ), sequence ) );
      }

      const DiscreteFunctionSpaceType &space () const { return space_; }

      const GhostMapperType &ghostMapper () const { assert( ghostMapper_ ); return ghostMapper_->first; }
      const ParallelMapperType &parallelMapper () const { assert( parallelMapper_ ); return parallelMapper_->first; }

      PetscInt ghostIndex ( const typename BlockMapperType::GlobalKeyType &index ) const
      {
        assert( static_cast< std::size_t >( index ) < size() );
        return ghostMapper().mapping()[ index ];
      }

      PetscInt parallelIndex ( const typename BlockMapperType::GlobalKeyType &index ) const
      {
        return parallelMapper().mapping()[ ghostIndex( index ) ];
      }

      std::size_t size () const { return ghostMapper().mapping().size(); }

      void update ()
      {
        const int sequence = space().sequence();

        assert( ghostMapper_ );
        update( *ghostMapper_, sequence );

        assert( parallelMapper_ );
        update( *parallelMapper_, sequence );
      }

    private:
      template< class Mapper >
      static void update ( std::pair< Mapper, int > &mapper, int sequence )
      {
        if( mapper.second != sequence )
        {
          mapper.first.update();
          mapper.second = sequence;
        }
      }

      const DiscreteFunctionSpaceType &space_;
      std::unique_ptr< std::pair< GhostMapperType, int >, typename GhostMapperProviderType::Deleter > ghostMapper_;
      std::unique_ptr< std::pair< ParallelMapperType, int >, typename ParallelMapperProviderType::Deleter > parallelMapper_;
    };
#endif // #if HAVE_PETSC

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_MAPPER_PETSC_HH
