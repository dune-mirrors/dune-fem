#ifndef DUNE_FEM_SOLVER_COMMUNICATION_OWNEROVERLAPCOPY_HH
#define DUNE_FEM_SOLVER_COMMUNICATION_OWNEROVERLAPCOPY_HH

#include <cassert>

#include <map>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/parallel/remoteindices.hh>

#include <dune/grid/common/datahandleif.hh>

#include <dune/istl/owneroverlapcopy.hh>

#include <dune/fem/space/mapper/parallel.hh>


namespace Dune
{

  namespace Fem
  {

    namespace ISTL
    {

      // OwnerOverlapCopyCommunication
      // -----------------------------

      template< class DiscreteFunctionSpace >
      using OwnerOverlapCopyCommunication = Dune::OwnerOverlapCopyCommunication< std::size_t, typename DiscreteFunctionSpace::BlockMapperType::GlobalKeyType >;



      // BuildRemoteIndicesDataHandle
      // ----------------------------

      template< class Mapper, class GlobalLookup >
      struct BuildRemoteIndicesDataHandle
        : public Dune::CommDataHandleIF< BuildRemoteIndicesDataHandle< Mapper, GlobalLookup >, int >
      {
        typedef typename GlobalLookup::GlobalIndex GlobalIndexType;
        typedef typename GlobalLookup::LocalIndex::Attribute AttributeType;

        BuildRemoteIndicesDataHandle ( int rank, const Mapper &mapper, const GlobalLookup &globalLookup )
          : rank_( rank ), mapper_( mapper ), globalLookup_( globalLookup )
        {}

        bool contains ( int dim, int codim ) const { return mapper_.contains( codim ); }
        bool fixedSize( int dim, int codim ) const { return true; }

        template< class Buffer, class Entity >
        void gather ( Buffer &buffer, const Entity &entity ) const
        {
          buffer.write( rank_ );
          int attribute = -1;
          mapper_.mapEachEntityDof( entity, [ this, &attribute ] ( int, auto index ) {
              auto *pair = globalLookup_.pair( index );
              assert( pair && ((attribute == -1) || (attribute == pair->local().attribute())) );
              attribute = pair->local().attribute();
            } );
          buffer.write( static_cast< int >( attribute ) );
        }

        template< class Buffer, class Entity >
        void scatter ( Buffer &buffer, const Entity &entity, std::size_t n )
        {
          int rank, attribute;
          buffer.read( rank );
          buffer.read( attribute );
          assert( (attribute != -1) || (mapper_.numEntityDofs( entity ) == 0) );
          mapper_.mapEachEntityDof( entity, [ this, rank, attribute ] ( int, auto index ) {
              auto *pair = globalLookup_.pair( index );
              assert( pair );
              remotes[ rank ].emplace_back( static_cast< AttributeType >( attribute ), pair );
            } );
        }

        template< class Entity >
        std::size_t size ( const Entity &entity ) const
        {
          return 2;
        }

        std::map< int, std::vector< Dune::RemoteIndex< GlobalIndexType, AttributeType > > > remotes;

      private:
        int rank_;
        const Mapper &mapper_;
        const GlobalLookup &globalLookup_;
      };



      template< class DiscreteFunctionSpace, class GlobalId, class LocalId >
      void buildCommunication ( const DiscreteFunctionSpace &dfSpace,
                                Dune::SolverCategory::Category solverCategory,
                                std::shared_ptr< Dune::OwnerOverlapCopyCommunication< GlobalId, LocalId > > &communication )
      {
        typedef typename DiscreteFunctionSpace::GridPartType GridPartType;
        typedef typename DiscreteFunctionSpace::BlockMapperType LocalMapperType;

        typedef typename Dune::OwnerOverlapCopyCommunication< GlobalId, LocalId >::GlobalLookupIndexSet GlobalLookupType;

        typedef typename GlobalLookupType::LocalIndex LocalIndexType;

        communication.reset( new Dune::OwnerOverlapCopyCommunication< GlobalId, LocalId >( solverCategory ) );

        const GridPartType &gridPart = dfSpace.gridPart();
        LocalMapperType &localMapper = dfSpace.blockMapper();

        // create global index mapping
        Dune::Fem::ParallelDofMapper< GridPartType, LocalMapperType, GlobalId > globalMapper( gridPart, localMapper );

        // construct local attributes
        std::vector< typename LocalIndexType::Attribute > attribute( localMapper.size(), Dune::OwnerOverlapCopyAttributeSet::owner );
        for( const auto &auxiliary : dfSpace.auxiliaryDofs() )
          attribute[ auxiliary ] = Dune::OwnerOverlapCopyAttributeSet::copy;

        // build parallel index set
        communication->indexSet().beginResize();
        for( LocalId i = 0, size = localMapper.size(); i < size; ++i )
          communication->indexSet().add( globalMapper.mapping()[ i ], LocalIndexType( i, attribute[ i ] ) );
        communication->indexSet().endResize();

        // build remote indices
        communication->buildGlobalLookup();
        BuildRemoteIndicesDataHandle< LocalMapperType, GlobalLookupType > buildRemoteIndicesDataHandle( gridPart.comm().rank(), localMapper, communication->globalLookup() );
        gridPart.communicate( buildRemoteIndicesDataHandle, Dune::All_All_Interface, Dune::ForwardCommunication );
        communication->freeGlobalLookup();

        communication->remoteIndices().setIndexSets( communication->indexSet(), communication->indexSet(), communication->communicator() );
        if( !buildRemoteIndicesDataHandle.remotes.empty() )
        {
          for( auto &remote : buildRemoteIndicesDataHandle.remotes )
          {
            std::sort( remote.second.begin(), remote.second.end(), [] ( const auto &a, const auto &b ) { return (a.localIndexPair().global() < b.localIndexPair().global()); } );
            auto modifier = communication->remoteIndices().template getModifier< false, true >( remote.first );
            for( const auto &remoteIndex : remote.second )
              modifier.insert( remoteIndex );
          }
        }
        else
          communication->remoteIndices().template getModifier< false, true >( 0 );
      }



      // SupportsAMG for OwnerOverlapCopyCommunication
      // ---------------------------------------------

      template< class T >
      struct SupportsAMG;

      template< class GlobalId, class LocalId >
      struct SupportsAMG< Dune::OwnerOverlapCopyCommunication< GlobalId, LocalId > >
        : public std::true_type
      {};

    } // namespace ISTL

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SOLVER_COMMUNICATION_OWNEROVERLAPCOPY_HH
