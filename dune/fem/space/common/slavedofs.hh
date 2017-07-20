#ifndef DUNE_FEM_SLAVEDOFS_HH
#define DUNE_FEM_SLAVEDOFS_HH

#include <cassert>

#include <iostream>
#include <memory>
#include <set>
#include <map>
#include <limits>
#include <algorithm>

#include <dune/common/exceptions.hh>
#include <dune/common/genericiterator.hh>
#include <dune/common/ftraits.hh>
#include <dune/common/typetraits.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/datahandleif.hh>

#include <dune/fem/storage/singletonlist.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/space/common/commindexmap.hh>

namespace Dune
{

  namespace Fem
  {

  /** @addtogroup Communication Communication
      @{
  **/

    template< class Space, class Mapper >
    class SlaveDofs
    {
      typedef SlaveDofs< Space, Mapper > ThisType;

    public:
      class SingletonKey;

    private:
      class LinkBuilder;

    public:
      //! type of discrete function space
      typedef Space SpaceType;

      //! for convenience
      typedef SpaceType DiscreteFunctionSpaceType ;

      //! type of used mapper
      typedef Mapper MapperType;

    protected:
      typedef Fem :: CommunicationIndexMap IndexMapType;

      const MapperType &mapper_;

      const int myRank_;
      const int mySize_;

      // type of communication indices
      IndexMapType slaves_;
      std::set<int> slaveSet_;

      //! know grid sequence number
      int sequence_;

    public:
      struct ConstIterator
      {
        ConstIterator () = default;
        ConstIterator ( const IndexMapType &slaves, int index ) : slaves_( &slaves ), index_( index ) {}

        const int &operator* () const { return (*slaves_)[ index_ ]; }
        const int *operator-> () const { return &(*slaves_)[ index_ ]; }

        bool operator== ( const ConstIterator &other ) const { return (index_ == other.index_); }
        bool operator!= ( const ConstIterator &other ) const { return (index_ != other.index_); }

        ConstIterator &operator++ () { ++index_; return *this; }
        ConstIterator operator++ ( int ) { ConstIterator copy( *this ); ++(*this); return copy; }

      private:
        const IndexMapType *slaves_ = nullptr;
        int index_ = 0;
      };

      //! constructor taking space
      SlaveDofs ( const SingletonKey &key )
      : mapper_( key.mapper() ),
        myRank_( key.gridPart().comm().rank() ),
        mySize_( key.gridPart().comm().size() ),
        slaves_(),
        sequence_( -1 )
      {}

      SlaveDofs ( const SlaveDofs& ) = delete;

      //! return dof number of salve with index
      int operator [] ( const int index ) const
      {
        return slaves_[ index ];
      }

      //! return number of slave dofs
      int size () const
      {
        return slaves_.size();
      }

      ConstIterator begin () const { return ConstIterator( slaves_, 0 ); }
      ConstIterator end () const { assert( size() > 0 ); return ConstIterator( slaves_, size()-1 ); }

      //! return true if index is contained, meaning is a slave dof
      bool isSlave( const int index ) const
      {
        typedef GenericIterator<const IndexMapType, const int> IteratorType;

        return std::binary_search(IteratorType(slaves_, 0),
                                  IteratorType(slaves_, size()-1),
                                  index);
      }

      //! insert index
      void insert( const int index )
      {
        slaveSet_.insert( index );
      }

      //! initialize
      void initialize ()
      {
        sequence_ = -1;
        slaveSet_.clear();
        slaves_.clear();
      }

      //! finalize
      void finalize (const SpaceType &space)
      {
        // insert slaves
        slaves_.set( slaveSet_ );

        // remove memory
        slaveSet_.clear();

        // store actual sequence number
        sequence_ = space.sequence();
      }

      //! check if grid has changed and rebuild cache if necessary
      void rebuild ()
      {
        static_assert( AlwaysFalse< ThisType >::value, "don't call rebuild() on slavedof class - use the rebuild method taking a space or use the method slavedofs on the space directly" );
      }

      void rebuild (const SpaceType &space)
      {
        // check whether grid has changed.
        if( sequence_ != space.sequence() )
        {
          initialize();
          buildMaps(space);
          finalize(space);
        }
      }

    protected:
      void buildMaps (const SpaceType &space)
      {
        // build linkage and index maps
        if( !space.continuous() )
          buildDiscontinuousMaps(space);
        else
          buildCommunicatedMaps(space);
      }

      void buildDiscontinuousMaps (const SpaceType &space)
      {
        // for discontinuous spaces we don't have to communicate
        const auto idxpitype = SpaceType::GridPartType :: indexSetPartitionType;
        const auto endit = space.gridPart().template end< 0, idxpitype >();
        for( auto it = space.gridPart().template begin< 0, idxpitype >(); it != endit; ++it )
        {
          const auto& entity = *it;
          if( entity.partitionType() != Dune::InteriorEntity )
            mapper_.mapEachEntityDof( entity, [this]( const int, const auto& value ){this->insert( value );} );
        }
        // insert overall size at the end
        insert( mapper_.size() );
      }

      void buildCommunicatedMaps (const SpaceType &space)
      {
        // we have to skip communication when parallel program is executed only on one processor
        // otherwise YaspGrid and Lagrange polorder=2 fails :(
        if( space.gridPart().comm().size() > 1 )
        {
          try
          {
            LinkBuilder handle( *this, space , mapper_ );
            space.gridPart().communicate( handle, SpaceType::GridPartType::indexSetInterfaceType, ForwardCommunication );
          }
          catch( const Exception &e )
          {
            std::cerr << e << std::endl << "Exception thrown in: " << __FILE__ << " line:" << __LINE__ << std::endl;
            abort();
          }
        }
        // insert overall size at the end
        insert( mapper_.size() );
      }
    };



    //! Key for CommManager singleton list
    template< class Space, class Mapper >
    class SlaveDofs< Space, Mapper > :: SingletonKey
    {
    public:
      typedef Space SpaceType;
      typedef typename SpaceType::GridPartType GridPartType;
      typedef Mapper MapperType;

    protected:
      const GridPartType &gridPart_;
      const MapperType *const mapper_;

    public:
      //! constructor taking space
      SingletonKey ( const typename SpaceType::GridPartType &gridPart, const MapperType &mapper )
      : gridPart_( gridPart ), mapper_( &mapper )
      {}

      //! copy constructor
      SingletonKey ( const SingletonKey &other )
      : gridPart_( other.gridPart_ ), mapper_( other.mapper_ )
      {}

      //! returns true if indexSet pointer and numDofs are equal
      bool operator== ( const SingletonKey &other ) const
      {
        return (&gridPart_ == &other.gridPart_) && (mapper_ == other.mapper_);
      }

      //! return reference to index set
      const GridPartType &gridPart () const
      {
        return gridPart_;
      }

      //! return reference to index set
      const MapperType &mapper () const
      {
        return *mapper_;
      }
    };



    template< class Space, class Mapper >
    class SlaveDofs< Space,Mapper > :: LinkBuilder
    : public CommDataHandleIF< LinkBuilder, int >
    {
    public:
      typedef Space SpaceType;
      typedef Mapper MapperType;

      enum { nCodim = SpaceType :: GridType :: dimension + 1 };

      typedef int DataType;

      const int myRank_;
      const int mySize_;

      typedef SlaveDofs< Space,Mapper > IndexMapType;
      IndexMapType &slaves_;

      const SpaceType &space_;
      const MapperType &mapper_;

      LinkBuilder( IndexMapType &slaves, const SpaceType &space, const MapperType& mapper )
      : myRank_( space.gridPart().comm().rank() ), mySize_( space.gridPart().comm().size() ),
        slaves_( slaves ), space_( space ), mapper_( mapper )
      {}

      bool contains ( int dim, int codim ) const
      {
        return mapper_.contains( codim );
      }

      bool fixedsize ( int dim, int codim ) const
      {
        return false;
      }

      //! read buffer and apply operation
      template< class MessageBuffer, class Entity >
      void gather ( MessageBuffer &buffer, const Entity &entity ) const
      {
        // for sending ranks write rank
        if( sendRank( entity ) )
          buffer.write( myRank_ );
      }

      //! read buffer and apply operation
      //! scatter is called for one every entity
      //! several times depending on how much data
      //! was gathered
      template< class MessageBuffer, class EntityType >
      void scatter ( MessageBuffer &buffer, const EntityType &entity, size_t n )
      {
        if( n == 1 )
        {
          int rank;
          buffer.read( rank );

          assert( (rank >= 0) && (rank < mySize_) );

          // if entity in not interiorBorder insert anyway
          if ( rank < myRank_ || ! sendRank( entity ) )
            mapper_.mapEachEntityDof( entity, [this]( const int , const auto& value ){slaves_.insert( value );} );
        }
      }

      //! return local dof size to be communicated
      template< class Entity >
      size_t size ( const Entity &entity ) const
      {
        return (sendRank( entity )) ? 1 : 0;
      }

    protected:
      template <class Entity>
      bool sendRank(const Entity& entity) const
      {
        const PartitionType ptype = entity.partitionType();
        return (ptype == InteriorEntity) || (ptype == BorderEntity);
      }
    };

  //@}

  } // end namespace Fem

} // end namespace Dune
#endif // #ifndef DUNE_FEM_SLAVEDOFS_HH
