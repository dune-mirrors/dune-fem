#ifndef DUNE_FEM_SLAVEDOFS_HH
#define DUNE_FEM_SLAVEDOFS_HH

#include <cassert>
#include <cstddef>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <set>

#include <dune/common/exceptions.hh>
#include <dune/common/genericiterator.hh>
#include <dune/common/ftraits.hh>
#include <dune/common/typetraits.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/datahandleif.hh>

#include <dune/fem/storage/singletonlist.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/space/common/commindexmap.hh>
#include <dune/fem/storage/envelope.hh>

namespace Dune
{

  namespace Fem
  {

  /** @addtogroup Communication Communication
      @{
  **/

    template< class GridPart, class Mapper >
    class SlaveDofs
    {
      typedef SlaveDofs< GridPart, Mapper > ThisType;

      class LinkBuilder;

    public:
      /** \brief type of grid part **/
      typedef GridPart GridPartType;

      //! type of used mapper
      typedef Mapper MapperType;

    protected:
      typedef Fem :: CommunicationIndexMap IndexMapType;

      const GridPartType &gridPart_;
      const MapperType &mapper_;

      // type of communication indices
      IndexMapType slaves_;

    public:
      struct ConstIterator
        : public std::iterator< std::forward_iterator_tag, const int, int >
      {
        ConstIterator () = default;
        ConstIterator ( const IndexMapType &slaves, int index ) : slaves_( &slaves ), index_( index ) {}

        const int &operator* () const { return (*slaves_)[ index_ ]; }
        const int *operator-> () const { return &(*slaves_)[ index_ ]; }

        const int &operator[] ( int n ) const noexcept { return (*slaves_)[ index_ + n ]; }

        bool operator== ( const ConstIterator &other ) const { return (index_ == other.index_); }
        bool operator!= ( const ConstIterator &other ) const { return (index_ != other.index_); }

        ConstIterator &operator++ () { ++index_; return *this; }
        ConstIterator operator++ ( int ) { ConstIterator copy( *this ); ++(*this); return copy; }

        ThisType &operator-- () noexcept { --index_; return *this; }
        ThisType operator-- ( int ) noexcept { ThisType copy( *this ); --(*this); return copy; }

        ThisType &operator+= ( int n ) noexcept { index_ += n; return *this; }
        ThisType &operator-= ( int n ) noexcept { index_ -= n; return *this; }

        ThisType operator+ ( int n ) const noexcept { return ThisType( index_ + n ); }
        ThisType operator- ( int n ) const noexcept { return ThisType( index_ - n ); }

        friend ThisType operator+ ( int n, const ThisType &i ) noexcept { return i + n; }

        int operator- ( const ThisType &other ) const noexcept { return (index_ - other.index_); }

        bool operator< ( const ThisType &other ) const noexcept { return (index_ < other.index_); }
        bool operator<= ( const ThisType &other ) const noexcept { return (index_ <= other.index_); }
        bool operator>= ( const ThisType &other ) const noexcept { return (index_ >= other.index_); }
        bool operator> ( const ThisType &other ) const noexcept { return (index_ > other.index_); }

      private:
        const IndexMapType *slaves_ = nullptr;
        int index_ = 0;
      };

      SlaveDofs ( const GridPartType &gridPart, const MapperType &mapper )
        : gridPart_( gridPart ), mapper_( mapper )
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
      bool isSlave ( int index ) const { return std::binary_search( begin(), end(), index ); }

      void rebuild ()
      {
        std::set< int > slaveSet;
        buildMaps( slaveSet );
        slaves_.clear();
        slaves_.set( slaveSet );
      }

      const GridPartType &gridPart () const { return gridPart_; }

    protected:
      void buildMaps ( std::set< int > &slaveSet )
      {
        // build linkage and index maps
        for( int codim = 1; codim <= GridPartType::dimension; ++codim )
        {
          if( mapper_.contains( codim ) )
            return buildCommunicatedMaps( slaveSet );
        }
        return buildDiscontinuousMaps( slaveSet );
      }

      void buildDiscontinuousMaps ( std::set< int > &slaveSet )
      {
        // if DoFs are only attached to codimension 0, we do not have to communicate
        const auto idxpitype = GridPartType::indexSetPartitionType;
        for( auto it = gridPart().template begin< 0, idxpitype >(), end = gridPart().template end< 0, idxpitype >(); it != end; ++it )
        {
          const auto& entity = *it;
          if( entity.partitionType() != Dune::InteriorEntity )
            mapper_.mapEachEntityDof( entity, [ &slaveSet ] ( int, int value ) { slaveSet.insert( value ); } );
        }
        // insert overall size at the end
        slaveSet.insert( mapper_.size() );
      }

      void buildCommunicatedMaps ( std::set< int > &slaveSet )
      {
        // we have to skip communication when parallel program is executed only on one processor
        // otherwise YaspGrid and Lagrange polorder=2 fails :(
        if( gridPart().comm().size() > 1 )
        {
          try
          {
            LinkBuilder handle( slaveSet, gridPart(), mapper_ );
            gridPart().communicate( handle, GridPartType::indexSetInterfaceType, ForwardCommunication );
          }
          catch( const Exception &e )
          {
            std::cerr << e << std::endl << "Exception thrown in: " << __FILE__ << " line:" << __LINE__ << std::endl;
            std::abort();
          }
        }
        // insert overall size at the end
        slaveSet.insert( mapper_.size() );
      }
    };



    // SlaveDofs::LinkBuilder
    // ----------------------

    template< class GridPart, class Mapper >
    class SlaveDofs< GridPart, Mapper >::LinkBuilder
      : public CommDataHandleIF< LinkBuilder, int >
    {
    public:
      LinkBuilder( std::set< int > &slaveSet, const GridPartType &gridPart, const MapperType &mapper )
        : myRank_( gridPart.comm().rank() ), mySize_( gridPart.comm().size() ),
          slaveSet_( slaveSet ), mapper_( mapper )
      {}

      bool contains ( int dim, int codim ) const { return mapper_.contains( codim ); }
      bool fixedSize ( int dim, int codim ) const { return false; }

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
      void scatter ( MessageBuffer &buffer, const EntityType &entity, std::size_t n )
      {
        if( n == 1 )
        {
          int rank;
          buffer.read( rank );

          assert( (rank >= 0) && (rank < mySize_) );

          // if entity in not interiorBorder insert anyway
          if ( rank < myRank_ || ! sendRank( entity ) )
            mapper_.mapEachEntityDof( entity, [this]( const int , const auto& value ){slaveSet_.insert( value );} );
        }
      }

      //! return local dof size to be communicated
      template< class Entity >
      std::size_t size ( const Entity &entity ) const
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

    private:
      int myRank_, mySize_;
      std::set< int > &slaveSet_;
      const MapperType &mapper_;
    };



    // MasterDofs
    // ----------

    template< class SlaveDofs >
    struct MasterDofs;

    template< class GridPart, class Mapper >
    struct MasterDofs< SlaveDofs< GridPart, Mapper > >
    {
      typedef SlaveDofs< GridPart, Mapper > SlaveDofsType;

      struct ConstIterator
        : public std::iterator< std::forward_iterator_tag, int, std::ptrdiff_t, Envelope< int >, int >
      {
        ConstIterator () = default;

        ConstIterator ( int index, int slave )
          : index_( index ), slave_( slave )
        {}

        ConstIterator ( const SlaveDofsType &slaveDofs, int index, int slave )
          : slaveDofs_( &slaveDofs ), index_( index ), slave_( slave )
        {
          skipSlaves();
        }

        int operator* () const { return index_; }
        Envelope< int > operator-> () const { return Envelope< int >( index_ ); }

        bool operator== ( const ConstIterator &other ) const { return (index_ == other.index_); }
        bool operator!= ( const ConstIterator &other ) const { return (index_ != other.index_); }

        ConstIterator &operator++ () { ++index_; skipSlaves(); return *this; }
        ConstIterator operator++ ( int ) { ConstIterator copy( *this ); ++(*this); return copy; }

        const SlaveDofsType &slaveDofs () const { assert( slaveDofs_ ); return *slaveDofs_; }

      private:
        void skipSlaves ()
        {
          assert( slave_ < slaveDofs().size() );
          for( ; (index_ == slaveDofs()[ slave_ ]) && (++slave_ != slaveDofs().size()); ++index_ )
            continue;
        }

        const SlaveDofsType *slaveDofs_ = nullptr;
        int index_ = 0, slave_ = 0;
      };

      explicit MasterDofs ( const SlaveDofsType &slaveDofs )
        : slaveDofs_( slaveDofs )
      {}

      ConstIterator begin () const { return ConstIterator( slaveDofs_, 0, 0 ); }
      ConstIterator end () const { return ConstIterator( slaveDofs_[ slaveDofs_.size()-1 ], slaveDofs_.size() ); }

      int size () const { return slaveDofs_[ slaveDofs_.size()-1 ] - (slaveDofs_.size()-1); }

    private:
      const SlaveDofsType &slaveDofs_;
    };



    // masterDofs
    // ----------

    template< class SlaveDofs >
    inline static MasterDofs< SlaveDofs > masterDofs ( const SlaveDofs &slaveDofs )
    {
      return MasterDofs< SlaveDofs >( slaveDofs );
    }

  //@}

  } // end namespace Fem

} // end namespace Dune
#endif // #ifndef DUNE_FEM_SLAVEDOFS_HH
