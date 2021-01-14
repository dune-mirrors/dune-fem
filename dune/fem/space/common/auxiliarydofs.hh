#ifndef DUNE_FEM_SPACE_COMMON_AUXILIARYDOFS_HH
#define DUNE_FEM_SPACE_COMMON_AUXILIARYDOFS_HH

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

    /** \brief @ingroup Communication
     *
     *  In parallel computations the dofs of a discrete function are made up by
     *  all primary dofs. For technical reasons some dofs exists on multiply
     *  processes but are only primary on exactly one process. Dofs on processes
     *  that are not primary are called auxiliary.
     */
    template< class GridPart, class Mapper >
    class AuxiliaryDofs
    {
      typedef AuxiliaryDofs< GridPart, Mapper > ThisType;

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
      IndexMapType auxiliarys_;

    public:
      struct ConstIterator
        : public std::iterator< std::forward_iterator_tag, const int, int >
      {
        ConstIterator () = default;
        ConstIterator ( const IndexMapType &auxiliarys, int index ) : auxiliarys_( &auxiliarys ), index_( index ) {}

        const int &operator* () const { return (*auxiliarys_)[ index_ ]; }
        const int *operator-> () const { return &(*auxiliarys_)[ index_ ]; }

        const int &operator[] ( int n ) const noexcept { return (*auxiliarys_)[ index_ + n ]; }

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
        const IndexMapType *auxiliarys_ = nullptr;
        int index_ = 0;
      };

      AuxiliaryDofs ( const GridPartType &gridPart, const MapperType &mapper )
        : gridPart_( gridPart ), mapper_( mapper )
      {}

      AuxiliaryDofs ( const AuxiliaryDofs& ) = delete;

      //! return dof number of auxiliary for index
      int operator [] ( const int index ) const
      {
        return auxiliarys_[ index ];
      }

      //! return number of auxiliary dofs
      int size () const
      {
        return auxiliarys_.size();
      }

      ConstIterator begin () const { return ConstIterator( auxiliarys_, 0 ); }
      ConstIterator end () const { assert( size() > 0 ); return ConstIterator( auxiliarys_, size()-1 ); }

      //! return true if index is contained, meaning it is a auxiliary dof
      bool contains ( int index ) const { return std::binary_search( begin(), end(), index ); }

      [[deprecated("Use contains instead")]]
      bool isSlave ( int index ) const { return contains( index ); }

      void rebuild ()
      {
        std::set< int > auxiliarySet;
        buildMaps( auxiliarySet );
        auxiliarys_.clear();
        auxiliarys_.set( auxiliarySet );
      }

      const GridPartType &gridPart () const { return gridPart_; }

    protected:
      void buildMaps ( std::set< int > &auxiliarySet )
      {
        // build linkage and index maps
        for( int codim = 1; codim <= GridPartType::dimension; ++codim )
        {
          if( mapper_.contains( codim ) )
            return buildCommunicatedMaps( auxiliarySet );
        }
        return buildDiscontinuousMaps( auxiliarySet );
      }

      void buildDiscontinuousMaps ( std::set< int > &auxiliarySet )
      {
        // if DoFs are only attached to codimension 0, we do not have to communicate
        const auto idxpitype = GridPartType::indexSetPartitionType;
        for( auto it = gridPart().template begin< 0, idxpitype >(), end = gridPart().template end< 0, idxpitype >(); it != end; ++it )
        {
          const auto& entity = *it;
          if( entity.partitionType() != Dune::InteriorEntity )
            mapper_.mapEachEntityDof( entity, [ &auxiliarySet ] ( int, int value ) { auxiliarySet.insert( value ); } );
        }
        // insert overall size at the end
        auxiliarySet.insert( mapper_.size() );
      }

      void buildCommunicatedMaps ( std::set< int > &auxiliarySet )
      {
        // we have to skip communication when parallel program is executed only on one processor
        // otherwise YaspGrid and Lagrange polorder=2 fails :(
        if( gridPart().comm().size() > 1 )
        {
          try
          {
            LinkBuilder handle( auxiliarySet, gridPart(), mapper_ );
            gridPart().communicate( handle, GridPartType::indexSetInterfaceType, ForwardCommunication );
          }
          catch( const Exception &e )
          {
            std::cerr << e << std::endl << "Exception thrown in: " << __FILE__ << " line:" << __LINE__ << std::endl;
            std::abort();
          }
        }
        // insert overall size at the end
        auxiliarySet.insert( mapper_.size() );
      }
    };



    // AuxiliaryDofs::LinkBuilder
    // ----------------------

    template< class GridPart, class Mapper >
    class AuxiliaryDofs< GridPart, Mapper >::LinkBuilder
      : public CommDataHandleIF< LinkBuilder, int >
    {
    public:
      LinkBuilder( std::set< int > &auxiliarySet, const GridPartType &gridPart, const MapperType &mapper )
        : myRank_( gridPart.comm().rank() ), mySize_( gridPart.comm().size() ),
          auxiliarySet_( auxiliarySet ), mapper_( mapper )
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
            mapper_.mapEachEntityDof( entity, [this]( const int , const auto& value ){auxiliarySet_.insert( value );} );
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
      std::set< int > &auxiliarySet_;
      const MapperType &mapper_;
    };



    // PrimaryDofs
    // -----------

    /** \brief @ingroup Communication
     *
     *  In parallel computations the dofs of a discrete function are made up by
     *  all primary dofs. For technical reasons some dofs exists on multiply
     *  processes but are only primary on exactly one process.
     */
    template< class AuxiliaryDofs >
    struct PrimaryDofs;

    template< class GridPart, class Mapper >
    struct PrimaryDofs< AuxiliaryDofs< GridPart, Mapper > >
    {
      typedef AuxiliaryDofs< GridPart, Mapper > AuxiliaryDofsType;

      struct ConstIterator
        : public std::iterator< std::forward_iterator_tag, int, std::ptrdiff_t, Envelope< int >, int >
      {
        ConstIterator () = default;

        ConstIterator ( int index, int auxiliary )
          : index_( index ), auxiliary_( auxiliary )
        {}

        ConstIterator ( const AuxiliaryDofsType &auxiliaryDofs, int index, int auxiliary )
          : auxiliaryDofs_( &auxiliaryDofs ), index_( index ), auxiliary_( auxiliary )
        {
          skipAuxiliarys();
        }

        int operator* () const { return index_; }
        Envelope< int > operator-> () const { return Envelope< int >( index_ ); }

        bool operator== ( const ConstIterator &other ) const { return (index_ == other.index_); }
        bool operator!= ( const ConstIterator &other ) const { return (index_ != other.index_); }

        ConstIterator &operator++ () { ++index_; skipAuxiliarys(); return *this; }
        ConstIterator operator++ ( int ) { ConstIterator copy( *this ); ++(*this); return copy; }

        const AuxiliaryDofsType &auxiliaryDofs () const { assert( auxiliaryDofs_ ); return *auxiliaryDofs_; }

        bool contains( const int index ) const { return ! auxiliaryDofs().contains( index ); }

      private:
        void skipAuxiliarys ()
        {
          const int aSize = auxiliaryDofs().size();
          assert( auxiliary_ < aSize );
          for( ; (index_ == auxiliaryDofs()[ auxiliary_ ]) && (++auxiliary_ != aSize); ++index_ )
            continue;
        }

        const AuxiliaryDofsType *auxiliaryDofs_ = nullptr;
        int index_ = 0, auxiliary_ = 0;
      };

      explicit PrimaryDofs ( const AuxiliaryDofsType &auxiliaryDofs )
        : auxiliaryDofs_( auxiliaryDofs )
      {}

      ConstIterator begin () const { return ConstIterator( auxiliaryDofs_, 0, 0 ); }
      ConstIterator end () const { return ConstIterator( auxiliaryDofs_[ auxiliaryDofs_.size()-1 ], auxiliaryDofs_.size() ); }

      int size () const { return auxiliaryDofs_[ auxiliaryDofs_.size()-1 ] - (auxiliaryDofs_.size()-1); }

    private:
      const AuxiliaryDofsType &auxiliaryDofs_;
    };



    // primaryDofs
    // -----------

    template< class AuxiliaryDofs >
    inline static PrimaryDofs< AuxiliaryDofs > primaryDofs ( const AuxiliaryDofs &auxiliaryDofs )
    {
      return PrimaryDofs< AuxiliaryDofs >( auxiliaryDofs );
    }

    template< class AuxiliaryDofs >
    [[deprecated("Use primaryDofs instead!" )]]
    inline static PrimaryDofs< AuxiliaryDofs > masterDofs ( const AuxiliaryDofs &auxiliaryDofs )
    {
      return PrimaryDofs< AuxiliaryDofs >( auxiliaryDofs );
    }

  //@}

  } // end namespace Fem

} // end namespace Dune
#endif // #ifndef DUNE_FEM_SPACE_COMMON_AUXILIARYDOFS_HH
