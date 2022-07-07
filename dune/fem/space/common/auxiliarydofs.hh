#ifndef DUNE_FEM_SPACE_COMMON_AUXILIARYDOFS_HH
#define DUNE_FEM_SPACE_COMMON_AUXILIARYDOFS_HH

#include <cassert>
#include <cstddef>

#include <algorithm>
#include <iostream>
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
      typedef typename IndexMapType :: IndexType IndexType;

      struct ConstIterator
      {
        typedef std::forward_iterator_tag iterator_category;
        typedef const IndexType value_type;
        typedef IndexType difference_type;
        typedef IndexType* pointer;
        typedef IndexType& reference;

        ConstIterator () = default;
        ConstIterator ( const IndexMapType &auxiliarys, IndexType index ) : auxiliarys_( &auxiliarys ), index_( index ) {}

        const IndexType &operator* () const { return (*auxiliarys_)[ index_ ]; }
        const IndexType *operator-> () const { return &(*auxiliarys_)[ index_ ]; }

        const IndexType &operator[] ( IndexType n ) const noexcept { return (*auxiliarys_)[ index_ + n ]; }

        bool operator== ( const ConstIterator &other ) const { return (index_ == other.index_); }
        bool operator!= ( const ConstIterator &other ) const { return (index_ != other.index_); }

        ConstIterator &operator++ () { ++index_; return *this; }
        ConstIterator operator++ ( int ) { ConstIterator copy( *this ); ++(*this); return copy; }

        ThisType &operator-- () noexcept { --index_; return *this; }
        ThisType operator-- ( int ) noexcept { ThisType copy( *this ); --(*this); return copy; }

        ThisType &operator+= ( IndexType n ) noexcept { index_ += n; return *this; }
        ThisType &operator-= ( IndexType n ) noexcept { index_ -= n; return *this; }

        ThisType operator+ ( IndexType n ) const noexcept { return ThisType( index_ + n ); }
        ThisType operator- ( IndexType n ) const noexcept { return ThisType( index_ - n ); }

        friend ThisType operator+ ( IndexType n, const ThisType &i ) noexcept { return i + n; }

        IndexType operator- ( const ThisType &other ) const noexcept { return (index_ - other.index_); }

        bool operator< ( const ThisType &other ) const noexcept { return (index_ < other.index_); }
        bool operator<= ( const ThisType &other ) const noexcept { return (index_ <= other.index_); }
        bool operator>= ( const ThisType &other ) const noexcept { return (index_ >= other.index_); }
        bool operator> ( const ThisType &other ) const noexcept { return (index_ > other.index_); }

      private:
        const IndexMapType *auxiliarys_ = nullptr;
        IndexType index_ = 0;
      };

      AuxiliaryDofs ( const GridPartType &gridPart, const MapperType &mapper )
        : gridPart_( gridPart ), mapper_( mapper )
      {}

      AuxiliaryDofs ( const AuxiliaryDofs& ) = delete;

      //! return dof number of auxiliary for index
      IndexType operator [] ( const IndexType index ) const
      {
        return auxiliarys_[ index ];
      }

      //! return number of auxiliary dofs
      IndexType size () const
      {
        return auxiliarys_.size();
      }

      //! return number of primaryDofs
      IndexType primarySize () const
      {
        const IndexType last = auxiliarys_.size() - 1;
        // last entry is the overall size, thus substract size - 1 from the
        // overall number to obtain the number of primary dofs
        return auxiliarys_[ last ] - last ;
      }

      ConstIterator begin () const { return ConstIterator( auxiliarys_, 0 ); }
      ConstIterator end () const { assert( size() > 0 ); return ConstIterator( auxiliarys_, size()-1 ); }

      //! return true if index is contained, meaning it is a auxiliary dof
      bool contains ( IndexType index ) const { return std::binary_search( begin(), end(), index ); }

      [[deprecated("Use contains instead")]]
      bool isSlave ( IndexType index ) const { return contains( index ); }

      void rebuild ()
      {
        std::set< IndexType > auxiliarySet;
        buildMaps( auxiliarySet );
        auxiliarys_.clear();
        auxiliarys_.set( auxiliarySet );
      }

      const GridPartType &gridPart () const { return gridPart_; }

    protected:
      void buildMaps ( std::set< IndexType > &auxiliarySet )
      {
        // build linkage and index maps
        for( int codim = 1; codim <= GridPartType::dimension; ++codim )
        {
          if( mapper_.contains( codim ) )
            return buildCommunicatedMaps( auxiliarySet );
        }
        return buildDiscontinuousMaps( auxiliarySet );
      }

      void buildDiscontinuousMaps ( std::set< IndexType > &auxiliarySet )
      {
        // if DoFs are only attached to codimension 0, we do not have to communicate
        const auto idxpitype = GridPartType::indexSetPartitionType;
        for( auto it = gridPart().template begin< 0, idxpitype >(), end = gridPart().template end< 0, idxpitype >(); it != end; ++it )
        {
          const auto& entity = *it;
          if( entity.partitionType() != Dune::InteriorEntity )
            mapper_.mapEachEntityDof( entity, [ &auxiliarySet ] ( IndexType, IndexType value ) { auxiliarySet.insert( value ); } );
        }
        // insert overall size at the end
        auxiliarySet.insert( mapper_.size() );
      }

      void buildCommunicatedMaps ( std::set< IndexType > &auxiliarySet )
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
      : public CommDataHandleIF< LinkBuilder, int > // int is data type to be communicated
    {
    public:
      LinkBuilder( std::set< IndexType > &auxiliarySet, const GridPartType &gridPart, const MapperType &mapper )
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

    protected:
      int myRank_, mySize_;
      std::set< IndexType > &auxiliarySet_;
      const MapperType &mapper_;
    };


    /** \brief @ingroup Communication
     *
     *  Apply action encoded in Functor f to all auxiliary dofs.
     *  \param[in]  auxiliaryDofs  AuxiliaryDofs instance (from space)
     *  \param[in]  f  a Functor or lambda offering operator()( const size_t dof )
     */
    template <class AuxiliaryDofs, class F>
    static void forEachAuxiliaryDof( const AuxiliaryDofs& auxiliaryDofs, F&& f )
    {
      // don't delete the last since this is the overall Size
      const size_t auxiliarySize = auxiliaryDofs.size() - 1;
      for(size_t auxiliary = 0; auxiliary<auxiliarySize; ++auxiliary)
      {
        // apply action to auxiliary dof
        f( auxiliaryDofs[auxiliary] );
      }
    }

    /** \brief @ingroup Communication
     *
     *  Apply action encoded in Functor f to all primary dofs.
     *  \param[in]  auxiliaryDofs  AuxiliaryDofs instance (from space)
     *  \param[in]  f  a Functor or lambda offering operator()( const size_t dof )
     */
    template <class AuxiliaryDofs, class F>
    static void forEachPrimaryDof( const AuxiliaryDofs& auxiliaryDofs, F&& f )
    {
      const size_t numAuxiliarys = auxiliaryDofs.size();
      for( size_t auxiliary = 0, dof = 0 ; auxiliary < numAuxiliarys; ++auxiliary, ++dof  )
      {
        const size_t nextAuxiliary = auxiliaryDofs[ auxiliary ];
        for(; dof < nextAuxiliary; ++dof )
        {
          // apply action to primary dof
          f( dof );
        }
      }
    }


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
      typedef typename AuxiliaryDofsType :: IndexType IndexType;

      struct ConstIterator
      {
        typedef std::forward_iterator_tag iterator_category;
        typedef IndexType value_type;
        typedef std::ptrdiff_t difference_type;
        typedef Envelope< IndexType > pointer;
        typedef IndexType reference;

        ConstIterator () = default;

        ConstIterator ( IndexType index, IndexType auxiliary )
          : index_( index ), auxiliary_( auxiliary )
        {}

        ConstIterator ( const AuxiliaryDofsType &auxiliaryDofs, IndexType index, IndexType auxiliary )
          : auxiliaryDofs_( &auxiliaryDofs ), index_( index ), auxiliary_( auxiliary )
        {
          skipAuxiliarys();
        }

        IndexType operator* () const { return index_; }
        Envelope< IndexType > operator-> () const { return Envelope< IndexType >( index_ ); }

        bool operator== ( const ConstIterator &other ) const { return (index_ == other.index_); }
        bool operator!= ( const ConstIterator &other ) const { return (index_ != other.index_); }

        ConstIterator &operator++ () { ++index_; skipAuxiliarys(); return *this; }
        ConstIterator operator++ ( int ) { ConstIterator copy( *this ); ++(*this); return copy; }

        const AuxiliaryDofsType &auxiliaryDofs () const { assert( auxiliaryDofs_ ); return *auxiliaryDofs_; }

        bool contains( const IndexType index ) const { return ! auxiliaryDofs().contains( index ); }

      private:
        void skipAuxiliarys ()
        {
          const IndexType aSize = auxiliaryDofs().size();
          assert( auxiliary_ < aSize );
          for( ; (index_ == auxiliaryDofs()[ auxiliary_ ]) && (++auxiliary_ != aSize); ++index_ )
            continue;
        }

        const AuxiliaryDofsType *auxiliaryDofs_ = nullptr;
        IndexType index_ = 0, auxiliary_ = 0;
      };

      [[deprecated("Use forEachPrimaryDof instead!")]]
      explicit PrimaryDofs ( const AuxiliaryDofsType &auxiliaryDofs )
        : auxiliaryDofs_( auxiliaryDofs )
      {}

      ConstIterator begin () const { return ConstIterator( auxiliaryDofs_, 0, 0 ); }
      ConstIterator end () const { return ConstIterator( auxiliaryDofs_[ auxiliaryDofs_.size()-1 ], auxiliaryDofs_.size() ); }

      IndexType size () const { return auxiliaryDofs_[ auxiliaryDofs_.size()-1 ] - (auxiliaryDofs_.size()-1); }

    private:
      const AuxiliaryDofsType &auxiliaryDofs_;
    };



    // primaryDofs
    // -----------

    template< class AuxiliaryDofs >
    [[deprecated("Use forEachPrimaryDof instead!" )]]
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
