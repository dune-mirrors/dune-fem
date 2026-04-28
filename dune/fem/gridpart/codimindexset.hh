#ifndef DUNE_FEM_CODIMINDEXSET_HH
#define DUNE_FEM_CODIMINDEXSET_HH

#include <algorithm>
#include <set>

#include <dune/grid/utility/persistentcontainer.hh>
#include <dune/grid/utility/persistentcontainervector.hh>
#include <dune/grid/utility/persistentcontainerwrapper.hh>
#include <dune/grid/utility/persistentcontainermap.hh>

#include <dune/fem/io/streams/streams.hh>
#include <dune/fem/storage/dynamicarray.hh>
#include <dune/fem/space/common/commindexmap.hh>

namespace Dune
{

  namespace Fem
  {

    //***********************************************************************
    //
    //  Index Set for one codimension
    //  --CodimIndexSet
    //
    //***********************************************************************
    template <class GridImp>
    class CodimIndexSet
    {
    protected:
      typedef GridImp GridType;
      typedef CodimIndexSet < GridType >  ThisType;

      // reference to grid
      const GridType& grid_;

      typedef unsigned char INDEXSTATE;
      const INDEXSTATE stateUnused_ = 0; // unused indices
      const INDEXSTATE stateUsed_   = 1; // used indices
      const INDEXSTATE stateNew_    = 2; // new indices

    public:
      // type of exported index (local to one core)
      typedef typename Fem::CommunicationIndexMap::IndexType IndexType;

      // indices in this status have not been initialized
      static IndexType invalidIndex() { return IndexType(-1); }

    protected:
      // array type for indices
      typedef DynamicArray< IndexType  > IndexArrayType;
      typedef DynamicArray< INDEXSTATE > IndexStateArrayType;

      // use the imporved PersistentContainer
      typedef PersistentContainer< GridType, IndexType > IndexContainerType;

      // the mapping of the global to leaf index
      IndexContainerType  leafIndex_;
      IndexStateArrayType indexState_;

      // stack for holes
      IndexArrayType holes_;

      // Array that only remeber the occuring
      // holes (for compress of data)
      IndexArrayType oldIdx_;
      IndexArrayType newIdx_;

      // last size of set before compress (needed in parallel runs)
      IndexType lastSize_;

      // codim for which index is provided
      const int myCodim_;

      // actual number of holes
      IndexType numberHoles_;

    public:

      //! Constructor taking memory factor (default = 1.1)
      CodimIndexSet (const GridType& grid,
                     const int codim,
                     const double memoryFactor = 1.1)
        : grid_( grid )
        , leafIndex_( grid, codim, invalidIndex() )
        , indexState_( 0 )
        , holes_(0)
        , oldIdx_(0)
        , newIdx_(0)
        , lastSize_ (0)
        , myCodim_( codim )
        , numberHoles_(0)
      {
        setMemoryFactor(memoryFactor);
      }

      //! set memory overestimation factor
      void setMemoryFactor(const double memoryFactor)
      {
        indexState_.setMemoryFactor( memoryFactor );
        holes_.setMemoryFactor(memoryFactor);
        oldIdx_.setMemoryFactor(memoryFactor);
        newIdx_.setMemoryFactor(memoryFactor);
      }

      //! reallocate the vectors
      void resize () { leafIndex_.resize( invalidIndex() ); }

      //! prepare for setup (nothing to do here)
      void prepareCompress () {}

    public:
      //! clear set
      void clear()
      {
        // set all values to invalidIndex
        leafIndex_.fill( invalidIndex() );
        // free all indices
        indexState_.resize( 0 );
      }

      //! set all entries to unused
      void resetUsed ()
      {
        std::fill( indexState_.begin(), indexState_.end(), stateUnused_ );
      }

      bool consecutive ()
      {
        std::set< IndexType > found ;
        // Something's wrong here: This method _must_ always return true
        typedef typename IndexContainerType::Iterator Iterator;
        bool consecutive = true;
        const Iterator end = leafIndex_.end();
        for( Iterator it = leafIndex_.begin(); it != end; ++it )
        {
          if( *it != invalidIndex() )
          {
            if( found.find( *it ) != found.end() )
            {
              std::cout << "index " << *it << " exists twice " << std::endl;
            }
            assert( found.find( *it ) == found.end() );
            found.insert( *it );
            //consecutive &= (*it < IndexType( indexState_.size() ));
          }
          //else
          //  consecutive = false;
          consecutive &= (*it < IndexType( indexState_.size() ));
        }
        return consecutive;
      }

      //! set all entries to unused
      void checkConsecutive () { assert( consecutive() ); }

      //! clear holes, i.e. set number of holes to zero
      void clearHoles()
      {
        // set number of holes to zero
        numberHoles_ = 0;
        // remember actual size
        lastSize_ = indexState_.size();
      }

      //! make to index numbers consecutive
      //! return true, if at least one hole was closed
      bool compress ()
      {
        const IndexType sizeOfVecs = indexState_.size();
        holes_.resize( sizeOfVecs );

        // true if a least one dof must be copied
        bool haveToCopy = false;

        // mark holes
        IndexType actHole = 0;
        for( IndexType index = 0; index < sizeOfVecs; ++index )
        {
          // create vector with all holes
          if( indexState_[ index ] == stateUnused_ )
            holes_[ actHole++ ] = index;
        }

        // the new size is the actual size minus the holes
        const IndexType actSize = sizeOfVecs - actHole;

        // resize hole storing vectors
        oldIdx_.resize(actHole);
        newIdx_.resize(actHole);

        // only compress if number of holes > 0
        if(actHole > 0)
        {
          // for int this is -1 or a large number for unsigned types
          const IndexType invIndex = invalidIndex();

          // close holes
          IndexType holes = 0; // number of real holes
          typedef typename IndexContainerType::Iterator Iterator;
          const Iterator end = leafIndex_.end();
          for( Iterator it = leafIndex_.begin(); it != end; ++it )
          {
            IndexType& index = *it;
            if( index == invIndex )
            {
              continue ;
            }
            else if( indexState_[ index ] == stateUnused_ )
            {
              index = invIndex;
              continue ;
            }

            // a index that is used but larger then actual size
            // has to move to a hole
            // if used index lies behind size, then index has to move
            // to one of the holes
            if( index >= actSize )
            {
              //std::cout << "Check index " << index << std::endl;
              // serach next hole that is smaler than actual size
              --actHole;
              // if actHole < 0 then error, because we have index larger then
              // actual size
              assert(actHole >= 0);
              while ( holes_[actHole] >= actSize )
              {
                assert(actHole > 0);
                --actHole;

                if( actHole == invIndex ) break;
              }

              assert(actHole >= 0);

#if HAVE_MPI
              // only for none-ghost elements hole storage is applied
              // this is because ghost indices might have in introduced
              // after the resize was done.
              if( indexState_[ index ] == stateUsed_ )
#endif
              {
                // remember old and new index
                oldIdx_[holes] = index;
                newIdx_[holes] = holes_[actHole];
                ++holes;
              }

              index = holes_[actHole];

              // means that dof manager has to copy the mem
              haveToCopy = true;
            }
          }

          // this call only sets the size of the vectors
          oldIdx_.resize(holes);
          newIdx_.resize(holes);

          // mark holes as new
          // note: This needs to be done after reassignment, so that their
          //       original entry will still see them as stateUnused_.
          for( IndexType hole = 0; hole < holes; ++hole )
            indexState_[ newIdx_[ hole ] ] = stateNew_;

        } // end if actHole > 0

        // store number of actual holes
        numberHoles_ = oldIdx_.size();

        // adjust size of container to correct value
        leafIndex_.resize( invalidIndex() );

        // resize vector of index states
        indexState_.resize( actSize );

#ifndef NDEBUG
        for( IndexType i=0; i<actSize; ++i )
          assert( indexState_[ i ] == stateUsed_ ||
                  indexState_[ i ] == stateUnused_ ||
                  indexState_[ i ] == stateNew_ );

        checkConsecutive();
#endif
        return haveToCopy;
      }

      //! return how much extra memory is needed for restriction
      IndexType additionalSizeEstimate () const { return indexState_.size(); }

      //! return size of grid entities per level and codim
      IndexType size () const { return indexState_.size(); }

      //! return size of grid entities per level and codim
      IndexType realSize () const
      {
        return leafIndex_.size();
      }

      //! return leaf index for given entity
      //- --index
      template <class EntityType>
      IndexType index ( const EntityType& entity ) const
      {
        assert( myCodim_ == EntityType :: codimension );
        assert( checkValidIndex( leafIndex_[ entity ] ) );
        return leafIndex_[ entity ];
      }

      //! return leaf index for given entity
      template <class EntityType>
      IndexType subIndex ( const EntityType& entity,
                           const int subNumber ) const
      {
        return subIndex( entity, subNumber,
                         std::integral_constant<bool, EntityType::codimension == 0 > () );
      }

      //! return leaf index for given entity
      template <class EntityType>
      IndexType subIndex ( const EntityType& entity,
                           const int subNumber,
                           const std::integral_constant<bool,false> codim0 ) const
      {
        DUNE_THROW(NotImplemented,"CodimIndexSet::subIndex: not implemented for entities with codim > 0" );
        return IndexType( -1 );
      }

      //! return leaf index for given entity
      template <class EntityType>
      IndexType subIndex ( const EntityType& entity,
                           const int subNumber,
                           const std::integral_constant<bool,true> codim0 ) const
      {
        assert( checkValidIndex( leafIndex_( entity, subNumber ) ) );
        return leafIndex_( entity, subNumber );
      }

      //! return state of index for given hierarchic number
      template <class EntityType>
      bool exists ( const EntityType& entity ) const
      {
        assert( myCodim_ == EntityType :: codimension );
        const IndexType &index = leafIndex_[ entity ];
        // if index is invalid (-1) it does not exist
        if (index==invalidIndex()) return false;
        assert( index < IndexType( indexState_.size() ) );
        return (indexState_[ index ] != stateUnused_);
      }

      template <class EntityType>
      bool exists ( const EntityType& entity ,
                    const int subNumber ) const
      {
        assert( 0 == EntityType :: codimension );
        const IndexType &index = leafIndex_( entity, subNumber );
        // if index is invalid (-1) it does not exist
        if (index==invalidIndex()) return false;
        assert( index < IndexType( indexState_.size() ) );
        return (indexState_[ index ] != stateUnused_);
      }

      //! return number of holes
      IndexType numberOfHoles () const
      {
        return numberHoles_;
      }

      //! return old index, for dof manager only
      IndexType oldIndex ( IndexType elNum ) const
      {
        assert( numberHoles_ == IndexType(oldIdx_.size()) );
        return oldIdx_[elNum];
      }

      //! return new index, for dof manager only returns index
      IndexType newIndex ( IndexType elNum) const
      {
        assert( numberHoles_ == IndexType(newIdx_.size()) );
        return newIdx_[elNum];
      }

      // insert element and create index for element number
      template <class EntityType>
      void insert (const EntityType& entity )
      {
        assert( myCodim_ == EntityType :: codimension );
        insertIdx( leafIndex_[ entity ] );
      }

      // insert element and create index for element number
      template <class EntityType>
      void insertSubEntity (const EntityType& entity,
                            const int subNumber)
      {
        assert( 0 == EntityType :: codimension );
        insertIdx( leafIndex_( entity, subNumber ) );
      }

      // insert element as ghost and create index for element number
      template <class EntityType>
      void insertGhost (const EntityType& entity )
      {
        assert( myCodim_ == EntityType :: codimension );
        // insert index
        IndexType &leafIdx = leafIndex_[ entity ];
        insertIdx( leafIdx );

        // if index is also larger than lastSize
        // mark as new to skip old-new index lists
        if( leafIdx >= lastSize_ )
          indexState_[ leafIdx ] = stateNew_;
      }

      // insert element and create index for element number
      template <class EntityType>
      void markForRemoval( const EntityType& entity )
      {
        assert( myCodim_ == EntityType :: codimension );
        const IndexType &index = leafIndex_[ entity ];
        if( index != invalidIndex() )
          indexState_[ index ] = stateUnused_;
      }

      // insert element as ghost and create index for element number
      template <class EntityType>
      bool validIndex (const EntityType& entity ) const
      {
        assert( myCodim_ == EntityType :: codimension );
        return (leafIndex_[ entity ] >= 0);
      }

      void print( std::ostream& out ) const
      {
        typedef typename IndexContainerType::ConstIterator Iterator;
        const Iterator end = leafIndex_.end();
        for( Iterator it = leafIndex_.begin(); it != end; ++it )
        {
          const IndexType &leafIdx = *it;
          if( leafIdx < indexState_.size() )
          {
            out << "idx: " << leafIdx << "  stat: " << indexState_[ leafIdx ] << std::endl;
          }
        }
      }

    protected:
      // return true if the index idx is valid
      bool checkValidIndex( const IndexType& idx ) const
      {
        assert( idx != invalidIndex() );
        assert( idx  < size() );
        return (idx != invalidIndex() ) && ( idx < size() );
      }

      // insert element and create index for element number
      void insertIdx ( IndexType &index )
      {
        if( index == invalidIndex() )
        {
          index = indexState_.size();
          indexState_.resize( index+1 );
        }
        assert( index < IndexType( indexState_.size() ) );
        indexState_[ index ] = stateUsed_;
      }

    public:
      // write to stream
      template <class StreamTraits>
      bool write(OutStreamInterface< StreamTraits >& out) const
      {
        // store current index set size
        // don't write something like  out << indexState_.size()
        // since on read you then don't know exactly what
        // type has been written, it must be the same types
        const uint32_t indexSize = indexState_.size();
        out << indexSize;

        // for consistency checking, write size as 64bit integer
        const uint64_t mysize = leafIndex_.size();
        out << mysize ;

        // backup indices
        typedef typename IndexContainerType::ConstIterator ConstIterator;
        const ConstIterator end = leafIndex_.end();
        for( ConstIterator it = leafIndex_.begin(); it != end; ++it )
          out << *it;

        return true;
      }

      // read from stream
      template <class StreamTraits>
      bool read(InStreamInterface< StreamTraits >& in)
      {
        // read current index set size
        uint32_t size = 0;
        in >> size;

        // mark all indices used
        indexState_.resize( size );
        std::fill( indexState_.begin(), indexState_.end(), stateUsed_ );

        // for consistency checking
        uint64_t storedSize = 0;
        in >> storedSize ;

        uint64_t leafsize = leafIndex_.size();
        // the stored size can be larger (visualization of parallel grids in serial)
        if( storedSize < leafsize )
        {
          DUNE_THROW(InvalidStateException,"CodimIndexSet: size consistency check failed during restore!");
        }

        // restore indices
        typedef typename IndexContainerType::Iterator Iterator;
        const Iterator end = leafIndex_.end();
        uint64_t count = 0 ;
        for( Iterator it = leafIndex_.begin(); it != end; ++it, ++count )
          in >> *it;

        // also read indices that were stored but are not needed on read
        if( count < storedSize )
        {
          IndexType value ;
          const uint64_t leftOver = storedSize - count ;
          for( uint64_t i = 0; i < leftOver; ++i )
            in >> value ;
        }

        return true;
      }

    }; // end of CodimIndexSet

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_CODIMINDEXSET_HH
