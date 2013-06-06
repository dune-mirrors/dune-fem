#ifndef DUNE_FEM_CODIMINDEXSET_HH
#define DUNE_FEM_CODIMINDEXSET_HH

#include <algorithm>

//- Dune includes 
#include <dune/common/misc.hh>
#include <dune/fem/space/common/arrays.hh>
#include <dune/fem/gridpart/defaultindexsets.hh>

#include <dune/fem/io/streams/xdrstreams.hh>
#include <dune/grid/utility/persistentcontainer.hh>

#include <dune/grid/utility/persistentcontainervector.hh>
#include <dune/grid/utility/persistentcontainerwrapper.hh>
#include <dune/grid/utility/persistentcontainermap.hh>

#ifdef ENABLE_ADAPTIVELEAFINDEXSET_FOR_YASPGRID
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/sgrid.hh>

namespace Dune 
{

  // PersistentContainer for YaspGrid
  // --------------------------------

  template< int dim, class Data >
  class PersistentContainer< YaspGrid< dim >, Data >
  : public PersistentContainerVector< YaspGrid< dim >, 
                                      typename YaspGrid< dim >::LeafIndexSet,
                                      std::vector<Data> >
  {
    typedef YaspGrid< dim > Grid;
    typedef PersistentContainerVector< Grid, typename Grid::LeafIndexSet, std::vector<Data> > BaseType;

  public:
    //! Constructor filling the container with values using the default constructor 
    //! Depending on the implementation this could be achieved without allocating memory
    PersistentContainer ( const Grid &grid, const int codim, const Data& value = Data() )
    : BaseType( grid, codim, grid.leafIndexSet(), 1.0, value )
    {}
  };

  // PersistentContainer for SGrid
  // -------------------------------

  template< int dim, int dimworld, class ctype, class Data >
  class PersistentContainer< SGrid< dim, dimworld, ctype >, Data >
  : public PersistentContainerVector< SGrid< dim, dimworld, ctype >, 
                                      typename SGrid< dim, dimworld, ctype >::LeafIndexSet,
                                      std::vector<Data> >
  {
    typedef SGrid< dim, dimworld, ctype > Grid ;
    typedef PersistentContainerVector< Grid, typename Grid::LeafIndexSet, std::vector<Data> > BaseType;

  public:
    //! Constructor filling the container with values using the default constructor 
    //! Depending on the implementation this could be achieved without allocating memory
    PersistentContainer ( const Grid &grid, const int codim, const Data& value = Data() )
    : BaseType( grid, codim, grid.leafIndexSet(), 1.0, value )
    {}
  };

}
#endif


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

    private:
      enum INDEXSTATE { UNUSED = 0,  // unused indices
                        USED   = 1,  // used indices
                        NEW    = 2 };//  new indices 

      // reference to grid 
      const GridType& grid_;

    public:
      // type of exported index 
      typedef int IndexType ;

    protected:  
      // indices in this status have not been initialized 
      static IndexType invalidIndex() { return -1; }

      // array type for indices 
      typedef MutableArray< IndexType > IndexArrayType;
      typedef MutableArray< INDEXSTATE > IndexStateArrayType;

      class IndexPersistentContainer 
      : public PersistentContainer< GridType, IndexType >
      {
        typedef PersistentContainer< GridType, IndexType > BaseType;

        // classes to make protected members public 
        template <class G, class T> 
        struct PublicPersistentContainerWrapper : public PersistentContainerWrapper< G, T >
        {
          using PersistentContainerWrapper< G, T > :: hostContainer_;
        };

        template <class G, class I, class V> 
        struct PublicPersistentContainerVector : public PersistentContainerVector< G, I, V >
        {
          using PersistentContainerVector< G, I, V > :: indexSet;
        };

        template <class G, class I, class M> 
        struct PublicPersistentContainerMap : public PersistentContainerMap< G, I, M >
        {
          using PersistentContainerMap< G, I, M > :: idSet;
        };

      public:
        using BaseType :: size ;
        using BaseType :: resize ;
        using BaseType :: codimension ;

        typedef typename BaseType :: Value Value;
        typedef typename BaseType :: Size  Size ;

        //! constructor needed by CodimIndexSet 
        IndexPersistentContainer( const GridType& grid, const int codim, const Value& value )
          : BaseType( grid, codim, value ) 
        {}

        //! only do a resize if the current size is smaller then the needed size 
        void enlargeOnly( const Value& value = Value() ) 
        {
          // call corrected implementation 
          enlargeImpl( *this, value ); 
        }

      protected:  
        // specialization for PersistentContainerWrapper that calles the other methods 
        template <class G, class T> 
        void enlargeImpl( PersistentContainerWrapper< G, T >& container, const Value& value ) 
        {
          enlargeImpl( ((PublicPersistentContainerWrapper< G, T > &) container).hostContainer_, value );
        }

        // enlarge implementation for persistent containers based on vectors 
        template < class G, class IndexSet, class Vector >
        void enlargeImpl( PersistentContainerVector< G, IndexSet, Vector >& container, const Value& value ) 
        {
          // get size of index set 
          const Size indexSetSize = 
            ((PublicPersistentContainerVector< G, IndexSet, Vector >& ) container).indexSet().size( codimension() );
          // is current size is to small then do a resize, otherwise do nothing
          if( size() < indexSetSize ) 
            resize( value ); 
        }

        // enlarge implementation for persistent containers based on maps
        template < class G, class IdSet, class Map >
        void enlargeImpl( PersistentContainerMap< G, IdSet, Map >& container, const Value& value ) 
        {
          // this needs a revision 
          PersistentContainerMap< G, IdSet, Map > checkSize( container );
          checkSize.resize( value );
          // is current size is to small then do a resize, otherwise do nothing
          if( size() < checkSize.size() ) 
            resize( value );
        }
      };

      // use the imporved PersistentContainer
      typedef IndexPersistentContainer IndexContainerType ;

      // the mapping of the global to leaf index 
      IndexContainerType leafIndex_;
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
      void resize ()
      {
        // enlarge index container, do not shrink, because the old indices are still needed during compress 
        leafIndex_.enlargeOnly( invalidIndex() );
      }

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
        std::fill( indexState_.begin(), indexState_.end(), UNUSED );
      }

      bool consecutive ()
      {
        // Something's wrong here: This method _must_ always return true
        typedef typename IndexContainerType::Iterator Iterator;
        bool consecutive = true;
        const Iterator end = leafIndex_.end();
        for( Iterator it = leafIndex_.begin(); it != end; ++it )
          consecutive &= (*it < IndexType( indexState_.size() ));
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
        const int sizeOfVecs = indexState_.size();
        holes_.resize( sizeOfVecs );

        // true if a least one dof must be copied 
        bool haveToCopy = false;
        
        // mark holes 
        int actHole = 0;
        for( int index = 0; index < sizeOfVecs; ++index )
        {
          // create vector with all holes 
          if( indexState_[ index ] == UNUSED )
            holes_[ actHole++ ] = index;
        }

        // the new size is the actual size minus the holes 
        int actSize = sizeOfVecs - actHole;

        // resize hole storing vectors 
        oldIdx_.resize(actHole);
        newIdx_.resize(actHole);

        // only compress if number of holes > 0    
        if(actHole > 0)
        {
          // close holes 
          int holes = 0; // number of real holes 
          typedef typename IndexContainerType::Iterator Iterator;
          const Iterator end = leafIndex_.end();
          for( Iterator it = leafIndex_.begin(); it != end; ++it )
          {
            // a index that is used but larger then actual size 
            // has to move to a hole 
            if( (*it == invalidIndex()) || (indexState_[ *it ] == UNUSED) )
            {
              // all unused indices are reset to invalidIndex
              *it = invalidIndex();
            }
            else 
            {
              // if used index lies behind size, then index has to move 
              // to one of the holes 
              if( *it >= actSize )
              {
                // serach next hole that is smaler than actual size 
                --actHole;
                // if actHole < 0 then error, because we have index larger then
                // actual size 
                assert(actHole >= 0);
                while ( holes_[actHole] >= actSize )
                {
                  --actHole;
                  if(actHole < 0) break;
                }

                assert(actHole >= 0);

#if HAVE_MPI 
                // only for none-ghost elements hole storage is applied
                // this is because ghost indices might have in introduced 
                // after the resize was done. 
                if( indexState_[ *it ] == USED )
#endif
                {
                  // remember old and new index 
                  oldIdx_[holes] = *it;
                  newIdx_[holes] = holes_[actHole];
                  ++holes;
                }
                
                *it = holes_[actHole];

                // means that dof manager has to copy the mem
                indexState_[ *it ] = NEW;
                haveToCopy = true;
              }
            }
          }

          // this call only sets the size of the vectors 
          oldIdx_.resize(holes);
          newIdx_.resize(holes);
        } // end if actHole > 0  
       
        // store number of actual holes 
        numberHoles_ = oldIdx_.size();

        // adjust size of container to correct value 
        leafIndex_.resize( invalidIndex() );

        // resize vector of index states
        indexState_.resize( actSize );

#ifndef NDEBUG
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
        assert( 0 == EntityType :: codimension );
        assert( checkValidIndex( leafIndex_( entity, subNumber ) ) );
        return leafIndex_( entity, subNumber );
      }
      
      //! return state of index for given hierarchic number  
      template <class EntityType> 
      bool exists ( const EntityType& entity ) const
      {
        assert( myCodim_ == EntityType :: codimension );
        const IndexType &index = leafIndex_[ entity ];
        return (indexState_[ index ] != UNUSED);
      }
     
      template <class EntityType> 
      bool exists ( const EntityType& entity ,
                    const int subNumber ) const 
      {
        assert( 0 == EntityType :: codimension );
        const IndexType &index = leafIndex_( entity, subNumber );
        return (indexState_[ index ] != UNUSED);
      }
     
      //! return number of holes 
      IndexType numberOfHoles () const
      {
        return numberHoles_;
      }

      //! return old index, for dof manager only 
      IndexType oldIndex (int elNum ) const
      {
        assert( numberHoles_ == IndexType(oldIdx_.size()) );
        return oldIdx_[elNum]; 
      }

      //! return new index, for dof manager only returns index 
      IndexType newIndex (int elNum) const
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
          indexState_[ leafIdx ] = NEW;
      }

      // insert element and create index for element number 
      template <class EntityType> 
      void markForRemoval( const EntityType& entity )
      {
        assert( myCodim_ == EntityType :: codimension );
        const IndexType &index = leafIndex_[ entity ];
        if( index != invalidIndex() )
          indexState_[ index ] = UNUSED;
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
          out << "idx: " << leafIdx << "  stat: " << indexState_[ leafIdx ] << std::endl;
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
        indexState_[ index ] = USED;
      }

    public:  
      // write to stream 
      template <class StreamTraits> 
      bool write(OutStreamInterface< StreamTraits >& out) const
      {
        // store current index set size 
        out << indexState_.size();
        
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
        int size;
        in >> size;

        // mark all indices used
        indexState_.resize( size );
        std::fill( indexState_.begin(), indexState_.end(), USED );
        
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
