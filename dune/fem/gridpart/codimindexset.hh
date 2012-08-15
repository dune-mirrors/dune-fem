#ifndef DUNE_FEM_CODIMINDEXSET_HH
#define DUNE_FEM_CODIMINDEXSET_HH

#include <algorithm>

//- Dune includes 
#include <dune/common/misc.hh>
#include <dune/fem/space/common/arrays.hh>
#include <dune/fem/gridpart/defaultindexsets.hh>

#include <dune/fem/io/streams/xdrstreams.hh>
#include <dune/grid/utility/persistentcontainer.hh>

#ifdef ENABLE_ADAPTIVELEAFINDEXSET_FOR_YASPGRID
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/sgrid.hh>

namespace Dune 
{

  // PersistentContainer for YaspGrid
  // --------------------------------

  template< int dim, class Data, class Allocator >
  class PersistentContainer< YaspGrid< dim >, Data, Allocator >
  : public PersistentContainerVector< YaspGrid< dim >, 
                                      typename YaspGrid< dim >::LeafIndexSet,
                                      std::vector<Data,Allocator> >
  {
    typedef YaspGrid< dim > Grid;
    typedef PersistentContainerVector< Grid, typename Grid::LeafIndexSet, std::vector<Data,Allocator> > BaseType;

  public:
    //! Constructor filling the container with values using the default constructor 
    //! Depending on the implementation this could be achieved without allocating memory
    PersistentContainer ( const Grid &grid, const int codim, const Allocator &allocator = Allocator() )
    : BaseType( grid, codim, grid.leafIndexSet(), 1.0, allocator )
    {}
  };

  // PersistentContainer for SGrid
  // -------------------------------

  template< int dim, int dimworld, class ctype, class Data, class Allocator >
  class PersistentContainer< SGrid< dim, dimworld, ctype >, Data, Allocator >
  : public PersistentContainerVector< SGrid< dim, dimworld, ctype >, 
                                      typename SGrid< dim, dimworld, ctype >::LeafIndexSet,
                                      std::vector<Data,Allocator> >
  {
    typedef SGrid< dim, dimworld, ctype > Grid ;
    typedef PersistentContainerVector< Grid, typename Grid::LeafIndexSet, std::vector<Data,Allocator> > BaseType;

  public:
    //! Constructor filling the container with values using the default constructor 
    //! Depending on the implementation this could be achieved without allocating memory
    PersistentContainer ( const Grid &grid, const int codim, const Allocator &allocator = Allocator() )
    : BaseType( grid, codim, grid.leafIndexSet(), 1.0, allocator )
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

      // Index pair that is stored, derive from pair to overload 
      // default constructor which sets the correct default data 
      struct IndexPair : public std::pair< IndexType, INDEXSTATE > 
      {
        typedef std::pair< IndexType, INDEXSTATE > BaseType;
        // default constructor 
        IndexPair() 
          : BaseType( ThisType::invalidIndex(), UNUSED ) {}
      }; 


      // array type for indices 
      typedef MutableArray< IndexType > IndexArrayType;

      typedef PersistentContainer< GridImp, IndexPair > IndexContainerType;

      // the mapping of the global to leaf index 
      IndexContainerType leafIndex_;

      // stack for holes 
      IndexArrayType holes_; 
     
      // Array that only remeber the occuring 
      // holes (for compress of data)
      IndexArrayType oldIdx_; 
      IndexArrayType newIdx_; 
     
      // next index to give away 
      IndexType nextFreeIndex_;

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
        , leafIndex_(grid, codim)
        , holes_(0)
        , oldIdx_(0)
        , newIdx_(0)
        , nextFreeIndex_ (0)
        , lastSize_ (0)
        , myCodim_( codim ) 
        , numberHoles_(0)
      {
        setMemoryFactor(memoryFactor);
      }

      //! set memory overestimation factor 
      void setMemoryFactor(const double memoryFactor)
      {
        holes_.setMemoryFactor(memoryFactor);
        oldIdx_.setMemoryFactor(memoryFactor);
        newIdx_.setMemoryFactor(memoryFactor);
      }

      //! reallocate the vectors
      void resize ()
      {
        leafIndex_.reserve();
      }

      //! prepare for setup (nothing to do here)
      void prepareCompress ()
      {
      }

    public:  
      //! clear set 
      void clear() 
      {
        // set all values to invalidIndex  
        std::fill( leafIndex_.begin(), leafIndex_.end(), IndexPair() );
        // reset next free index 
        nextFreeIndex_ = 0;
      }

      //! set all entries to unused 
      void resetUsed() 
      {
        typedef typename IndexContainerType :: Iterator Iterator;
        const Iterator endit = leafIndex_.end();
        for( Iterator it = leafIndex_.begin(); it != endit; ++it )
        {
          (*it).second = UNUSED;
        }
      }

      //! set all entries to unused 
      void checkConsecutive() 
      {
        typedef typename IndexContainerType :: Iterator Iterator;
        const Iterator endit = leafIndex_.end();
        for( Iterator it = leafIndex_.begin(); it != endit; ++it )
        {
          const int idx = (*it).first; 
          assert( idx < nextFreeIndex_ );
        }
      }

      //! clear holes, i.e. set number of holes to zero 
      void clearHoles() 
      {
        // set number of holes to zero 
        numberHoles_ = 0;
        // remember actual size 
        lastSize_ = nextFreeIndex_;
      }

      //! make to index numbers consecutive 
      //! return true, if at least one hole was closed 
      bool compress ()
      {
        const int sizeOfVecs = leafIndex_.size();
        holes_.resize( sizeOfVecs );

        // true if a least one dof must be copied 
        bool haveToCopy = false;
        
        // mark holes 
        int actHole = 0;
        int newActSize = 0;
        typedef typename IndexContainerType :: Iterator Iterator;
        const Iterator endit = leafIndex_.end();
        for( Iterator it = leafIndex_.begin(); it != endit; ++it )
        {
          const IndexPair& leafIdx = *it;
          if( leafIdx.first >= 0 )
          {
            // create vector with all holes 
            if( leafIdx.second == UNUSED )
            {
              holes_[actHole] = leafIdx.first;
              ++actHole;
            }

            // count the size of the leaf indices 
            ++newActSize;
          }
        }

        assert( newActSize >= actHole );
        // the new size is the actual size minus the holes 
        int actSize = newActSize - actHole;

        // resize hole storing vectors 
        oldIdx_.resize(actHole);
        newIdx_.resize(actHole);

        // only compress if number of holes > 0    
        if(actHole > 0)
        {
          // close holes 
          //
          // NOTE: here the holes closing should be done in 
          // the opposite way. future work. 
          int holes = 0; // number of real holes 
          //size_t i = 0;
          const Iterator endit = leafIndex_.end();
          for( Iterator it = leafIndex_.begin(); it != endit; ++it )
          {
            IndexPair& leafIdx = *it;
            // a index that is used but larger then actual size 
            // has to move to a hole 
            if( leafIdx.second == UNUSED) 
            {
              // all unused indices are reset to invalidIndex  
              leafIdx.first = invalidIndex() ;
            }
            else 
            {
              // if used index lies behind size, then index has to move 
              // to one of the holes 
              if(leafIdx.first >= actSize)
              {
                // serach next hole that is smaler than actual size 
                actHole--;
                // if actHole < 0 then error, because we have index larger then
                // actual size 
                assert(actHole >= 0);
                while ( holes_[actHole] >= actSize )
                {
                  actHole--;
                  if(actHole < 0) break;
                }

                assert(actHole >= 0);

#if HAVE_MPI 
                // only for none-ghost elements hole storage is applied
                // this is because ghost indices might have in introduced 
                // after the resize was done. 
                if( leafIdx.second == USED ) 
#endif
                {
                  // remember old and new index 
                  oldIdx_[holes] = leafIdx.first; 
                  newIdx_[holes] = holes_[actHole];
                  ++holes;
                }
                
                leafIdx.first = holes_[actHole];

                // means that dof manager has to copy the mem
                leafIdx.second = NEW;
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

        // adjust size 
        leafIndex_.update();
        
        // the next index that can be given away is equal to size
        nextFreeIndex_ = actSize;

#ifndef NDEBUG
        checkConsecutive();
#endif

        return haveToCopy;
      }

      //! return how much extra memory is needed for restriction 
      IndexType additionalSizeEstimate () const { return nextFreeIndex_; }

      //! return size of grid entities per level and codim 
      IndexType size () const
      {
        return nextFreeIndex_;
      }
      
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
        assert( checkValidIndex( leafIndex_[ entity ].first ) );
        return leafIndex_[ entity ].first;
      }
      
      //! return leaf index for given entity   
      template <class EntityType>
      IndexType subIndex ( const EntityType& entity,
                           const int subNumber ) const 
      {
        assert( 0 == EntityType :: codimension );
        assert( checkValidIndex( leafIndex_( entity, subNumber ).first ) );
        return leafIndex_( entity, subNumber ).first;
      }
      
      //! return state of index for given hierarchic number  
      template <class EntityType> 
      bool exists ( const EntityType& entity ) const
      {
        assert( myCodim_ == EntityType :: codimension );
        return leafIndex_[ entity ].second != UNUSED;
      }
     
      template <class EntityType> 
      bool exists ( const EntityType& entity ,
                    const int subNumber ) const 
      {
        assert( 0 == EntityType :: codimension );
        return leafIndex_( entity, subNumber).second != UNUSED;
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
        IndexPair& leafIdx = leafIndex_[ entity ];
        insertIdx( leafIdx );

        // if index is also larger than lastSize
        // mark as new to skip old-new index lists 
        if( leafIdx.first >= lastSize_ ) 
        {
          leafIdx.second = NEW;
        }
      }

      // insert element and create index for element number 
      template <class EntityType> 
      void markForRemoval( const EntityType& entity )
      {
        assert( myCodim_ == EntityType :: codimension );
        leafIndex_[ entity ].second = UNUSED;
      }

      // insert element as ghost and create index for element number 
      template <class EntityType> 
      bool validIndex (const EntityType& entity ) const
      {
        assert( myCodim_ == EntityType :: codimension );
        return leafIndex_[ entity ].first >= 0; 
      }

      void print( std::ostream& out ) const 
      {
        typedef typename IndexContainerType :: ConstIterator Iterator;
        const Iterator endit = leafIndex_.end();
        for( Iterator it = leafIndex_.begin(); it != endit; ++it )
        {
          const IndexPair& leafIdx = *it;
          out << "idx: " << leafIdx.first << "  stat: " << leafIdx.second << std::endl;
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
      void insertIdx ( IndexPair& leafIdx )
      {
        if( leafIdx.first == invalidIndex() )
          leafIdx.first = nextFreeIndex_ ++ ;

        leafIdx.second = USED;
      }

    public:  
      // write to stream 
      template <class StreamTraits> 
      bool write(OutStreamInterface< StreamTraits >& out) const
      {
        // store current index set size 
        out << nextFreeIndex_ ;
        
        // for consistency checking, write size as 64bit integer
        const uint64_t mysize = leafIndex_.size();
        out << mysize ;

        // backup indices 
        typedef typename IndexContainerType :: ConstIterator ConstIterator;
        const ConstIterator endit = leafIndex_.end();
        for( ConstIterator it = leafIndex_.begin(); it != endit; ++it )
        {
          out << (*it).first ;
        }

        return true;
      }
      
      // read from stream 
      template <class StreamTraits> 
      bool read(InStreamInterface< StreamTraits >& in)
      {
        // read current index set size 
        in >> nextFreeIndex_ ;
        
        // for consistency checking 
        uint64_t storedSize = 0;
        in >> storedSize ;

        uint64_t leafsize = leafIndex_.size() ;
        // the stored size can be larger (visualization of parallel grids in serial)
        if( storedSize < leafsize ) 
        {
          DUNE_THROW(InvalidStateException,"CodimIndexSet: size consistency check failed during restore!"); 
        }

        // restore indices  
        typedef typename IndexContainerType :: Iterator Iterator;
        const Iterator endit = leafIndex_.end();
        uint64_t count = 0 ;
        for( Iterator it = leafIndex_.begin(); it != endit; ++it, ++count )
        {
          in >> (*it).first ;
        }

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
