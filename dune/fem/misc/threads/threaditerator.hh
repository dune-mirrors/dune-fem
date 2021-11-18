#ifndef DUNE_FEM_THREADITERATOR_HH
#define DUNE_FEM_THREADITERATOR_HH

#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/fem/gridpart/filter/domainfilter.hh>
#include <dune/fem/misc/threads/threaditeratorstorage.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/storage/dynamicarray.hh>

namespace Dune
{

  namespace Fem
  {

    /** \brief Thread iterators */
    template <class GridPart, PartitionIteratorType ptype = InteriorBorder_Partition >
    class ThreadIterator
    {
      ThreadIterator( const ThreadIterator& );
      ThreadIterator& operator= ( const ThreadIterator& );
    public:
      // partition type of iterators used
      static const PartitionIteratorType pitype = ptype ;

      typedef GridPart GridPartType;
      typedef typename GridPartType :: GridType  GridType;
      typedef typename GridPartType :: template Codim< 0 > :: template Partition< pitype > :: IteratorType       IteratorType ;
      typedef typename GridPartType :: template Codim< 0 > :: EntityType         EntityType ;
      typedef typename GridPartType :: IndexSetType IndexSetType ;
      typedef DofManager< GridType > DofManagerType;

      typedef DomainFilter<GridPartType> FilterType;


    protected:
      const GridPartType& gridPart_;
      const DofManagerType& dofManager_;
      const IndexSetType& indexSet_;

      int sequence_;
      int numThreads_;

      std::vector< IteratorType > iterators_;
      DynamicArray< int > threadNum_;
      std::vector< std::vector< int > > threadId_;
      std::vector< std::unique_ptr< FilterType > > filters_;

      // if true, thread 0 does only communication and no computation
      const bool communicationThread_;
      const bool verbose_ ;

    public:
      //! contructor creating thread iterators
      explicit ThreadIterator( const GridPartType& gridPart, const ParameterReader &parameter = Parameter::container() )
        : gridPart_( gridPart )
        , dofManager_( DofManagerType :: instance( gridPart_.grid() ) )
        , indexSet_( gridPart_.indexSet() )
        , sequence_( -1 )
        , numThreads_( MPIManager::numThreads() )
        , iterators_( MPIManager::maxThreads() + 1 , gridPart_.template end< 0, pitype >() )
        , threadId_( MPIManager::maxThreads() )
        , filters_( MPIManager::maxThreads() )
        , communicationThread_( parameter.getValue<bool>("fem.threads.communicationthread", false)
                    &&  Fem :: MPIManager :: maxThreads() > 1 ) // only possible if maxThreads > 1
        , verbose_( Parameter::verbose() &&
                    parameter.getValue<bool>("fem.threads.verbose", false ) )
      {
        threadNum_.setMemoryFactor( 1.1 );
        for(int thread=0; thread < Fem :: MPIManager :: maxThreads(); ++thread )
        {
          filters_[ thread ].reset( new FilterType( gridPart_, threadNum_, thread ) );
        }
        update();
      }

      //! return filter for given thread
      const FilterType& filter( const unsigned int thread ) const
      {
        assert( thread < filters_.size() );
        return *(filters_[ thread ]);
      }

      //! update internal list of iterators
      void update()
      {
        const int sequence = dofManager_.sequence();
        // if grid got updated also update iterators
        if( sequence_ != sequence || numThreads_ != MPIManager :: numThreads() )
        {
          if( ! MPIManager :: singleThreadMode() )
          {
            std::cerr << "Don't call ThreadIterator::update in a parallel environment!" << std::endl;
            assert( false );
            abort();
          }

          // update currently used thread numbers
          numThreads_  = MPIManager :: numThreads() ;
          const size_t numThreads = numThreads_;

          // get end iterator
          const IteratorType endit = gridPart_.template end< 0, pitype >();

          // pass default value to resize to initialize all iterators
          iterators_.resize( numThreads+1, endit );

          IteratorType it = gridPart_.template begin< 0, pitype >();
          if( it == endit )
          {
            // set all iterators to end iterators
            for( size_t thread = 0; thread <= numThreads; ++thread )
              iterators_[ thread ] = endit ;

            // free memory here
            threadNum_.resize( 0 );

            // update sequence number
            sequence_ = sequence;
            return ;
          }

          // thread 0 starts at begin
          iterators_[ 0 ] = it ;

          // get size for index set (this only works well when pitype == All_Partition)
          // otherwise element have to be counted
          const size_t iterSize = countElements( it, endit );
          const size_t size = indexSet_.size( 0 );

          // resize threads storage
          threadNum_.resize( size );
          // set all values to default value
          for(size_t i = 0; i<size; ++i) threadNum_[ i ] = -1;

          // here use iterator to count
          size_t checkSize = 0;
          const size_t roundOff = (iterSize % numThreads);
          const size_t counterBase = ((size_t) iterSize / numThreads );

          // just for diagnostics
          std::vector< int > nElems( numThreads, 0 );

          for( size_t thread = 1; thread <= numThreads; ++thread )
          {
            size_t i = 0;
            const size_t counter = counterBase + (( (thread-1) < roundOff ) ? 1 : 0);
            nElems[ thread-1 ] = counter ;
            checkSize += counter ;
            //std::cout << counter << " for thread " << thread-1 << std::endl;
            while( (i < counter) && (it != endit) )
            {
              const EntityType &entity = *it;
              assert( std::size_t( indexSet_.index( entity ) ) < std::size_t( threadNum_.size() ) );
              threadNum_[ indexSet_.index( entity ) ] = thread - 1;
              ++i;
              ++it;
            }
            iterators_[ thread ] = it ;
          }
          iterators_[ numThreads ] = endit ;

          if( checkSize != iterSize )
          {
            assert( checkSize == iterSize );
            DUNE_THROW(InvalidStateException,"Partitioning inconsistent!");
          }

          // update sequence number
          sequence_ = sequence;

          if( verbose_ )
          {
            std::cout << "ThreadIterator: sequence = " << sequence_ << " size = " << checkSize << std::endl;
            const size_t counterSize = nElems.size();
            for(size_t i = 0; i<counterSize; ++i )
              std::cout << "ThreadIterator: T[" << i << "] = " << nElems[ i ] << std::endl;
          }

          checkConsistency( iterSize );

          //for(size_t i = 0; i<size; ++i )
          //  std::cout << threadNum_[ i ] << std::endl;
        }
      }

      //! return begin iterator for current thread
      IteratorType begin() const
      {
        if( MPIManager :: singleThreadMode() )
        {
          return gridPart_.template begin< 0, pitype >();
        }
        // in multi thread mode return iterators for each thread
        else
        {
          assert( MPIManager :: thread() < numThreads_ );
          return iterators_[ MPIManager :: thread() ];
        }
      }
      IteratorType begin(int thread) const
      {
        return iterators_[ thread ];
      }

      //! return end iterator for current thread
      IteratorType end() const
      {
        if( MPIManager :: singleThreadMode() )
        {
          return gridPart_.template end< 0, pitype >();
        }
        // in multi thread mode return iterators for each thread
        else
        {
          assert( MPIManager :: thread() < numThreads_ );
          return iterators_[ MPIManager :: thread() + 1 ];
        }
      }
      IteratorType end(int thread) const
      {
        return iterators_[ thread + 1 ];
      }

      //! return thread number this entity belongs to
      int index( const EntityType& entity ) const
      {
        return indexSet_.index( entity );
      }

      int threadParallel( const EntityType& entity ) const
      {
        assert( std::size_t( threadNum_.size() ) > std::size_t( indexSet_.index( entity ) ) );
        // NOTE: this number can also be negative for ghost elements or elements
        // that do not belong to the set covered by the space iterators
        return threadNum_[ indexSet_.index( entity ) ];
      }
      //! return thread number this entity belongs to
      int thread( const EntityType& entity ) const
      {
        if( MPIManager::singleThreadMode() )
          return 0;
        else
          return threadParallel(entity);
      }

      //! set ratio between master thread and other threads in comp time
      void setMasterRatio( const double ratio )
      {
      }

    protected:
      template < class Iterator >
      size_t countElements( const Iterator& begin, const Iterator& end ) const
      {
        size_t count = 0;
        for( Iterator it = begin; it != end; ++ it )
          ++count ;
        return count ;
      }

      // check that we have a non-overlapping iterator decomposition
      void checkConsistency( const size_t totalElements )
      {
#ifndef NDEBUG
        const int numThreads = MPIManager :: numThreads() ;
        std::set< int > indices ;
        for( int thread = 0; thread < numThreads; ++ thread )
        {
          const IteratorType end = iterators_[ thread+1 ];
          for( IteratorType it = iterators_[ thread ]; it != end; ++it )
          {
            const int idx = gridPart_.indexSet().index( *it );
            assert( indices.find( idx ) == indices.end() ) ;
            indices.insert( idx );
          }
        }
        assert( indices.size() == totalElements );
#endif
      }
    };

    /** \brief Storage of thread iterators */
    template <class GridPart, PartitionIteratorType pitype = InteriorBorder_Partition >
    class ThreadIteratorStorage
      : public ThreadIteratorStorageBase< ThreadIterator< GridPart, pitype > >
    {
      typedef ThreadIteratorStorageBase< ThreadIterator< GridPart, pitype > > BaseType ;
    public:
      ThreadIteratorStorage( const GridPart& gridPart )
        : BaseType( gridPart )
      {}
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_THREADITERATOR_HH
