#ifndef DUNE_FEM_THREADITERATOR_HH
#define DUNE_FEM_THREADITERATOR_HH

#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/fem/gridpart/filter/domainfilter.hh>
#include <dune/fem/misc/threads/threadmanager.hh>
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

#ifdef USE_SMP_PARALLEL
      int sequence_;
      std::vector< IteratorType > iterators_;
      DynamicArray< int > threadNum_;
      std::vector< std::vector< int > > threadId_;
      std::vector< FilterType* > filters_;
#endif

      // if true, thread 0 does only communication and no computation
      const bool communicationThread_;
      const bool verbose_ ;

    public:
      //! contructor creating thread iterators
      explicit ThreadIterator( const GridPartType& gridPart, const ParameterReader &parameter = Parameter::container() )
        : gridPart_( gridPart )
        , dofManager_( DofManagerType :: instance( gridPart_.grid() ) )
        , indexSet_( gridPart_.indexSet() )
#ifdef USE_SMP_PARALLEL
        , sequence_( -1 )
        , iterators_( ThreadManager::maxThreads() + 1 , gridPart_.template end< 0, pitype >() )
        , threadId_( ThreadManager::maxThreads() )
#endif
        , communicationThread_( parameter.getValue<bool>("fem.threads.communicationthread", false)
                    &&  Fem :: ThreadManager :: maxThreads() > 1 ) // only possible if maxThreads > 1
        , verbose_( Parameter::verbose() &&
                    parameter.getValue<bool>("fem.threads.verbose", false ) )
      {
#ifdef USE_SMP_PARALLEL
        threadNum_.setMemoryFactor( 1.1 );
        filters_.resize( Fem :: ThreadManager :: maxThreads(), (FilterType *) 0 );
        for(int thread=0; thread < Fem :: ThreadManager :: maxThreads(); ++thread )
        {
          filters_[ thread ] = new FilterType( gridPart_, threadNum_, thread );
        }
#endif
        update();
      }

#ifdef USE_SMP_PARALLEL
      ~ThreadIterator()
      {
        for(size_t i = 0; i<filters_.size(); ++i )
        {
          delete filters_[ i ];
        }
      }

      //! return filter for given thread
      const FilterType& filter( const unsigned int thread ) const
      {
        assert( thread < filters_.size() );
        return *(filters_[ thread ]);
      }
#endif

      //! update internal list of iterators
      void update()
      {
#ifdef USE_SMP_PARALLEL
        const int sequence = gridPart_.sequence();
        // if grid got updated also update iterators
        if( sequence_ != sequence )
        {
          if( ! ThreadManager :: singleThreadMode() )
          {
            std::cerr << "Don't call ThreadIterator::update in a parallel environment!" << std::endl;
            assert( false );
            abort();
          }

          const size_t maxThreads = ThreadManager :: maxThreads() ;

          // get end iterator
          const IteratorType endit = gridPart_.template end< 0, pitype >();
          IteratorType it = gridPart_.template begin< 0, pitype >();
          if( it == endit )
          {
            // set all iterators to end iterators
            for( size_t thread = 0; thread <= maxThreads; ++thread )
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
          const size_t roundOff = (iterSize % maxThreads);
          const size_t counterBase = ((size_t) iterSize / maxThreads );

          // just for diagnostics
          std::vector< int > nElems( maxThreads, 0 );

          for( size_t thread = 1; thread <= maxThreads; ++thread )
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
          iterators_[ maxThreads ] = endit ;

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
#endif
      }

      //! return begin iterator for current thread
      IteratorType begin() const
      {
#ifdef USE_SMP_PARALLEL
        return iterators_[ ThreadManager :: thread() ];
#else
        return gridPart_.template begin< 0, pitype >();
#endif
      }

      //! return end iterator for current thread
      IteratorType end() const
      {
#ifdef USE_SMP_PARALLEL
        return iterators_[ ThreadManager :: thread() + 1 ];
#else
        return gridPart_.template end< 0, pitype >();
#endif
      }

      //! return thread number this entity belongs to
      int index( const EntityType& entity ) const
      {
        return indexSet_.index( entity );
      }

      //! return thread number this entity belongs to
      int thread( const EntityType& entity ) const
      {
#ifdef USE_SMP_PARALLEL
        assert( std::size_t( threadNum_.size() ) > std::size_t( indexSet_.index( entity ) ) );
        // NOTE: this number can also be negative for ghost elements or elements
        // that do not belong to the set covered by the space iterators
        return threadNum_[ indexSet_.index( entity ) ];
#else
        return 0;
#endif
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

#ifdef USE_SMP_PARALLEL
      // check that we have a non-overlapping iterator decomposition
      void checkConsistency( const size_t totalElements )
      {
#ifndef NDEBUG
        const int maxThreads = ThreadManager :: maxThreads() ;
        std::set< int > indices ;
        for( int thread = 0; thread < maxThreads; ++ thread )
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
#endif
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
