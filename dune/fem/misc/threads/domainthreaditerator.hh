#ifndef DUNE_FEM_DOMAINTHREADITERATOR_HH
#define DUNE_FEM_DOMAINTHREADITERATOR_HH

#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/gridpart/filteredgridpart.hh>
#include <dune/fem/gridpart/filter/domainfilter.hh>
#include <dune/fem/misc/threads/threaditeratorstorage.hh>

#ifdef USE_SMP_PARALLEL
#if HAVE_DUNE_ALUGRID
#define USE_THREADPARTITIONER
#include <dune/fem/misc/threads/threadpartitioner.hh>
#endif
#endif

namespace Dune {

  namespace Fem {

    /** \brief Thread iterators using domain decomposition */
    template <class GridPart>
    class DomainDecomposedIterator
    {
    public:
      typedef GridPart  GridPartType;
      typedef typename GridPartType :: GridType      GridType;
      typedef typename GridPartType :: IndexSetType  IndexSetType;

      typedef DomainFilter<GridPartType> FilterType;
#ifdef USE_THREADPARTITIONER
      typedef FilteredGridPart< GridPartType, FilterType, false > FilteredGridPartType ;

      typedef typename FilteredGridPartType :: template Codim< 0 > :: IteratorType
        IteratorType;
#else
      typedef typename GridPartType :: template Codim< 0 > :: IteratorType IteratorType ;
#endif

      typedef typename IteratorType :: Entity EntityType ;

      typedef DofManager< GridType > DofManagerType;

      static const PartitionIteratorType pitype = GridPartType :: indexSetPartitionType ;

#ifdef USE_THREADPARTITIONER
      typedef ThreadPartitioner< GridPartType >           ThreadPartitionerType;
      typedef typename ThreadPartitionerType :: Method    PartitioningMethodType ;
#endif // #ifdef USE_SMP_PARALLEL

    protected:
      const GridPartType& gridPart_;
      const DofManagerType& dofManager_;
      const IndexSetType& indexSet_;

#ifdef USE_THREADPARTITIONER
      int sequence_;
      int numThreads_;
      std::vector< std::unique_ptr< FilteredGridPartType > > filteredGridParts_;

      typedef typename FilterType :: DomainArrayType ThreadArrayType;
      ThreadArrayType threadNum_;
#endif
      // ratio of computing time needed by the master thread compared to the other threads
      double masterRatio_ ;

#ifdef USE_THREADPARTITIONER
      const PartitioningMethodType method_;
#endif // #ifdef USE_SMP_PARALLEL

      // if true, thread 0 does only communication and no computation
      const bool communicationThread_;
      const bool verbose_ ;

    protected:
#ifdef USE_THREADPARTITIONER
      static PartitioningMethodType getMethod ( const ParameterReader &parameter )
      {
        // default is recursive
        const std::string methodNames[] = { "recursive", "kway", "sfc" };
        return (PartitioningMethodType ) parameter.getEnum("fem.threads.partitioningmethod", methodNames, 0 );
      }
#endif // #ifdef USE_SMP_PARALLEL

    public:
      //! contructor creating thread iterators
      explicit DomainDecomposedIterator( const GridPartType& gridPart, const ParameterReader &parameter = Parameter::container() )
        : gridPart_( gridPart ),
          dofManager_( DofManagerType :: instance( gridPart_.grid() ) ),
          indexSet_( gridPart_.indexSet() )
#ifdef USE_THREADPARTITIONER
        , sequence_( -1 )
        , numThreads_( Fem :: MPIManager :: numThreads() )
        , filteredGridParts_( Fem :: MPIManager :: maxThreads() )
#endif
        , masterRatio_( 1.0 )
#ifdef USE_THREADPARTITIONER
        , method_( getMethod( parameter ) )
#endif // #ifdef USE_SMP_PARALLEL
        , communicationThread_( parameter.getValue<bool>("fem.threads.communicationthread", false)
                    &&  Fem :: MPIManager :: maxThreads() > 1 ) // only possible if maxThreads > 1
        , verbose_( Parameter::verbose() &&
                    parameter.getValue<bool>("fem.threads.verbose", false ) )
      {
#ifdef USE_THREADPARTITIONER
        for(int thread=0; thread < Fem :: MPIManager :: maxThreads(); ++thread )
        {
          // thread is the thread number of this filter
          filteredGridParts_[ thread ].reset(
              new FilteredGridPartType( const_cast<GridPartType &> (gridPart_),
                                        FilterType( gridPart_, threadNum_, thread ) ) );
        }

        threadNum_.setMemoryFactor( 1.1 );
#endif
        update();
      }

#ifdef USE_THREADPARTITIONER
      //! return filter for given thread
      const FilterType& filter( const int thread ) const
      {
        return filteredGridParts_[ thread ]->filter();
      }
#endif // USE_SMP_PARALLEL

      //! update internal list of iterators
      void update()
      {
#ifdef USE_THREADPARTITIONER
        const int sequence = dofManager_.sequence() ;
        // if grid got updated also update iterators
        if( sequence_ != sequence || numThreads_ != MPIManager :: numThreads() )
        {
          if( ! MPIManager :: singleThreadMode() )
          {
            std::cerr << "Don't call ThreadIterator::update in a parallel environment!" << std::endl;
            assert( false );
            abort();
          }

          const int commThread = communicationThread_ ? 1 : 0;
          // get number of partitions possible
          const size_t partitions = MPIManager :: numThreads() - commThread ;

          // create partitioner
          ThreadPartitionerType db( gridPart_, partitions, masterRatio_ );
          // do partitioning
          db.serialPartition( method_ );

          // get end iterator
          typedef typename GridPartType :: template Codim< 0 > :: IteratorType GPIteratorType;
          const GPIteratorType endit = gridPart_.template end< 0 >();

          // get size for index set
          const size_t size = indexSet_.size( 0 );

          // resize threads storage
          threadNum_.resize( size );
          // set all values to default value
          for(size_t i = 0; i<size; ++i) threadNum_[ i ] = -1;

          {
            // just for diagnostics
            std::vector< int > counter( partitions+commThread , 0 );

            int numInteriorElems = 0;
            for(GPIteratorType it = gridPart_.template begin< 0 >();
                it != endit; ++it, ++numInteriorElems )
            {
              const EntityType& entity  = * it;
              const int rank = db.getRank( entity ) + commThread ;
              assert( rank >= 0 );
              //std::cout << "Got rank = " << rank << "\n";
              threadNum_[ indexSet_.index( entity ) ] = rank ;
              ++counter[ rank ];
            }

            // update sequence number
            sequence_ = sequence;

            // update numThreads_
            numThreads_ = MPIManager :: numThreads();

            if( verbose_ )
            {
              std::cout << "DomainDecomposedIterator: sequence = " << sequence_ << " size = " << numInteriorElems << std::endl;
              const size_t counterSize = counter.size();
              for(size_t i = 0; i<counterSize; ++i )
                std::cout << "DomainDecomposedIterator: T[" << i << "] = " << counter[ i ] << std::endl;
            }

#ifndef NDEBUG
            // check that all interior elements have got a valid thread number
            for(GPIteratorType it = gridPart_.template begin< 0 >(); it != endit; ++it )
            {
              assert( threadNum_[ indexSet_.index( *it ) ] >= 0 );
            }
#endif
            //for(size_t i = 0; i<size; ++i )
            //{
            //  //std::cout << threadNum_[ i ] << std::endl;
            //}
          }
        }
#endif
      }

      //! return begin iterator for current thread
      IteratorType begin() const
      {
#ifdef USE_THREADPARTITIONER
        if( ! MPIManager :: singleThreadMode () )
        {
          const int thread = MPIManager :: thread() ;
          if( communicationThread_ && thread == 0 )
            return filteredGridParts_[ thread ]->template end< 0 > ();
          else
            return filteredGridParts_[ thread ]->template begin< 0 > ();
        }
        else
#endif
        {
          return gridPart_.template begin< 0 > ();
        }
      }

      //! return end iterator for current thread
      IteratorType end() const
      {
#ifdef USE_THREADPARTITIONER
        if( ! MPIManager :: singleThreadMode () )
        {
          return filteredGridParts_[ MPIManager :: thread() ]->template end< 0 > ();
        }
        else
#endif
        {
          return gridPart_.template end< 0 > ();
        }
      }

      //! return thread number this entity belongs to
      int index(const EntityType& entity ) const
      {
        return indexSet_.index( entity );
      }

      //! return thread number this entity belongs to
      int thread(const EntityType& entity ) const
      {
#ifdef USE_THREADPARTITIONER
        assert( (int)threadNum_.size() > (int) indexSet_.index( entity ) );
        // NOTE: this number can also be negative for ghost elements or elements
        // that do not belong to the set covered by the space iterators
        return threadNum_[ indexSet_.index( entity ) ];
#else
        return 0;
#endif
      }

      void setMasterRatio( const double ratio )
      {
        masterRatio_ = 0.5 * (ratio + masterRatio_);
      }
    };


    /** \brief Storage of thread iterators using domain decomposition */
    template <class GridPart>
    class DomainDecomposedIteratorStorage
      : public ThreadIteratorStorageBase< DomainDecomposedIterator< GridPart > >
    {
      typedef ThreadIteratorStorageBase< DomainDecomposedIterator< GridPart > > BaseType ;
    public:
      DomainDecomposedIteratorStorage( const GridPart& gridPart )
        : BaseType( gridPart )
      {}
    };

  } // end namespace Fem

} // end namespace Dune

#endif // #ifndef DUNE_FEM_DG_DOMAINTHREADITERATOR_HH
