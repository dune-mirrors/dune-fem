#ifndef DUNE_FEM_THREADITERATOR_HH
#define DUNE_FEM_THREADITERATOR_HH

#include <vector>

#include <dune/fem/misc/threadmanager.hh>

namespace Dune {
  namespace Fem {

    /** \brief Thread iterator */
    template <class DiscreteFunctionSpace>  
    class ThreadIterator
    {
      typedef DiscreteFunctionSpace SpaceType;
    public:  
      typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;
      typedef typename SpaceType :: IteratorType IteratorType ;
      typedef typename IteratorType :: Entity EntityType ;
      typedef typename SpaceType :: GridType :: template Codim<0> ::
        EntityPointer EntityPointer ;
      typedef typename SpaceType :: IndexSetType IndexSetType ;
    protected:  
      const SpaceType& space_ ;
      const IndexSetType& indexSet_;

#ifdef _OPENMP
      int sequence_;
      std::vector< IteratorType > iterators_;
      std::vector< int > threadNum_;
      std::vector< std::vector< int > > threadId_; 
#endif
    public:  
      explicit ThreadIterator( const SpaceType& spc )
        : space_( spc )
        , indexSet_( space_.indexSet() )
#ifdef _OPENMP
        , sequence_( -1 )  
        , iterators_( ThreadManager::maxThreads() + 1 , space_.end() )
        , threadId_( ThreadManager::maxThreads() )
#endif
      {
        update();
      }

      void update() 
      {
#ifdef _OPENMP
        if( sequence_ != space_.sequence() )
        {
          if( omp_get_num_threads() > 1 ) 
          {
            std::cerr << "Don't call ThreadIterator::update in a parallel environment!" << std::endl;
            assert( false );
            abort();
          }

          const int maxThreads = omp_get_max_threads() ;
          IteratorType it = space_.begin(); 
          const IteratorType endit = space_.end();
          if( it == endit ) 
          {
            // update sequence number 
            sequence_ = space_.sequence();
            return ;
          }
          // thread 0 starts at begin 
          iterators_[ 0 ] = it ;

          int size = 0;

          for( IteratorType countit = space_.begin(); countit != endit; ++countit, ++size ) {}
          // resize threads 
          threadNum_.resize( size );

          // here ruse iterator to count 
          const int counter = ((int) size / maxThreads );
          for( int thread = 1; thread <= maxThreads; ++thread ) 
          {
            int i = 0; 
            while( (i < counter) && (it != endit) )
            {
              threadNum_[ indexSet_.index( *it ) ] = thread - 1;
              ++i;
              ++it;
            }
            iterators_[ thread ] = it ;
          }
          iterators_[ maxThreads ] = endit ;

          // update sequence number 
          sequence_ = space_.sequence();
          //for(int i = 0; i<size; ++i ) 
          //  std::cout << threadNum_[ i ] << std::endl;
        }
#endif
      }

      IteratorType begin() const 
      {
#ifdef _OPENMP
        return iterators_[ omp_get_thread_num() ];
#else 
       return space_.begin();
#endif
      }

      IteratorType end() const 
      {
#ifdef _OPENMP
        return iterators_[ omp_get_thread_num() + 1 ];
#else 
       return space_.end();
#endif
      }

      int thread(const EntityType& entity ) const 
      {
#ifdef _OPENMP
        assert( threadNum_.size() > indexSet_.index( entity ) );
        return threadNum_[ indexSet_.index( entity ) ];
#else 
        return 0;
#endif
      }
    };


    /** \brief Thread iterator */
    template <class DiscreteFunctionSpace>  
    class ThreadIteratorPointer
    {
      typedef DiscreteFunctionSpace SpaceType;
    public:  
      typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;
      typedef typename SpaceType :: GridType :: template Codim<0> :: Entity EntityType ;
      typedef typename SpaceType :: GridType :: template Codim<0> ::
        EntityPointer EntityPointer ;
      typedef typename SpaceType :: IndexSetType IndexSetType ;
    protected:  
      const SpaceType& space_ ;
      const IndexSetType& indexSet_;

      class Iterator
      {
        const std::vector< EntityPointer >& vec_;
        size_t pos_ ;
      public:
        typedef typename SpaceType :: GridType :: template Codim<0> :: Entity EntityType ;
        typedef EntityType Entity ;
        Iterator( const std::vector< EntityPointer >& vec, const size_t pos)
          : vec_( vec ), pos_( pos )
        {}

        const EntityType& operator* () const 
        {
          assert( pos_ < vec_.size() );
          return *( vec_[ pos_ ] );
        }

        void operator ++ () { ++pos_; }
        bool operator != (const Iterator & other ) const 
        {
          return pos_ != other.pos_ ;
        }
      };

      int sequence_;
      std::vector< std::vector< EntityPointer > > pointers_;
      std::vector< int > threadNum_;
    public:  
      typedef Iterator IteratorType ;
      explicit ThreadIteratorPointer( const SpaceType& spc )
        : space_( spc )
        , indexSet_( space_.indexSet() )
#ifdef _OPENMP
        , sequence_( -1 )  
        , pointers_( omp_get_max_threads() + 1 )
#endif
      {
        update();
      }

      void update() 
      {
#ifdef _OPENMP
        if( sequence_ != space_.sequence() )
        {
          if( omp_get_num_threads() > 1 ) 
          {
            std::cerr << "Don't call ThreadIterator::update in a parallel environment!" << std::endl;
            assert( false );
            abort();
          }

          const int maxThreads = omp_get_max_threads() ;
          typedef typename SpaceType :: IteratorType SpcIteratorType ;
          SpcIteratorType it = space_.begin(); 
          const SpcIteratorType endit = space_.end();
          if( it == endit ) 
          {
            // update sequence number 
            sequence_ = space_.sequence();
            return ;
          }

          int size = 0;
          for( SpcIteratorType countit = space_.begin(); countit != endit; ++countit, ++size ) {}
          // resize threads 
          threadNum_.resize( size );

          // here ruse iterator to count 
          const int counter = ((int) size / maxThreads );
          for( int thread = 0; thread <= maxThreads; ++thread ) 
          {
            int i = 0; 
            pointers_[ thread ].resize( counter, it );
            while( (i < counter) && (it != endit) )
            {
              pointers_[ thread ][ i ] = it ;
              threadNum_[ indexSet_.index( *it ) ] = thread ;
              ++i;
              ++it;
            }
          }
          // update sequence number 
          sequence_ = space_.sequence();
        }
#endif
      }

      Iterator begin() const 
      {
#ifdef _OPENMP
        return Iterator( pointers_[ omp_get_thread_num() ], 0 );
#else 
        return Iterator( pointers_[ 0 ], 0 );
#endif
      }

      Iterator end() const 
      {
#ifdef _OPENMP
        return Iterator( pointers_[ omp_get_thread_num() ], pointers_[ omp_get_thread_num() ].size() );
#else 
        return Iterator( pointers_[ 0 ], pointers_[ 0 ].size() );
#endif
      }

      int thread(const EntityType& entity ) const 
      {
        assert( ((int) threadNum_.size()) > indexSet_.index( entity ) );
        return threadNum_[ indexSet_.index( entity ) ];
      }
    };
  }
}

#endif