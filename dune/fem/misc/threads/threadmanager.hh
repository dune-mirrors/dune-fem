#ifndef DUNE_FEM_OMPMANAGER_HH
#define DUNE_FEM_OMPMANAGER_HH

#include <cassert>
#include <cstdlib>

#if defined USE_PTHREADS && HAVE_PTHREAD == 0 
#warning "pthreads were not found!"
#undef USE_PTHREADS 
#endif 

#if defined _OPENMP || defined USE_PTHREADS
#ifndef USE_SMP_PARALLEL
#define USE_SMP_PARALLEL
#endif
#endif

#ifdef USE_PTHREADS 
#include <pthread.h>
#include <map>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif


namespace Dune 
{

  namespace Fem 
  {

#ifdef _OPENMP
    struct ThreadManager 
    {
      /** return maximal number of threads possbile in the current run 
          \note can be set by the OMP_NUM_THREADS environment variable 
                from outside of the code */
      static inline int maxThreads() 
      {
        return omp_get_max_threads();
      }

      //! return number of current threads 
      static inline int currentThreads() 
      {
        return omp_get_num_threads();
      }

      //! return thread number 
      static inline int thread() 
      {
        return omp_get_thread_num();
      }

      //! return true if the current thread is the master thread (i.e. thread 0)
      static inline bool isMaster() 
      {
        return thread() == 0 ;
      }

      //! set maximal number of threads available during run  
      static inline void setMaxNumberThreads( const int numThreads ) 
      {
        omp_set_num_threads( numThreads );
      }

      //! returns true if program is operating on one thread currently
      static inline bool singleThreadMode() 
      {
        return currentThreads() == 1 ;
      }
    }; // end class ThreadManager 
#elif defined USE_PTHREADS 
#warning "ThreadManager: using pthreads"
    struct ThreadManager 
    {
    private:  
      struct Manager ;
      typedef int thread_number_t(Manager&);
      typedef std::map< const pthread_t, int > ThreadIdMap ;
      struct Manager 
      {
        const pthread_t master_;
        ThreadIdMap idmap_;
        thread_number_t *threadNum_;
        int maxThreads_;
        int activeThreads_;

        static inline Manager& instance() 
        {
          static Manager mg ;
          return mg ;
        }

        inline void setThreadNumber(  const pthread_t& theadId, const int threadNum )
        {
          assert( isMaster() );
          // thread number 0 is reserved for the master thread 
          assert( threadNum > 0 );
          idmap_[ theadId ] = threadNum ;
        }

        inline void setMaxThreads( const int maxThreads ) 
        {
          assert( isMaster() );
          maxThreads_ = maxThreads;
        }

        inline void singleThreadMode() 
        {
          // make sure we are in single thread mode 
          assert( isMaster() );
          // set current thread to one 
          activeThreads_ = 1;
          // set pointer to thread number function 
          threadNum_ = & Manager :: singleThreadNumber ;
        }

        inline void multiThreadMode(const int nThreads ) 
        {
          // make sure we are in single thread mode 
          assert( isMaster() );
          assert( threadNum_ == & Manager :: singleThreadNumber  );
          assert( nThreads <= maxThreads_ );
          assert( nThreads > 0 );

          // set number of current threads 
          activeThreads_ = nThreads;

          // set pointer to thread number function 
          threadNum_ = & Manager :: multiThreadNumber ;
        }

        inline int currentThreads() const { return activeThreads_; }
        inline int maxThreads() const { return maxThreads_; }
        inline int thread() { return threadNum_( *this ); }

        bool isMaster() const { return master_ == pthread_self(); }

      private:  
        //! default thread number 
        static inline int singleThreadNumber( Manager& ) { return 0; }

        //! thread number in multi thread mode 
        static inline int multiThreadNumber( Manager& mg ) 
        {
          // get thread id 
          const pthread_t self = pthread_self();
          // the master thread has always the thread number 0 
          if( mg.master_ == self ) 
          {
            return 0 ;
          }
          else 
          {
            // otherwise find mapped thread number 
            ThreadIdMap :: iterator num = mg.idmap_.find( self ) ;
            assert( num != mg.idmap_.end() );
            return (*num).second ;
          }
        }

        //! default constructor 
        Manager() 
          : master_( pthread_self() ), // get id of master thread 
            idmap_(), 
            threadNum_( &Manager :: singleThreadNumber ),
            maxThreads_( 1 ), activeThreads_( 1 )
        {
        }
      };

      static inline Manager& manager() 
      {
        return Manager :: instance();
      }
    public:  
      ///////////////////////////////////////////////////////
      //  begin of pthread specific interface 
      ///////////////////////////////////////////////////////
      //! initialize single thread mode 
      static inline void initSingleThreadMode() 
      {
        manager().singleThreadMode();
      }

      //! initialize multi thread mode 
      static inline void initMultiThreadMode( const int nThreads ) 
      {
        manager().multiThreadMode( nThreads );
      }

      //! set thread number for given thead id 
      static inline void setThreadNumber( const pthread_t& threadId, const int threadNum ) 
      {
        manager().setThreadNumber( threadId, threadNum );
      }

      ///////////////////////////////////////////////////////
      //  INTERFACE 
      ///////////////////////////////////////////////////////
      //! return maximal number of threads possbile in the current run 
      static inline int maxThreads() 
      {
        return manager().maxThreads(); 
      }

      //! return number of current threads 
      static inline int currentThreads() 
      {
        return manager().currentThreads(); 
      }

      //! return thread number 
      static inline int thread() 
      {
        return manager().thread(); 
      }

      //! set maximal number of threads available during run  
      static inline void setMaxNumberThreads( const int numThreads ) 
      {
        manager().setMaxThreads( numThreads );
      }

      //! return true if the current thread is the master thread (i.e. thread 0)
      static inline bool isMaster() 
      {
        return manager().isMaster() ;
      }

      //! returns true if program is operating on one thread currently
      static inline bool singleThreadMode() 
      {
        return currentThreads() == 1 ;
      }
    }; // end class ThreadManager 
#else 
    struct ThreadManager 
    {
      //! return maximal number of threads possbile in the current run 
      static inline int maxThreads() { return 1;  }

      //! return number of current threads 
      static inline int currentThreads() { return 1; }

      //! return thread number 
      static inline int thread() {  return 0; }

      //! return true if the current thread is the master thread (i.e. thread 0)
      static inline bool isMaster() { return true ; }

      //! set maximal number of threads available during run  
      static inline void setMaxNumberThreads( const int numThreads ) { } 

      //! returns true if program is operating on one thread currently
      static inline bool singleThreadMode() { return true ; }
    }; // end class ThreadManager 
#endif

  } // namespace Fem 

} // namespace Dune 

#endif // #ifndef DUNE_FEM_OMPMANAGER_HH
