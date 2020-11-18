#ifndef DUNE_FEM_OMPMANAGER_HH
#define DUNE_FEM_OMPMANAGER_HH

#include <cassert>
#include <cstdlib>

#ifdef USE_PTHREADS
#if HAVE_PTHREAD == 0
#warning "pthreads were not found!"
#undef USE_PTHREADS
#endif
#endif

#if defined _OPENMP || defined(USE_PTHREADS)
#ifndef USE_SMP_PARALLEL
#define USE_SMP_PARALLEL
#endif
#endif

#if HAVE_PTHREAD
#include <pthread.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include <dune/common/visibility.hh>

namespace Dune
{

  namespace Fem
  {

    /** \class ThreadManager
     *  \ingroup Utility
     *  \brief The ThreadManager wrapps basic shared memory functionality
     *         provided by OpenMP or pthreads such as thread id, number of threads, etc.
     *
     *  \note All methods are static members.
     */
    struct EmptyThreadManager
    {
      //! true if pthreads are used
      static constexpr bool pthreads = false ;

      //! \brief initialize single thread mode (when in multithread mode)
      static inline void initSingleThreadMode() {}

      //! \brief initialize multi thread mode (when in single thread mode)
      static inline void initMultiThreadMode( const int nThreads ) {}

      //! \brief set max number of threads and thread number for this thread
      static inline void initThread( const int maxThreads, const int threadNum ) {}

      //! \brief return maximal number of threads possbile in the current run
      static inline int maxThreads() { return 1;  }

      //! \brief return number of current threads
      static inline int currentThreads() { return 1; }

      //! \brief return thread number
      static inline int thread() {  return 0; }

      //! \brief return true if the current thread is the master thread (i.e. thread 0)
      static inline bool isMaster() { return true ; }

      //! \brief set maximal number of threads available during run
      static inline void setMaxNumberThreads( const int numThreads ) { }

      //! \brief returns true if program is operating on one thread currently
      static inline bool singleThreadMode() { return true ; }
    }; // end class ThreadManager

#ifdef _OPENMP
    struct OpenMPThreadManager : public EmptyThreadManager
    {
      //! true if pthreads are used
      static constexpr bool pthreads = false ;

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
#endif


#if HAVE_PTHREAD
    struct PThreadsManager
    {
      //! true if pthreads are used
      static constexpr bool pthreads = true ;

    private:
      struct Manager
      {
        DUNE_EXPORT static inline Manager& instance()
        {
          static thread_local Manager mg ;
          return mg ;
        }
        DUNE_EXPORT int &maxThreads_()
        {
          static int maxThreads = 1;
          return maxThreads;
        }

        inline void initThread( const int maxThreads, const int threadNum )
        {
          // thread number 0 is reserved for the master thread
          maxThreads_() = maxThreads;
          threadNum_  = threadNum ;
        }

        inline void singleThreadMode()
        {
          activeThreads_ = 1;
        }

        inline void multiThreadMode(const int nThreads )
        {
          activeThreads_ = nThreads;
        }

        inline int maxThreads() { return maxThreads_(); }
        inline int currentThreads() const { return activeThreads_; }
        inline int thread()
        {
          assert( threadNum_ >= 0 );
          return threadNum_;
        }

      private:
        int threadNum_;
        int activeThreads_;

        Manager()
          : threadNum_( 0 ), activeThreads_( 1 )
        {}
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

      //! set max number of threads and thread number for this thread
      static inline void initThread( const int maxThreads, const int threadNum )
      {
        manager().initThread( maxThreads, threadNum );
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
        // this call also initiates the master thread
        manager().initThread( numThreads, 0 );
      }

      //! return true if the current thread is the master thread (i.e. thread 0)
      static inline bool isMaster()
      {
        return thread() == 0;
      }

      //! returns true if program is operating on one thread currently
      static inline bool singleThreadMode()
      {
        return currentThreads() == 1 ;
      }
    }; // end class ThreadManager (pthreads)
#endif

#ifdef _OPENMP
// in debug mode show which threading model is used
#ifndef NDEBUG
#warning "ThreadManager: using OpenMP"
#endif
    using ThreadManager = OpenMPThreadManager;
#elif defined(USE_PTHREADS)
// in debug mode show which threading model is used
#ifndef NDEBUG
#warning "ThreadManager: using pthreads"
#endif
    using ThreadManager = PThreadsManager;
#else
    using ThreadManager = EmptyThreadManager;
#endif

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_OMPMANAGER_HH
