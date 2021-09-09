#ifndef DUNE_FEM_OMPMANAGER_HH
#define DUNE_FEM_OMPMANAGER_HH

#include <atomic>
#include <cassert>
#include <cstdlib>
#include <thread>

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

//- dune-common includes
#include <dune/common/exceptions.hh>
#include <dune/common/visibility.hh>

namespace Dune
{

  namespace Fem
  {
    /*! \brief Exception thrown when a code segment that is supposed to be only accessed in
     *         single thread mode is accessed in multi thread mode. For example,
     *         creation of quadratures or basis function caching cannot work in
     *         multi thread mode.
     *  \note  This exception is not derived from Dune::Exception because there static variables are used.
     */
    class SingleThreadModeError : public std::exception
    {
    public:
      void message(const std::string &msg){}
      const char* what() const noexcept override
      {
        return "SingleThreadModeError";
      }
    };

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

      /** \brief initialize ThreadManager */
      static inline void initialize() {}

      //! \brief initialize single thread mode (when in multithread mode)
      static inline void initSingleThreadMode() {}

      [[deprecated("Use initMultiThreadMode();")]]
      static inline void initMultiThreadMode( const int ) {}

      //! \brief initialize multi thread mode (when in single thread mode)
      static inline void initMultiThreadMode() {}

      //! \brief set max number of threads and thread number for this thread
      static inline void initThread( const int maxThreads, const int threadNum ) {}

      //! \brief return maximal number of threads possible in the current run
      static inline int maxThreads() { return 1;  }

      //! \brief return number of current threads
      static inline int numThreads() { return 1; }

      //! \brief return thread number
      static inline int thread() {  return 0; }

      //! \brief return true if the current thread is the master thread (i.e. thread 0)
      static inline bool isMaster() { return true ; }

      //! \brief set number of threads available during next run
      static inline void setNumThreads( const int numThreads ) { }

      //! \brief set maximal number of threads available during entire run
      static inline void setMaxNumberThreads( const int maxThreads ) { }

      //! \brief returns true if program is operating on one thread currently
      static inline bool singleThreadMode() { return true ; }
    }; // end class ThreadManager

    namespace detail {
      /** \brief read environment variables DUNE_NUM_THREADS and OMP_NUM_THREADS
       * (in that order) to obtain the maximal available number of threads.
       */
      static inline const int getEnvNumberThreads ()
      {
        int maxThreads = 1;
#ifdef USE_SMP_PARALLEL
        // use environment variable (for both openmp or pthreads) if set
        {
          const char* mThreads = getenv("DUNE_NUM_THREADS");
          if( mThreads )
          {
            maxThreads = std::max( int(1), atoi( mThreads ) );
            return maxThreads;
          }
        }
        {
          const char* mThreads = getenv("OMP_NUM_THREADS");
          if( mThreads )
          {
            maxThreads = std::max( int(1), atoi( mThreads ) );
            return maxThreads;
          }
        }
#endif
        return maxThreads;
      }
    } // end namespace detail

#ifdef _OPENMP
    struct OpenMPThreadManager : public EmptyThreadManager
    {
    private:
      static int& maxThreads_()
      {
        // initialize with the system default, change with env vars
        static int m = std::max(1u, std::thread::hardware_concurrency());;
        return m;
      }
    public:
      //! true if pthreads are used
      static constexpr bool pthreads = false ;

      /** \brief initialize ThreadManager with maximal number of threads available.
       *  \note  Can be set by the DUNE_NUM_THREADS or OMP_NUM_THREADS environment variable
       *         from outside of the code.
       */
      static inline void initialize()
      {
        // call maxThreads to initialize the variable
        maxThreads_();

        // set max threads (if set by DUNE_NUM_THREADS )
        omp_set_num_threads( detail::getEnvNumberThreads() );
      }

      /** return maximal number of threads possible during operation */
      static inline int maxThreads()
      {
        return maxThreads_();
      }

      //! return number of threads set by user
      static inline int numThreads()
      {
        return omp_get_max_threads();
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

      //! set number of threads available during next parallel section
      static inline void setNumThreads( const int numThreads )
      {
        omp_set_num_threads( numThreads );
      }

      //! set maximal number of threads usable
      static inline void setMaxNumberThreads( const int numThreads )
      {
        setNumThreads( numThreads );
      }

      [[deprecated("Use initMultiThreadMode();")]]
      static inline void initMultiThreadMode( const int nThreads )
      {
      }

      static inline void initMultiThreadMode()
      {
        // nothing to do here, this is done automatically in the parallel section
      }

      //! returns true if program is operating on one thread currently
      static inline bool singleThreadMode()
      {
        return omp_get_num_threads() == 1 ;
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

      public:
        inline void initThread( const int maxThreads, const int threadNum )
        {
          maxThreads_ = maxThreads ;
          threadNum_  = threadNum ;
        }

        inline void setNumThreads( const int nThreads )
        {
          assert( nThreads <= maxThreads_ );
          numThreads_ = nThreads;
        }

        inline void initSingleThreadMode()
        {
          activeThreads_ = 1;
        }

        inline void initMultiThreadMode( )
        {
          activeThreads_ = numThreads_;
        }

        bool singleThreadMode() const { return activeThreads_ == 1; }

        inline int maxThreads() const { return maxThreads_; }
        inline int numThreads() const { return numThreads_; }
        inline int thread()
        {
          assert( threadNum_ >= 0 );
          return threadNum_;
        }

      private:
        int maxThreads_;
        int numThreads_;
        int activeThreads_;
        int threadNum_;

        Manager()
          : maxThreads_( std::max(1u, std::thread::hardware_concurrency() ) ),
            numThreads_( detail::getEnvNumberThreads() ), activeThreads_( 1 ), threadNum_( 0 )
        {}
      };

      static inline Manager& manager()
      {
        return Manager :: instance();
      }

    public:
      //! initialize thread manager
      static inline void initialize ()
      {
        // this call also initiates the master thread
        // (other threads are set in ThreadPool)
        initThread( maxThreads(), 0 );
      }

      ///////////////////////////////////////////////////////
      //  begin of pthread specific interface
      ///////////////////////////////////////////////////////
      //! initialize single thread mode
      static inline void initSingleThreadMode()
      {
        manager().initSingleThreadMode();
      }

      //! initialize multi thread mode
      static inline void initMultiThreadMode( )
      {
        manager().initMultiThreadMode( );
      }

      //! set max number of threads and thread number for this thread
      static inline void initThread( const int maxThreads, const int threadNum )
      {
        manager().initThread( maxThreads, threadNum );
      }

      //! initialize multi thread mode
      [[deprecated("Use initMultiThreadMode()")]]
      static inline void initMultiThreadMode( const int )
      {
        initMultiThreadMode();
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
      static inline int numThreads()
      {
        return manager().numThreads();
      }

      //! return thread number
      static inline int thread()
      {
        return manager().thread();
      }

      //! set maximal number of threads available during run
      static inline void setNumThreads( const int nThreads )
      {
        manager().setNumThreads( nThreads );
      }

      //! set maximal number of threads available during run
      static inline void setMaxNumberThreads( const int maxThreads )
      {
        manager().initThread( maxThreads, thread() );
      }

      //! return true if the current thread is the master thread (i.e. thread 0)
      static inline bool isMaster()
      {
        return thread() == 0;
      }

      //! returns true if program is operating on one thread currently
      static inline bool singleThreadMode()
      {
        return manager().singleThreadMode();
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
