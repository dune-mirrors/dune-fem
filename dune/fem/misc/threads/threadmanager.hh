#ifndef DUNE_FEM_OMPMANAGER_HH
#define DUNE_FEM_OMPMANAGER_HH

#include <atomic>
#include <cassert>
#include <cstdlib>
#include <thread>
#include <iostream>

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
#include <dune/fem/storage/singleton.hh>

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
#ifndef NDEBUG
      // for performance reasons we only copy messages when in debug mode
      std::string msg_;
      void message(const std::string &msg) { msg_ = msg; }
      const char* what() const noexcept override { return msg_.c_str(); }
#else
      void message(const std::string &msg) {}
      const char* what() const noexcept override
      {
        return "SingleThreadModeError: remove -DNDEBUG to obtain a more detailed message!";
      }
#endif
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
      static inline const unsigned int getEnvNumberThreads (unsigned int defaultValue)
      {
#ifdef USE_SMP_PARALLEL
        unsigned int maxThreads = defaultValue;
        // use environment variable (for both openmp or pthreads) if set
        const char* mThreads = std::getenv("DUNE_NUM_THREADS");
        if( mThreads )
          maxThreads = std::max( int(1), atoi( mThreads ) );
        else
        {
          const char* mThreads = std::getenv("OMP_NUM_THREADS");
          if( mThreads )
            maxThreads = std::max( int(1), atoi( mThreads ) );
        }
#else
        unsigned int maxThreads = 1;
#endif
        return maxThreads;
      }
    } // end namespace detail

/*
#ifdef _OPENMP
    class OpenMPThreadManager : public EmptyThreadManager
    {
      unsigned int maxThreads__;
      unsigned int numThreads__;
      inline void setNumThreads_( const unsigned int nThreads )
      {
        assert( nThreads <= maxThreads__ );
        if (nThreads > maxThreads__)
          DUNE_THROW( InvalidStateException, "requested number of threads exceeds allowed maximum of "+
                         std::to_string(maxThreads__)+
                         " which is fixed at simulation start. Set 'DUNE_NUM_THREADS' environment variable to increase the maximum");
        numThreads__ = nThreads;
      }
      inline void setMaxThreads_( const unsigned int maxThreads )
      {
        maxThreads__ = maxThreads;
      }
      inline int maxThreads_() const { return maxThreads__; }
      inline int numThreads_() const { return numThreads__; }

    public:

      static inline OpenMPThreadManager& manager()
      {
        return Singleton<OpenMPThreadManager> :: instance();
      }

      //! true if pthreads are used
      static constexpr bool pthreads = false ;

      OpenMPThreadManager()
      : maxThreads__( std::max(1u,
                              detail::getEnvNumberThreads( std::thread::hardware_concurrency() )
                             ) ),
        numThreads__( detail::getEnvNumberThreads(1) )
      {}

      //! return maximal number of threads possbile in the current run
      static inline int maxThreads()
      {
        return manager().maxThreads_();
      }

      //! return number of current threads
      static inline int numThreads()
      {
        return manager().numThreads_();
      }

      //! return thread number
      static inline int thread()
      {
        return omp_get_thread_num();
      }

      //! set maximal number of threads available during run
      static inline void setNumThreads( const int nThreads )
      {
        manager().setNumThreads_( nThreads );
      }

      //! set maximal number of threads available during run
      static inline void setMaxNumberThreads( const int maxThreads )
      {
        manager().setMaxThreads_( maxThreads );
      }

      //! return true if the current thread is the master thread (i.e. thread 0)
      static inline bool isMaster()
      {
        return thread() == 0;
      }

      [[deprecated("Use initMultiThreadMode();")]]
      static inline void initMultiThreadMode( const int ) {}

      //! \brief initialize multi thread mode (when in single thread mode)
      static inline void initMultiThreadMode()
      {
        // nothing to do here, this is done automatically in the parallel section
      }

      //! returns true if program is operating on one thread currently
      static inline bool singleThreadMode()
      {
        return omp_get_num_threads() == 1 ;
      }
    }; // end class ThreadManager (pthreads)
#endif
*/

    template <bool usePThreads>
    struct PThreadsManager
    {
      //! true if pthreads are used
      static constexpr bool pthreads = usePThreads ;

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
          if ( nThreads > maxThreads_ )
            DUNE_THROW( InvalidStateException, "requested number of threads exceeds allowed maximum of "+
                           std::to_string(maxThreads_)+
                           " which is fixed at simulation start. Set 'DUNE_NUM_THREADS' environment variable to increase the maximum");
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
          : maxThreads_( std::max(1u,
                                  detail::getEnvNumberThreads( std::thread::hardware_concurrency() )
                                 ) ),
            numThreads_( detail::getEnvNumberThreads(1) ), activeThreads_( 1 ), threadNum_( 0 )
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
    }; // end class ThreadManager (pthreads and OpenMP)

#ifdef _OPENMP
// in debug mode show which threading model is used
#ifndef NDEBUG
#warning "ThreadManager: using OpenMP"
#endif
    //using ThreadManager = OpenMPThreadManager;
    using ThreadManager = PThreadsManager< false >;
#elif defined(USE_PTHREADS)
// in debug mode show which threading model is used
#ifndef NDEBUG
#warning "ThreadManager: using pthreads"
#endif
    using ThreadManager = PThreadsManager< true >;
#else
    using ThreadManager = EmptyThreadManager;
#endif

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_OMPMANAGER_HH
