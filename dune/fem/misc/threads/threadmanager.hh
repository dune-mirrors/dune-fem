#ifndef DUNE_FEM_OMPMANAGER_HH
#define DUNE_FEM_OMPMANAGER_HH

#error DEPRECATED

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
#include <dune/fem/misc/mpimanager.hh>

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

    /** \class MPIManager
     *  \ingroup Utility
     *  \brief The MPIManager wrapps basic shared memory functionality
     *         provided by OpenMP or pthreads such as thread id, number of threads, etc.
     *
     *  \note All methods are static members.
     */
    struct EmptyMPIManager
    {
      //! true if pthreads are used
      static constexpr bool pthreads = false ;

      /** \brief initialize MPIManager */
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
    }; // end class MPIManager


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
          // return Singleton< Manager > :: instance();
          static thread_local Manager mg; // = Singleton< Manager > :: instance();
          return mg ;
        }
        // friend class Singleton< Manager >;

      public:
        inline void initThread( const int threadNum )
        {
          maxThreads_ = MPIManager::maxThreads() ;
          threadNum_  = threadNum ;
          std::cout << "Manager " << this << " init thread " << threadNum_ << std::endl;
        }

        inline void setNumThreads( )
        {
          uint nThreads = MPIManager::useThreads();
          assert( nThreads <= maxThreads_ );
          if ( nThreads > maxThreads_ )
            DUNE_THROW( InvalidStateException, "requested number of threads exceeds allowed maximum of "+
                           std::to_string(maxThreads_)+
                           " which is fixed at simulation start. Set 'DUNE_NUM_THREADS' environment variable to increase the maximum");
          numThreads_ = nThreads;
          std::cout << "Manager " << this << " set " << numThreads_ << std::endl;
        }

        inline int maxThreads() const { return maxThreads_; }
        inline int numThreads() const
        {
          std::cout << "Manager " << this << " numThreads=" << numThreads_ << std::endl;
          return numThreads_;
        }
        inline int thread()
        {
          assert( threadNum_ >= 0 );
          // std::cout << "Manager " << this << " threadNum=" << threadNum_ << std::endl;
          return threadNum_;
        }

      private:
        int maxThreads_;
        int numThreads_;
        int threadNum_;

        Manager()
          : maxThreads_( MPIManager::maxThreads() )
          , numThreads_( MPIManager::useThreads() )
          , threadNum_ ( 0 )
        {
          std::cout << "Manager::Manager " << this << " = " << maxThreads_
                    << " " << numThreads_ << std::endl;
        }
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
        initThread( -1, 0 );
      }

      ///////////////////////////////////////////////////////
      //  begin of pthread specific interface
      ///////////////////////////////////////////////////////
      //! initialize single thread mode
      static inline void initSingleThreadMode()
      {
        MPIManager::initSingleThreadMode();
      }

      //! initialize multi thread mode
      static inline void initMultiThreadMode( )
      {
        MPIManager::initMultiThreadMode( );
      }

      //! set max number of threads and thread number for this thread
      static inline void initThread( const int maxThreads, const int threadNum )
      {
        manager().initThread( threadNum );
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
        MPIManager::setUseThreads( nThreads );
        manager().setNumThreads( );
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
        return MPIManager::singleThreadMode();
      }
    }; // end class MPIManager (pthreads and OpenMP)

#ifdef _OPENMP
// in debug mode show which threading model is used
#ifndef NDEBUG
#warning "MPIManager: using OpenMP"
#endif
    using MPIManager = PThreadsManager< false >;
#elif defined(USE_PTHREADS)
// in debug mode show which threading model is used
#ifndef NDEBUG
#warning "MPIManager: using pthreads"
#endif
    using MPIManager = PThreadsManager< true >;
#else
    using MPIManager = EmptyMPIManager;
#endif

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_OMPMANAGER_HH
