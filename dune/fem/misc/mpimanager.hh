#ifndef DUNE_FEM_MPIMANAGER_HH
#define DUNE_FEM_MPIMANAGER_HH

#if defined _OPENMP || defined(USE_PTHREADS)
#ifndef USE_SMP_PARALLEL
#define USE_SMP_PARALLEL
#endif
#endif

#include <memory>
#include <condition_variable>
#include <thread>
#include <chrono>
#include <functional>
#include <shared_mutex>
#include <atomic>

#include <dune/common/parallel/mpicommunication.hh>
#include <dune/common/parallel/mpihelper.hh>

#if HAVE_PETSC
#include <dune/fem/misc/petsc/petsccommon.hh>
#endif

#include <dune/fem/storage/singleton.hh>

#ifdef _OPENMP
#include <omp.h>
#endif

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

      class ThreadPool
      {
        static const bool useStdThreads = true ;
#ifndef _OPENMP
        static_assert( useStdThreads, "useStdThreads is disabled but OpenMP has not been found!");
#endif

        // maximum number of threads spawned
        int maxThreads_;
        // number of threads to be used in next parallel region
        int numThreads_;
        int activeThreads_;

        std::vector<std::thread> threads_;
        std::unordered_map<std::thread::id,int> numbers_; // still used for possible debugging can be removed if thread_local thread number works
        std::condition_variable_any waitA_;
        std::shared_mutex lockA_;
        std::condition_variable_any waitB_;
        std::shared_mutex lockB_;

        // function to run
        std::function<void(void)> run_;
        // stop thread
        bool finalized_;

#if 1   // this doesn't work as expected
        // store a static thread local variable for the thread number
        static int& threadNumber_()
        {
          static thread_local int number = -1;
          return number;
        }
#endif
        // method executed by each thread
        void wait(int t)
        {
          // set thread number (static thread local)
          ThreadPool::threadNumber_() = t;

          std::shared_lock<std::shared_mutex> lkA(lockA_);
          std::shared_lock<std::shared_mutex> lkB(lockB_);

          while (!finalized_)
          {
            // wait until a new task has been set or until threads are to be finalized
            // unlock 'A' and wait until master thread either changed run_ // or finalizes
            // reaquire the (shared) lock after that
            while (!run_ && !finalized_)
              waitA_.wait(lkA);
            // check if to finalize
            if (finalized_) break;
            ThreadPool::threadNumber_() = t;
            // run the code is required - note that both shared locks are
            // held so the main thread has to wait to uniquely acquire
            // lock 'B' until 'run_' was finished by all threads
            if (t<numThreads())
              run_();
            // this is the same 'waiting' done above but on the 'B' lock. In this case
            // we wait until 'run_' has been cleared again by the main thread
            // which can only happen after all threads have enter the
            // 'wait' which releases the 'B' lock.
            // This is needed to make sure a thread doesn't execute the same 'run_' twice
            while (run_)
              waitB_.wait(lkB);
          }
        }
<<<<<<< HEAD

        public:
        ThreadPool()
        : maxThreads_( std::max(1u, detail::getEnvNumberThreads( std::thread::hardware_concurrency() )) )
        , numThreads_( detail::getEnvNumberThreads(1) )
        , activeThreads_(1)
        , threads_()
        , run_(nullptr)
        , finalized_(false)
        {
          // spawn max number of threads to use
          ThreadPool::threadNumber_() = 0;
#ifdef USE_SMP_PARALLEL
          if( useStdThreads )
          {
            numbers_[std::this_thread::get_id()] = 0;
            for (int t=1;t<maxThreads_;++t)
            {
              threads_.push_back( std::thread( [this,t]() { wait(t); } ) );
              numbers_[threads_[t-1].get_id()] = t;
            }
          }
#endif
        }
        ~ThreadPool()
        {
#ifdef USE_SMP_PARALLEL
          if( useStdThreads )
          {
            std::unique_lock<std::shared_mutex> lkA(lockA_);
            finalized_ = true;
          }
          waitA_.notify_all();
          // join all threads
          std::for_each(threads_.begin(),threads_.end(), std::mem_fn(&std::thread::join));
#endif
        }

        template<typename F, typename... Args>
        void runConsec(F&& f, Args&&... args)
        {
          if(! useStdThreads )
          {
            runOpenMP(f, args...);
            return ;
          }

          if ( numThreads_==1 )
            f(args...);
          else
          {
            initMultiThreadMode();
            bool caughtException(false);
            run_ = [&]() {
                  try { f(args...); }
                  catch (const SingleThreadModeError& e )
                  { caughtException = true; }
            };
            for (int t=0;t<numThreads();++t)
            {
              ThreadPool::threadNumber_() = t;
              run_(args...);
            }
            ThreadPool::threadNumber_() = 0;
            initSingleThreadMode();
            if( caughtException )
              DUNE_THROW(SingleThreadModeError, "ThreadPool::run: single thread mode violation occurred!");
          }
        }

=======
>>>>>>> thread local variables simply don't work with clang so use std::unordered_map for thread numbers
        template<typename F, typename... Args>
        void runOpenMP(F&& f, Args&&... args)
        {
#ifdef _OPENMP
          const int nThreads = numThreads();
          if( nThreads == 1 )
          {
            f(args...);
            return ;
          }

          std::atomic< bool > singleThreadModeError( false );

          initMultiThreadMode();
#pragma omp parallel num_threads(nThreads)
          {
            // set thread number to thread_local variable
            threadNumber_() = omp_get_thread_num();
            // execute given code in parallel
            try
            {
              f(args...);
            }
            catch (const Dune::Fem::SingleThreadModeError& e)
            {
//#ifndef NDEBUG
//            std::cout << "thread[" << ThreadManager::thread() << "] " << e.what() << std::endl;
//#endif
              singleThreadModeError = true ;
            }

          } // end parallel region

          // enter single thread mode again
          initSingleThreadMode();

          // only throw one exception to the outside world
          if( singleThreadModeError )
          {
            DUNE_THROW(SingleThreadModeError, "ThreadPool::run: single thread mode violation occurred!");
          }
#endif
        }

        public:
        ThreadPool()
        : maxThreads_( std::max(1u, detail::getEnvNumberThreads( std::thread::hardware_concurrency() )) )
        , numThreads_( detail::getEnvNumberThreads(1) )
        , activeThreads_(1)
        , threads_()
        , run_(nullptr)
        , finalized_(false)
        {
          // spawn max number of threads to use
          ThreadPool::threadNumber_() = 0;
#ifdef USE_SMP_PARALLEL
          if( useStdThreads )
          {
            numbers_[std::this_thread::get_id()] = 0;
            for (int t=1;t<maxThreads_;++t)
            {
              threads_.push_back( std::thread( [this,t]() { wait(t); } ) );
              numbers_[threads_[t-1].get_id()] = t;
            }
          }
#endif
        }
        ~ThreadPool()
        {
#ifdef USE_SMP_PARALLEL
          if( useStdThreads )
          {
            // all threads should be in the 'start' waiting phase - notify of change of 'finalize variable
            {
              std::unique_lock<std::shared_mutex> lk(lockA_);
              finalized_ = true;
            }
            waitA_.notify_all();
            // join all threads
            std::for_each(threads_.begin(),threads_.end(), std::mem_fn(&std::thread::join));
          }
#endif
        }

        template<typename F, typename... Args>
        void run(F&& f, Args&&... args)
        {
          if(! useStdThreads )
          {
            runOpenMP(f, args...);
            return ;
          }

          if ( numThreads_==1 )
            f(args...);
          else
          {
            // see explanation in 'wait' function
            initMultiThreadMode();
            std::atomic<bool> caughtException(false);
            {
              // acquire lock and set 'run_' - can only be done if all
              // threads are waiting at top of while loop
              std::lock_guard<std::shared_mutex> lkA(lockA_);
              run_ = [&]() {
                    try { f(args...); }
                    catch (const SingleThreadModeError& e )
                    { caughtException = true; }
              };
            }
            // notify all threads of new task - those will all acquire the lock (shared)
            waitA_.notify_all();
            // execute task on master thread
            ThreadPool::threadNumber_() = 0;
            run_(args...);
            {
              // try to acquire lock in non shared mode - this is only possible if all threads have
              // finished the current task and are waiting at bottom of loop
              std::lock_guard<std::shared_mutex> lkB(lockB_);
              run_ = nullptr;
            }
            // notify all threads that task has been completed
            // this moves all threads back to beginning of while loop freeing 'A'
            waitB_.notify_all();

            initSingleThreadMode();
            if( caughtException )
              DUNE_THROW(SingleThreadModeError, "ThreadPool::run: single thread mode violation occurred!");
          }
        }

        int numThreads() { return numThreads_; }
        int maxThreads() { return maxThreads_; }
#if 0
        int threadNumber()
        {
          auto t = ThreadPool::threadNumber_();
          assert( t>=0 );
          return t;
        }
#else
        int threadNumber()
        {
#ifdef _OPENMP
          if(! useStdThreads )
            return omp_get_thread_num();
          else
#endif
          return numbers_.at(std::this_thread::get_id());
        }
#endif
        void initSingleThreadMode() { activeThreads_ = 1; }
        void initMultiThreadMode() { activeThreads_ = numThreads_; }
        bool singleThreadMode() { return activeThreads_ == 1; }
        void setNumThreads( int use ) { assert( singleThreadMode() ); numThreads_ = use; }
        bool isMaster() { return threadNumber() == 0; }
      };

    } // end namespace detail


    struct MPIManager
    {
      typedef Dune::CollectiveCommunication< MPIHelper::MPICommunicator >
        CollectiveCommunication;
    private:
      static MPIManager &instance ()
      {
        return Singleton< MPIManager > :: instance();
      }

      static bool mpiFinalized ()
      {
        bool finalized = false ;
#if HAVE_MPI
        // check that MPI was not already finalized
        {
          int wasFinalized = -1;
          MPI_Finalized( &wasFinalized );
          finalized = bool( wasFinalized );
        }
#endif // #if HAVE_MPI
        return finalized ;
      }

    public:
      //! destructor calling finalize if this has not been done
      ~MPIManager()
      {
        _finalize();
      }

      void _finalize()
      {
        if( ! mpiFinalized() )
        {
#if HAVE_PETSC
          if( petscWasInitializedHere_ )
            ::Dune::Petsc::finalize();
#endif
          // if MPI_Init was called here and finalize has not been
          // called yet, then this is the place to call it
          if( wasInitializedHere_ )
          {
#if HAVE_MPI
            MPI_Finalize();
#endif
          }
        }
      }

      static void finalize()
      {
        instance()._finalize();
      }

      static void initialize ( int &argc, char **&argv );

      static const CollectiveCommunication &comm ()
      {
        const std::unique_ptr< CollectiveCommunication > &comm = instance().comm_;
        if( !comm )
          DUNE_THROW( InvalidStateException, "MPIManager has not been initialized." );
        return *comm;
      }

      static int rank ()
      {
        return comm().rank();
      }

      static int size ()
      {
        return comm().size();
      }

      //! \brief initialize single thread mode (when in multithread mode)
      static inline void initSingleThreadMode() { instance().pool_.initSingleThreadMode(); }

      //! \brief initialize multi thread mode (when in single thread mode)
      static inline void initMultiThreadMode() { instance().pool_.initMultiThreadMode(); }

      //! \brief return maximal number of threads possible in the current run
      static int maxThreads() { return instance().pool_.maxThreads(); }

      //! \brief return number of current threads
      static int numThreads() { return instance().pool_.numThreads(); }

      //! \brief return thread number
      static int thread() { return instance().pool_.threadNumber(); }

      //! \brief return true if the current thread is the master thread (i.e. thread 0)
      static bool isMaster() { return instance().pool_.isMaster(); }

      //! \brief set number of threads available during next run
      static void setNumThreads( int use ) { instance().pool_.setNumThreads(use); }

      //! \brief returns true if program is operating on one thread currently
      static bool singleThreadMode() { return instance().pool_.singleThreadMode(); }

      //! \brief run functor f with given arguments args in threaded mode
      template<typename F, typename... Args>
      static void run(F&& f, Args&&... args) { instance().pool_.run(f,args...); }

    private:
      MPIHelper *helper_ = nullptr;
      std::unique_ptr< CollectiveCommunication > comm_;
      bool wasInitializedHere_ = false ;
#if HAVE_PETSC
      bool petscWasInitializedHere_ = false ;
#endif
      detail::ThreadPool pool_;
    };

    using ThreadManager = MPIManager;
    using ThreadPool = MPIManager;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_MPIMANAGER_HH
