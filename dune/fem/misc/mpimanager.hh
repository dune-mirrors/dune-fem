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
        // maximum number of threads spawned
        int maxThreads_;
        // number of threads to be used in next parallel region
        int numThreads_;
        int activeThreads_;

        std::vector<std::thread> threads_;
        std::unordered_map<std::thread::id,int> numbers_; // still used for possible debugging can be removed if thread_local thread number works
        std::condition_variable_any wait_;
        std::shared_mutex lock_;
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
          ThreadPool::threadNumber_() = t;
          // set thread number (static thread local)
          while (!finalized_)
          {
            // wait until a new task has been set or until threads are to be finalized
            {
              std::shared_lock<std::shared_mutex> lk(lock_);
              auto condition = [this]() { return run_ || finalized_; };
              while (!condition()) // spurious wakeup can happen...
                wait_.wait(lk, condition);
              if (finalized_) break;
              ThreadPool::threadNumber_() = t;
              if (t<numThreads()) run_();
            }
            // /* debug: check that thread local variable for thread number is working
            if (numbers_[std::this_thread::get_id()] != t || threadNumber() != t)
              std::cout << "error with thread numbering: should be " << t
                        << " map says " << numbers_[std::this_thread::get_id()]
                        << " and thread local variable " << threadNumber()
                        << std::endl;
            assert(numbers_[std::this_thread::get_id()] == t);
            assert(ThreadPool::threadNumber_() == t);
            assert(threadNumber() == t);
            // */
            {
              std::shared_lock<std::shared_mutex> lk(lock_);
              auto condition = [this]() { return !run_; };
              while (!condition())
                wait_.wait(lk, condition);
            }
          }
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
          numbers_[std::this_thread::get_id()] = 0;
          for (int t=1;t<maxThreads_;++t)
          {
            threads_.push_back( std::thread( [this,t]() { wait(t); } ) );
            numbers_[threads_[t-1].get_id()] = t;
          }
        }
        ~ThreadPool()
        {
          // all threads should be in the 'start' waiting phase - notify of change of 'finalize variable
          {
            std::unique_lock<std::shared_mutex> lk(lock_);
            finalized_ = true;
          }
          wait_.notify_all();
          // join all threads
          std::for_each(threads_.begin(),threads_.end(), std::mem_fn(&std::thread::join));
        }

        template<typename F, typename... Args>
        void run(F&& f, Args&&... args)
        {
          if ( numThreads_==1 )
            f(args...);
          else
          {
            initMultiThreadMode();
            std::atomic<bool> caughtException(false);
            // notify all threads to new task
            {
              std::unique_lock<std::shared_mutex> lk(lock_);
              run_ = [&]() {
                    try { f(args...); }
                    catch (const SingleThreadModeError& e )
                    { caughtException = true; }
              };
            }
            wait_.notify_all();
            // execture task on master thread
            ThreadPool::threadNumber_() = 0;
            run_(args...);
            // notify all threads that task has been completed
            {
              std::unique_lock<std::shared_mutex> lk(lock_);
              run_ = nullptr;
            }
            wait_.notify_all();
            initSingleThreadMode();
            if( caughtException )
              DUNE_THROW(SingleThreadModeError, "ThreadPool::run: single thread mode violation occurred!");
          }
        }

        int numThreads() { return numThreads_; }
        int maxThreads() { return maxThreads_; }
#if 1
        int threadNumber()
        {
          auto t = ThreadPool::threadNumber_();
          assert( t>=0 );
          return t;
        }
#else
        int threadNumber() { return numbers_.at(std::this_thread::get_id()); }
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

      static inline void initSingleThreadMode() { instance().pool_.initSingleThreadMode(); }
      static inline void initMultiThreadMode() { instance().pool_.initMultiThreadMode(); }
      static bool singleThreadMode() { return instance().pool_.singleThreadMode(); }
      static int numThreads() { return instance().pool_.numThreads(); }
      static int maxThreads() { return instance().pool_.maxThreads(); }
      static void setNumThreads( int use ) { instance().pool_.setNumThreads(use); }
      static int thread() { return instance().pool_.threadNumber(); }
      template<typename F, typename... Args>
      static void run(F&& f, Args&&... args) { instance().pool_.run(f,args...); }
      static bool isMaster() { return instance().pool_.isMaster(); }

    private:
      MPIHelper *helper_ = nullptr;
      std::unique_ptr< CollectiveCommunication > comm_;
      bool wasInitializedHere_ = false ;
#if HAVE_PETSC
      bool petscWasInitializedHere_ = false ;
#endif
      detail::ThreadPool pool_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_MPIMANAGER_HH
