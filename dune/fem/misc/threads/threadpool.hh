#ifndef DUNE_FEM_MISC_THREADS_THREADPOOL_HH
#define DUNE_FEM_MISC_THREADS_THREADPOOL_HH

#include <cassert>
#include <mutex>
#include <vector>
#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/fem/misc/threads/threadmanager.hh>
#include <dune/fem/storage/singleton.hh>


namespace Dune
{

namespace Fem
{

  class ThreadPool
  {
    class ObjectIF
    {
    protected:
      ObjectIF() {}
    public:
      virtual ~ObjectIF() {}
      virtual bool run() = 0;
    };

    template <class Object>
    class ObjectWrapper : public ObjectIF
    {
      Object& obj_;
      std::mutex* mutex_;
    public:
      ObjectWrapper( Object& obj, std::mutex* mtx = nullptr )
        : obj_( obj ), mutex_( mtx )
      {}

      bool run ()
      {
        // will be set to true if exception is caught
        bool singleThreadModeError = false ;

        // if mutex was set lock here
        if( mutex_ )
          mutex_->lock();

        try
        {
          // call obj to execute code
          obj_();
        }
        catch ( const Dune::Fem::SingleThreadModeError& e )
        {
          // indicate that this exception was caught be returning true
          singleThreadModeError = true ;
        }

        // if mutex exists then unlock
        if( mutex_ )
          mutex_->unlock();

        return singleThreadModeError;
      }
    };

#ifdef USE_PTHREADS
    ////////////////////////////////////////////
    // class ThreadPoolObject
    ////////////////////////////////////////////
    class ThreadPoolObject
    {
      ObjectIF* objPtr_;
      pthread_barrier_t* barrierBegin_ ;
      pthread_barrier_t* barrierEnd_ ;
      pthread_t threadId_ ;
      int maxThreads_ ;
      int threadNumber_ ;

      std::atomic< bool > singleThreadModeError_;

      bool isNotMainThread () const { return threadNumber_ > 0; }

    public:
      // constructor creating thread with given thread number
      ThreadPoolObject(pthread_barrier_t* barrierBegin,
                       pthread_barrier_t* barrierEnd,
                       const int maxThreads,
                       const int threadNumber )
        : objPtr_( nullptr ),
          barrierBegin_ ( barrierBegin ),
          barrierEnd_ ( barrierEnd ),
          threadId_( 0 ),
          maxThreads_( maxThreads ),
          threadNumber_( threadNumber ),
          singleThreadModeError_( false )
      {
        assert( threadNumber > 0 );
      }

      // constructor creating master thread
      explicit ThreadPoolObject(pthread_barrier_t* barrierBegin,
                                pthread_barrier_t* barrierEnd,
                                const int maxThreads)
        : objPtr_( nullptr ),
          barrierBegin_ ( barrierBegin ),
          barrierEnd_ ( barrierEnd ),
          threadId_( pthread_self() ),
          maxThreads_( maxThreads ),
          threadNumber_( 0 ),
          singleThreadModeError_( false )
      {
      }

      // copy constructor
      ThreadPoolObject(const ThreadPoolObject& other)
        : objPtr_( other.objPtr_ ),
          barrierBegin_( other.barrierBegin_ ),
          barrierEnd_( other.barrierEnd_ ),
          threadId_( other.threadId_ ),
          maxThreads_( other.maxThreads_ ),
          threadNumber_( other.threadNumber_ ),
          singleThreadModeError_( false )
      {}

      // assignment operator
      ThreadPoolObject& operator = ( const ThreadPoolObject& other)
      {
        objPtr_       = other.objPtr_ ;
        barrierBegin_ = other.barrierBegin_ ;
        barrierEnd_   = other.barrierEnd_ ;
        threadId_     = other.threadId_;
        maxThreads_   = other.maxThreads_;
        threadNumber_ = other.threadNumber_;
        singleThreadModeError_ = false;
        return *this;
      }

      // Create the thread and start work
      void start( ObjectIF* obj )
      {
        // init object
        objPtr_ = obj;

        if( isNotMainThread() )
        {
          // if thread has not been initialized
          if( threadId_ == 0 )
          {
            // create a thread that can be joined later
            pthread_create(&threadId_, 0, &ThreadPoolObject::startThread, (void *) this);
          }
        }
        else
        {
          // on master thread there is no need to start an extra thread
          run();
        }
      }

      //! return 1 of thread is stoped, 0 otherwise
      int stoped() const
      {
        return ( objPtr_ == nullptr ) ? 1 : 0 ;
      }

      //! return true if SingleThreadModeError was caught during run
      int singleThreadModeError() const
      {
        return singleThreadModeError_ ? 1 : 0;
      }

      // do the work
      void run()
      {
        assert( barrierBegin_ );
        assert( barrierEnd_ );

        // wait for all threads
        pthread_barrier_wait( barrierBegin_ );

        singleThreadModeError_ = false ;

        // we should be in single thread mode here
        assert( ThreadManager::singleThreadMode() );

        // when object pointer is set call run, else terminate
        if( objPtr_ )
        {
          // update thread_local currentThreads number
          ThreadManager::initMultiThreadMode();

          singleThreadModeError_ = objPtr_->run();
        }
        else
        {
          // this is terminating the threads
          return ;
        }

        // work finished, set objPtr to zero
        objPtr_ = nullptr ;

        // wait for all threads
        pthread_barrier_wait( barrierEnd_ );

        // reset to single thread
        ThreadManager::initSingleThreadMode();

        // when thread is not master then
        // just call run and wait at barrier
        if( isNotMainThread() )
        {
          run();
        }
      }

      //! destroy thread by calling pthread_join
      void destroy()
      {
        if( isNotMainThread() )
          pthread_join(threadId_, 0);
      }

    private:
      // This is the static class function that serves as a
      // C style function pointer for the pthread_create call
      static void* startThread(void *obj)
      {
        // set maxThreads and threadNumber for secondary thread
        ThreadManager :: initThread( ((ThreadPoolObject *) obj)->maxThreads_, ((ThreadPoolObject *) obj)->threadNumber_ );

        // do the work
        ((ThreadPoolObject *) obj)->run();

        return 0;
      }
    }; // end ThreadPoolObject
    ////////////////////////////////////////////////////
    //  end ThreadPoolObject
    ////////////////////////////////////////////////////

    std::vector< ThreadPoolObject > threads_;
    pthread_barrier_t waitBegin_ ;
    pthread_barrier_t waitEnd_ ;
    const int maxThreads_ ;

  private:
    // prohibit copying
    ThreadPool( const ThreadPool& ) = delete;
  public:
    // default constructor
    ThreadPool()
      : threads_()
      , waitBegin_()
      , waitEnd_()
      , maxThreads_( ThreadManager :: maxThreads() )
    {
      // initialize barrier
      pthread_barrier_init( &waitBegin_, 0, maxThreads_ );

      // initialize barrier
      pthread_barrier_init( &waitEnd_, 0, maxThreads_ );

      // initialize slave threads
      for(int i=1; i<maxThreads_; ++i)
      {
        // create thread handles for pthreads
        threads_.push_back( ThreadPoolObject( &waitBegin_, &waitEnd_, maxThreads_, i ) );
      }

      // insert master thread at last because this thread creates
      // all other threads before it start its calculations
      threads_.push_back( ThreadPoolObject( &waitBegin_, &waitEnd_, maxThreads_ ) );
    } // end constructor

  protected:
    //! start all threads to do the job
    bool startThreads( ObjectIF* obj = 0 )
    {
      // start threads, this will call the runThread method
      // and call initMultiThreadMode on ThreadManager
      for(int i=0; i<maxThreads_; ++i)
      {
        threads_[ i ].start( obj );
      }

      // wait until all threads are done
      int count = 0;
      while( count < maxThreads_ )
      {
        count = 0;
        // join threads
        for(int i=0; i<maxThreads_; ++i)
        {
          count += threads_[ i ].stoped() ;
        }
      }

      // activate initSingleThreadMode is called on each thread after obj was run

      int singleError = 0 ;
      // check whether a SingleThreadModeError occurred in one of the threads
      for(int i=0; i<maxThreads_; ++i)
      {
        singleError += int(threads_[ i ].singleThreadModeError());
      }

      return singleError > 0 ;
    }

    //! run all threads
    template <class Functor>
    bool runThreads( Functor& functor, std::mutex* mtx = nullptr )
    {
      // create object wrapper
      ObjectWrapper< Functor > objPtr( functor, mtx );

      // start parallel execution
      return startThreads( &objPtr ) ;
    }

    // return instance of ThreadPool
    static ThreadPool& instance()
    {
      return Singleton< ThreadPool > :: instance();
    }

  public:
    //! destructor deleting threads
    ~ThreadPool()
    {
      // start threads with null object which will terminate each sub thread
      startThreads () ;

      // call thread join
      for(int i=0; i<maxThreads_; ++i)
      {
        threads_[ i ].destroy();
      }

      // destroy barrier
      pthread_barrier_destroy( &waitEnd_ );
      // destroy barrier
      pthread_barrier_destroy( &waitBegin_ );
    }

#endif // end HAVE_PTHREAD

  public:
    /** \brief run functor in multithreaded environment
     *  \param[in] functor  callable object to execute code that
     *                      should be parallel (i.e. a lambda), the signature
     *                      is:
     *                      \code
     *                      void operator() () { };
     *                      \code
     */
    template <class Functor>
    static void run ( Functor& functor )
    {
      // this routine should not be called in multiThreadMode, since
      // this routine is actually starting the multiThreadMode
      if( ! ThreadManager :: singleThreadMode() )
        DUNE_THROW(InvalidStateException,"ThreadPool::run called from thread parallel region!");

#ifdef USE_PTHREADS
      if( ThreadManager :: pthreads )
      {
        bool singleThreadError = instance().runThreads( functor );
        if( singleThreadError )
        {
          DUNE_THROW(SingleThreadModeError, "ThreadPool::run: single thread mode violation occurred!" );
        }
      }
      else
#endif
      {
        std::atomic< bool > singleThreadModeError( false );

        // OpenMP parallel region
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
          // execute code in parallel
          try
          {
            // execute code in parallel
            functor();
          }
          catch (const Dune::Fem::SingleThreadModeError& e)
          {
            singleThreadModeError = true ;
          }
        } // end parallel region

        // only throw one exception to the outside world
        if( singleThreadModeError )
        {
          DUNE_THROW(SingleThreadModeError, "ThreadPool::run: single thread mode violation occurred!");
        }
      }
    }

    /** \brief run functor in multithreaded environment but locked with mutex
     *  \param[in] functor  callable object to execute code that
     *                      should be parallel (i.e. a lambda), the signature
     *                      is:
     *                      \code
     *                      void operator() () { };
     *                      \code
     */
    template <class Functor>
    static void runLocked ( Functor& functor )
    {
      // this routine should not be called in multiThreadMode, since
      // this routine is actually starting the multiThreadMode
      if( ! ThreadManager :: singleThreadMode() )
        DUNE_THROW(InvalidStateException,"ThreadPool::runLocked called from thread parallel region!");

      // run threads in blocking mode
      std::mutex mtx;

#ifdef USE_PTHREADS
      if( ThreadManager :: pthreads )
      {
        // pthread version
        bool singleThreadError = instance().runThreads( functor, &mtx );
        if( singleThreadError )
        {
          DUNE_THROW(SingleThreadModeError, "ThreadPool::runLocked: single thread mode violation occurred!" );
        }
      }
      else
#endif
      {
        // create object wrapper with mutex which executed a locked run
        ObjectWrapper< Functor > obj( functor, &mtx );

        std::atomic< bool > singleThreadModeError ( false );
        // OpenMP parallel region
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
          try
          {
            // execute code in parallel
            obj.run();
          }
          catch (const Dune::Fem::SingleThreadModeError& e)
          {
            singleThreadModeError = true;
          }
        } // end parallel region

        // only throw one exception to the outside world
        if( singleThreadModeError )
        {
          DUNE_THROW(SingleThreadModeError, "ThreadPool::runLocked: single thread mode violation occurred!");
        }
      }
    }

  };

} // end namespace Fem

} // end namespace Dune
#endif
