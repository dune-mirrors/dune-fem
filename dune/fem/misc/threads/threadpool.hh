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
      virtual void run() = 0;
    };

    template <class Object>
    class ObjectWrapper : public ObjectIF
    {
      Object& obj_;
      std::mutex* mutex_;
    public:
      ObjectWrapper( Object& obj, std::mutex* mtx = nullptr )
        : obj_( obj ), mutex_( mtx ) {}

      void run ()
      {
        // if mutex was setm lock here
        if( mutex_ )
          mutex_->lock();

        // call obj to execute code
        obj_();

        // if mutex exists then unlock
        if( mutex_ )
          mutex_->unlock();
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

      bool isSlave () const { return threadNumber_ > 0; }

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
          threadNumber_( threadNumber )
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
          threadNumber_( 0 )
      {
      }

      // copy constructor
      ThreadPoolObject(const ThreadPoolObject& other)
        : objPtr_( other.objPtr_ ),
          barrierBegin_( other.barrierBegin_ ),
          barrierEnd_( other.barrierEnd_ ),
          threadId_( other.threadId_ ),
          maxThreads_( other.maxThreads_ ),
          threadNumber_( other.threadNumber_ )
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
        return *this;
      }

      // Create the thread and start work
      void start( ObjectIF* obj )
      {
        // init object
        objPtr_ = obj;

        if( isSlave() )
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

      // do the work
      void run()
      {
        assert( barrierBegin_ );
        assert( barrierEnd_ );

        // wait for all threads
        pthread_barrier_wait( barrierBegin_ );

        // when object pointer is set call run, else terminate
        if( objPtr_ )
        {
          objPtr_->run();
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

        // when thread is not master then
        // just call run and wait at barrier
        if( isSlave() )
        {
          run();
        }
      }

      //! destroy thread by calling pthread_join
      void destroy()
      {
        if( isSlave() )
          pthread_join(threadId_, 0);
      }

    private:
      // This is the static class function that serves as a
      // C style function pointer for the pthread_create call
      static void* startThread(void *obj)
      {
        // set maxThreads and threadNumber for slave thread
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
    void startThreads( ObjectIF* obj = 0 )
    {
      // set number of active threads
      ThreadManager :: initMultiThreadMode( maxThreads_ );

      // start threads, this will call the runThread method
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

      // activate initSingleThreadMode again
      Fem :: ThreadManager :: initSingleThreadMode();
    }

    //! run all threads
    template <class Functor>
    void runThreads( Functor& functor, std::mutex* mtx = nullptr )
    {
      // create object wrapper
      ObjectWrapper< Functor > objPtr( functor, mtx );

      // start parallel execution
      startThreads( &objPtr ) ;
    }

    // return instance of ThreadPool
    static ThreadPool& instance()
    {
      return Singleton< ThreadPool > :: instance();
      //static Singleton< std::unique_ptr< ThreadPool > handle( new ThreadPool() );
      //return *handle;
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
        DUNE_THROW(InvalidStateException,"ThreadPool :: run called from thread parallel region!");

#ifdef USE_PTHREADS
      if( ThreadManager :: pthreads )
      {
        // pthread version
        instance().runThreads( functor );
      }
      else
#endif
      {
        // OpenMP parallel region
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
          // execute code in parallel
          functor();
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
        DUNE_THROW(InvalidStateException,"ThreadPool :: run called from thread parallel region!");

      // run threads in blocking mode
      std::mutex mtx;

#ifdef USE_PTHREADS
      if( ThreadManager :: pthreads )
      {
        // pthread version
        instance().runThreads( functor, &mtx );
      }
      else
#endif
      {
        // create object wrapper with mutex which executed a locked run
        ObjectWrapper< Functor > obj( functor, &mtx );

        // OpenMP parallel region
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
          {
            // execute code in parallel
            obj.run();
          }
        }
      }
    }

  };

} // end namespace Fem

} // end namespace Dune
#endif
