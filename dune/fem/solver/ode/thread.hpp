#ifndef THREAD_HPP
#define THREAD_HPP


#include <pthread.h>
#include <cassert>
#include <iostream>


extern "C"{ 
  void* _run(void *ptr);
}



class Thread
{
public:
  Thread(int num_of_locks=0);
  virtual ~Thread();
  void start();
  void stop();
  void join();
  void lock(int i=0);
  void unlock(int i=0);

protected:
  virtual void run() = 0;

private:
  friend void* _run(void *ptr);

  const int num_of_locks;
  pthread_mutex_t *mutex;
  pthread_t thread;
};





inline
void* _run(void *ptr) 
{
  Thread *thread = (Thread *)ptr;
  thread->run();
  return NULL;
}


inline
Thread::Thread(int num_of_locks) : num_of_locks(num_of_locks)
{
  mutex = new pthread_mutex_t[num_of_locks];
  assert(mutex);
  for(int i=0; i<num_of_locks; i++) pthread_mutex_init(mutex+i, NULL);
}


inline
Thread::~Thread()
{
  delete[] mutex;
}


inline
void Thread::start()
{
  pthread_create( &thread, NULL, _run, (void *)this );
}


inline
void Thread::stop()
{
  pthread_cancel(thread);
}


inline
void Thread::join()
{
  pthread_join( thread, NULL );
}


inline
void Thread::lock(int i)
{
  assert(i>=0 && i<num_of_locks);
  pthread_mutex_lock(mutex+i);
}


inline
void Thread::unlock(int i)
{
  assert(i>=0 && i<num_of_locks);
  pthread_mutex_unlock(mutex+i);
}


#endif
