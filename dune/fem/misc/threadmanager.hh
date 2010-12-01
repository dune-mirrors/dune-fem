#ifndef DUNE_FEM_OMPMANAGER_HH
#define DUNE_FEM_OMPMANAGER_HH

#ifdef _OPENMP
#include <omp.h>
#endif

namespace Dune 
{
  namespace Fem {

  struct ThreadManager 
  {
    //! return maximal number of threads possbile in the current run 
    static int maxThreads() 
    {
#ifdef _OPENMP
      return omp_get_max_threads();
#else 
      return 1;
#endif
    }

    //! return number of current threads 
    static int currentThreads() 
    {
#ifdef _OPENMP
      return omp_get_num_threads();
#else 
      return 1;
#endif
    }

    //! return thread number 
    static int thread() 
    {
#ifdef _OPENMP
      return omp_get_thread_num();
#else 
      return 0;
#endif
    }

    //! make all threads wait until all reached the barrier
    static void barrier () 
    {
#ifdef _OPENMP
      {
#pragma omp barrier 
      }
#endif
    }

    //! return true if the current thread is the master thread (i.e. thread 0)
    static bool isMaster() 
    {
      return thread() == 0 ;
    }

  }; // end class ThreadManager 

  } // namespace Fem 

} // end namespace Dune 
#endif
