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
    static inline int maxThreads() 
    {
#ifdef _OPENMP
      return omp_get_max_threads();
#else 
      return 1;
#endif
    }

    //! return number of current threads 
    static inline int currentThreads() 
    {
#ifdef _OPENMP
      return omp_get_num_threads();
#else 
      return 1;
#endif
    }

    //! return thread number 
    static inline int thread() 
    {
#ifdef _OPENMP
      return omp_get_thread_num();
#else 
      return 0;
#endif
    }

    //! make all threads wait until all reached the barrier
    static inline void barrier () 
    {
#ifdef _OPENMP
      {
#pragma omp barrier 
      }
#endif
    }

    //! return true if the current thread is the master thread (i.e. thread 0)
    static inline bool isMaster() 
    {
      return thread() == 0 ;
    }

    //! set maximal number of threads available during run  
    static inline void setMaxNumberThreads( const int numThreads ) 
    {
#ifdef _OPENMP
      omp_set_num_threads( numThreads );
#endif
    }

    //! returns true if program is operating on one thread currently
    static inline bool singleThreadMode() 
    {
      return currentThreads() == 1 ;
    }

  }; // end class ThreadManager 

  } // namespace Fem 

} // end namespace Dune 
#endif
