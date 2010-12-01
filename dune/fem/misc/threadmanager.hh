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
    static int maxThreads() 
    {
#ifdef _OPENMP
      return omp_get_max_threads();
#else 
      return 1;
#endif
    }

    static int currentThreads() 
    {
#ifdef _OPENMP
      return omp_get_num_threads();
#else 
      return 1;
#endif
    }

    static int thread() 
    {
#ifdef _OPENMP
      return omp_get_thread_num();
#else 
      return 0;
#endif
    }

    static void barrier () 
    {
#ifdef _OPENMP
      {
#pragma omp barrier 
      }
#endif
    }

  }; // end class ThreadManager 

  } // namespace Fem 

} // end namespace Dune 
#endif
