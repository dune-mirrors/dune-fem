#ifndef DUNE_FEM_THREADSAFEVALUES_HH
#define DUNE_FEM_THREADSAFEVALUES_HH

#include <vector>
#include <dune/fem/misc/threads/threadmanager.hh>

namespace Dune {

  namespace Fem {

    template <class T>
    class ThreadSafeValue
    {
#ifdef USE_SMP_PARALLEL
      std::vector< T > value_;
#else 
      T value_;
#endif
    public:
      typedef T ValueType ;

      ThreadSafeValue( const ValueType& init ) 
        : value_( 
#ifdef USE_SMP_PARALLEL
            ThreadManager::maxThreads(), 
#endif
            init )
      {}

      ThreadSafeValue() 
        : value_( 
#ifdef USE_SMP_PARALLEL
            ThreadManager::maxThreads()
#endif
             )
      {}

      size_t size() const { return ThreadManager::maxThreads(); }

      ValueType& operator * () { return this->operator[]( ThreadManager::thread() ); }
      const ValueType& operator * () const { return this->operator[]( ThreadManager::thread() ); }

      ValueType& operator [] ( const int thread ) { 
        return value_
#ifdef USE_SMP_PARALLEL
          [ thread ]
#endif    
          ;
      }

      const ValueType& operator [] ( const int thread ) const { 
        return value_
#ifdef USE_SMP_PARALLEL
          [ thread ]
#endif    
          ;
      }
    };

  } // end namespace Fem

} // end namespace Dune 


#endif
