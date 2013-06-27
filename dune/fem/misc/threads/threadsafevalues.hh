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
      //! type of value to be thread safe
      typedef T ValueType ;

      //! \brief constructor initializing values for all threads given a init value
      ThreadSafeValue( const ValueType& init ) 
        : value_( 
#ifdef USE_SMP_PARALLEL
            ThreadManager::maxThreads(), 
#endif
            init )
      {}

      //! \brief default constructor 
      ThreadSafeValue() 
        : value_( 
#ifdef USE_SMP_PARALLEL
            ThreadManager::maxThreads()
#endif
             )
      {}

      //! \brief return number of threads 
      size_t size() const { return ThreadManager::maxThreads(); }

      //! \brief return reference to thread private value 
      ValueType& operator * () { return this->operator[]( ThreadManager::thread() ); }
      //! \brief return reference to thread private value 
      const ValueType& operator * () const { return this->operator[]( ThreadManager::thread() ); }

      //! \brief return reference to private value for given thread number 
      ValueType& operator [] ( const unsigned int thread ) { 
        assert( thread < size() );
        return value_
#ifdef USE_SMP_PARALLEL
          [ thread ]
#endif    
          ;
      }

      //! \brief return reference to private value for given thread number 
      const ValueType& operator [] ( const unsigned int thread ) const { 
        assert( thread < size() );
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
