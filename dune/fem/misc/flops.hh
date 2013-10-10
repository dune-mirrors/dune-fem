#ifndef DUNE_FEM_FLOPS_HH
#define DUNE_FEM_FLOPS_HH

#if HAVE_PAPI
#include <papi.h>
#endif

//- system includes 
#include <iostream>
#include <vector>
#include <cassert>

//- dune-fem includes
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/misc/threads/threadmanager.hh>
#include <dune/fem/misc/threads/threadsafevalue.hh>

namespace Dune {

  namespace Fem {
    
    class FlopCounter
    {
      typedef std::vector< float > values_t ;
      ThreadSafeValue< values_t > values_;
      bool masterStarted_;

      void evaluateCounters( float& realTime, 
                             float& procTime, 
                             float& mFlops ) 
      {
#if HAVE_PAPI
        long long flops ; // we are only interrested in MFLOPs
        int retval = PAPI_flops(&realTime, &procTime, &flops, &mFlops);
        if( retval < PAPI_OK ) 
        {
          std::cerr << "ERROR: PAPI_FP_OPS event is not available, check papi_avail!" << std::endl;
        }
#endif
      }

      FlopCounter () 
        : values_( values_t(3, float(0.0)) ),
          masterStarted_( false )
      {
      }

      void startCounter( const bool startCounters ) 
      {
        const bool reallyStart = (ThreadManager :: isMaster() && startCounters) || 
                                 masterStarted_ ; 
        if( reallyStart ) 
        {
#if HAVE_PAPI
          PAPI_thread_init((unsigned long(*)(void))(pthread_self));
          PAPI_register_thread();
#endif
          float realtime, proctime, mflops;
          evaluateCounters( realtime, proctime, mflops );

          // notify that master thread started counters 
          if( ThreadManager :: isMaster() ) 
            masterStarted_ = true ;
        }
      }

      void stopCounter() 
      {
        if( masterStarted_ ) 
        {
          // get reference to thread local value 
          values_t& values = *values_;
          evaluateCounters( values[ 0 ], values[ 1 ], values[ 2 ] );
          std::cout << "counters[ " << ThreadManager::thread() << "]: " << values[ 2 ] << std::endl;
          std::cout << "Counter: "<< this << std::endl;
        }
      }

      void printCounter( std::ostream& out ) const
      {
        const int threads = ThreadManager :: maxThreads ();
        for( int i=0; i<threads; ++i ) 
        {
          std::cout << "print" << values_[ i ][ 2 ] << std::endl;
        }

        values_t values( values_[ 0 ] );
        // tkae maximum for times and sum flops for all threads 
        for( int i=1; i<threads; ++i ) 
        {
          values[ 0 ] = std::max( values[ 0 ], values_[ i ][ 0 ] );
          values[ 1 ] = std::max( values[ 1 ], values_[ i ][ 1 ] );
          std::cout << "Other threads = " << values_[ i ][ 2 ] << std::endl;
          values[ 2 ] += values_[ i ][ 2 ];
        }

        values_t max( values );
        values_t min( values );
        values_t sum( values );

        typedef MPIManager :: CollectiveCommunication CollectiveCommunication;
        const CollectiveCommunication& comm = MPIManager :: comm();

        const int size = max.size();
        // compute max, min, and sum of flop values 
        comm.max( &max[ 0 ], size ); 
        comm.min( &min[ 0 ], size ); 
        comm.sum( &sum[ 0 ], size ); 

        if( comm.rank() == 0 ) 
        {
          printValues( out, "FlopCounter::sum: ", sum );
          printValues( out, "FlopCounter::max: ", max );
          printValues( out, "FlopCounter::min: ", min );
        }
      }

      static FlopCounter& instance() 
      {
        static FlopCounter counter;
        return counter;
      }

    public:
      static void start( const bool startCounters = true ) 
      {
        instance().startCounter( startCounters );
      }

      static void stop( ) 
      {
        instance().stopCounter();
      }

      static void print( std::ostream& out ) 
      {
        assert( ThreadManager :: singleThreadMode () );
        instance().printCounter( out );
      }

    protected:
      void printValues( std::ostream& out, const std::string name, const values_t& values ) const
      {
        out << name << " ";
        for( unsigned int i=0; i<values.size(); ++i )
        {
          out << values[ i ] << "  ";
        }
        out << std::endl;
      }
    };

  } // namespace Fem
} // namespace Dune
#endif
