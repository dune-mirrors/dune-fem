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
    
    // FlopCounter
    // -----------

    /** 
     * @brief A class wrapper for the function PAPI_flops 
     *        from the package PAPI. The results are 
     *        CPU time, real and process local and the number 
     *        of floating point operations in MFLOP/s.
     **/
    class FlopCounter
    {
      typedef std::vector< float > values_t ;
      ThreadSafeValue< values_t >  values_;
      ThreadSafeValue< int >       stopped_;

      // call PAPI_flops for given values
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

      // constructor
      FlopCounter () 
        : values_( values_t(3, float(0.0)) ),
          stopped_( 0 )
      {
      }

      // initialize counters 
      void startCounter()
      {
        if( ! ThreadManager :: singleThreadMode() ) 
        {
#if HAVE_PAPI
          PAPI_thread_init((unsigned long(*)(void))(ThreadManager :: thread));
          PAPI_register_thread();
#endif
        }
        float realtime, proctime, mflops;
        evaluateCounters( realtime, proctime, mflops );
        // mark as not stopped 
        *stopped_ = 0;
      }

      // stop counters and store values
      void stopCounter() 
      {
        if( *stopped_ == 0 ) 
        {
          // get reference to thread local value 
          values_t& values = *values_;
          evaluateCounters( values[ 0 ], values[ 1 ], values[ 2 ] );

          // mark thread as stopped 
          *stopped_ = 1 ;
        }
      }

      // print values to given ostream, all values are gathered to 
      // the master rank
      void printCounter( std::ostream& out ) const
      {
        // make sure this method is called in single thread mode only 
        assert( ThreadManager :: singleThreadMode () );

        int allStopped = 0 ;
        const int threads = ThreadManager :: maxThreads ();
        for( int i=0; i<threads; ++i ) 
        {
          allStopped += stopped_[ i ];
        }

        // make sure all other thread have been stopped, otherwise 
        // the results wont be coorect 
        if( allStopped != threads ) 
          DUNE_THROW(InvalidStateException,"Not all thread have been stopped");

        values_t values( values_[ 0 ] );
        // tkae maximum for times and sum flops for all threads 
        for( int i=1; i<threads; ++i ) 
        {
          values[ 0 ] = std::max( values[ 0 ], values_[ i ][ 0 ] );
          values[ 1 ] = std::max( values[ 1 ], values_[ i ][ 1 ] );
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
          out << "FlopCounter::typ:   real  proc  mflops " << std::endl;
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
      /** \brief Start counters.
       *
       * \note  Call this method for the master thread before all other 
       *        threads call this method if used in a multi-thread environment.
       */
      static void start( )
      {
        instance().startCounter();
      }

      /** \brief stop counters */
      static void stop( ) 
      {
        instance().stopCounter();
      }

      /** \brief print values to given ostream, all values are gathered to 
       *         the master rank before printing 
       */         
      static void print( std::ostream& out ) 
      {
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
