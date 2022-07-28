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
#include <dune/fem/misc/threads/threadsafevalue.hh>
#include <dune/fem/storage/singleton.hh>

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
      typedef std::vector< float >  values_t ;
      ThreadSafeValue< values_t  >  values_;
      ThreadSafeValue< long long >  flop_;
      ThreadSafeValue< int >        stopped_;

      // call PAPI_flops for given values
      void evaluateCounters( float& realTime,
                             float& procTime,
                             float& mFlops,
                             long long& flop )
      {
#if HAVE_PAPI
        int retval = PAPI_flops(&realTime, &procTime, &flop, &mFlops);
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

      static unsigned long threadId ()
      {
        return MPIManager :: thread();
      }

      // initialize counters
      void startCounter()
      {
        if( ! MPIManager :: singleThreadMode() )
        {
#if HAVE_PAPI
          PAPI_thread_init( threadId );
          PAPI_register_thread();
#endif
        }
        float realtime, proctime, mflops;
        long long flop ;
        evaluateCounters( realtime, proctime, mflops, flop );
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
          long long& flop = *flop_;
          evaluateCounters( values[ 0 ], values[ 1 ], values[ 2 ], flop );

          // mark thread as stopped
          *stopped_ = 1 ;
        }
      }

      // print values to given ostream, all values are gathered to
      // the master rank
      void printCounter( std::ostream& out ) const
      {
        // make sure this method is called in single thread mode only
        assert( MPIManager :: singleThreadMode () );

        int allStopped = 0 ;
        const int threads = MPIManager :: maxThreads ();
        for( int i=0; i<threads; ++i )
        {
          allStopped += stopped_[ i ];
        }

        // make sure all other thread have been stopped, otherwise
        // the results wont be coorect
        if( allStopped != threads )
          DUNE_THROW(InvalidStateException,"Not all thread have been stopped");

        typedef std::vector< double > result_t ;
        result_t values( 5, 0.0 );

        for( int i=0; i<3; ++i )
          values[ i ] = values_[ 0 ][ i ];
        values[ 3 ] = flop_[ 0 ];

        // tkae maximum for times and sum flops for all threads
        for( int i=1; i<threads; ++i )
        {
          values[ 0 ]  = std::max( values[ 0 ], double(values_[ i ][ 0 ]) );
          values[ 1 ]  = std::max( values[ 1 ], double(values_[ i ][ 1 ]) );
          values[ 2 ] += values_[ i ][ 2 ];
          values[ 3 ] += flop_[ i ];
        }
        // convert to GFLOP
        values[ 3 ] /= 1.0e9 ;
        // compute mflops ourselfs
        values[ 4 ] = values[ 3 ] / values[ 0 ];

        result_t max( values );
        result_t min( values );
        result_t sum( values );

        typedef MPIManager :: Communication Communication;
        const Communication& comm = MPIManager :: comm();

        const int size = max.size();
        // compute max, min, and sum of flop values
        comm.max( &max[ 0 ], size );
        comm.min( &min[ 0 ], size );
        comm.sum( &sum[ 0 ], size );

        if( comm.rank() == 0 )
        {
          out << "FlopCounter::typ:   real  proc  mflops flop  flop/real " << std::endl;
          printValues( out, "FlopCounter::sum: ", sum );
          printValues( out, "FlopCounter::max: ", max );
          printValues( out, "FlopCounter::min: ", min );
        }
      }

      friend class Dune::Fem::Singleton< FlopCounter >;

      static FlopCounter& instance()
      {
        return Dune::Fem::Singleton< FlopCounter >::instance();
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
      template <class vec_t>
      void printValues( std::ostream& out, const std::string name, const vec_t& values ) const
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
