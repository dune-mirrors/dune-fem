#ifndef DUNE_FEM_FLOPS_HH
#define DUNE_FEM_FLOPS_HH

#if HAVE_PAPI
#include <papi.h>
#endif

//- system includes 
#include <iostream>
#include <vector>

//- dune-fem includes
#include <dune/fem/misc/mpimanager.hh>

namespace Dune {

  namespace Fem {
    
    class FlopCounter
    {
      typedef std::vector< float > values_t ;

      values_t values_;
      long long flops_;

      void evaluateCounters( float& realTime,  float& procTime, 
                             long long& flops, float& mFlops ) 
      {
#if HAVE_PAPI
        int retval = PAPI_flops(&realTime, &procTime, &flops, &mFlops);
        if( retval < PAPI_OK ) 
        {
          DUNE_THROW(InvalidStateException,"PAPI_FP_OPS event is not available");
        }
#endif
      }

    public:
      FlopCounter () 
      : values_( 3, float(0.0 ) ),
        flops_( 0 )
      {
        // reset flop counters
        // start();
      }

      void start() 
      {
        float realtime, proctime, mflops;
        long long flops;
        evaluateCounters( realtime, proctime, flops, mflops );
      }

      void stop() 
      {
        evaluateCounters( values_[ 0 ], values_[ 1 ], flops_, values_[ 2 ] );
      }

      void print( std::ostream& out ) const
      {
        values_t max( values_ );
        values_t min( values_ );
        values_t sum( values_ );

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
