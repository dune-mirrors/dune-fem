#include <config.h>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/misc/femtimer.hh>

namespace Dune
{

  namespace Fem
  {

    Timer< true >::Timer ()
    : timesS_(),
      timers_(),
      output_(),
      stepCount_(0),
      changed_(true)
    {
      push_time();
    }


    unsigned int Timer< true >::add ( const std::string &name, int nr )
    {
      unsigned int id;
      const unsigned int numTimers = timers_.size();
      for( id = 0; id < numTimers; ++id )
      {
        if( timers_[ id ].times.empty() )
          break;
      }

      if( id == numTimers )
        timers_.push_back( TimerInfo( name, nr ) );
      else
        timers_[ id ] = TimerInfo( name, nr );

      reset_timer( id );
      changed_ = true;
      return id;
    }


    void Timer< true >::remove ( unsigned int id )
    {
      timers_[ id ] = TimerInfo( "", 0 );
      changed_ = true;
    }


    void Timer< true >::remove ()
    {
      timers_.clear();
      changed_ = true;
    }


    void Timer< true >::print_timer ( std::ostream &out, int id )
    {
      const TimerInfo &info = timers_[ id ];
      const unsigned int numTimers = timers_.size();
      if( numTimers > 1 )
        out << "(" << info.name << ":";

      out << info.times[ 0 ];
      for( unsigned int i = 1; i < info.times.size(); ++i )
        out << "," << info.times[ i ] / info.times[ 0 ];

      if( numTimers > 1 )
        out << ") ";
    }


    void Timer< true >::print_timer ( std::ostream &out, const std::string &msg )
    {
      out << msg << " : ";
      for( unsigned int i = 0; i < timers_.size(); ++i )
        print_timer( out, i );
      out << std::endl;
    }

    void Timer< true >::printToFile ()
    {
      for( unsigned int i = 0; i< timers_.size(); ++i )
      {
        const TimerInfo &info = timers_[ i ];
        if( info.times.empty() )
          continue;

        output_ << std::setw( 6 ) << inMS( info.times[ 0 ] ) << "ms";
        const unsigned int numTimes = info.times.size();
        if( numTimes > 1 )
        {
          output_ << " ( ";
          for( unsigned int nr = 1; nr < numTimes; ++nr )
            output_ << std::setw( 3 ) << inProz( info.times[ nr ], info.times[ 0 ] ) << "% ";
          output_ << ") ";
        }
      }
      output_ << std::endl;
    }


    void Timer< true >::printToFile ( const std::string &fileName, int step )
    {
      if (!output_.is_open()) {
        output_.open(fileName.c_str());
        if( !output_ )
          DUNE_THROW( IOError, "FemTimer: Unable to open '" << fileName << "' for writing." );
        changed_=true;
      }
      if( changed_ )
      {
        for( unsigned int i = 0; i < timers_.size(); ++i )
        {
          const TimerInfo &info = timers_[ i ];
          if( !info.times.empty() )
            output_ << std::setw( 12 + (info.times.size()-1)*5 ) << info.name;
        }
        output_ << std::endl;
        stepCount_ = 0;
        changed_ = false;
      }
      if( stepCount_ % step == 0 )
        printToFile();
      ++stepCount_;
    }


    void Timer< true >::printToFile ( const TimeProviderBase &tp, const std::string &fileName, int step )
    {
      if (!output_.is_open()) {
        output_.open(fileName.c_str());
        if( !output_ )
          DUNE_THROW( IOError, "FemTimer: Unable to open '" << fileName << "' for writing." );
        changed_=true;
      }
      if( changed_ )
      {
        output_ << std::endl << std::endl;
        output_ << std::setw( 12 ) << "Time" << " ";
        output_ << std::setw( 12 ) << "dt" << " ";
        for( unsigned int i = 0; i < timers_.size(); ++i )
        {
          const TimerInfo &info = timers_[ i ];
          if( !info.times.empty() )
            output_ << std::setw( 12+(info.times.size()-1)*5 ) << info.name;
        }
        output_ << std::endl;
        stepCount_ = 0;
        changed_ = false;
      }
      if( stepCount_ % step == 0 )
      {
        output_ << std::setw( 10 ) << std::scientific << tp.time() << " ";
        output_ << std::setw( 10 ) << std::scientific << tp.deltaT() << " ";
        printToFile();
      }
      ++stepCount_;
    }

  } // namespace Fem

} // namespace Dune
