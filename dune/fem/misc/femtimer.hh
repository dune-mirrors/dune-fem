#ifndef DUNE_FEMTIMER_HH
#define DUNE_FEMTIMER_HH

#include <stack>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <limits>

//#include <dune/common/exceptions.hh>
#include <dune/common/timer.hh>

#include <dune/fem/solver/timeprovider.hh>

namespace Dune
{

  namespace Fem
  {

    // Timer
    // -----

    template< bool enable >
    struct Timer;


    template<>
    struct Timer< false >
    {
      typedef enum { max, sum } operation;

      static unsigned int addTo ( const std::string &name, int nr = 0 ) { return 0; }
      static void removeFrom ( unsigned int id ) {}
      static void removeAll () {}

      static void start ( int id, int nr = 0 ) {}
      static double stop ( int id, int nr = 0, operation op = sum ) { return -1; }
      static double stop ( int id, operation op ) { return -1; }

      static void reset() {}
      static void reset( int id ) {}
      static void reset( int id, int nr ) {}

      static void print ( std::ostream &out, int id ) {}
      static void print ( std::ostream &out, const std::string &msg = "" ) {}

      static void printFile ( const std::string &fileName, int step = 1 ) {}
      static void printFile ( const TimeProviderBase &tp,
                              const std::string &fileName, int step = 1 ) {}
    };


    template<>
    struct Timer< true >
    {
      typedef enum { max, sum } operation;

    private:
      struct TimerInfo
      {
        std::vector< double > startTimes, times;
        std::string name;

        TimerInfo ( const std::string &n, const unsigned int nr )
        : startTimes( nr ), times( nr ), name( n )
        {}
      };

      Timer ();
      ~Timer ();

      void push_time() { timesS_.push( timer_.elapsed() ); }

      double pop_time()
      {
        const double elapsed = timer_.elapsed() - timesS_.top();
        timesS_.pop();
        return elapsed;
      }

      unsigned int add ( const std::string &name, int nr );
      void remove ( unsigned int id );
      void remove ();

      void start_timer( int id, int nr )
      {
        timers_[ id ].startTimes[ nr ] = timer_.elapsed();
        assert( timers_[ id ].startTimes[ 0 ] >= double( 0 ) );
      }

      double stop_timer ( int id, int nr, operation op )
      {
        TimerInfo &info = timers_[ id ];
        assert( (info.startTimes[ nr ] >= double( 0 )) && (info.startTimes[ 0 ] >= double( 0 )) );
        double elapsed = timer_.elapsed() - info.startTimes[ nr ];
        info.startTimes[ nr ] = double( -1 );
        switch( op )
        {
        case sum:
          info.times[ nr ] += elapsed;
          break;
        case max:
          info.times[ nr ] = std::max( info.times[ nr ], elapsed );
          break;
        }
        return elapsed;
      }

      void reset_timer ( int id, int nr )
      {
        timers_[ id ].times[ nr ] = double( 0 );
        timers_[ id ].startTimes[ nr ] = double( -1 );
      }

      void reset_timer ( int id )
      {
        for( unsigned int i = 0; i < timers_[ id ].times.size(); ++i )
          reset_timer( id, i );
      }

      void reset_timer ()
      {
        for( unsigned int i = 0; i < timers_.size(); ++i )
          reset_timer( i );
      }

      void print_timer ( std::ostream &out, int id );
      void print_timer ( std::ostream &out, const std::string &msg );

      size_t inMS ( const double t )
      {
        return (size_t (t * 1e3));
      }

      size_t inProz ( const double p, double rel )
      {
        size_t ret = (size_t)((p / rel) * 100.);
        return std :: min( ret, size_t(100) );
      }

      void printToFile ();
      void printToFile ( const std::string &fileName, int step );
      void printToFile ( const TimeProviderBase &tp, const std::string &fileName, int step );

      static Timer &instance ()
      {
        static Timer instance_;
        return instance_;
      }

    public:
      //! push a new timer to the stack
      static void start () { instance().push_time(); }

      //! retrieve a timer from the stack
      static double stop () { return instance().pop_time(); }

      //! add a new timer with description 
      //! @param name description for output
      //! @param nr   number of subtimers to store
      //! @return     id used to identify this timer in all following calls
      static unsigned int addTo ( const std::string &name, int nr = 0 )
      {
        return instance().add(name,nr+1);
      }

      //! remove a timer with given id
      static void removeFrom ( unsigned int id ) { instance().remove( id ); }

      //! remove all timers 
      static void removeAll () { instance().remove(); }

      //! start a given timer (or subtimer)
      //! @param id the id returned by the method addTo
      //! @param nr the number of the subtimer
      static void start ( int id, int nr = 0 ) { instance().start_timer( id, nr ); }

      //! stop a given timer (or subtimer)
      //! @param id the id of the timer as returned by the method addTo
      //! @param nr the number of the subtimer
      //! @param op the operation to perform (sum or max)
      //! @return the total (or max) time used by this timer since the last reset
      static double stop ( int id, int nr = 0, operation op = sum )
      {
        return instance().stop_timer( id, nr, op );
      }

      //! stop a given timer (with subtimer 0)
      //! @param id the id of the timer as returned by the method addTo
      //! @param op the operation to perform (sum or max)
      //! @return the total (or max) time used by this timer since the last reset
      static double stop ( int id, operation op )
      {
        return instance().stop_timer( id, 0, op );
      }

      //! reset all timers to zero
      static void reset () { instance().reset_timer(); }

      //! reset a given timer with all its subtimers
      static void reset ( int id ) { instance().reset_timer( id ); }

      //! rest a given subtimer
      static void reset ( int id, int nr ) { instance().reset_timer( id, nr ); }

      //! print the values of a given timer (plus subtimers) to a stream
      static void print ( std::ostream &out, int id ) { instance().print_timer( out, id ); }

      //! print the values of all timers to a stream
      static void print ( std::ostream &out, const std::string &msg = "" )
      {
        instance().print_timer(out,msg);
      }

      //! print the values of all timers to a file
      //! if the file is open a new line is appended
      //! @param fileName name of the file
      //! @param step only add a line to the file each step calls of this method
      static void printFile ( const std::string &fileName, int step = 1 )
      {
        instance().printToFile(fileName,step);
      }

      //! print the values of all timers to a file
      //! if the file is open a new line is appended 
      //! information taken from a time provider is also added to the file
      //! @param tp the time provider
      //! @param fileName name of the file
      //! @param step only add a line to the file each step calls of this method
      static void printFile ( const TimeProviderBase &tp,
                              const std::string &fileName, int step = 1 )
      {
        instance().printToFile(tp,fileName,step);
      }

    private:
      Dune::Timer timer_;
      std::stack< double > timesS_;
      std::vector< TimerInfo > timers_;
      std::ofstream output_;
      int stepCount_;
      bool changed_;
    };

  }



   /** \class   FemTimer
     *  \ingroup HelperClasses
     *  \brief   class with singleton instance managing 
     *           timming for parts of program.
     *
     *  The simplest way of timing one line of
     *  code is to enclose it with the
     *  \c TIMEDEXECUTION
     *  macro; the return value corresponds to the
     *  elapsed time, e.g., instead of writting
     *  \code
     *  double error = calcError(u,uh);
     *  \endcode
     *  use
     *  \code
     *  double used = TIMEDEXECUTION(
     *    double error = calcError(u,uh);
     *  )
     *  \endcode
     *
     *  A more general usage is through the 
     *  FemTimer class.
     *  The singleton instance can either be accessed
     *  through FemTimer::instance or through the
     *  reference femTimer.
     *  Note that the following usage is only
     *  available if \c FEMTIMER is defined 
     *  otherwise all methods of the class
     *  FemTimer with the exception of \c start
     *  and \c stop used for the \c TIMEDEXECUTION
     *  macro are empty.
     *  
     *  For computing the execution time of any part
     *  of a code, first get a unique id from the
     *  FemTimer by calling 
     *  \code
     *  id = femTimer.addTo(name,subMarkers);
     *  \endcode
     *  where \c name is a string used for output
     *  and subMarkers is an integer value greaten
     *  or equal to one, which can be used to 
     *  time parts of the program. This can for example
     *  be done in the constructor of an operator.
     *  Remember to return the id to the FemTimer
     *  by calling 
     *  \code
     *  femTimer.removeFrom(id);
     *  \endcode
     *
     *  To start and stop the time keeping for a given program
     *  write
     *  \code
     *  femTimer.start(id);
     *  ...
     *  femTimer.stop(id);
     *  \endcode
     *  Execution time is summed up over all calls
     *  to start and stop. It is possible to pass an operation
     *  argument which changes this behavior; \c sum and \c max are
     *  implemented Using 
     *  \code
     *  femTimer.reset(id);
     *  \endcode
     *  the main timer (and all its subtimings)
     *  are set back to zero. Calling \c reset
     *  without an argument resets all stored
     *  timers.
     *
     *  The use of sub timers works as shown in
     *  the following example:
     *  \code
     *  unsigned int id = femTimer.addTo("test",2);
     *  ...
     *  femTimer.start(id);
     *  ...
     *  femTimer.start(id,1);
     *  f1(); // call to a first function
     *  femTimer.end(id,1);
     *  ...
     *  femTimer.start(id,2);
     *  f2(); // call to a second function
     *  femTimer.end(id,2);
     *  ...
     *  femTimer.end(id);
     *  \endcode
     *  Using \c femTimer.print(out,"test"); 
     *  the result of all timings is printed to
     *  an \c ostream. Subtimings are given
     *  relative to the main timing, i.e.,
     *  first the main time is printed and following
     *  that the relative time used for each subpart
     *  of the algorithm. 
     *  In the same manner the timing information
     *  can be stored in a file using
     *  \c printFile(filename).
     *  The first call opens the file and prints
     *  the string identifying each timing;
     *  each successive call prints one line containing
     *  all the timing information, again given
     *  first the main timing followed by the
     *  relative time used in each sub timing.
     */
#ifdef FEMTIMER
    typedef Fem::Timer< true > FemTimer;
#else
    typedef Fem::Timer< false > FemTimer;
#endif



#define TIMEDEXECUTION(command) \
    (femTimer.start(),command,femTimer.stop())
   /** \class   ExecutionTimer
     *  \ingroup HelperClasses
     *  \brief  class with a start and stop method for
     *          timing parts of a program.
     **/
  class ExecutionTimer {
    public:
    ExecutionTimer() : total_(0) {
    }
    void start() {
      start_=time_.elapsed();
    }
    void end() {
      total_=start_-time_.elapsed();
    }
    double read() {
      return total_;
    }
    void reset() {
      total_=0;
    }
    double total_;
    double start_;
    Timer time_;
  };

}

#endif // #ifndef DUME_FEMTIMER_HH
