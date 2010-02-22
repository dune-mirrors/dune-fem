#ifndef FEMTIMER_HH
#define FEMTIMER_HH

#include <stack>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <limits>

#include <dune/common/exceptions.hh>
#include <dune/common/timer.hh>

#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/misc/mpimanager.hh>

namespace Dune {
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
class FemTimer {
  public:
  typedef enum {max,sum} operation;
  private:
  FemTimer() : timesS_(), 
               startTimesV_(),
               timesV_(0), timesVName_(0),
               output_(),
               stepCount_(0),
               changed_(true)
  {
    push_time();
  }

  ~FemTimer ()
  {
    double totalTime = pop_time();

    if( output_.is_open() )
    {
      output_ << "#  ******** TOTAL RUNTIME: " << totalTime
              << "   ******** " << std::endl;
      output_.close();
    }

    const MPIManager::CollectiveCommunication &comm = MPIManager::comm();
    if( comm.rank() == 0 )
    {
      double *totalTimes = new double[ comm.size() ];
      comm.gather( &totalTime, totalTimes, 1, 0 );
      double avgTime = 0.0;
      double minTime = std::numeric_limits< double >::max();
      double maxTime = std::numeric_limits< double >::min();
      for( int i = 0; i < comm.size(); ++i )
      {
        avgTime += totalTimes[ i ];
        minTime = std::min( minTime, totalTimes[ i ] );
        maxTime = std::max( maxTime, totalTimes[ i ] );
      }
      avgTime /= comm.size();
      delete[] totalTimes;

      std::cerr << "#  ******** TOTAL RUNTIME: average = " << avgTime
                << ", minimum = " << minTime << ", maximum = " << maxTime
                << "   ******** " << std::endl;
    }
    else comm.gather( &totalTime, (double *)0, 1, 0 );
  }

  void push_time() {
    timesS_.push(timer_.elapsed());
  }
  double pop_time() {
    double ret = timer_.elapsed()-timesS_.top();
    timesS_.pop();
    return ret;
  }
  // ********************************************************
  unsigned int add(const std::string& name,int nr) {
    int id;
    for (id=0;id<int(timesV_.size());++id) {
      if (timesV_[id].size()==0) 
        break;
    }
    if (id==int(timesV_.size())) {
      startTimesV_.push_back(std::vector<double>(nr));
      timesV_.push_back(std::vector<double>(nr));
      timesVName_.push_back(name);
    } else {
      startTimesV_[id] = std::vector<double>(nr);
      timesV_[id] = std::vector<double>(nr);
      timesVName_[id] = name;
    }
    reset_timer(id);
    changed_=true;
    return id;
  }
  void remove(unsigned int id) {
    timesV_[id].clear();
    startTimesV_[id].clear();
    timesVName_[id] = "";
    changed_=true;
  }
  void remove() {
    timesV_.clear();
    startTimesV_.clear();
    timesVName_.clear();
    changed_=true;
  }
  // ******************************************
  void start_timer(int id,int nr)
  {
    startTimesV_[id][nr] = timer_.elapsed();
    assert( startTimesV_[ id ][ 0 ] >= 0. );
  }
  double stop_timer(int id,int nr,operation op) {
    assert( (startTimesV_[ id ][ nr ] >= 0.) && (startTimesV_[ id ][ 0 ] >= 0.) );
    double ret = timer_.elapsed() - startTimesV_[id][nr];
    startTimesV_[ id ][ nr ] = -1.;
    switch (op) {
    case sum: timesV_[id][nr] += ret; 
              break;
    case max: timesV_[id][nr] = std::max(ret,timesV_[id][nr]); 
              break;
    }
    return ret;
  }
  void reset_timer(int id,int nr) {
    timesV_[id][nr] = 0.;
    startTimesV_[ id ][ nr ] = -1.;
  }
  void reset_timer(int id) {
    for (unsigned int i=0;i<timesV_[id].size();++i)
      reset_timer(id,i);
  }
  void reset_timer() {
    for (unsigned int i=0;i<timesV_.size();++i)
      reset_timer(i);
  }
  // *****************************************************
  void print_timer(std::ostream& out,int id) {
    if (timesV_.size()>1) {
      out << "(" << timesVName_[id] << ":";
    }
    out << timesV_[id][0];
    for (unsigned int i=1;i<timesV_[id].size();++i)
      out << "," << timesV_[id][i]/timesV_[id][0];
    if (timesV_.size()>1) {
      out << ") ";
    }
  }
  void print_timer(std::ostream& out,const std::string& msg) {
    out << msg << " : ";
    for (unsigned int i=0;i<timesV_.size();++i)
      print_timer(out,i);
    out << std::endl;
  }
  inline size_t inMS(const double t) 
  {
    return (size_t (t * 1e3));
  }
  inline size_t inProz(const double p,double rel) 
  {
    size_t ret = (size_t)((p / rel) * 100.);
    return std :: min( ret, size_t(100) );
  }
  void printToFile() {
    for (unsigned int i=0;i<timesV_.size();++i) {
      if (timesV_[i].size()>0) {
        output_ << std::setw(6) << inMS(timesV_[i][0]);
        if (timesV_[i].size()>1) {
          output_ << " ( ";
          for (unsigned int nr=1;nr<timesV_[i].size();++nr) {
            output_ << std::setw(3) << inProz(timesV_[i][nr],timesV_[i][0])
                    << "% ";
          }
          output_ << ") ";
        }
      }
    }
    output_ << std::endl;
  }
  void printToFile(const std::string& fileName,
                   int step) {
    if (!output_.is_open()) {
      output_.open(fileName.c_str());
      if( !output_ )
        DUNE_THROW( IOError, "FemTimer: Unable to open '" << fileName << "' for writing." );
      changed_=true;
    }
    if (changed_) {
      for (unsigned int i=0;i<timesV_.size();++i) {
        if (timesV_[i].size()>0) 
          output_ << std::setw(10+(timesV_[i].size()-1)*5) 
                  << timesVName_[i];
      }
      output_ << std::endl;
      stepCount_=0;
      changed_=false;
    }
    if (stepCount_%step==0) {
      printToFile();
    }
    stepCount_++;
  }
  void printToFile(const TimeProviderBase& tp,
                   const std::string& fileName,
                   int step) {
    if (!output_.is_open()) {
      output_.open(fileName.c_str());
      if( !output_ )
        DUNE_THROW( IOError, "FemTimer: Unable to open '" << fileName << "' for writing." );
      changed_=true;
    }
    if (changed_) {
      output_ << std::endl << std::endl;
      output_ << std::setw(12) << "Time" << " ";
      output_ << std::setw(12) << "dt" << " ";
      for (unsigned int i=0;i<timesV_.size();++i) {
        if (timesV_[i].size()>0) 
          output_ << std::setw(10+(timesV_[i].size()-1)*5) 
                  << timesVName_[i];
      }
      output_ << std::endl;
      stepCount_=0;
      changed_=false;
    }
    if (stepCount_%step==0) {
      output_ << std::setw(10) << std::scientific << tp.time() << " ";
      output_ << std::setw(10) << std::scientific << tp.deltaT() << " ";
      printToFile();
    }
    stepCount_++;
  }
  Timer timer_;
  std::stack<double> timesS_;
  std::vector<std::vector<double> > startTimesV_;
  std::vector<std::vector<double> > timesV_;
  std::vector<std::string> timesVName_;
  std::ofstream output_;
  int stepCount_;
  bool changed_;
  // **************************************
  // **************************************
  // **************************************
  public:
  static FemTimer& instance() {
    static FemTimer instance_;
    return instance_;
  }
  //! push a new timer to the stack
  static void start() {
    instance().push_time();
  }
  //! retrieve a timer from the stack
  static double stop() {
    return instance().pop_time();
  }
  //! add a new timer with description 
  //! @param name description for output
  //! @param nr   number of subtimers to store
  //! @return     id used to identify this timer in all following calls
  static unsigned int addTo(const std::string& name, int nr=0) {
#ifdef FEMTIMER
    return instance().add(name,nr+1);
#else
    return 0;
#endif
  }
  //! remove a timer with given id
  static void removeFrom(unsigned int id) {
#ifdef FEMTIMER
    instance().remove(id);
#endif
  }
  //! remove all timers 
  static void removeAll() {
#ifdef FEMTIMER
    instance().remove();
#endif
  }
  //! start a given timer (or subtimer)
  //! @param id the id returned by the method addTo
  //! @param nr the number of the subtimer
  static void start(int id,int nr=0) {
#ifdef FEMTIMER
    instance().start_timer(id,nr);
#endif
  }
  //! stop a given timer (or subtimer)
  //! @param id the id of the timer as returned by the method addTo
  //! @param nr the number of the subtimer
  //! @param op the operation to perform (sum or max)
  //! @return the total (or max) time used by this timer since the last reset
  static double stop(int id,int nr=0,operation op=sum) {
#ifdef FEMTIMER
    return instance().stop_timer(id,nr,op);
#else
    return -1;
#endif
  }
  //! stop a given timer (with subtimer 0)
  //! @param id the id of the timer as returned by the method addTo
  //! @param op the operation to perform (sum or max)
  //! @return the total (or max) time used by this timer since the last reset
  static double stop(int id,operation op) {
#ifdef FEMTIMER
    return instance().stop_timer(id,0,op);
#else
    return -1;
#endif
  }
  //! reset all timers to zero
  static void reset() {
#ifdef FEMTIMER
    instance().reset_timer();
#endif
  }
  //! reset a given timer with all its subtimers
  static void reset(int id) {
#ifdef FEMTIMER
    instance().reset_timer(id);
#endif
  }
  //! rest a given subtimer
  static void reset(int id,int nr) {
#ifdef FEMTIMER
    instance().reset_timer(id,nr);
#endif
  }
  //! print the values of a given timer (plus subtimers) to a stream
  static void print(std::ostream& out,int id) {
#ifdef FEMTIMER
    instance().print_timer(out,id);
#endif
  }
  //! print the values of all timers to a stream
  static void print(std::ostream& out,const std::string msg="") {
#ifdef FEMTIMER
    instance().print_timer(out,msg);
#endif
  }
  //! print the values of all timers to a file
  //! if the file is open a new line is appended
  //! @param fileName name of the file
  //! @param step only add a line to the file each step calls of this method
  static void printFile(const std::string& fileName, int step=1) {
#ifdef FEMTIMER
    instance().printToFile(fileName,step);
#endif
  }
  //! print the values of all timers to a file
  //! if the file is open a new line is appended 
  //! information taken from a time provider is also added to the file
  //! @param tp the time provider
  //! @param fileName name of the file
  //! @param step only add a line to the file each step calls of this method
  static void printFile(const TimeProviderBase& tp,
                        const std::string& fileName, int step=1) {
#ifdef FEMTIMER
    instance().printToFile(tp,fileName,step);
#endif
  }
};

// This construction is incompatible with the MPIHelper singleton
#if 0
namespace {
  FemTimer& femTimer = FemTimer::instance();
};
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
#endif
