#ifndef DUNE_TIMEUTILITY_HH
#define DUNE_TIMEUTILITY_HH

//- system includes 
#include <limits>

//- Dune includes 
#include <dune/fem/io/file/asciiparser.hh>

namespace Dune {

  //! class for time and time step estimate handling
  class TimeProvider 
  {
  public:
    //! constructor taking initial time, default is zero
    TimeProvider(double startTime = 0.0) : 
      time_(startTime),
      dt_(-1.0),
      dtEstimate_(0.0),
      cfl_(1.0),
      timeStep_(0)
    {
      resetTimeStepEstimate();
    }
    
    TimeProvider(double startTime , double cfl ) : 
      time_(startTime),
      dt_(-1.0),
      dtEstimate_(0.0),
      cfl_(cfl),
      timeStep_(0)
    {
      resetTimeStepEstimate();
    }
    
    //! destructor 
    virtual ~TimeProvider() {}

    bool notInitialized() const { return (dt_ < 0.0); }

    //! return internal time 
    double time() const { return time_; }
    
    //! set time step estimate to minimum of given value and
    //! internal time step estiamte 
    void provideTimeStepEstimate(double dtEstimate) {
      dtEstimate_ = std::min(dtEstimate_, dtEstimate);
    }
    
    //! set internal time to given time 
    void setTime(double time) { time_ = time; }

    //! set internal time to given time 
    //! and time step counter to given counter 
    void setTime(double time, int timeStep )  
    { 
      time_ = time; 
      timeStep_ = timeStep;
    }

    //! increase internal time by internal timeStep 
    void augmentTime() 
    { 
      time_ += deltaT(); 
      ++timeStep_;
    }
    
    //! set time step estimate to big value 
    void resetTimeStepEstimate() 
    {
      // reset estimate 
      dtEstimate_ = std::numeric_limits<double>::max();
    }
    
    //! return time step estimate 
    double timeStepEstimate() const {
      return dtEstimate_;
    }

    //! return cfl number 
    double cfl () const { return cfl_; } 
    
    //! set internal cfl to minimum of given value and
    //! internal clf value  
    void provideCflEstimate(double cfl) 
    {
      cfl_ = std::min(cfl_, cfl );
    }
    
    //! return time step estimate times cfl number 
    double deltaT () const 
    {
      assert( (dt_ * cfl_) > 0.0 );
      return dt_ * cfl_;
    }
    
    //! syncronize time step, i.e. set timeStep to values of current
    //! estimate and reset
    void syncTimeStep() 
    {
      // save current time step 
      dt_ = dtEstimate_;
      // reset estimate 
      resetTimeStepEstimate();
    }

    int timeStep () const 
    {
      return timeStep_;
    }
  private:
    //! do not copy this class 
    TimeProvider(const TimeProvider&);
    TimeProvider& operator=(const TimeProvider&);
    
  protected:
    double time_;
    double dt_;
    double dtEstimate_;
    double cfl_;
    int timeStep_;
  };
  
  //! improved class for 
  //! for time and time step estimate handling
  //! alos reads parameter from parameter file 
  class ImprovedTimeProvider : public TimeProvider 
  {
  public:
    //! constructor taking initial time, default is zero
    ImprovedTimeProvider(const std::string paramFile = "") 
      : TimeProvider() 
    {
      this->time_ = readStartTime(paramFile);
      this->cfl_  = readCFL(paramFile);
    }

  private:
    //! do not copy this class 
    ImprovedTimeProvider(const ImprovedTimeProvider&);
    ImprovedTimeProvider& operator=(const ImprovedTimeProvider&);

    // read parameter start time from given file
    double readStartTime(const std::string& file) const 
    {
      double startTime = 0.0;
      readParameter(file,"StartTime",startTime);
      return startTime;
    }
    
    // read parameter CFL from given file
    double readCFL(const std::string& file) const 
    {
      double cfl = 1.0;
      readParameter(file,"CFL",cfl);
      return cfl;
    }
  };

  template <class CommunicatorType>
  class ParallelTimeProvider
  {
  public:
    ParallelTimeProvider(const CommunicatorType& comm,
                         TimeProvider& tp)
      : comm_(comm), tp_(tp)
    {}
    
    bool notInitialized() const { return tp_.notInitialized(); }
    
    //! return internal time 
    double time() const { return tp_.time(); }
    
    //! set time step estimate to minimum of given value and
    //! internal time step estiamte 
    void provideTimeStepEstimate(double dtEstimate) 
    {
      tp_.provideTimeStepEstimate(dtEstimate);
    }
    
    //! set internal time to given time 
    void setTime(double time) { tp_.setTime(time); }

    //! set time step estimate to big value 
    void resetTimeStepEstimate() {
      tp_.resetTimeStepEstimate(); 
    }
    
    //! return time step estimate 
    double timeStepEstimate() const 
    {
      return tp_.timeStepEstimate();
    }

    //! return cfl number 
    double cfl () const { return tp_.cfl(); } 
    
    //! return time step estimate times cfl 
    double deltaT () const { return tp_.deltaT(); }
    
    //! return number of current time step
    int timeStep() const 
    {
      return tp_.timeStep();
    }

    //! syncronize time step between processors 
    void syncTimeStep() 
    {
      // get time step estimate 
      double dt = timeStepEstimate(); 
      // do min over all processors 
      dt = comm_.min( dt );
      // set time step estimate 
      tp_.provideTimeStepEstimate(dt);
      // set timeStep and reset estimate 
      tp_.syncTimeStep();
    }

    //! augment time and return new value 
    double augmentTime() 
    {
      // increase time 
      tp_.augmentTime();
      // return new time 
      return time();
    }

  private:
    const CommunicatorType& comm_;
    TimeProvider& tp_;
  };
  
} // end namespace Dune

#endif
