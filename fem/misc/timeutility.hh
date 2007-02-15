#ifndef DUNE_TIMEUTILITY_HH
#define DUNE_TIMEUTILITY_HH

#include <limits>

namespace Dune {

  //! class for time and time step estimate handling
  class TimeProvider 
  {
  public:
    //! return internal time 
    double time() const { return time_; }
    
    //! set time step estimate to minimum of given value and
    //! internal time step estiamte 
    void provideTimeStepEstimate(double dtEstimate) {
      dtEstimate_ = std::min(dtEstimate_, dtEstimate);
    }
  public:
    //! constructor taking initial time, default is zero
    TimeProvider(double startTime = 0.0) : 
      time_(startTime),
      dtEstimate_()
    {
      resetTimeStepEstimate();
    }
    
    //! destructor 
    virtual ~TimeProvider() {}

    //! set internal time to given time 
    void setTime(double time) { time_ = time; }

    //! increase internal time by dt 
    void augmentTime(double dt) { time_ += dt; }
    
    //! set time step estimate to big value 
    void resetTimeStepEstimate() {
      dtEstimate_ = std::numeric_limits<double>::max();
    }
    
    //! return time step estimate 
    double timeStepEstimate() const {
      return dtEstimate_;
    }
    
  private:
    //! do not copy this class 
    TimeProvider(const TimeProvider&);
    TimeProvider& operator=(const TimeProvider&);
    
  private:
    double time_;
    double dtEstimate_;
  };
  
} // end namespace Dune

#endif
