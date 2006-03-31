#ifndef DUNE_TIMEUTILITY_HH
#define DUNE_TIMEUTILITY_HH

#include <limits>

namespace Dune {

  class TimeProvider {
  public:
    double time() const { return time_; }
    
    void provideTimeStepEstimate(double dtEstimate) {
      dtEstimate_ = std::min(dtEstimate_, dtEstimate);
    }
  public:
    TimeProvider(double startTime = 0.0) : 
      time_(startTime),
      dtEstimate_()
    {
      resetTimeStepEstimate();
    }
    
    virtual ~TimeProvider() {}

    void setTime(double time) { time_ = time; }
    void augmentTime(double dt) { time_ += dt; }
    
    void resetTimeStepEstimate() {
      dtEstimate_ = std::numeric_limits<double>::max();
    }
    
    double timeStepEstimate() const {
      return dtEstimate_;
    }
  private:
    TimeProvider(const TimeProvider&);
    TimeProvider& operator=(const TimeProvider&);
    
  private:
    double time_;
    double dtEstimate_;
  };
  
} // end namespace Dune

#endif
