#ifndef DUNE_TIMEPROVIDER_HH
#define DUNE_TIMEPROVIDER_HH

//- system includes 
#include <limits>
#include <cassert>

//- Dune includes 
#include <dune/fem/io/file/asciiparser.hh>

namespace Dune {

/** @ingroup ODESolver 

    \brief A TimeProvider is a class for managing the global simulation time of
    time-dependent simulations. 

    \remarks 
    There exsist three implementations of time provider, two for serial
    runs (i.e. TimeProvider and ImprovedTimeProvider) 
    and a wrapper using these two for parallel runs
    (i.e. ParallelTimeProvider).
  */
  class TimeProvider 
  {
  public:
    /** \brief constructor taking initial time, default is zero
        \param[in] startTime initial time 
     */
    TimeProvider(double startTime = 0.0) : 
      time_(startTime),
      savedTime_(time_),
      dt_(-1.0),
      dtEstimate_(0.0),
      cfl_(1.0),
      timeStep_(0),
      lock_(true)
    {
      resetTimeStepEstimate();
    }
    
    /** \brief constructor taking start time and cfl number 
        \param[in] startTime initial time 
        \param[in] cfl initial cfl number 
    */
    TimeProvider(double startTime , double cfl ) : 
      time_(startTime),
      savedTime_(time_),
      dt_(-1.0),
      dtEstimate_(0.0),
      cfl_(cfl),
      timeStep_(0),
      lock_(true)
    {
      resetTimeStepEstimate();
    }
    
    //! destructor 
    ~TimeProvider() {}

    /** \brief return current time 
        \return current time 
    */
    double time() const { return time_; }

    /** \brief stores current global time */
    void unlock() 
    {
      lock_ = false;
      savedTime_ = time_;
    }
    
    /** \brief restores saved global time */
    void lock () 
    {
      lock_ = true;
      time_ = savedTime_;
    }

    /** \brief set time step estimate to minimum of given value and
               internal time step estiamte 
         \param[in] dtEstimate time step size estimate 
    */
    void provideTimeStepEstimate(const double dtEstimate) {
      dtEstimate_ = std::min(dtEstimate_, dtEstimate);
    }
    
    /** \brief set internal time to given time 
         \param[in] time new time 
    */
    void setTime(const double time) 
    { 
      assert( ! lock_ );
      time_ = time; 
    }

    /** \brief set internal time to given time 
         and time step counter to given counter 
         \param[in] time new time 
         \param[in] timeStep new time step counter 
    */
    void setTime(const double time, const int timeStep )  
    { 
      assert( ! lock_ );
      time_ = time; 
      timeStep_ = timeStep;
    }

    /** \brief augment time , i.e. \f$t = t + \triangle t\f$ and
        increase time step counter  */
    void augmentTime()
    { 
      assert( lock_ );
      time_ += deltaT(); 
      ++timeStep_;
    }
    
    /** \brief reset set time step estimate 
        by setting ti to big value 
    */
    void resetTimeStepEstimate() 
    {
      // reset estimate 
      dtEstimate_ = std::numeric_limits<double>::max();
    }
    
    /** \brief  return time step size estimate  
        \return time step size estimate  
    */
    double timeStepEstimate() const {
      return dtEstimate_;
    }

    /** \brief  return cfl number 
        \return cfl number 
    */
    double cfl () const { return cfl_; } 
    
    /** \brief set internal cfl to given value
        \param[in] cfl new cfl number 
    */
    void setCfl(const double cfl) 
    {
      cfl_ = cfl;
    }
    
    /** \brief set internal cfl to minimum of given value and internal
        cfl number 
        \param[in] cfl cfl estimate 
    */
    void provideCflEstimate(const double cfl) 
    {
      cfl_ = std::min(cfl_, cfl );
    }
    
    /** \brief  return time step size times cfl number
        \return \f$\triangle t \cdot CFL\f$ 
    */
    double deltaT () const 
    {
      assert( (dt_ * cfl_) > 0.0 );
      return dt_ * cfl_;
    }

    /** \brief  sets time step size 
    */
    void setDeltaT (double dt)
    {
      assert(dt > 0.0 );
      dt_ = dt/cfl_;
    }
    
    /** \brief  syncronize time step, i.e. set timeStep to values of current 
        estimate and reset estimate 
    */
    void syncTimeStep() 
    {
      // save current time step 
      dt_ = dtEstimate_;
      // reset estimate 
      resetTimeStepEstimate();
    }

    /** \brief return current time step counter 
        \return current time step counter 
    */
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
    double savedTime_;
    double dt_;
    double dtEstimate_;
    double cfl_;
    int timeStep_;
    bool lock_;
  };
  
  /** \brief improved class for 
      for time and time step estimate handling
      also reads parameter from parameter file 
  */
  class ImprovedTimeProvider : public TimeProvider 
  {
  public:
    /** \brief constructor 
        \param[in] paramFile parameter file to read start time and cfl number
                   keywords are "StartTime" and "CFL" 
        \param[in] rank rank of process, output is only shown for process 0 (default value is 0) 
    */
    ImprovedTimeProvider(const std::string paramFile = "",
                         const int rank = 0) 
      : TimeProvider() 
    {
      this->time_ = readStartTime(paramFile, rank);
      this->cfl_  = readCFL(paramFile, rank);
    }

  private:
    //! do not copy this class 
    ImprovedTimeProvider(const ImprovedTimeProvider&);
    ImprovedTimeProvider& operator=(const ImprovedTimeProvider&);

    // read parameter start time from given file
    double readStartTime(const std::string& file, const int rank) const 
    {
      double startTime = 0.0;
      readParameter(file,"StartTime",startTime, (rank == 0) );
      return startTime;
    }
    
    // read parameter CFL from given file
    double readCFL(const std::string& file, const int rank) const 
    {
      double cfl = 1.0;
      readParameter(file,"CFL",cfl, (rank == 0) );
      return cfl;
    }
  };

  /** \brief TimeProvider that is working in a parallel program */
  template <class CommunicatorType>
  class ParallelTimeProvider
  {
  public:
    /** \brief constructor communicator and serial time provider
        \param[in] comm communicator 
        \param[in] tp serial time provider 
    */
    ParallelTimeProvider(const CommunicatorType& comm,
                         TimeProvider& tp)
      : comm_(comm), tp_(tp)
    {}
    
    //! return internal time 
    double time() const { return tp_.time(); }
    
    /** @copydoc TimeProvider::provideTimeStepEstimate */
    void provideTimeStepEstimate(const double dtEstimate) 
    {
      tp_.provideTimeStepEstimate(dtEstimate);
    }

    /** @copydoc TimeProvider::lock */
    void lock() { tp_.lock(); }
    
    /** @copydoc TimeProvider::unlock */
    void unlock() { tp_.unlock(); }
    
    /** @copydoc TimeProvider::setTime */
    void setTime(const double time) { tp_.setTime(time); }

    /** @copydoc TimeProvider::resetTimeStepEstimate */
    void resetTimeStepEstimate() 
    {
      tp_.resetTimeStepEstimate(); 
    }
    
    /** @copydoc TimeProvider::timeStepEstimate */
    double timeStepEstimate() const 
    {
      return tp_.timeStepEstimate();
    }

    /** @copydoc TimeProvider::cfl */
    double cfl () const { return tp_.cfl(); } 
    
    /** @copydoc TimeProvider::deltaT */
    double deltaT () const { return tp_.deltaT(); }
    
    /** @copydoc TimeProvider::timeStep */
    int timeStep() const 
    {
      return tp_.timeStep();
    }

    /** @copydoc TimeProvider::syncTimeStep 
        
        \note Here a global communication to minimize the time step size
        over all time step sizes from all processors is done. 
    */ 
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

    /** @copydoc TimeProvider::augmentTime */
    double augmentTime() 
    {
      // increase time 
      tp_.augmentTime();
      // return new time 
      return time();
    }

  private:
    //! communicator  
    const CommunicatorType& comm_;
    //! serial time provider 
    TimeProvider& tp_;
  };
  
} // end namespace Dune
#endif
