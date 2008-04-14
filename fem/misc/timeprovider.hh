#ifndef DUNE_TIMEPROVIDER_HH
#define DUNE_TIMEPROVIDER_HH

//- system includes 
#include <limits>
#include <cassert>

//- Dune includes 
#include <dune/common/exceptions.hh>
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

    \note 
    An example time loop could look as follows:

    @code

    // create time provider 
    TimeProvider tp ( startTime, clf );

    // creates time step estimate 
    odeSolver.initialize( U );

    for( tp.init(); tp.time() < endTime; tp.next() )
    {
      // do calculation 
    }

    @endcode
  */
  class TimeProvider 
  {
  public:
    /** \brief constructor taking initial time, default is zero
        \param[in] startTime initial time 
     */
    TimeProvider(const double startTime = 0.0) : 
      time_(startTime),
      dt_(-1.0),
      dtEstimate_(0.0),
      cfl_(1.0),
      timeStep_(0),
      synced_(false)
    {
      resetTimeStepEstimate();
    }
    
    /** \brief constructor taking start time and cfl number 
        \param[in] startTime initial time 
        \param[in] cfl initial cfl number 
    */
    TimeProvider(const double startTime , const double cfl ) : 
      time_(startTime),
      dt_(-1.0),
      dtEstimate_(0.0),
      cfl_(cfl),
      timeStep_(0),
      synced_(false)
    {
      resetTimeStepEstimate();
    }
    
    //! destructor 
    ~TimeProvider() {}

    /** \brief init dt with given estimate */
    void init() 
    {
      syncTimeStep();
    }
    
    /** \brief goto next time step */
    void next() 
    {
      // if already syncronized do nothing 
      if( syncronized() ) return ;
      
      // increase time by delta t 
      augmentTime(); 
      
      // set new delta t 
      syncTimeStep();
    }

    /** \brief return current time 
        \return current time 
    */
    double time() const { return time_; }
    
    /** \brief set time step estimate to minimum of given value and
               internal time step estiamte 
         \param[in] dtEstimate time step size estimate 
    */
    void provideTimeStepEstimate(const double dtEstimate) {
      dtEstimate_ = std::min(dtEstimate_, dtEstimate);
      synced_ = false ;
    }

    /** \brief count current time step a not valid */
    void invalidateTimeStep() 
    {
      DUNE_THROW(InvalidStateException,"TimeStep invalid!");
    }
    
    /** \brief restore time and timestep from outside 
         (i.e. from former calculation)  
         \param[in] time new time 
         \param[in] timeStep new time step counter 
    */
    void restore(const double time, const int timeStep )  
    { 
      time_ = time; 
      timeStep_ = timeStep;
    }
    
    /** \brief restore time and timestep from outside 
         (i.e. from former calculation)  
         \param[in] time new time 
         \param[in] timeStep new time step counter 
    */
    void setTime(const double time, const int timeStep ) DUNE_DEPRECATED
    { 
      restore(time, timeStep);
    }

    /** \brief augment time , i.e. \f$t = t + \triangle t\f$ and
        increase time step counter  */
    double augmentTime()
    { 
      // if already syncronized do nothing 
      if( syncronized() ) return time_ ;

      time_ += deltaT(); 
      ++timeStep_;
      return time_;
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

    /** \brief sets time step size to size of dt 
        \param dt new time step size 
    */
    void setDeltaT (const double dt)
    {
      resetTimeStepEstimate();
      provideTimeStepEstimate(dt);
      syncTimeStep();
    }
    
    /** \brief  syncronize time step, i.e. set timeStep to values of current 
        estimate and reset estimate 
    */
    void syncTimeStep() 
    {
      // if already syncronized do nothing 
      if( syncronized() ) return ;

      // save current time step 
      dt_ = dtEstimate_;
      // reset estimate 
      resetTimeStepEstimate();

      // now up to date 
      synced_ = true ;
    }

    /** \brief return current time step counter 
        \return current time step counter 
    */
    int timeStep () const {  return timeStep_;  }

    /** \brief returns true if TimeProvider is in syncronized state,
        i.e. after next has been called. 
    */
    bool syncronized () const { return synced_; }

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
    bool synced_;
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
    
    /** @copydoc TimeProvider::time */ 
    double time() const { return tp_.time(); }
    
    /** @copydoc TimeProvider::provideTimeStepEstimate */
    void provideTimeStepEstimate(const double dtEstimate) 
    {
      tp_.provideTimeStepEstimate(dtEstimate);
    }

    /** @copydoc TimeProvider::resetTimeStepEstimate */
    void resetTimeStepEstimate() 
    {
      tp_.resetTimeStepEstimate(); 
    }

    /** @copydoc TimeProvider::init */
    void init() 
    {
      syncTimeStep();
    }
    
    /** @copydoc TimeProvider::next */
    void next() 
    {
      // if already synced do nothing 
      if( tp_.syncronized() ) return ;
      
      // increase time by delta t 
      augmentTime(); 
      // sync new time step size 
      syncTimeStep();
    }
    
    /** @copydoc TimeProvider::invalidateTimeStep */
    void invalidateTimeStep() 
    {
      tp_.invalidateTimeStep();
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
    int timeStep() const { return tp_.timeStep(); }

    /** @copydoc TimeProvider::syncTimeStep 
        
        \note Here a global communication to minimize the time step size
        over all time step sizes from all processors is done. 
    */ 
    void syncTimeStep() 
    {
      // do nothing if already synced 
      if( tp_.syncronized() ) return ;

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
      // increase time (if not syncronized )
      return tp_.augmentTime();
    }

  protected:
    //! communicator  
    const CommunicatorType& comm_;
    //! serial time provider 
    TimeProvider& tp_;
  };
  
} // end namespace Dune
#endif
