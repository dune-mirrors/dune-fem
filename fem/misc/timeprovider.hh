#ifndef DUNE_TIMEPROVIDER_HH
#define DUNE_TIMEPROVIDER_HH

//- system includes 
#include <limits>
#include <cassert>

//- Dune includes 
#include <dune/common/exceptions.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/file/asciiparser.hh>

namespace Dune
{

  /** \class   TimeProvider
   *  \ingroup ODESolver
   *  \brief   manager for global simulation time of time-dependent solutions
   *
   *  When calculating possibly multiple time-dependent solutions, it is often
   *  necessary to use the same time in all calculations. This means that we
   *  have to use the same time step for all our calculations. A TimeProvider
   *  keeps track of this information in a simple and unified way.
   *
   *  An example time loop could look as follows:
   *  \code
   *  // create time provider
   *  TimeProvider tp( startTime );
   *
   *  // set the initial time step estimate
   *  odeSolver.initialize( U );
   *
   *  // time loop
   *  for( tp.init(); tp.time() < endTime; tp.next() )
   *  {
   *    // do calculation
   *  }
   *  \endcode
   *
   *  Within the time loop, both tp.time() and tp.deltaT() are fixed and cannot
   *  be altered. Within the time loop, the user provides (upper) estimates for
   *  the next time step. This is usually implicitly done by the ODE solvers.
   *  The minimum of all those estimates is taken as the basis for the next
   *  time step.
   *
   *  Obviously, we need to provide an initial estimate. In the above example,
   *  this is done by the initialize method of the ODE solver. On tp.init(),
   *  the first time step (deltaT) is set based on the estimate.
   *
   *  In order to allow the user to incluence the calculation of the next time
   *  step from the estimate, the time provider also maintains a CFL constant
   *  (which is constant during the entire simulation). The time stap is then
   *  calculated as follows:
   *  \f[
   *  \mathrm{deltaT} = \mathrm{cfl} * \mathrm{timeStepEstimate}.
   *  \f]
   *
   *  \remark There exist two implementations of a time provider, one for
   *          serial runs (TimeProvider) and a wrapper for parallel runs
   *          (ParallelTimeProvider).
   */
  class TimeProvider 
  {
    typedef TimeProvider ThisType;

  public:
    /** \brief default constructor */
    inline TimeProvider()
    : time_( Parameter :: getValue( "fem.timeprovider.starttime",
                                    (double)0.0 ) ),
      dt_( -1.0 ),
      dtEstimate_( 0.0 ),
      cfl_( Parameter :: getValidValue( "fem.timeprovider.cfl", (double)1.0,
                                        ValidateGreater< double >( 0.0 ) ) ),
      timeStep_( 0 ),
      synced_( false )
    {
      resetTimeStepEstimate();
    }

    /** \brief constructor taking start time
     *
     *  \param[in]  startTime  initial time
     */
    inline explicit TimeProvider( const double startTime )
    : time_( startTime ),
      dt_( -1.0 ),
      dtEstimate_( 0.0 ),
      cfl_( Parameter :: getValidValue( "fem.timeprovider.timestep.factor", (double)1.0,
                                        ValidateGreater< double >( 0.0 ) ) ),
      timeStep_( 0 ),
      synced_( false )
    {
      resetTimeStepEstimate();
    }
    
    /** \brief constructor taking start time and CFL constant
     *
     *  \param[in]  startTime  initial time
     *  \param[in]  cfl        CFL constant
     */
    TimeProvider( const double startTime, const double cfl ) DUNE_DEPRECATED
    : time_( startTime ),
      dt_( -1.0 ),
      dtEstimate_( 0.0 ),
      cfl_( cfl ),
      timeStep_( 0 ),
      synced_( false )
    {
      resetTimeStepEstimate();
    }
    
    //! destructor 
    ~TimeProvider() {}

    /** \brief init dt with given estimate
     *
     *  \param[in]  maxTimeStep  maximum allowd time step (default to
     *                           numeric_limits< double > :: max())
     */
    void init( double maxTimeStep = std :: numeric_limits< double > :: max() ) 
    {
      provideTimeStepEstimate( maxTimeStep );
      syncTimeStep();
    }
    
    /** \brief goto next time step
     * 
     *  \param[in]  maxTimeStep  maximum allowed time step (defaults to
     *                           numeric_limits< double > :: max())
     */
    void next ( double maxTimeStep = std :: numeric_limits< double > :: max() ) 
    {
      provideTimeStepEstimate( maxTimeStep );

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
    TimeProvider( const ThisType & );
    ThisType &operator=( const ThisType & );
    
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

      \deprecated
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
                         const int rank = 0) DUNE_DEPRECATED
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
