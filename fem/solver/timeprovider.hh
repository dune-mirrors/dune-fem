#ifndef DUNE_FEM_TIMEPROVIDER_HH
#define DUNE_FEM_TIMEPROVIDER_HH

//- system includes 
#include <limits>
#include <cassert>

//- Dune includes 
#include <dune/common/exceptions.hh>

#include <dune/fem/misc/commhelper.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/file/asciiparser.hh>

namespace Dune
{

  /** \class   TimeProviderBase
   *  \ingroup ODESolver
   *  \brief   gereral base for time providers
   */
  class TimeProviderBase
  {
    typedef TimeProviderBase ThisType;

  protected:
    double time_;
    int timeStep_;
    double dt_;
    bool valid_;
    double dtEstimate_;

  public:
    inline TimeProviderBase ()
    : time_( Parameter :: getValue( "fem.timeprovider.starttime",
                                    (double)0.0 ) ),
      timeStep_( 0 ),
      valid_( false )
    {
      initTimeStepEstimate();
    }

    inline explicit TimeProviderBase ( const double startTime )
    : time_( startTime ),
      timeStep_( 0 ),
      valid_( false )
    {
      initTimeStepEstimate();
    }

  private:
    TimeProviderBase ( const ThisType & );
    ThisType &operator= ( const ThisType & );

  public:
    /** \brief obtain the current time
     *
     *  \returns the current time
     */
    inline double time () const
    {
      return time_;
    }
    
    /** \brief obtain number of the current time step
     *
     *  \return the current time step counter
     */
    inline int timeStep () const
    {
      assert( timeStepValid() );
      return timeStep_;
    }
 
    /** \brief obtain the size of the current time step
     *
     *  \returns the size of the current time step
     */
    inline double deltaT () const
    {
      assert( timeStepValid() );
      return dt_;
    }

    inline bool timeStepValid () const
    {
      return valid_;
    }
   
    /** \brief set time step estimate to minimum of given value and
               internal time step estiamte 
         \param[in] dtEstimate time step size estimate 
    */
    inline void provideTimeStepEstimate ( const double dtEstimate )
    {
      dtEstimate_ = std :: min( dtEstimate_, dtEstimate );
    }

    /** \brief count current time step a not valid */
    inline void invalidateTimeStep ()
    {
      valid_ = false;
    }

  protected:
    inline void advance ()
    {
      time_ += deltaT();
      ++timeStep_;
    }

    inline void initTimeStepEstimate ()
    {
      dtEstimate_ = std :: numeric_limits< double > :: max();
    }
  };



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
  template< class CommProvider = DefaultCollectiveCommunicationType >
  class TimeProvider;



  template< class C >
  class TimeProvider< CollectiveCommunication< C > >
  : public TimeProviderBase
  {
    typedef TimeProvider< CollectiveCommunication< C > > ThisType;
    typedef TimeProviderBase BaseType;

  public:
    typedef CollectiveCommunication< C > CollectiveCommunicationType;

  protected:
    typedef CollectiveCommunicationHelper< CollectiveCommunicationType >
      CollectiveCommHelperType;
    
  protected:
    const CollectiveCommunicationType &comm_;
    const double cfl_;

    using BaseType :: dt_;
    using BaseType :: dtEstimate_;
    using BaseType :: valid_;
    using BaseType :: timeStep_;

  public:
    /** \brief default constructor
     *
     *  \param[in]  comm  collective communication (optional)
     */
    inline explicit
    TimeProvider ( const CollectiveCommunicationType &comm
                     = CollectiveCommHelperType :: defaultCommunication() )
    : BaseType(),
      comm_( comm ),
      cfl_( Parameter :: getValidValue( "fem.timeprovider.factor", (double)1.0,
                                        ValidateGreater< double >( 0.0 ) ) )
    {}

    /** \brief constructor taking start time
     *
     *  \param[in]  startTime  initial time
     *  \param[in]  comm       collective communication (optional)
     */
    inline explicit
    TimeProvider ( const double startTime,
                   const CollectiveCommunicationType &comm
                     = CollectiveCommHelperType :: defaultCommunication() )
    : BaseType( startTime ),
      comm_( comm ),
      cfl_( Parameter :: getValidValue( "fem.timeprovider.timestep.factor", (double)1.0,
                                        ValidateGreater< double >( 0.0 ) ) )
    {}
    
    /** \brief constructor taking start time and CFL constant
     *
     *  \param[in]  startTime  initial time
     *  \param[in]  cfl        CFL constant
     *  \param[in]  comm       collective communication (optional)
     */
    inline
    TimeProvider ( const double startTime,
                   const double cfl,
                   const CollectiveCommunicationType &comm
                     = CollectiveCommHelperType :: defaultCommunication() )
    : BaseType( startTime ),
      comm_( comm ),
      cfl_( cfl )
    {}
    
  private:
    TimeProvider ( const ThisType & );
    ThisType &operator= ( const ThisType & );

  public:
    /** \brief init dt with given estimate
     *
     *  \param[in]  maxTimeStep  maximum allowd time step (default to
     *                           numeric_limits< double > :: max())
     */
    void init( double maxTimeStep = std :: numeric_limits< double > :: max() ) 
    {
      provideTimeStepEstimate( maxTimeStep );
      initTimeStep();
    }
    
    /** \brief goto next time step
     * 
     *  \param[in]  maxTimeStep  maximum allowed time step (defaults to
     *                           numeric_limits< double > :: max())
     */
    void next ( double maxTimeStep = std :: numeric_limits< double > :: max() ) 
    {
      provideTimeStepEstimate( maxTimeStep );

      advance();
      initTimeStep();
    }

  protected:
    using BaseType :: advance;
    using BaseType :: initTimeStepEstimate;

    inline void initTimeStep ()
    {
      dt_ = cfl_ * dtEstimate_;
      dt_ = comm_.min( dt_ );
      assert( dt_ > 0.0 );
      valid_ = true;

      initTimeStepEstimate();
    }


    // old methods, possibly deprecated in future
    // ------------------------------------------
    
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
    double augmentTime() DUNE_DEPRECATED
    { 
      advance();
      return time_;
    }
    
    /** \brief reset set time step estimate 
        by setting ti to big value 
    */
    void resetTimeStepEstimate() DUNE_DEPRECATED
    {
      initTimeStepEstimate();
    }
    
    /** \brief  return time step size estimate  
        \return time step size estimate  
    */
    double timeStepEstimate() const DUNE_DEPRECATED
    {
      return dtEstimate_;
    }

    /** \brief  return cfl number 
        \return cfl number 
    */
    double cfl () const DUNE_DEPRECATED
    {
      return cfl_;
    }
    
    /** \brief set internal cfl to given value
        \param[in] cfl new cfl number 
    */
    void setCfl(const double cfl) DUNE_DEPRECATED DUNE_DEPRECATED
    {
      //cfl_ = cfl;
    }
    
    /** \brief set internal cfl to minimum of given value and internal
        cfl number 
        \param[in] cfl cfl estimate 
    */
    void provideCflEstimate(const double cfl) DUNE_DEPRECATED
    {
      //cfl_ = std::min(cfl_, cfl );
    }
    
    /** \brief sets time step size to size of dt 
        \param dt new time step size 
    */
    void setDeltaT (const double dt) DUNE_DEPRECATED
    {
      //resetTimeStepEstimate();
      //provideTimeStepEstimate(dt);
      //syncTimeStep();
    }
    
    /** \brief  syncronize time step, i.e. set timeStep to values of current 
        estimate and reset estimate 
    */
    void syncTimeStep() DUNE_DEPRECATED 
    {
      initTimeStep();
    }

    /** \brief returns true if TimeProvider is in syncronized state,
        i.e. after next has been called. 
    */
    bool syncronized () const DUNE_DEPRECATED
    {
      return false;
    }
  };



  template< class CommProvider >
  class TimeProvider
  : public TimeProvider
    < typename CommProvider :: Traits :: CollectiveCommunication >
  {
    typedef TimeProvider< CommProvider > ThisType;
    typedef TimeProvider
      < typename CommProvider :: Traits :: CollectiveCommunication >
      BaseType;

  public:
    typedef typename CommProvider :: Traits :: CollectiveCommunication
      CollectiveCommunicationType;

  public:
    inline explicit TimeProvider ( const CommProvider &comm )
    : BaseType( comm.comm() )
    {}

    inline TimeProvider ( const double startTime,
                          const CommProvider &comm )
    : BaseType( startTime, comm.comm() )
    {}
    
    inline TimeProvider ( const double startTime,
                          const double cfl,
                          const CommProvider &comm )
    : BaseType( startTime, cfl, comm.comm() )
    {}
  };

} // end namespace Dune

#endif
