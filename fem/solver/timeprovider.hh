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
#include <dune/fem/io/file/persistencemanager.hh>
#include <dune/fem/misc/femtuples.hh>

namespace Dune
{

  /** \class   TimeProviderBase
   *  \ingroup ODESolver
   *  \brief   gereral base for time providers 
   *
   *  This class consists of the methods required for example in the
   *  ODE Solvers, e.g., provideTimeStepEstimate and 
   *  provideTimeStepUpperBound.
   *  InvalidateTimeStep can be used to mark this time step as invalid.
   *  Furthermore, method for accessing the simulation time, the
   *  time step counter and the time step size are provided.
   *
   *  The derived class TimeProvider provides the additional method
   *  required for implementing a time loop.
   *  
   */
  class TimeProviderBase : public AutoPersistentObject
  {
    typedef TimeProviderBase ThisType;

  protected:
    double time_;
    int timeStep_;
    double dt_;
    bool valid_;
    bool dtValid_;
    double dtEstimate_;
    double dtUpperBound_;

  public:
    inline TimeProviderBase ()
    : time_( Parameter :: getValue( "fem.timeprovider.starttime",
                                    (double)0.0 ) ),
      timeStep_( 0 ),
      valid_( false ),
      dtValid_( false )
    {
      initTimeStepEstimate();
    }

    inline explicit TimeProviderBase ( const double startTime )
    : time_( startTime ),
      timeStep_( 0 ),
      valid_( false ),
      dtValid_( false )
    {
      initTimeStepEstimate();
    }

    virtual ~TimeProviderBase() {}

    void backup() const {
      Tuple<const double&,const int&,const double&,const bool&,const double&>
        values(time_,timeStep_,dt_,valid_,dtEstimate_);
      PersistenceManager::backupValue("timeprovider",values);
    }
    void restore() {
      Tuple<double&,int&,double&,bool&,double&>
        values(time_,timeStep_,dt_,valid_,dtEstimate_);
      PersistenceManager::restoreValue("timeprovider",values);
      dtValid_=true;
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

    /** \brief set time step estimate to minimum of given value and
               internal time step estiamte 
         \param[in] dtEstimate time step size estimate 
    */
    inline void provideTimeStepEstimate ( const double dtEstimate )
    {
      dtEstimate_ = std :: min( dtEstimate_, dtEstimate );
      dtValid_ = true;
    }
    /** \brief set upper bound for time step to minimum of given value and
               internal bound
         \param[in] upperBound time step size estimate 
    */
    inline void provideTimeStepUpperBound ( const double upperBound )
    {
      dtUpperBound_ = std :: min( dtUpperBound_, upperBound );
      dtValid_ = true;
    }
    
    /** \brief count current time step a not valid */
    inline void invalidateTimeStep ()
    {
      valid_ = false;
    }

    /** \brief return if this time step should be used */
    inline bool timeStepValid () const
    {
      return valid_;
    }
   

  protected:
    inline void advance ()
    {
      if (timeStepValid()) {
        time_ += deltaT();
        ++timeStep_;
      }
    }

    inline void initTimeStepEstimate ()
    {
      dtEstimate_ = std :: numeric_limits< double > :: max();
      dtUpperBound_ = std :: numeric_limits< double > :: max();
      dtValid_ = false;
    }
  };



  /** 
   *  \ingroup ODESolver
   *  \brief   manager for global simulation time of time-dependent solutions
   *
   *  When calculating possibly multiple time-dependent solutions, it is often
   *  necessary to use the same time in all calculations. This means that we
   *  have to use the same time step for all our calculations. A TimeProvider
   *  keeps track of this information in a simple and unified way.
   *
   *  An example of a time loop could look as follows:
   *  \code
   *  // create time provider
   *  TimeProvider tp( startTime );
   *
   *  SpaceOperator spaceOperator;
   *  typedef SpaceOperator::DestinationType DestinationType;
   *  OdeSolver<DestinationType> odeSolver(spaceOperator,tp,order);
   *
   *  DestinationType U;
   *  initialize(U);
   *
   *  // set the initial time step estimate
   *  odeSolver.initialize( U );
   *
   *  // time loop
   *  for( tp.init(); tp.time() < endTime; tp.next() )
   *  {
   *    // do calculation
   *    odeSolver.solve(U);
   *  }
   *  \endcode
   *
   *  Within the time loop, both tp.time() and tp.deltaT() are fixed and cannot
   *  be altered and an the next time step should be fixed in the loop,
   *  e.g., in the method solve of the ode solver an upper estimate
   *  for the next time step is provided; if more than one time
   *  step restriction has to be imposed, the minimum is taken for
   *  the next time step.
   *  By calling the method provideTimeStepEstimate(maxDt) in the body of the
   *  loop an upper estimate for the next time step can be supplied;
   *  to fix the next time step (ignoring the estimates) an optinal
   *  argument can be passed to the next method on the
   *  Dune::TimeProvider.
   *
   *  Obviously, we need to provide an initial estimate. In the above example,
   *  this is done by the initialize method of the ODE solver. In tp.init(),
   *  the first time step (deltaT) is set based on the estimate and 
   *  this value can also be fixed independent of the estimate through
   *  an optinal argument. The following loop would fix the time step
   *  to 1e-3
   *  \code
   *  for( tp.init(1e-3); tp.time() < endTime; tp.next(1e-3) )
   *  {
   *    // do calculation
   *    odeSolver.solve(U);
   *  }
   *  \endcode
   *
   *  In order to allow the user to incfluence the calculation of the next time
   *  step from the estimate, the time provider also maintains an additional
   *  factor (which is constant during the entire simulation). 
   *  Therefore the acctual time step used, is calculated as follows:
   *  \f[
   *  \mathrm{deltaT} = \mathrm{factor} * \mathrm{timeStepEstimate}.
   *  \f]
   *  Therefore in the above example 1e-3 might not be the acctual
   *  time step depending on the value of the factor in the
   *  TimeProvider.
   *  The default value for this factor is equal to one but can be changed
   *  either during the construction of the Dune::TimeProvider or
   *  by using the parameter \c fem.timeprovider.factor.
   *  A further parameter read by the Dune::TimeProvider is
   *  fem.timeprovider.starttime defining the starting time of
   *    the simulation (default is zero).
   *
   *  The most general implementation is given in the class
   *  Dune::TimeProvider< CollectiveCommunication< C > >  which
   *  takes a Dune::CollectiveCommunication instance in the 
   *  constructor which is used in parallel computations is
   *  syncronize the time step. It defaults to 
   *  Dune::CollectiveCommHelperType :: defaultCommunication()
   *  and also works for seriell runs where the template argument
   *  does not have to be prescribed.
   *  If the communication manager from a given grid is to be used
   *  the class Dune::GridTimeProvider using the GridType as
   *  template argument can be used instead, with the same
   *  functionality.
   *
     \parametername \c fem.timeprovider.factor \n
                    multiplication factor to use for each time step;
                    defaults to 1.
     \parametername \c fem.timeprovider.starttime \n
                    time used for initializing the starting time
                    defaults to zero.
   */
  template< class CommProvider = DefaultCollectiveCommunicationType >
  class TimeProvider {
  };
  
  /** \ingroup ODESolver
   *  \brief   the basic Dune::TimeProvider implementation.
   *
   *  This implementation of a timeprovider takes a CollectiveCommunicate 
   *  for parallel runs which default to a default communicator
   *  which also works for seriellel simulations.
   *
   */
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
    using BaseType :: dtUpperBound_;
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
      cfl_( Parameter :: getValidValue( "fem.timeprovider.factor", (double)1.0,
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

    virtual~TimeProvider() {}
    
  private:
    TimeProvider ( const ThisType & );
    ThisType &operator= ( const ThisType & );

  public:
    /** \brief init dt with time step estimate
     */
    void init() 
    {
      initTimeStep(dtEstimate_);
    }
    /** \brief init dt with provided time step
     *
     *  \param[in]  timeStep  value of the first time step (is multiplied with
     *                        factor)
     */
    void init( double timeStep ) 
    {
      initTimeStep(timeStep);
    }
    
    /** \brief goto next time step
     *
     * Sets the size of the next time step to the current time step estimate
     * and sets the estimate to infinity.
     */
    void next ( ) 
    {
      assert(this->dtValid_);
      advance();
      initTimeStep(dtEstimate_);
    }
    /** \brief goto next time step
     * 
     * Sets the size of the next time step to the provided time step value
     * and sets the estimate to infinity.
     * 
     *  \param[in]  timeStep  value of the next time step (is multiplied with
     *                        factor)
     */
    void next ( double timeStep ) 
    {
      advance();
      initTimeStep(timeStep);
    }

    /** \brief  return the global factor number 
        \return time step factor 
    */
    double factor () const
    {
      return cfl_;
    }

  protected:
    using BaseType :: advance;
    using BaseType :: initTimeStepEstimate;

    inline void initTimeStep (double dtEstimate)
    {
      dt_ = std::min(cfl_ * dtEstimate,dtUpperBound_);
      dt_ = comm_.min( dt_ );
      assert( dt_ > 0.0 );
      valid_ = true;

      initTimeStepEstimate();
    }

  public:
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

    virtual void backup() const {
      BaseType::backup();
    }
    virtual void restore() {
      BaseType::restore();
      const_cast<double&>(cfl_) 
        =  Parameter :: getValidValue<double>( "fem.timeprovider.factor",
                           ValidateGreater< double >( 0.0 ) ) ;
    }
  };

  /** \class   GridTimeProvider
   *  \ingroup ODESolver
   *  \brief   the same functionality as the Dune::TimeProvider.
   *
   *  This implementation of a timeprovider takes the CollectiveCommunicate 
   *  from a Dune::Grid instance.
   */
  template< class Grid >
  class GridTimeProvider
  : public TimeProvider
    < typename Grid :: Traits :: CollectiveCommunication >
  {
    typedef GridTimeProvider< Grid > ThisType;
    typedef TimeProvider
      < typename Grid :: Traits :: CollectiveCommunication >
      BaseType;

  public:
    typedef typename Grid :: Traits :: CollectiveCommunication
      CollectiveCommunicationType;

  public:
    inline explicit GridTimeProvider ( const Grid &grid )
    : BaseType( grid.comm() )
    {}

    inline GridTimeProvider ( const double startTime,
                              const Grid &grid )
    : BaseType( startTime, grid.comm() )
    {}
    
    inline GridTimeProvider ( const double startTime,
                              const double cfl,
                              const Grid &grid )
    : BaseType( startTime, cfl, grid.comm() )
    {}
    
    virtual ~GridTimeProvider() {}
  };

} // end namespace Dune

#endif
