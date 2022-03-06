#ifndef DUNE_FEM_TIMEPROVIDER_HH
#define DUNE_FEM_TIMEPROVIDER_HH

#include <cassert>
#include <limits>
#include <tuple>

#include <dune/fem/io/file/persistencemanager.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/space/common/dofmanager.hh>

#include <dune/fem/misc/mpimanager.hh>

namespace Dune
{

  namespace Fem
  {

    /** \class   TimeProviderBase
     *  \ingroup ODESolver
     *  \brief   general base for time providers
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

    public:
      inline TimeProviderBase ( const ParameterReader &parameter = Parameter::container() )
      : time_( parameter.getValue( "fem.timeprovider.starttime",
                                      static_cast<double>(0.0) ) ),
        timeStep_( 0 ),
        dt_( 0.0 ),
        invdt_( HUGE_VAL ),
        valid_( false ),
        dtEstimateValid_( false ),
        parameter_( parameter )
      {
        initTimeStepEstimate();
      }

      inline explicit TimeProviderBase ( const double startTime, const ParameterReader &parameter = Parameter::container() )
      : time_( startTime ),
        timeStep_( 0 ),
        dt_( 0.0 ),
        invdt_( HUGE_VAL ),
        valid_( false ),
        dtEstimateValid_( false ),
        parameter_( parameter )
      {
        initTimeStepEstimate();
      }

      inline virtual ~TimeProviderBase()
      {}

      void backup() const
      {
        std::tuple<const double&,const int&,const double&,const bool&,const double&>
          values(time_,timeStep_,dt_,valid_,dtEstimate_);
        PersistenceManager::backupValue("timeprovider",values);
      }

      void restore()
      {
        std::tuple<double&,int&,double&,bool&,double&>
          values(time_,timeStep_,dt_,valid_,dtEstimate_);
        PersistenceManager::restoreValue("timeprovider",values);
        dtEstimateValid_ = true;
        invdt_ = 1.0 / dt_;
      }


      TimeProviderBase ( const ThisType & ) = delete;

      ThisType &operator= ( const ThisType & ) = delete;


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

      /** \brief obtain the size of the inverse of the current time step
       *
       *  \returns the size of the inverse of the current time step
       */
      inline double inverseDeltaT () const
      {
        assert( timeStepValid() );
        return invdt_;
      }

      /** \brief obtain current estimate on time step
       *
       *  \returns the current estimate for the time step
       */
      inline double timeStepEstimate () const
      {
        return dtEstimate_;
      }

      /** \brief set time step estimate to minimum of given value and
                 internal time step estiamte
           \param[in] dtEstimate time step size estimate
      */
      inline void provideTimeStepEstimate ( const double dtEstimate )
      {
        dtEstimate_ = std::min( dtEstimate_, dtEstimate );
        dtEstimateValid_ = true;
      }
      /** \brief set upper bound for time step to minimum of given value and
                 internal bound
           \param[in] upperBound time step size estimate
      */
      inline void provideTimeStepUpperBound ( const double upperBound )
      {
        dtUpperBound_ = std::min( dtUpperBound_, upperBound );
        dtEstimateValid_ = true;
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
      double time_;
      int timeStep_;
      double dt_;
      double invdt_;
      bool valid_;
      bool dtEstimateValid_;
      double dtEstimate_;
      double dtUpperBound_;
      ParameterReader parameter_;

      void advance ()
      {
        if( timeStepValid() )
        {
          time_ += deltaT();
          ++timeStep_;
        }
      }

      void initTimeStepEstimate ()
      {
        dtEstimate_ = std::numeric_limits< double >::max();
        dtUpperBound_ = std::numeric_limits< double >::max();
        dtEstimateValid_ = false;
      }
    };


    /** \class   FixedStepTimerProvider
     *  \ingroup ODESolver
     *  \brief   simple time provider for constant time steps
     *
     *  An example of a time loop could look as follows:
     *  \code
     *  // create time provider
     *  FixedStepTimeProvider tp( startTime, timeStepSize );
     *
     *  // time loop
     *  for( ; tp.time() < endTime; tp.next() )
     *  {
     *    // do stuff
     *  }
     *  \endcode
     *
     */
    template< class Communication = typename MPIManager::Communication >
    class FixedStepTimeProvider
    : public TimeProviderBase
    {
      typedef FixedStepTimeProvider< Communication > ThisType;
      typedef TimeProviderBase BaseType;

    public:
      typedef Communication CommunicationType;

      /** \brief constructor
       *
       *  \param[in]  startTime     initial time
       *  \param[in]  timeStepSize  time step size
       *  \param[in]  comm          collective communication (default Dune::Fem::MPIManager::comm())
       */
      explicit FixedStepTimeProvider ( const double startTime, const double timeStepSize,
                                       const CommunicationType &comm,
                                       const ParameterReader &parameter = Parameter::container() )
        : BaseType( startTime, parameter ), comm_( comm )
      {
        dt_ = timeStepSize;
        initTimeStep();
      }

      explicit FixedStepTimeProvider ( const double startTime, const double timeStepSize,
                                       const ParameterReader &parameter = Parameter::container() )
        : BaseType( startTime, parameter ), comm_( MPIManager::comm() )
      {
        dt_ = timeStepSize;
        initTimeStep();
      }

      /** \brief constructor
       *
       *  \param[in]  comm  collective communication (default Dune::Fem::MPIManager::comm())
       *
       *  The initial time need to be provided using the parameter fem.timeprovider.starttime while
       *  the time step size need to be provided using the parameter fem.timeprovider.fixedtimestep.
       */
      explicit FixedStepTimeProvider ( const ParameterReader &parameter = Parameter::container() )
        : BaseType( parameter.getValue< double>( "fem.timeprovider.starttime", 0.0 ), parameter ), comm_( MPIManager::comm() )
      {
        dt_ = parameter.getValidValue< double >("fem.timeprovider.fixedtimestep", [] ( double v ) { return v > 0.0;} );
        initTimeStep();
      }

      explicit FixedStepTimeProvider ( const CommunicationType &comm,
                                       const ParameterReader &parameter = Parameter::container() )
        : BaseType( parameter.getValue< double>( "fem.timeprovider.starttime", 0.0 ), parameter ), comm_( comm )
      {
        dt_ = parameter.getValidValue< double >("fem.timeprovider.fixedtimestep", [] ( double v ) { return v > 0.0;} );
        initTimeStep();
      }
      virtual ~FixedStepTimeProvider () {}

      FixedStepTimeProvider ( const ThisType & ) = delete;
      FixedStepTimeProvider ( ThisType && ) = delete;
      ThisType &operator= ( const ThisType & ) = delete;
      ThisType &operator= ( ThisType && ) = delete;

      /** \brief goto next time step */
      void next ()
      {
        if( !timeStepValid() )
          DUNE_THROW( InvalidStateException, "Invalid Time Step in FixedStepTimeProvider" );
        advance();
        initTimeStep();
      }

    protected:
      using BaseType::advance;
      using BaseType::initTimeStepEstimate;
      using BaseType::time_;
      using BaseType::dt_;
      using BaseType::valid_;

      inline void initTimeStep ()
      {
        valid_ = true;
        initTimeStepEstimate();
      }

      const CommunicationType &comm_;
    };


    /**
       \ingroup ODESolver
       \brief   manager for global simulation time of time-dependent solutions

       When calculating possibly multiple time-dependent solutions, it is often
       necessary to use the same time in all calculations. This means that we
       have to use the same time step for all our calculations. A TimeProvider
       keeps track of this information in a simple and unified way.

       An example of a time loop could look as follows:
       \code
       // create time provider
       TimeProvider tp( startTime );

       SpaceOperator spaceOperator;
       typedef SpaceOperator::DestinationType DestinationType;
       OdeSolver<DestinationType> odeSolver(spaceOperator,tp,order);

       DestinationType U;
       initialize(U);

       // set the initial time step estimate
       odeSolver.initialize( U );

       // time loop
       for( tp.init(); tp.time() < endTime; tp.next() )
       {
         // do calculation
         odeSolver.solve(U);
       }
       \endcode

       Within the time loop, both tp.time() and tp.deltaT() are fixed and cannot
       be altered and an the next time step should be fixed in the loop,
       e.g., in the method solve of the ode solver an upper estimate
       for the next time step is provided; if more than one time
       step restriction has to be imposed, the minimum is taken for
       the next time step.
       By calling the method provideTimeStepEstimate(maxDt) in the body of the
       loop an upper estimate for the next time step can be supplied;
       to fix the next time step (ignoring the estimates) an optinal
       argument can be passed to the next method on the
       Dune::Fem::TimeProvider.

       Obviously, we need to provide an initial estimate. In the above example,
       this is done by the initialize method of the ODE solver. In tp.init(),
       the first time step (deltaT) is set based on the estimate and
       this value can also be fixed independent of the estimate through
       an optional argument. The following loop would fix the time step
       to 1e-3
       \code
       for( tp.init(1e-3); tp.time() < endTime; tp.next(1e-3) )
       {
         // do calculation
         odeSolver.solve(U);
       }
       \endcode

       In order to allow the user to influence the calculation of the next time
       step from the estimate, the time provider also maintains an additional
       factor (which is constant during the entire simulation).
       Therefore the actual time step used, is calculated as follows:
       \f[
       \mathrm{deltaT} = \mathrm{factor} * \mathrm{timeStepEstimate}.
       \f]
       Therefore in the above example 1e-3 might not be the acctual
       time step depending on the value of the factor in the
       TimeProvider.
       The default value for this factor is equal to one but can be changed
       either during the construction of the Dune::TimeProvider or
       by using the parameter \c fem.timeprovider.factor.
       A further parameter read by the Dune::TimeProvider is
       fem.timeprovider.starttime defining the starting time of
       the simulation (default is zero).

       The most general implementation is given in the class
       Dune::Fem::TimeProvider< Communication< C > >  which
       takes a Dune::Communication instance in the
       constructor which is used in parallel computations is
       syncronize the time step. It defaults to
       Dune::Fem::MPIManager::comm() and also works for serial runs.

       If the communication manager from a given grid is to be used
       the class Dune::Fem::GridTimeProvider using the GridType as
       template argument can be used instead, with the same
       functionality.

       \parametername \c fem.timeprovider.factor \n
                      multiplication factor to use for each time step;
                      defaults to 1.
       \parametername \c fem.timeprovider.starttime \n
                      time used for initializing the starting time
                      defaults to zero.
       \parametername \c fem.timeprovider.updatestep \n
                      only do the update of the time step size
                      every 'updatestep' to avoid the
                      expensive communication to achieve this
                      (for testing only);
                      defaults to 1
     */
    template< class Communication = typename MPIManager::Communication >
    class TimeProvider
    : public TimeProviderBase
    {
      typedef TimeProvider< Communication > ThisType;
      typedef TimeProviderBase BaseType;

    public:
      typedef Communication CommunicationType;

    protected:

      using BaseType::parameter_;

      inline double getCflFactor() const
      {
        return parameter_.getValidValue( "fem.timeprovider.factor", static_cast<double>(1.0),
            [] ( double val ) { return val > 0.0; } );
      }

      inline int getUpdateStep () const
      {
        return parameter_.getValidValue( "fem.timeprovider.updatestep", static_cast<int>(1),
            [] ( int step ) { return step > 0; } );
      }

    public:
      /** \brief default constructor
       *
       *  \param[in]  comm  collective communication (default Dune::Fem::MPIManager::comm())
       */
      explicit TimeProvider ( const ParameterReader &parameter = Parameter::container() )
      : BaseType( parameter ),
        comm_( MPIManager::comm() ),
        cfl_( getCflFactor() ),
        updateStep_( getUpdateStep() ),
        counter_( updateStep_ )
      {}

      explicit TimeProvider ( const CommunicationType &comm, const ParameterReader &parameter = Parameter::container() )
      : BaseType( parameter ),
        comm_( comm ),
        cfl_( getCflFactor() ),
        updateStep_( getUpdateStep() ),
        counter_( updateStep_ )
      {}

      /** \brief constructor taking start time
       *
       *  \param[in]  startTime  initial time
       *  \param[in]  comm       collective communication (default Dune::Fem::MPIManager::comm())
       */
      explicit TimeProvider ( const double startTime,
                              const CommunicationType &comm = MPIManager::comm() )
      : BaseType( startTime ),
        comm_( comm ),
        cfl_( getCflFactor() ),
        updateStep_( getUpdateStep() ),
        counter_( updateStep_ )
      {}

      /** \brief constructor taking start time and CFL constant
       *
       *  \param[in]  startTime  initial time
       *  \param[in]  cfl        CFL constant
       *  \param[in]  comm       collective communication (default Dune::Fem::MPIManager::comm())
       */
      TimeProvider ( const double startTime,
                     const double cfl,
                     const CommunicationType &comm = MPIManager::comm() )
      : BaseType( startTime ),
        comm_( comm ),
        cfl_( cfl ),
        updateStep_( 1 ),
        counter_( updateStep_ )
      {}

      virtual ~TimeProvider()
      {}

      TimeProvider ( const ThisType & ) = delete;

      ThisType &operator= ( const ThisType & ) = delete;

      /** \brief init dt with time step estimate
       */
      inline void init ()
      {
        initTimeStep( dtEstimate_ );
      }

      /** \brief init dt with provided time step
       *
       *  \param[in]  timeStep  value of the first time step (is multiplied with
       *                        factor)
       */
      inline void init ( const double timeStep )
      {
        initTimeStep( timeStep );
      }

      /** \brief goto next time step
       *
       * Sets the size of the next time step to the current time step estimate
       * and sets the estimate to infinity.
       */
      inline void next ()
      {
        assert( this->dtEstimateValid_ );
        advance();
        initTimeStep( dtEstimate_ );
      }

      /** \brief goto next time step
       *
       * Sets the size of the next time step to the provided time step value
       * and sets the estimate to infinity.
       *
       *  \param[in]  timeStep  value of the next time step (is multiplied with
       *                        factor)
       */
      inline void next ( const double timeStep )
      {
        advance();
        initTimeStep(timeStep);
      }

      /** \brief  return the global factor number
          \return time step factor
      */
      inline double factor () const
      {
        return cfl_;
      }

    protected:
      using BaseType::advance;
      using BaseType::initTimeStepEstimate;

      void initTimeStep ( const double dtEstimate )
      {
        // increase counter
        ++counter_ ;

        if( counter_ >= updateStep_ )
        {
          // set timestep estimate
          dt_ = std::min(cfl_ * dtEstimate,dtUpperBound_);
          dt_ = comm_.min( dt_ );
          invdt_ = 1.0 / dt_;
          valid_ = (dt_ > 0.0);
          // reset counter
          counter_ = 0;
        }

        initTimeStepEstimate();
      }

    public:
      /** \brief restore time and timestep from outside
           (i.e. from former calculation)
           \param[in] time new time
           \param[in] timeStep new time step counter
      */
      inline void restore ( const double time, const int timeStep )
      {
        time_ = time;
        timeStep_ = timeStep;
      }

      inline virtual void backup () const
      {
        BaseType::backup();
      }

      inline virtual void restore ()
      {
        BaseType::restore();
        const_cast< double & >( cfl_ ) = getCflFactor();
      }

    protected:
      using BaseType::dt_;
      using BaseType::invdt_;
      using BaseType::dtEstimate_;
      using BaseType::dtUpperBound_;
      using BaseType::valid_;
      using BaseType::timeStep_;

      const CommunicationType& comm_;
      const double cfl_;
      const int updateStep_;
      int counter_;
    };



    /** \class   GridTimeProvider
     *  \ingroup ODESolver
     *  \brief   the same functionality as the Dune::TimeProvider.
     *
     *  This implementation of a timeprovider takes the Communication
     *  from a Dune::Grid instance.
     */
    template< class Grid >
    class GridTimeProvider
    : public TimeProvider< typename Grid::Traits::Communication >
    {
      typedef GridTimeProvider< Grid > ThisType;
      typedef TimeProvider< typename Grid::Traits::Communication > BaseType;

      // type of DofManager for sequence number
      typedef DofManager < Grid > DofManagerType ;

    public:
      typedef typename Grid::Traits::Communication CommunicationType;

      explicit GridTimeProvider ( const Grid &grid )
      : BaseType( grid.comm() ),
        dm_( DofManagerType ::instance( grid ) ),
        sequence_( -1 )
      {}

      GridTimeProvider ( const double startTime,
                         const Grid &grid )
      : BaseType( startTime, grid.comm() ),
        dm_( DofManagerType ::instance( grid ) ),
        sequence_( -1 )
      {}

      GridTimeProvider ( const double startTime,
                         const double cfl,
                         const Grid &grid )
      : BaseType( startTime, cfl, grid.comm() ),
        dm_( DofManagerType ::instance( grid ) ),
        sequence_( -1 )
      {}

      virtual ~GridTimeProvider() {}

    protected:
      using BaseType :: counter_ ;
      using BaseType :: updateStep_ ;

      // this initTimeStep method also check the sequence number
      // in case the grid has changed due to adaptivity
      void initTimeStep ( const double dtEstimate )
      {
        const int currentSequence = dm_.sequence();
        // check sequence number
        if( sequence_ != currentSequence )
        {
          // if sequence number changed, update in any case
          counter_  = updateStep_ ;
          sequence_ = currentSequence ;
        }

        // call initTimeStep on base class
        BaseType :: initTimeStep( dtEstimate );
      }

      const DofManagerType& dm_;
      int sequence_ ;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_TIMEPROVIDER_HH
