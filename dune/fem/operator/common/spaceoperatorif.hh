#ifndef DUNE_FEM_SPACEOPERATORIF_HH
#define DUNE_FEM_SPACEOPERATORIF_HH

//- system includes
#include <cassert>
#include <cstdlib>
#include <limits>
#include <utility>

//-Dune fem includes
#include <dune/fem/operator/common/automaticdifferenceoperator.hh>
#include <dune/fem/operator/common/objpointer.hh>
#include <dune/fem/function/adaptivefunction.hh>

namespace Dune
{

  namespace Fem
  {
    /** \class   SpaceOperatorInterface
      * \ingroup OperatorCommon
      *  \brief  interface for time evolution operators
      *
      * The SpaceOperatorInterface defines an interface for operators
      * \f$L: X \longrightarrow X\f$ from a discrete function space \f$X\f$ into
      * itself.
      * This interface is used to implement operators working with the ODE solvers.
      *
      * \tparam  DiscreteFunction  type of discretefunction modelling the elements
      *                            of \f$X\f$.
      *
      * \interfaceclass
      */
    template< class DiscreteFunction,
              class JacobianOperator = Fem::AutomaticDifferenceOperator< DiscreteFunction >  >
    class SpaceOperatorInterface
    : public JacobianOperator
    {
      typedef SpaceOperatorInterface< DiscreteFunction, JacobianOperator > ThisType;
      typedef Fem::Operator< DiscreteFunction > BaseType;

    public:
      //! type of argument and destination
      typedef DiscreteFunction DestinationType;

      //! type of discrete function space
      typedef typename DestinationType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      //! convenience typedef for space type
      typedef DiscreteFunctionSpaceType SpaceType;

      //! destructor
      virtual ~SpaceOperatorInterface() {}

      using BaseType::operator ();

      //! return reference to space (needed by ode solvers)
      virtual const DiscreteFunctionSpaceType& space() const = 0;

      /** \brief return size of discrete function space, i.e. number of unknowns */
      virtual int size () const { return space().size(); }

      /** \brief call operator once to calculate initial time step size
          \param U0  initial data to compute initial time step size
       */
      virtual void initializeTimeStepSize ( const DestinationType& U0 ) const;

      /** \brief return true if limit method is implemented

          \return true if limit is implemented
       */
      virtual bool hasLimiter () const { return false ; }

      /** \brief set time for operators
          \param time current time of evaluation
      */
      virtual void setTime ( const double time ) {}

      /** \brief estimate maximum time step
       *
       *  For an explicit time discretization, the time step has to be limited.
       *  An estimate for the maximum time step of an explicit Euler scheme is
       *  returned by this function.
       *  Maximum time steps for higher order Runge Kutta schemes can be derived
       *  from this value.
       *  */
      virtual double timeStepEstimate () const
      {
        return std::numeric_limits< double >::max();
      }

      /** \brief limiter application operator
          \param arg   argument, u
          \param dest  destination, Limiter(u)

          \note: Default implementation is to copy arg into dest.
       */
      virtual void limit (const DestinationType& arg, DestinationType& dest) const
      {
        // if this method is not overloaded then hasLimiter should return false
        assert( ! hasLimiter () );
        // default operation is copy arg to dest
        dest.assign( arg );
      }

      /** \brief limiter application operator
          \param[inout] U argument and destination to apply Limiter(u), needs internal copying

          \note: Default implementation is to do nothing (hasLimiter == false)
       */
      virtual void applyLimiter ( DestinationType& U ) const
      {
        if( hasLimiter() )
        {
          if( ! uTmp_ )
            uTmp_.reset( new DestinationType( "SpaceOpIF::uTmp_", space() ) );

          assert( uTmp_ );
          uTmp_->assign( U );
          limit( *uTmp_, U );
        }
      }

      //! return reference to pass's local memory
      virtual const DestinationType* destination() const { return nullptr; }

    protected:
      mutable std::unique_ptr< DestinationType > uTmp_;
    };

    //! only for keeping the pointer
    template <class OperatorType>
    class SpaceOperatorStorage
    : public ObjPointerStorage
    {
      //! copying not allowed
      SpaceOperatorStorage(const SpaceOperatorStorage& org);
      SpaceOperatorStorage& operator = (const SpaceOperatorStorage& org);

    protected:
      // operator storage
      mutable OperatorType* op_;
      // model storage
      ObjPointerStorage* model_;

    public:
      //! constructor storing pointer
      SpaceOperatorStorage(OperatorType * op)
        : op_(op), model_(0)
      {}

      //! constructor storing pointer
      SpaceOperatorStorage(OperatorType * op, ObjPointerStorage* model)
        : op_(op), model_(model)
      {}

      //! destructor deletes operator
      ~SpaceOperatorStorage()
      {
        // delete operator before destructor of base class is called
        delete op_; op_ = 0;
        delete model_; model_ = 0;
      }

      //! return reference to pass
      OperatorType& pass() const
      {
        assert( op_ );
        return (*op_);
      }
    };

    //! only for keeping the pointer
    template <class OperatorType>
    class SpaceOperatorPtr
    : public SpaceOperatorStorage< OperatorType >,
      public SpaceOperatorInterface<typename OperatorType::DestinationType>
    {
      //! type of base class
      typedef SpaceOperatorStorage< OperatorType > BaseType;

			protected:
      // use pass method of base
      using BaseType :: pass;
			private:

      //! type of destination
      typedef typename OperatorType::DestinationType DestinationType;

      //! type of discrete function space
      typedef typename DestinationType :: DiscreteFunctionSpaceType SpaceType;

      //! copying not allowed
      SpaceOperatorPtr(const SpaceOperatorPtr& org);
      SpaceOperatorPtr& operator = (const SpaceOperatorPtr& org);

    public:
      //! constructor storing pointer
      SpaceOperatorPtr(OperatorType * op)
        : BaseType(op)
      {}

      //! constructor storing pointer
      SpaceOperatorPtr(OperatorType * op, ObjPointerStorage* model)
        : BaseType(op,model)
      {}

      //! destructor
      virtual ~SpaceOperatorPtr() {}

      //! application operator does nothing here
      virtual void operator () (const DestinationType& arg, DestinationType& dest) const
      {
        // this method should not be called
        assert(false);
        abort();
      }

      //! return reference to space
      const SpaceType& space() const { return pass().space(); }

      /** @copydoc SpaceOperatorInterface::setTime  */
      void setTime(const double time) { pass().setTime(time); }

      /** @copydoc SpaceOperatorInterface::timeStepEstimate */
      double timeStepEstimate () const { return pass().timeStepEstimate(); }

      //! return reference to pass's local memory
      const DestinationType* destination() const
      {
        pass().allocateLocalMemory();
        return & (pass().destination());
      }
    };

    //! apply wrapper
    template <class OperatorType>
    class SpaceOperatorWrapper
    : public SpaceOperatorPtr< OperatorType >
    {
      //! type of base class
      typedef SpaceOperatorPtr< OperatorType > BaseType;

      // use pass method of base
      using BaseType :: pass;

      //! copying not allowed
      SpaceOperatorWrapper(const SpaceOperatorWrapper& org);
      SpaceOperatorWrapper& operator = (const SpaceOperatorWrapper& org);
    public:
      //! type of Argument and Destination
      typedef typename OperatorType::DestinationType DestinationType;
      //! type of discrete function space
      typedef typename DestinationType :: DiscreteFunctionSpaceType SpaceType;

      //! constructor storing pointer
      SpaceOperatorWrapper(OperatorType * op)
        : BaseType(op)
      {}

      //! constructor storing pointer
      SpaceOperatorWrapper(OperatorType * op, ObjPointerStorage* model)
        : BaseType(op,model)
      {}

      //! call application operator of internal operator
      void operator () (const DestinationType& arg, DestinationType& dest) const
      {
        pass()(arg,dest);
      }
    };


    // Implementation of SpaceOperatorInterface
    // ----------------------------------------

    template< class DiscreteFunction, class JacobianOperator >
    inline void SpaceOperatorInterface< DiscreteFunction, JacobianOperator >
      ::initializeTimeStepSize ( const DestinationType &U0 ) const
    {
      // create temporary variable
      DestinationType tmp( U0 );
      // call operator
      (*this)( U0, tmp );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACEOPERATORIF_HH
