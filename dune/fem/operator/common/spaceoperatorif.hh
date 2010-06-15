#ifndef DUNE_SPACEOPERATORIF_HH
#define DUNE_SPACEOPERATORIF_HH

//- system includes 
#include <limits>
#include <utility>

//-Dune fem includes 
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/objpointer.hh>

namespace Dune {

/** @ingroup OperatorCommon
  \brief ODESpaceOperatorInterface for Operators that work with PARDG ODE solvers 
  of the type \f$L: X \longrightarrow X\f$ where \f$X\f$ is a discrete function space.
  
  \interfaceclass
*/
template <class DestinationImp> 
class PARDGSpaceOperatorInterface 
{
protected:
  // only allow derived class to call this constructor  
  PARDGSpaceOperatorInterface() {}

public:
  //! type of argument and destination 
  typedef DestinationImp DestinationType;
  
  //! destructor 
  virtual ~PARDGSpaceOperatorInterface() {}
  
  /** \brief return size of discrete function space, i.e. number of unknowns */
  virtual int size() const = 0 ;

  /** \brief call operator once to calculate initial time step size 
      \param U0  initial data to compute initial time step size 
   */
  virtual void initializeTimeStepSize( const DestinationType& U0 ) const = 0;

  /** \brief application operator to apply right hand side 
      \param u  argument, u 
      \param f  destination, f(u)
   */
  virtual void operator() (const double* u, double *f ) const = 0;

  /** \brief set time for operators 
      \param time current time of evaluation 
  */
  virtual void setTime(const double time) {}

  /** \brief returns maximal possible dt to assure stabil 
       explicit Runge Kutta ODE Solver. */
  virtual double timeStepEstimate () const 
  {
    return std::numeric_limits<double>::max();  
  }
};

/** @ingroup OperatorCommon
  \brief SpaceOperatorInterface for Operators of the type 
  \f$L: X \longrightarrow X\f$ where \f$X\f$ is a discrete function space.
  Use this interface to implement Operators working on DiscreteFunctions for ODE
  solvers.
  
  \interfaceclass
*/
template <class DestinationImp> 
class SpaceOperatorInterface 
: public Fem::Operator< DestinationImp, DestinationImp >,
  public PARDGSpaceOperatorInterface< DestinationImp >
{
  typedef Fem::Operator< DestinationImp, DestinationImp > BaseType;
  using BaseType :: operator ();

protected:
  // only allow derived class to call this constructor  
  SpaceOperatorInterface() {}

public:
  //! type of argument and destination 
  typedef DestinationImp DestinationType;
  
  //! type of discrete function space 
  typedef typename DestinationType :: DiscreteFunctionSpaceType   SpaceType;
  
  //! destructor 
  virtual ~SpaceOperatorInterface() {}

  //! return reference to space (needed by ode solvers)
  virtual const SpaceType& space() const = 0;

  /** \copydoc Dune::PARDGSpaceOperatorInterface::size() */
  virtual int size() const { return space().size(); }

  /** \copydoc Dune::PARDGSpaceOperatorInterface::operator(const double* u, double *f) */
  virtual void operator() (const double* u, double *f ) const 
  {
    // get space instance 
    const SpaceType& spc = space();

    // convert argument to discrete function 
    const DestinationType arg ("SpaceOperatorIF::ARG" , spc, u);

    // convert argument to discrete function 
    DestinationType dest("SpaceOperatorIF::DEST", spc, f);

    // call operator apply 
    this->operator ()( arg, dest );
  }

  /** \copydoc Dune::PARDGSpaceOperatorInterface::initializeTimeStepSize( const DestinationType& U0 ) */
  virtual void initializeTimeStepSize( const DestinationType& U0 ) const 
  {
    // create temporary variable 
    DestinationType tmp( U0 );
    // call operator 
    this->operator() ( U0, tmp );
  }

  //! return reference to pass's local memory  
  virtual const DestinationType* destination() const { return 0; }
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

  // use pass method of base 
  using BaseType :: pass;
  
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

} // end namespace Dune 
#endif
