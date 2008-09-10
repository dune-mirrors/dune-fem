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
  \brief SpaceOperatorInterface for Operators of the type 
  \f$L: X \longrightarrow X\f$ where \f$X\f$ is a discrete function space.
  
  \interfaceclass
*/
template <class DestinationImp>
class SpaceOperatorInterface 
: public Operator< typename DestinationImp :: RangeFieldType,
                   typename DestinationImp :: RangeFieldType,
                   DestinationImp,
                   DestinationImp>
{
protected:
  // only allow derived class to call this constructor  
  SpaceOperatorInterface() {}

public:
  //! type of argument and destination 
  typedef DestinationImp DestinationType;
  
  //! type of discrete function space 
  typedef typename DestinationType :: DiscreteFunctionSpaceType SpaceType;
  
  //! destructor 
  virtual ~SpaceOperatorInterface() {}

  //! apply operator 
  virtual void operator () (const DestinationType& arg, DestinationType& dest) const = 0;

  //! return reference to space (needed by ode solvers)
  virtual const SpaceType& space() const = 0;

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

  //! return reference to pass's local memory  
  virtual const DestinationType* destination() const { return 0; }

  template <class TimeProviderImp> 
  void DUNE_DEPRECATED timeProvider(TimeProviderImp* tp)
  {
    assert( tp );
    // deprecated method 
    this->setTime( tp->time() );
  }
};

//! only for kepping the pointer 
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

//! only for kepping the pointer 
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
