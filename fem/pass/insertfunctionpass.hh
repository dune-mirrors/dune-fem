#ifndef DUNE_INSERTFUNCTIONPASS_HH
#define DUNE_INSERTFUNCTIONPASS_HH

//- System includes
#include <string>

//- local includes 
#include <dune/fem/misc/femtuples.hh>
#include <dune/fem/function/common/discretefunction.hh>

#include "discretemodel.hh"
#include "pass.hh"

namespace Dune {

  //! Traits for InsertFunctionPass to create dummy discrete model
  template <class DiscreteFunctionImp>
  struct EmptyDiscreteModelTraits
  {
    typedef DiscreteFunctionImp DiscreteFunctionType;
    typedef EmptyDiscreteModelTraits<DiscreteFunctionType> ThisType;  
    typedef DiscreteModelDefault<ThisType> DiscreteModelType; 

    typedef DiscreteFunctionImp DestinationType;
    typedef typename DiscreteFunctionImp::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
    typedef typename DiscreteFunctionSpaceType::GridType GridType;
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
    typedef typename GridType::template Codim<0>::Entity EntityType;
    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
  };
  
  /**
   * @brief Base class for specific pass implementations.
     InsertFunctionPass simply inserts a discrete function from outside of the pass tree 
     into the current pass tree, for example when calculating the species 
     transport the velocity function comes from a different pass but has to
     be inserted into the species pass. 
   */
  template< class DiscreteFunctionImp , class PreviousPassImp , int passIdImp  = -1 >
  class InsertFunctionPass 
    : public Pass< DiscreteModelDefault< EmptyDiscreteModelTraits< DiscreteFunctionImp > >
                   , PreviousPassImp , passIdImp >
  {
    //! make sure DiscreteFunctionImp provides local functions 
    enum { hasLocalFunction = Conversion<DiscreteFunctionImp,HasLocalFunction>:: exists };
    CompileTimeChecker <hasLocalFunction> discrete_function_does_not_provide_local_functions;
    
  public:
    //- Typedefs and enums
    //! type of traits for this class  
    typedef EmptyDiscreteModelTraits<DiscreteFunctionImp> Traits;
    //! type of discrete model for this class 
    typedef DiscreteModelDefault<Traits> DiscreteModelImp;
    //! Base class
    typedef Pass<DiscreteModelImp, PreviousPassImp , passIdImp > BaseType;

    //! type of this class 
    typedef InsertFunctionPass<DiscreteFunctionImp,PreviousPassImp , passIdImp > ThisType;

    //! Repetition of template arguments
    typedef DiscreteModelImp DiscreteModelType;
    //! Repetition of template arguments
    typedef PreviousPassImp PreviousPassType;

    // Types from the base class
    typedef typename BaseType::TotalArgumentType ArgumentType;
    typedef typename BaseType::GlobalArgumentType GlobalArgumentType;

    //! The discrete function representing the return value of this pass
    typedef typename Traits::DestinationType DestinationType;
    //! The discrete function space belonging to DestinationType
    typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

  public:
    //- Public methods
    //! Constructor
    //! \param destination to be stored in this pass 
    //! \param pass Previous pass
    InsertFunctionPass(DestinationType & destination, 
                       PreviousPassType& pass)
      : BaseType(pass)
    {
      this->destination_ = &destination;
    }

    //- Public methods
    //! Constructor
    //! \param destination to be stored in this pass 
    //! \param pass Previous pass
    //! \param space to make constructor look like that from other passes
    InsertFunctionPass(DestinationType & destination, 
                       PreviousPassType& pass,
                       const DiscreteFunctionSpaceType& space) 
      : BaseType(pass)
    {
      this->destination_ = &destination;
    }

    //- Public methods
    //! Constructor
    //! \param pass Previous pass
    InsertFunctionPass(PreviousPassType& pass) 
      : BaseType(pass)
    {
      assert( this->destination_ == 0 );
    }

    //! Destructor
    virtual ~InsertFunctionPass() 
    {
      this->destination_ = 0;
    }

    //! Allocates the local memory of a pass, if needed.
    virtual void allocateLocalMemory() 
    {
      // do not allocate memory here 
    }

    //! return reference to space 
    const DiscreteFunctionSpaceType& space () const 
    {
      assert( this->destination_ );
      return this->destination_->space();
    }

    //! empty method here
    void operator () (const GlobalArgumentType& arg, DestinationType& dest) const
    {
    }

    //! empty method here
    void prepare(const ArgumentType& arg, DestinationType& dest) const
    {
    }

    //! empty method here
    void finalize(const ArgumentType& arg, DestinationType& dest) const
    {
    }

    //! set internal destination pointer to dest 
    void setDestination(DestinationType& dest) 
    {
      this->destination_ = &dest;
    }
    
  protected:
    void compute(const ArgumentType& arg, DestinationType& dest) const 
    {
      // do nothing here 
    }
  }; // end class InsertFunctionPass


} // end namespace Dune
#endif
