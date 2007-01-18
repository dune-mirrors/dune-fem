#ifndef DUNE_INSERTFUNCTIONPASS_HH
#define DUNE_INSERTFUNCTIONPASS_HH

//- System includes
#include <string>

//- local includes 
#include <dune/common/tuples.hh>
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
  template <class DiscreteFunctionImp, class PreviousPassImp>
  class InsertFunctionPass 
    : public LocalPass<DiscreteModelDefault<EmptyDiscreteModelTraits<DiscreteFunctionImp> > , 
             PreviousPassImp>
  {
  public:
    //- Typedefs and enums
    //! type of traits for this class  
    typedef EmptyDiscreteModelTraits<DiscreteFunctionImp> Traits;
    //! type of discrete model for this class 
    typedef DiscreteModelDefault<Traits> DiscreteModelImp;
    //! Base class
    typedef LocalPass<DiscreteModelImp, PreviousPassImp> BaseType;

    //! type of this class 
    typedef InsertFunctionPass<DiscreteFunctionImp,PreviousPassImp> ThisType;

    //! Repetition of template arguments
    typedef DiscreteModelImp DiscreteModelType;
    //! Repetition of template arguments
    typedef PreviousPassImp PreviousPassType;

    // Types from the base class
    typedef typename BaseType::Entity EntityType; 
    typedef typename BaseType::ArgumentType ArgumentType;
    typedef typename BaseType::GlobalArgumentType GlobalArgumentType;

    // Types from the traits
    typedef typename Traits::DestinationType DestinationType;
    typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  public:
    //- Public methods
    //! Constructor
    //! \param destination to be stored in this pass 
    //! \param pass Previous pass
    InsertFunctionPass(DestinationType & destination, 
                       PreviousPassType& pass) 
      : BaseType(pass,destination.space())
    {
      this->destination_ = &destination;
    }

    //- Public methods
    //! Constructor
    //! \param space to be stored in this base pass 
    //! \param pass Previous pass
    InsertFunctionPass(const DiscreteFunctionSpaceType& space, PreviousPassType& pass) 
      : BaseType(pass,space)
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

    //! empty method here
    void applyLocal(EntityType& en) const
    {
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
  }; // end class InsertFunctionPass


} // end namespace Dune
#endif
