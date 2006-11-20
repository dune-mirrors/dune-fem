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
    typedef typename GridType::template Codim<0>::Entity EntityType;
    typedef typename GridType::Traits::IntersectionIterator IntersectionIterator;
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
    InsertFunctionPass(DestinationType & dest, 
              PreviousPassType& pass) 
      : BaseType(pass,dest.space())
      , realDestination_(dest)
    {
      //std::cout << "create dummy pass " <<  realDestination_.name() <<  "\n";
    }

    //! Destructor
    virtual ~InsertFunctionPass() 
    {
      this->destination_ = 0;
    }

    //! Allocates the local memory of a pass, if needed.
    virtual void allocateLocalMemory() 
    {
      if(!this->destination_)
      {
        //std::cout << "Destination points to " << realDestination_.name() << "\n";
        this->destination_ = &realDestination_;
      }
      // make sure that destination_ point to inserted function 
      assert( this->destination_ == &realDestination_ );
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

  private:
    //! function to be stored from outside 
    DestinationType & realDestination_;
  }; // end class InsertFunctionPass


} // end namespace Dune
#endif
