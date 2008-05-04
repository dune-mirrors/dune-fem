#define DUNE_IDENTITY_HH
#ifndef DUNE_IDENTITY_HH

// Dune includes
#include <dune/fem/operator/common/operator.hh>

// Local includes
#include "inverseoperatorfactory.hh"


namespace Dune {
  // ! The identity operator is a class, that maps an argument on itself
  template <class DiscreteFunctionType>
  class Identity : 
    public Operator<typename DiscreteFunctionType::DomainFieldType,
                    typename DiscreteFunctionType::RangeFieldType,
                    DiscreteFunctionType,
                    DiscreteFunctionType> {
  public:
    //! Applying here just means copying
    virtual void operator() (const DiscreteFunctionType& arg,
                             DiscreteFunctionType& dest) const {
      dest.clear();
      dest.addScaled(arg, 1.0);
      // * Does this really a copy?
      //dest.assign(arg);
      //dest = arg;
    }
  };
  
  //! Solver inverting an identity (that was an easy one!)
  template <typename DiscreteFunctionType>
  class IdentitySolver :
    public Operator<typename DiscreteFunctionType::DomainFieldType,
                    typename DiscreteFunctionType::RangeFieldType,
                    DiscreteFunctionType,
                    DiscreteFunctionType> {
  public:
    typedef Mapping<typename DiscreteFunctionType::DomainFieldType,
                    typename DiscreteFunctionType::RangeFieldType,
                    DiscreteFunctionType,
                    DiscreteFunctionType> MappingType;

    //! Constructor
    IdentitySolver(const MappingType& dummy) : op_() {}

    //! Inverting an identity is the same as applying it
    void operator() (const DiscreteFunctionType& arg,
                     DiscreteFunctionType& dest) const {
      op_(arg, dest);
    }

  private:
    Identity<DiscreteFunctionType> op_;
  };

  //! Solver factory for an identity solver
  template <typename DiscreteFunctionType>
  class IdentitySolverFactory : 
    public InverseOperatorFactory<DiscreteFunctionType> {
  public:
    typedef typename InverseOperatorFactory<DiscreteFunctionType>::MappingType MappingType;
    
    IdentitySolverFactory() {}

    virtual MappingType* createOperator(const MappingType& op) {
      return new IdentitySolver<DiscreteFunctionType>(op);
    }
  private:
  };


} // end namespace Dune

#endif
