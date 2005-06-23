#ifndef __COMPOSITEDISCRETEFUNCTIONSPACE_HH__
#define __COMPOSITEDISCRETEFUNCTIONSPACE_HH__

#include <dune/grid/common/grid.hh>
#include <dune/fem/common/discretefunctionspace.hh>
#include <dune/fem/dofmanager.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/misc.hh>

using namespace Dune;

namespace Adi {
  //! \brief An assembly of scalar discrete function spaces for vectorial
  //! problems
  //! Wraps a single discrete function space and acts as if it was a vectorial
  //! one
  template <class ContainedFunctionSpaceImp, int dimRange>
  class CompositeDiscreteFunctionSpace :
    public DiscreteFunctionSpaceDefault<
    FunctionSpace<
    typename ContainedFunctionSpaceImp::DomainField,
    typename ContainedFunctionSpaceImp::RangeField,
    ContainedFunctionSpaceImp::DimDomain,
    dimRange
    >,
    typename ContainedFunctionSpaceImp::GridType,
    CompositeDiscreteFunctionSpace<ContainedFunctionSpaceImp, dimRange>,
    typename ContainedFunctionSpaceImp::BaseFunctionSetType
    >
  {
    typedef DiscreteFunctionSpaceDefault<
      FunctionSpace<
      typename ContainedFunctionSpaceImp::DomainField,
      typename ContainedFunctionSpaceImp::RangeField,
      ContainedFunctionSpaceImp::DimDomain,
      dimRange
      >,
      typename ContainedFunctionSpaceImp::GridType,
      CompositeDiscreteFunctionSpace<ContainedFunctionSpaceImp, dimRange>,
      typename ContainedFunctionSpaceImp::BaseFunctionSetType
      > BaseType;
  public:
    //- Typedefs
    //! Remember own type
    typedef CompositeDiscreteFunctionSpace<ContainedFunctionSpaceImp,
                                           dimRange> ThisType;

    //! The type of the contained discrete function space
    typedef ContainedFunctionSpaceImp ContainedFunctionSpaceType;
    //! The base function set type
    //! The same base function set as in the wrapped function space is used
    typedef typename ContainedFunctionSpaceImp::BaseFunctionSetType BaseFunctionSetType;

    //! The continuous function space
    typedef FunctionSpace<
      typename ContainedFunctionSpaceImp::DomainField,
      typename ContainedFunctionSpaceImp::RangeField,
      ContainedFunctionSpaceImp::DimDomain,
      dimRange
      > FunctionSpaceType;
    
    //! The domain type
    typedef typename FunctionSpaceType::Domain Domain;

    //! The range type
    typedef typename FunctionSpaceType::Range Range;

    //! The boundary manager
    typedef BoundaryManager<FunctionSpaceType> BoundaryManagerType;

    //- Enums and constants
    //! The modified dimension of the range vector space
    enum { DimRange = dimRange };
    
    //! This here IS the beast
    static const int id_ = 666;

    //- Methods
    //! Constructor
    CompositeDiscreteFunctionSpace(ContainedFunctionSpaceType& dfs,
                                   BoundaryManagerType bc) :
      BaseType(dfs.getGrid(), bc, id_, dfs.level()), 
      funcSpace_(dfs) 
    {
      /*
       CompileTimeChecker<
        SameType<
        typename ContainedFunctionSpaceImp::FunctionSpaceType,
        FunctionSpace<
        typename ContainedFunctionSpaceImp::FunctionSpaceType::DomainField,
        typename ContainedFunctionSpaceImp::FunctionSpaceType::RangeField,
        ContainedFunctionSpaceImp::FunctionSpaceType::DimDomain,
        1>
        >::value
        > only_scalar_vector_spaces_accepted;
      */
    }

    //! Destructor
    virtual ~CompositeDiscreteFunctionSpace() {}

    //! \brief Are base functions continuous?
    //! Forwards call to contained space
    bool continuous() const { return funcSpace_.continuous(); }

    //! Number of unknowns for this function space per range field element
    int size() const { return funcSpace_.size(); }

    //! Maximal polynom order
    int polynomOrder() const { return funcSpace_.polynomOrder(); }
    
    //! Local polynom order
    template <class EntityType>
    int localPolynomOrder(EntityType& en) const { 
      return funcSpace_.localPolynomOrder(en);
    }
    
    //! Local to global mapping
    template <class EntityType>
    int mapToGlobal(EntityType& en, int localNum) const {
      return funcSpace_.mapToGlobal(en, localNum);
    }
    
    //! Get base function set
    template <class EntityType>
    const BaseFunctionSetType& getBaseFunctionSet(EntityType& en) const {
      return funcSpace_.getBaseFunctionSet(en);
    }

    //! Get contained function space
    ContainedFunctionSpaceType& containedFunctionSpace() const {
      return funcSpace_;
    }

  private:
    //- Prohibited defaults
    //! Copy constructor
    CompositeDiscreteFunctionSpace(const ThisType&);
    
    //! Assignment opertor
    ThisType& operator=(const ThisType&);

    //- Member variables
    //! The wrapped discrete function space
    ContainedFunctionSpaceType& funcSpace_;
  };

} // end namespace Adi

#endif
