#ifndef DUNE_COMBINEDSPACE_HH
#define DUNE_COMBINEDSPACE_HH

//- System includes
#include <vector>

//- Dune includes
#include <dune/common/misc.hh>

//- Local includes
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/dofmapperinterface.hh>
#include <dune/fem/basefunctions/common/basefunctions.hh>
#include "subspace.hh"

#include "combineddofstorage.hh"

namespace Dune {
  
  // Forward declarations
  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  class CombinedSpace;
  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  class CombinedMapper;
  template <class BaseFunctionSetImp, int N, DofStoragePolicy policy>
  class CombinedBaseFunctionSet;

  //! Traits class for CombinedSpace
  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy = PointBased>
  struct CombinedSpaceTraits {
  private:
    typedef DiscreteFunctionSpaceImp ContainedDiscreteFunctionSpaceType;

    typedef typename ContainedDiscreteFunctionSpaceType::Traits 
    ContainedSpaceTraits;
    typedef typename ContainedSpaceTraits::FunctionSpaceType 
    ContainedFunctionSpaceType;
    typedef typename ContainedSpaceTraits::BaseFunctionSetType 
    ContainedBaseFunctionSetType;
    
    enum { ContainedDimRange = ContainedFunctionSpaceType::DimRange,
           ContainedDimDomain = ContainedFunctionSpaceType::DimDomain };
  public:
    typedef typename ContainedSpaceTraits::MapperType ContainedMapperType;

    typedef typename ContainedFunctionSpaceType::DomainFieldType 
    DomainFieldType;
    typedef typename ContainedFunctionSpaceType::RangeFieldType 
    RangeFieldType;
    typedef typename ContainedFunctionSpaceType::RangeType 
    ContainedRangeType;
    typedef typename ContainedFunctionSpaceType::JacobianRangeType
    ContainedJacobianRangeType;

   typedef CombinedSpace<
      DiscreteFunctionSpaceImp, N, policy> DiscreteFunctionSpaceType;
    typedef FunctionSpace<
      DomainFieldType, RangeFieldType, 
      ContainedDimDomain, ContainedDimRange*N > FunctionSpaceType;
    typedef CombinedBaseFunctionSet<
      DiscreteFunctionSpaceImp, N, policy> BaseFunctionSetType;
    typedef CombinedMapper<DiscreteFunctionSpaceImp, N, policy> MapperType;
   
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

    typedef typename ContainedSpaceTraits::GridType GridType;
    typedef typename ContainedSpaceTraits::GridPartType GridPartType;
    typedef typename ContainedSpaceTraits::IndexSetType IndexSetType;
    typedef typename ContainedSpaceTraits::IteratorType IteratorType;

    typedef CombinedDofConversionUtility<DiscreteFunctionSpaceImp,policy> DofConversionType;

    enum { DimRange = FunctionSpaceType::DimRange,
           DimDomain = FunctionSpaceType::DimDomain };
  public:
    //- Friends
    friend class CombinedSpace<DiscreteFunctionSpaceImp, N, policy>;
    friend class CombinedBaseFunctionSet<DiscreteFunctionSpaceImp, N, policy>;
    friend class CombinedMapper<DiscreteFunctionSpaceImp, N, policy>;
  };

  //! Class to combine N scalar spaces
  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy = PointBased>
  class CombinedSpace :
    public DiscreteFunctionSpaceDefault<
    CombinedSpaceTraits<DiscreteFunctionSpaceImp, N, policy> 
  > 
  {
  private:
    //- Private typedefs
    typedef DiscreteFunctionSpaceDefault<
    CombinedSpaceTraits<DiscreteFunctionSpaceImp, N, policy> 
    > BaseType;
  public:
    // polynomial Order is the same as for the single space 
    enum { polynomialOrder = DiscreteFunctionSpaceImp :: polynomialOrder };
    
    //- Public typedefs and enums
    typedef CombinedSpace<DiscreteFunctionSpaceImp, N, policy> ThisType;
    typedef CombinedSpaceTraits<DiscreteFunctionSpaceImp, N, policy> Traits;
    typedef DiscreteFunctionSpaceImp ContainedDiscreteFunctionSpaceType;
    
    typedef typename Traits::IteratorType IteratorType;

    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeFieldType RangeFieldType;
    typedef typename Traits::DomainFieldType DomainFieldType;

    typedef typename Traits::ContainedRangeType ContainedRangeType;
    typedef typename Traits::ContainedJacobianRangeType ContainedJacobianRangeType;

    typedef typename Traits::BaseFunctionSetType BaseFunctionSetType;
    typedef typename Traits::MapperType MapperType;
    typedef typename Traits::GridType GridType;
    typedef typename Traits::GridPartType GridPartType;
    typedef typename Traits::IndexSetType IndexSetType;

    typedef typename Traits::DofConversionType DofConversionType;

    CompileTimeChecker<(Traits::ContainedDimRange == 1)> use_CombinedSpace_only_with_scalar_spaces;
  public:
    //- Public methods
    //! constructor
    CombinedSpace(ContainedDiscreteFunctionSpaceType& spc);

    //! destructor
    ~CombinedSpace();

    //! type
    int type() const { return spaceId_; }

    //! continuous?
    bool continuous() const { return spc_.continuous(); }

    //! polynom order
    int polynomOrder() const { return spc_.polynomOrder(); }

    //! begin iterator
    IteratorType begin() const { return spc_.begin(); }

    //! end iterator
    IteratorType end() const { return spc_.end(); }

    //! total number of dofs
    int size() const { return mapper_.size(); }

    //! map a local dof number to a global one
    template <class EntityType>
    int mapToGlobal(EntityType& en, int local) const {
      return mapper_.mapToGlobal(en, local);
    }

    //! access to base function set
    template <class EntityType>
    const BaseFunctionSetType& getBaseFunctionSet(const EntityType& en) const 
    {
      GeometryIdentifier::IdentifierType id =
                GeometryIdentifier::fromGeometry(en.geometry()); 

      assert(id < (int) baseSetVec_.size());
      assert(id >= 0);

      assert( baseSetVec_[id] );
      return *baseSetVec_[id];
    }

    //! access to grid
    const GridType& grid() const { return spc_.grid(); }
    
    //! access to gridPart
    GridPartType& gridPart() { return spc_.gridPart(); }
    //! access to gridPart
    const GridPartType& gridPart() const { return spc_.gridPart(); }

    const IndexSetType& indexSet() const { return spc_.indexSet(); }

    //! access to mapper
    const MapperType& mapper() const { return mapper_; }

    //- Additional methods
    //! number of components
    int numComponents() const { return N; }

    //! policy of this space
    DofStoragePolicy myPolicy() const{ return DofConversionType::policy(); }
  private:
    //- Private typedefs
    typedef typename Traits::ContainedMapperType ContainedMapperType;
   
    //- Friend
    friend class SubSpace<ThisType>;
    
  private:
    //- Private methods
    CombinedSpace(const ThisType& other);

    const ContainedMapperType& containedMapper() const { 
      return mapper_.containedMapper(); 
    }

  private:
    //- Member data  
    ContainedDiscreteFunctionSpaceType& spc_;

    MapperType mapper_;
    std::vector<BaseFunctionSetType*> baseSetVec_;

    static const int spaceId_;

  }; // end class CombinedSpace  

  //! Wrapper class for base function sets. This class is used within 
  //! CombinedSpace
  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  class CombinedBaseFunctionSet : 
    public BaseFunctionSetDefault<
    CombinedSpaceTraits<DiscreteFunctionSpaceImp, N, policy>
  >
  {
  public:
    //- Typedefs and enums
    enum { numComponents = N };
    typedef CombinedBaseFunctionSet<
      DiscreteFunctionSpaceImp, N, policy> ThisType;
    typedef CombinedSpaceTraits<
      DiscreteFunctionSpaceImp, N, policy> Traits;

    typedef typename Traits::DiscreteFunctionSpaceType 
    DiscreteFunctionSpaceType;
    typedef typename Traits::ContainedBaseFunctionSetType ContainedBaseFunctionSetType;
    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::JacobianRangeType JacobianRangeType;
    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::ContainedRangeType ContainedRangeType;
    typedef typename Traits::ContainedJacobianRangeType ContainedJacobianRangeType;

    typedef typename Traits::RangeFieldType DofType;

    typedef DofConversionUtility<PointBased>
      DofConversionUtilityType;
  public:
    //- Public methods
    //! Constructor
    CombinedBaseFunctionSet(const ContainedBaseFunctionSetType& bfSet) :
      containedResult_(0.0), 
      baseFunctionSet_(bfSet),
      util_(N)
    {}

    //! Number of base functions
    //! The number of base functions equals the total number of degrees of
    //! freedom (dof), since the dofs are considered to be scalar and the 
    //! combined base functions to be vector valued
    int getNumberOfBaseFunctions() const DUNE_DEPRECATED {
      return baseFunctionSet_.getNumberOfBaseFunctions()*N;
    }

    //! Number of base functions
    //! The number of base functions equals the total number of degrees of
    //! freedom (dof), since the dofs are considered to be scalar and the 
    //! combined base functions to be vector valued
    int numBaseFunctions() const {
      return baseFunctionSet_.numBaseFunctions()*N;
    }

    //! evaluate base function
    template <int diffOrd>
    void evaluate (int baseFunct, 
                   const FieldVector<deriType, diffOrd> &diffVariable, 
                   const DomainType & x, RangeType & phi ) const;

    //! evaluate base function at quadrature point
    template <int diffOrd, class QuadratureType >
    void evaluate (int baseFunct, 
                   const FieldVector<deriType, diffOrd> &diffVariable, 
                   QuadratureType & quad, 
                   int quadPoint, RangeType & phi ) const;

    //- Additional methods
    //! Number of distinct (scalar) base functions
    int numDifferentBaseFunctions() const {
      return baseFunctionSet_.numBaseFunctions();
    }

   //! evaluate base function
    void evaluateScalar(int baseFunct, 
                        const DomainType& x, 
                        ContainedRangeType& phi) const;
    
    //! evaluate base function at quadrature point
    void jacobianScalar(int baseFunct, 
                        const DomainType& x,
                        ContainedJacobianRangeType& phi) const;
   
    DofType evaluateSingle(int baseFunct, 
                           const DomainType& xLocal,
                           const RangeType& factor) const;
    
    template <class QuadratureType> 
    DofType evaluateSingle(int baseFunct, 
                           const QuadratureType & quad, int qp, 
                           const RangeType& factor) const;
    
    template <class Entity, class QuadratureType>
    DofType evaluateGradientSingle(int baseFunct,
                                   Entity& en,
                                   const QuadratureType & quad, int qp, 
                                   const JacobianRangeType& factor) const;

    template <class Entity>
    DofType evaluateGradientSingle(int baseFunct,
                                   Entity& en,
                                   const DomainType& xLocal,
                                   const JacobianRangeType& factor) const;

  private:
    //- Private methods
    CombinedBaseFunctionSet(const ThisType& other);

    //int containedDof(int combinedDofNum) const;
    //int component(int combinedDofNum) const;
    void expand(int baseFunct, 
                const ContainedRangeType& arg, 
                RangeType& dest) const;

  private:
    //- Data members
    mutable ContainedRangeType containedResult_;
    mutable ContainedRangeType phi_;
    mutable ContainedJacobianRangeType grad_;
    mutable DomainType gradScaled_;
    const ContainedBaseFunctionSetType& baseFunctionSet_;
    const DofConversionUtilityType util_;
  }; // end class CombinedBaseFunctionSet

  //! Wrapper class for mappers. This class is to be used in conjunction with
  //! the CombinedSpace
  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  class CombinedMapper : public DofMapperDefault<
    CombinedMapper<DiscreteFunctionSpaceImp, N, policy> 
  > 
  {
  public:
    //- Friends
    friend class CombinedSpace<DiscreteFunctionSpaceImp, N, policy>;
    friend class SubSpace<CombinedSpace<DiscreteFunctionSpaceImp, N, policy> >;

  public:
    //- Typedefs and enums
    enum { numComponents = N };
    typedef CombinedMapper<DiscreteFunctionSpaceImp, N, policy> ThisType;

    typedef CombinedSpaceTraits<
      DiscreteFunctionSpaceImp, N, policy> SpaceTraits;
    typedef DiscreteFunctionSpaceImp ContainedDiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceImp::Traits::MapperType ContainedMapperType;

    typedef CombinedDofConversionUtility<DiscreteFunctionSpaceImp,PointBased> LocalDofConversionUtilityType; 
    typedef CombinedDofConversionUtility<DiscreteFunctionSpaceImp,policy> GlobalDofConversionUtilityType;
  public:
    //- Public methods
    //! Constructor
    CombinedMapper(const ContainedDiscreteFunctionSpaceType& spc,
                   const ContainedMapperType& mapper) :
      spc_(spc),
      mapper_(mapper),
      utilLocal_(spc_,N),
      utilGlobal_(spc_,N)
    {}

    //! Total number of degrees of freedom
    int size() const;

    //! Map a local degree of freedom on an entity to a global one
    template <class EntityType>
    int mapToGlobal(EntityType& en, int localNum) const;

    //- Method inherited from mapper interface
    //! if grid has changed determine new size 
    //! (to be called once per timestep, therefore virtual )
    int newSize() const { return mapper_.newSize()*N; }
  
    //! old size
    int oldSize() const { return mapper_.oldSize()*N; }

    //! return max number of local dofs per entity 
    int numberOfDofs () const DUNE_DEPRECATED { return mapper_.numberOfDofs()*N; }

    //! return max number of local dofs per entity 
    int numDofs () const { return mapper_.numDofs()*N; }

    //! returns true if index is new ( for dof compress )
    bool indexNew (int num) const { 
      assert(false); // * check correctness first: do I correctly map between combinedMapper and containedMapper indices with utilGlobal
      return mapper_.indexNew(utilGlobal_.containedDof(num));
    }
    
    //! return old index in dof array of given index ( for dof compress ) 
    int oldIndex (int num) const; 
    
    //! return new index in dof array 
    int newIndex (int num) const;

    //! return estimate for size that is addtional needed for restriction 
    //! of data
    int additionalSizeEstimate() const {
      return mapper_.additionalSizeEstimate()*N;
    }

  private:
    //- Private methods
    CombinedMapper(const ThisType& other);

    const ContainedMapperType& containedMapper() const {
      return mapper_;
    }

    static int chooseSize(int pointBased, int variableBased,
                   Int2Type<PointBased>) {
      return pointBased;
    }

    static int chooseSize(int pointBased, int variableBased, 
                   Int2Type<VariableBased>) {
      return variableBased;
    }

  private:
    //- Data members
    const ContainedDiscreteFunctionSpaceType& spc_;
    const ContainedMapperType& mapper_;

    const LocalDofConversionUtilityType utilLocal_;
    GlobalDofConversionUtilityType utilGlobal_;
  }; // end class CombinedMapper

} // end namespace Dune

  // include implementation
#include "combinedspace.cc"

#endif
