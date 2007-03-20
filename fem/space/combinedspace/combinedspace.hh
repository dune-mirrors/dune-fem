#ifndef DUNE_COMBINEDSPACE_HH
#define DUNE_COMBINEDSPACE_HH

//- System includes
#include <vector>

//- Dune includes
#include <dune/common/misc.hh>

//- Local includes
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/dofmapperinterface.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/common/basefunctioninterface.hh>
#include <dune/fem/space/basefunctions/basefunctionsets.hh>
#include <dune/fem/space/basefunctions/basefunctionstorage.hh>
#include "subspace.hh"

#include "combineddofstorage.hh"

namespace Dune {
  
  // Forward declarations
  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  class CombinedSpace;
  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  class CombinedMapper;

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
    enum { localBlockSize = N * ContainedSpaceTraits :: localBlockSize };
    
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

    // type of singleton factory 
    typedef VectorialBaseFunctionSet<FunctionSpace<double, double, ContainedDimDomain, N>, CachingStorage> BaseFunctionSetType;

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
    friend class CombinedMapper<DiscreteFunctionSpaceImp, N, policy>;
  };

  //! Class to combine N scalar spaces 
  //! Policies PointBased and VariableBased decide, how dof are stored in
  //! vectors, PointBased stores all local dofs consecutive, 
  //! VectorBased stores all dofs for one component consecutive 
  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy = PointBased>
  class CombinedSpace :
    public DiscreteFunctionSpaceDefault<
    CombinedSpaceTraits<DiscreteFunctionSpaceImp, N, policy> 
  > 
  {
  private:
    // CompileTimeChecker<policy==PointBased> OnlyWorksWithPointBasedPolicy;
    //- Private typedefs
    typedef DiscreteFunctionSpaceDefault<
    CombinedSpaceTraits<DiscreteFunctionSpaceImp, N, policy> 
    > BaseType;
  public:
    // polynomial Order is the same as for the single space 
    enum { CombinedFSpaceId = CombinedSpace_id };

    enum { polynomialOrder = DiscreteFunctionSpaceImp :: polynomialOrder };
    
    //- Public typedefs and enums
    typedef CombinedSpace<DiscreteFunctionSpaceImp, N, policy> ThisType;
    typedef CombinedSpaceTraits<DiscreteFunctionSpaceImp, N, policy> Traits;
    typedef DiscreteFunctionSpaceImp ContainedDiscreteFunctionSpaceType;
    typedef typename ContainedDiscreteFunctionSpaceType::FunctionSpaceType
    ContainedSpaceType;
    
    enum { localBlockSize = Traits :: localBlockSize };
    
    typedef DofManager<typename Traits::GridType> DofManagerType;
    typedef DofManagerFactory<DofManagerType> DofManagerFactoryType;

    typedef typename Traits::IteratorType IteratorType;

    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeFieldType RangeFieldType;
    typedef typename Traits::DomainFieldType DomainFieldType;

    typedef typename Traits::ContainedRangeType ContainedRangeType;
    typedef typename Traits::ContainedJacobianRangeType ContainedJacobianRangeType;

    typedef typename Traits::BaseFunctionSetType BaseFunctionSetType;
    typedef typename ContainedDiscreteFunctionSpaceType::ScalarFactoryType ScalarFactoryType;
    typedef BaseFunctionSetSingletonFactory<GeometryType,BaseFunctionSetType,
                ScalarFactoryType> SingletonFactoryType; 
    typedef SingletonList< GeometryType, BaseFunctionSetType,
            SingletonFactoryType > SingletonProviderType;
    typedef typename Traits::MapperType MapperType;
    typedef typename Traits::GridType GridType;
    typedef typename Traits::GridPartType GridPartType;
    typedef typename Traits::IndexSetType IndexSetType;

    typedef typename Traits::DofConversionType DofConversionType;
    typedef SubSpace<ThisType> SubSpaceType;

    enum { spaceId_ = 13 };
    
    CompileTimeChecker<(Traits::ContainedDimRange == 1)> use_CombinedSpace_only_with_scalar_spaces;
  public:
    //- Public methods
    //! constructor
    /*
    inline
    CombinedSpace(ContainedDiscreteFunctionSpaceType& spc);
    */
    //! constructor
    inline
    CombinedSpace(GridPartType& gridpart);

    //! destructor
    ~CombinedSpace();

    //! continuous?
    bool continuous() const { return spc_.continuous(); }

    //! polynom order
    int polynomOrder() const { return spc_.polynomOrder(); }

    //! polynom order
    int order() const { return spc_.order(); }

    //! begin iterator
    IteratorType begin() const { return spc_.begin(); }

    //! end iterator
    IteratorType end() const { return spc_.end(); }

    //! Return the identifier
    DFSpaceIdentifier type () const
    {
      return CombinedSpace_id;
    }

    //! total number of dofs
    int size() const { return mapper_.size(); }

    //! map a local dof number to a global one
    template <class EntityType>
    int mapToGlobal(EntityType& en, int local) const {
      return mapper_.mapToGlobal(en, local);
    }

    //! access to base function set
    template <class EntityType>
    const BaseFunctionSetType& baseFunctionSet(const EntityType& en) const 
    {
      GeometryIdentifier::IdentifierType id =
                GeometryIdentifier::fromGeometry(en.geometry()); 
      return this->baseFunctionSet( id );
    }

    //! access to base function set for given id 
    const BaseFunctionSetType& 
    baseFunctionSet(const GeometryIdentifier::IdentifierType id) const 
    {
      assert(id < (int) baseSetVec_.size());
      assert(id >= 0);

      assert( baseSetVec_[id] );
      return *baseSetVec_[id];
    }

    //! access to base function set
    template <class EntityType>
    const BaseFunctionSetType& getBaseFunctionSet(const EntityType& en) const 
    {
      return this->baseFunctionSet(en); 
    }

    //! access to grid
    const GridType& grid() const { return spc_.grid(); }
    
    //! access to gridPart
    GridPartType& gridPart() { return spc_.gridPart(); }
    //! access to gridPart
    const GridPartType& gridPart() const { return spc_.gridPart(); }

    //! \brief return corresponding index set  
    const IndexSetType& indexSet() const { return spc_.indexSet(); }

    //! access to mapper
    const MapperType& mapper() const { return mapper_; }

    //- Additional methods
    //! number of components
    int numComponents() const { return N; }

    //! return index in grid sequence 
    int sequence () const { return dm_.sequence(); }

    //! policy of this space
    DofStoragePolicy myPolicy() const{ return DofConversionType::policy(); }
 
    //! return subspace for ith component
    SubSpaceType& subSpace(int i) {
      assert( i >= 0 && i< N );
      assert( subSpaces_[i] );
      return *(subSpaces_[i]);
    }
  const ContainedDiscreteFunctionSpaceType& containedSpace() const
    {return spc_;}
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

  protected:
    //- Member data  
    ContainedDiscreteFunctionSpaceType spc_;

    MapperType mapper_;
    std::vector<BaseFunctionSetType*> baseSetVec_;
    std::vector<const BaseFunctionSetType*> baseSecVec_;
    FieldVector<SubSpaceType*,N> subSpaces_;
    const DofManagerType & dm_;

  }; // end class CombinedSpace  

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
    inline
    int size() const;

    //! Map a local degree of freedom on an entity to a global one
    template <class EntityType>
    inline
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

    //! return old index in dof array of given index ( for dof compress ) 
    inline
    int oldIndex (int num) const; 
    
    //! return new index in dof array 
    inline
    int newIndex (int num) const;

    //! return estimate for size that is addtional needed for restriction 
    //! of data
    int additionalSizeEstimate() const {
      return mapper_.additionalSizeEstimate()*N;
    }
    //! return true if compress will affect data  
    bool needsCompress () const 
    {
      return mapper_.needsCompress ();
    }
    //! return number of holes in the data 
    int numberOfHoles() const 
    {
      return mapper_.numberOfHoles()*N; 
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
