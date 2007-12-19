#ifndef DUNE_COMBINEDSPACE_HH
#define DUNE_COMBINEDSPACE_HH

//- System includes
#include <vector>
#include <map>

//- Dune includes
#include <dune/common/misc.hh>

//- Local includes
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/dofmapperinterface.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/common/basefunctioninterface.hh>
#include <dune/fem/space/basefunctions/basefunctionsets.hh>
#include <dune/fem/space/basefunctions/basefunctionstorage.hh>
#include <dune/fem/space/basefunctions/basefunctionproxy.hh>

#include "combineddofstorage.hh"

namespace Dune {
  /** @addtogroup CombinedSpace
      Class to combine N scalar discrete function spaces.  
      Policies PointBased and VariableBased decide, how dof are stored in
      vectors. PointBased stores all local dofs consecutive, 
      VectorBased stores all dofs for one component consecutive. 
   @{
  
  */
  
  // Forward declarations
  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  class CombinedSpace;

  template< class ContainedSpaceImp, int N, DofStoragePolicy policy >
  class CombinedMapper;

  template< class CombinedMapperTraits >
  class CombinedDofMapIterator;



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
    typedef VectorialBaseFunctionSet< FunctionSpaceType, CachingStorage >
      BaseFunctionSetImp;
    typedef VectorialBaseFunctionProxy<BaseFunctionSetImp> BaseFunctionSetType;

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

  
  /** @brief 
      Combined Space Function Space
      **/
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

    typedef typename Traits::BaseFunctionSetImp  BaseFunctionSetImp;
    typedef typename Traits::BaseFunctionSetType BaseFunctionSetType;
    typedef typename ContainedDiscreteFunctionSpaceType::ScalarFactoryType ScalarFactoryType;
    
    typedef BaseFunctionSetSingletonFactory<GeometryType,BaseFunctionSetImp,
                ScalarFactoryType> SingletonFactoryType; 
    typedef SingletonList< GeometryType, BaseFunctionSetImp,
            SingletonFactoryType > SingletonProviderType;
    typedef typename Traits::MapperType MapperType;
    typedef typename Traits::GridType GridType;
    typedef typename Traits::GridPartType GridPartType;
    typedef typename Traits::IndexSetType IndexSetType;

    typedef typename Traits::DofConversionType DofConversionType;
    typedef CombinedSubMapper<ThisType> SubMapperType;

    typedef typename ContainedDiscreteFunctionSpaceType :: 
    BlockMapperType BlockMapperType;
    enum { spaceId_ = 13 };
    
    CompileTimeChecker<(Traits::ContainedDimRange == 1)>
      use_CombinedSpace_only_with_scalar_spaces;
  public:
    //- Public methods
    //! constructor
    /*
    inline
    CombinedSpace(ContainedDiscreteFunctionSpaceType& spc);
    */
    //! constructor
    inline explicit CombinedSpace( GridPartType &gridpart );

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

    //! access to base function set
    template <class EntityType>
    const BaseFunctionSetType baseFunctionSet(const EntityType& en) const 
    {
      return this->baseFunctionSet( en.geometry().type() );
    }

    //! access to base function set for given id 
    const BaseFunctionSetType 
    baseFunctionSet(const GeometryType geomType) const 
    {
      assert(baseSetMap_.find( geomType ) != baseSetMap_.end());
      return BaseFunctionSetType(baseSetMap_[geomType]);
    }

    //! access to mapper
    MapperType& mapper() const { return mapper_; }

    //! access to mapper
    BlockMapperType& blockMapper() const { return spc_.blockMapper(); }

    //- Additional methods
    //! number of components
    int numComponents() const { return N; }

    //! return index in grid sequence 
    int sequence () const { return dm_.sequence(); }

    //! policy of this space
    DofStoragePolicy myPolicy() const{ return DofConversionType::policy(); }
 
    //! return reference to contained space  
    const ContainedDiscreteFunctionSpaceType& containedSpace() const  { return spc_; }

  private:
    //- Private typedefs
    typedef typename Traits::ContainedMapperType ContainedMapperType;
   
  private:
    //- Private methods
    CombinedSpace(const ThisType& other);

    ContainedMapperType& containedMapper() const { 
      return mapper_.containedMapper(); 
    }

  protected:
    //- Member data  
    ContainedDiscreteFunctionSpaceType spc_;

    mutable MapperType mapper_;
    typedef std::map< const GeometryType, BaseFunctionSetImp* > BaseFunctionMapType; 
    mutable BaseFunctionMapType baseSetMap_; 
    const DofManagerType & dm_;

  }; // end class CombinedSpace  



  template< class ContainedSpaceImp, int N, DofStoragePolicy policy >
  struct CombinedMapperTraits
  {
    typedef ContainedSpaceImp ContainedDiscreteFunctionSpaceType;

    enum { numComponents = N };

    typedef CombinedMapperTraits
      < ContainedDiscreteFunctionSpaceType, numComponents, policy >
      Traits;

    typedef typename ContainedDiscreteFunctionSpaceType :: MapperType
      ContainedMapperType;

    typedef typename ContainedDiscreteFunctionSpaceType :: GridPartType
      GridPartType;

    typedef typename GridPartType :: template Codim< 0 > :: IteratorType :: Entity
      EntityType;
    
    typedef CombinedDofConversionUtility
      < ContainedDiscreteFunctionSpaceType, policy >
      GlobalDofConversionUtilityType;
    
    typedef CombinedDofMapIterator< Traits > DofMapIteratorType;

    typedef CombinedMapper
      < ContainedDiscreteFunctionSpaceType, numComponents, policy >
      DofMapperType;
  };

  

  template< class CombinedMapperTraits >
  class CombinedDofMapIterator
  {
  public:
    typedef CombinedMapperTraits Traits;

    enum IteratorType { beginIterator, endIterator };

  private:
    typedef CombinedDofMapIterator< Traits > ThisType;

  public:
    typedef typename Traits :: EntityType EntityType;
    
    typedef typename Traits :: ContainedMapperType ContainedMapperType;

    typedef typename Traits :: GlobalDofConversionUtilityType
      GlobalDofConversionUtilityType;

    enum { numComponents = Traits :: numComponents };

  protected:
    const EntityType &entity_;
    const ContainedMapperType &mapper_;
    const int numScalarDofs_;
    const GlobalDofConversionUtilityType dofUtil_;

    int dof_, component_;
    int global_;
    
  public:
    inline CombinedDofMapIterator ( const IteratorType type,
                                    const EntityType &entity,
                                    const ContainedMapperType &mapper,
                                    const GlobalDofConversionUtilityType &dofUtil )
    : entity_( entity ),
      mapper_( mapper ),
      dofUtil_( dofUtil ),
      numScalarDofs_( mapper_.numDofs( entity ) ),
      dof_( type == beginIterator ? 0 : numScalarDofs_ ),
      component_( type == beginIterator ? 0 : numComponents ),
      global_( type == beginIterator ?  mapper_.mapToGlobal( entity_, dof_ ) : 0 )
    {}

    inline CombinedDofMapIterator ( const ThisType &other )
    : entity_( other.entity_ ),
      mapper_( other.mapper_ ),
      dofUtil_( other.dofUtil ),
      numScalarDofs_( other.numScalarDofs_ ),
      dof_( other.dof_ ),
      component_( other.component_ ),
      global_( other.global_ )
    {}

    inline ThisType &operator++ ()
    {
      ++component_;
      if( component_ >= numComponents )
      {
        ++dof_;
        if( dof_ < numScalarDofs_ )
          global_ = mapper_.mapToGlobal( entity_, dof_ );
        component_ = 0;
      }
      return *this;
    }

    inline bool operator== ( const ThisType &other ) const
    {
      return (dof_ == other.dof_) && (component_ == other.component_);
    }

    inline bool operator!= ( const ThisType &other ) const
    {
      return (dof_ != other.dof_) || (component_ != other.component_);
    }

    inline int local () const
    {
      assert( dof_ < numScalarDofs_ );
      return dof_ * numComponents + component_;
    }

    inline int global () const
    {
      assert( dof_ < numScalarDofs_ );
      return dofUtil_.combinedDof( global_, component_ );
    }

    inline int component () const
    {
      return component_;
    }

    inline int localScalar () const
    {
      return dof_;
    }

    inline int globalScalar () const
    {
      return global_;
    }
  };



  //! Wrapper class for mappers. This class is to be used in conjunction with
  //! the CombinedSpace
  template< class ContainedSpaceImp, int N, DofStoragePolicy policy >
  class CombinedMapper
  : public DofMapperDefault< CombinedMapperTraits< ContainedSpaceImp, N, policy > >
  {
  public:
    typedef CombinedMapperTraits< ContainedSpaceImp, N, policy > Traits;

    enum { numComponents = Traits :: numComponents };

    typedef typename Traits :: ContainedDiscreteFunctionSpaceType
      ContainedDiscreteFunctionSpaceType;

  private:
    typedef CombinedMapper
      < ContainedDiscreteFunctionSpaceType, numComponents, policy >
      ThisType;

    friend class CombinedSpace
      < ContainedDiscreteFunctionSpaceType, numComponents, policy>;

  public:
    typedef CombinedSpaceTraits
      < ContainedDiscreteFunctionSpaceType, numComponents, policy >
      SpaceTraits;

    typedef typename Traits :: ContainedMapperType ContainedMapperType;

    typedef typename Traits :: EntityType EntityType;

    typedef typename Traits :: DofMapIteratorType DofMapIteratorType;

    typedef typename Traits :: GlobalDofConversionUtilityType
      GlobalDofConversionUtilityType;
    
    typedef CombinedDofConversionUtility
      < ContainedDiscreteFunctionSpaceType, PointBased >
      LocalDofConversionUtilityType;

  public:
    //- Public methods
    //! Constructor
    CombinedMapper ( const ContainedDiscreteFunctionSpaceType &spc,
                     ContainedMapperType &mapper )
    : spc_( spc ),
      mapper_( mapper ),
      utilLocal_( spc_, numComponents ),
      utilGlobal_( spc_, numComponents ),
      oldSize_( spc_.size() ),
      size_( spc_.size() )
    {}

  private:
    // prohibit copying
    CombinedMapper ( const ThisType & );

  public:
    //! Total number of degrees of freedom
    inline int size () const;
    
    /** \copydoc Dune::DofMapperInterface::begin(const EntityType &entity) const */
    inline DofMapIteratorType begin ( const EntityType &entity ) const
    {
      return DofMapIteratorType
        ( DofMapIteratorType :: beginIterator, entity, mapper_, utilGlobal_ );
    }
    
    /** \copydoc Dune::DofMapperInterface::end(const EntityType &entity) const */
    inline DofMapIteratorType end ( const EntityType &entity ) const
    {
      return DofMapIteratorType
        ( DofMapIteratorType :: endIterator, entity, mapper_, utilGlobal_ );
    }

    //! Map a local degree of freedom on an entity to a global one
    inline int mapToGlobal( const EntityType &entity,
                            int localNum ) const;

    //- Method inherited from mapper interface
    //! if grid has changed determine new size 
    //! (to be called once per timestep, therefore virtual )
    int newSize() const
    {
      return mapper_.newSize() * numComponents;
    }
  
    /*
    //! old size
    int oldSize() const { return mapper_.oldSize()*N; }
    */

    /** \copydoc Dune::DofMapperInterface::numDofs() const */
    inline int numDofs () const
    {
      return mapper_.numDofs() * numComponents;
    }

    /** \copydoc Dune::DofMapperInterface::numDofs(const EntityType &entity) const */
    inline int numDofs ( const EntityType &entity ) const
    {
      return mapper_.numDofs( entity ) * numComponents;
    }

    //! return old index in dof array of given index ( for dof compress ) 
    inline int oldIndex ( const int hole, const int block ) const; 
    
    //! return new index in dof array 
    inline int newIndex ( const int hole, const int block ) const;

    //! return number of holes in the data 
    int numberOfHoles ( const int block ) const;
  
    //! returnn number of mem blocks 
    int numBlocks () const; 

    //! update offset information
    void update (); 
    
    //! return current old offset of block 
    int oldOffSet ( const int block ) const;

    //! return current offset of block 
    int offSet ( const int block ) const;

    //! return true if compress will affect data  
    bool needsCompress () const;

  protected:
    ContainedMapperType &containedMapper () const
    {
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

  protected:
    //- Data members
    const ContainedDiscreteFunctionSpaceType &spc_;
    mutable ContainedMapperType &mapper_;

    const LocalDofConversionUtilityType utilLocal_;
    GlobalDofConversionUtilityType utilGlobal_;
    int oldSize_, size_;
  }; // end class CombinedMapper
  
  /** @} **/  
  
} // end namespace Dune

// include implementation
#include "combinedspace.cc"
#include "combinedadaptmanager.hh"

#endif
