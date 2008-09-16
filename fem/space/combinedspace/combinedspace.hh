#ifndef DUNE_COMBINEDSPACE_HH
#define DUNE_COMBINEDSPACE_HH

//- System includes
#include <vector>
#include <map>

//- Dune includes
#include <dune/common/misc.hh>

//- Local includes
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/common/basefunctioninterface.hh>
#include <dune/fem/space/basefunctions/basefunctionsets.hh>
#include <dune/fem/space/basefunctions/basefunctionstorage.hh>
#include <dune/fem/space/basefunctions/basefunctionproxy.hh>

#include "combineddofstorage.hh"
#include "mapper.hh"

namespace Dune
{
  
  namespace CombinedSpaceHelper
  {

    template< class ContainedSpace, int N, DofStoragePolicy policy >
    struct BlockTraits;

    template< class ContainedSpace, int N >
    struct BlockTraits< ContainedSpace, N, PointBased >
    {
      enum { localBlockSize = N * ContainedSpace :: localBlockSize };

      typedef typename ContainedSpace :: Traits :: BlockMapperType
        ContainedBlockMapperType;

      typedef CombinedMapper< ContainedBlockMapperType, 1, PointBased >
        BlockMapperType;

      inline static ContainedBlockMapperType &
      containedBlockMapper( const ContainedSpace &space )
      {
        return space.blockMapper();
      }
    };

    template< class ContainedSpace, int N >
    struct BlockTraits< ContainedSpace, N, VariableBased >
    {
      enum { localBlockSize = 1 };
      
      typedef typename ContainedSpace :: Traits :: MapperType
        ContainedBlockMapperType;
      
      typedef CombinedMapper< ContainedBlockMapperType, N, VariableBased >
        BlockMapperType;

      inline static ContainedBlockMapperType &
      containedBlockMapper( const ContainedSpace &space )
      {
        return space.mapper();
      }

      /* This would be the better definition. The problem is that the DoFs
       * are not ordered blockwise
       */
#if 0
      enum { localBlockSize = ContainedSpace :: localBlockSize };

      typedef CombinedMapper
        < typename ContainedSpace :: Traits :: BlockMapperType, N, VariableBased >
        BlockMapperType;

      inline static ContainedBlockMapperType &
      containedBlockMapper( const ContainedSpace &space )
      {
        return space.blockMapper();
      }
#endif
    };

  }


  
  /** @addtogroup CombinedSpace
      Class to combine N scalar discrete function spaces.  
      Policies PointBased and VariableBased decide, how dof are stored in
      vectors. PointBased stores all local dofs consecutive, 
      VectorBased stores all dofs for one component consecutive. 
   @{
  
  */
  
  // Forward declarations
  template< class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy = PointBased >
  class CombinedSpace;



  //! Traits class for CombinedSpace
  template< class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy >
  struct CombinedSpaceTraits
  {
  private:
    typedef DiscreteFunctionSpaceImp ContainedDiscreteFunctionSpaceType;

    typedef typename ContainedDiscreteFunctionSpaceType::Traits 
    ContainedSpaceTraits;
    typedef typename ContainedSpaceTraits::FunctionSpaceType 
    ContainedFunctionSpaceType;
    typedef typename ContainedSpaceTraits::BaseFunctionSetType 
    ContainedBaseFunctionSetType;
    
    enum { ContainedDimRange = ContainedFunctionSpaceType::dimRange,
           ContainedDimDomain = ContainedFunctionSpaceType::dimDomain };

    typedef CombinedSpaceHelper :: BlockTraits< DiscreteFunctionSpaceImp, N, policy >
      BlockTraits;

  public:
    typedef typename ContainedSpaceTraits :: MapperType ContainedMapperType;

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

    typedef CombinedMapper< ContainedMapperType, N, policy > MapperType;
    
    enum { localBlockSize = BlockTraits :: localBlockSize };
    //enum { localBlockSize = N * ContainedSpaceTraits :: localBlockSize };
    typedef typename BlockTraits :: BlockMapperType BlockMapperType;
   
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

    typedef typename ContainedSpaceTraits::GridType GridType;
    typedef typename ContainedSpaceTraits::GridPartType GridPartType;
    typedef typename ContainedSpaceTraits::IndexSetType IndexSetType;
    typedef typename ContainedSpaceTraits::IteratorType IteratorType;

    typedef CombinedDofConversionUtility< ContainedMapperType, policy >
      DofConversionType;

    enum { dimRange = FunctionSpaceType :: dimRange,
           dimDomain = FunctionSpaceType :: dimDomain };

    //! \brief defines type of data handle for communication 
    template <class DiscreteFunctionImp,
              class OperationImp = DFCommunicationOperation :: Copy>
    struct CommDataHandle
    {
      //! \brief uses type of contained space 
      typedef typename ContainedDiscreteFunctionSpaceType:: template
        CommDataHandle<DiscreteFunctionImp,OperationImp> :: Type Type;
      //! \brief type of operation to perform on scatter 
      typedef typename ContainedDiscreteFunctionSpaceType:: template
        CommDataHandle<DiscreteFunctionImp,OperationImp> :: OperationType OperationType;
    };

  private:
    //- Friends
    friend class CombinedSpace< DiscreteFunctionSpaceImp, N, policy >;
    friend class CombinedMapper< ContainedMapperType, N, policy >;
  };


  

  

  /** @brief 
      Combined Space Function Space
      **/
  template< class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy >
  class CombinedSpace
  : public DiscreteFunctionSpaceDefault
    < CombinedSpaceTraits< DiscreteFunctionSpaceImp, N, policy > > 
  {
  public:
    typedef CombinedSpaceTraits< DiscreteFunctionSpaceImp, N, policy > Traits;
    
  private:
    typedef DiscreteFunctionSpaceDefault< Traits > BaseType;

  public:
    // polynomial Order is the same as for the single space 
    enum { CombinedFSpaceId = CombinedSpace_id };

    enum { polynomialOrder = DiscreteFunctionSpaceImp :: polynomialOrder };
    
    //- Public typedefs and enums
    typedef CombinedSpace<DiscreteFunctionSpaceImp, N, policy> ThisType;
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
    
    typedef typename Traits :: MapperType MapperType;
    typedef typename Traits :: BlockMapperType BlockMapperType;

    typedef typename Traits :: ContainedMapperType ContainedMapperType;
 
    typedef typename Traits::GridType GridType;
    typedef typename Traits::GridPartType GridPartType;
    typedef typename Traits::IndexSetType IndexSetType;

    typedef typename Traits::DofConversionType DofConversionType;
    typedef CombinedSubMapper<ThisType> SubMapperType;

    enum { spaceId_ = 13 };
    
    CompileTimeChecker<(Traits::ContainedDimRange == 1)>
      use_CombinedSpace_only_with_scalar_spaces;
  public:
    //! default communication interface 
    static const InterfaceType defaultInterface =
          ContainedDiscreteFunctionSpaceType :: defaultInterface;

    //! default communication direction 
    static const CommunicationDirection defaultDirection =
          ContainedDiscreteFunctionSpaceType :: defaultDirection;

    //- Public methods
    //! constructor
    inline explicit CombinedSpace( GridPartType &gridpart,
        const InterfaceType commInterface = defaultInterface ,
        const CommunicationDirection commDirection = defaultDirection );

  private:
    // prohibit copying
    CombinedSpace ( const ThisType & );

  public:
    //! destructor
    ~CombinedSpace();

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::contains(const int codim) const */
    inline bool contains ( const int codim ) const
    {
      return containedSpace().contains( codim );
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::continuous() const */
    inline bool continuous () const
    {
      return containedSpace().continuous();
    }

#if 0
    /** \copydoc Dune::DiscreteFunctionSpaceInterface::polynomOrder() const */
    inline int polynomOrder () const
    {
      return containedSpace().polynomOrder();
    }
#endif

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::order() const */
    inline int order () const
    {
      return containedSpace().order();
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::begin() const */
    inline IteratorType begin () const
    {
      return containedSpace().begin();
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::end() const */
    inline IteratorType end () const
    {
      return containedSpace().end();
    }

    //! Return the identifier
    inline DFSpaceIdentifier type () const
    {
      return CombinedSpace_id;
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::baseFunctionSet(const EntityType &entity) const */
    template< class EntityType >
    const BaseFunctionSetType baseFunctionSet ( const EntityType &entity ) const
    {
      return baseFunctionSet( entity.geometry().type() );
    }

    //! access to base function set for given id 
    const BaseFunctionSetType baseFunctionSet ( const GeometryType type ) const
    {
      assert( baseSetMap_.find( type ) != baseSetMap_.end() );
      return BaseFunctionSetType( baseSetMap_[ type ] );
    }

    //! access to mapper
    inline MapperType &mapper () const
    {
      return mapper_;
    }

    //! access to mapper
    inline BlockMapperType &blockMapper () const
    {
      return blockMapper_;
    }

    //- Additional methods
    //! number of components
    int numComponents() const { return N; }

    //! return index in grid sequence 
    int sequence () const { return dm_.sequence(); }

    //! policy of this space
    inline DofStoragePolicy myPolicy() const
    {
      return DofConversionType :: policy();
    }
 
    //! return reference to contained space  
    inline const ContainedDiscreteFunctionSpaceType &containedSpace () const
    {
      return containedSpace_;
    }

    //! return a reference to the contained space's mapper
    inline ContainedMapperType &containedMapper () const
    { 
      return containedSpace().mapper();
    }

  protected:
    //- Member data  
    ContainedDiscreteFunctionSpaceType containedSpace_;

    mutable MapperType mapper_;
    mutable BlockMapperType blockMapper_;

    typedef std::map< const GeometryType, BaseFunctionSetImp* > BaseFunctionMapType; 
    mutable BaseFunctionMapType baseSetMap_; 
    const DofManagerType & dm_;
  }; // end class CombinedSpace  

  /** @} **/  
  
} // end namespace Dune

// include implementation
#include "combinedspace.cc"
#include "combinedadaptmanager.hh"

#endif
