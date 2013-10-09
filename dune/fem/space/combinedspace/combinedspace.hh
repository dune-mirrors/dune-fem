#ifndef DUNE_FEM_COMBINEDSPACE_HH
#define DUNE_FEM_COMBINEDSPACE_HH


#warning "This file is broken and will be fixed in time."
//- System includes
#include <vector>
#include <map>

//- Dune includes
#include <dune/common/static_assert.hh>
#include <dune/fem/version.hh>

//- Local includes
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>
#include <dune/fem/space/basefunctions/basefunctionsets.hh>
#include <dune/fem/space/basefunctions/basefunctionstorage.hh>
#include <dune/fem/space/basefunctions/basefunctionproxy.hh>

#include "combineddofstorage.hh"

#include "mapper.hh"

namespace Dune
{

  namespace Fem 
  {

    template< class FunctionSpaceImp,
              class GridPartImp,
              int polOrder,
              template< class > class BaseFunctionStorageImp >
    class LagrangeDiscreteFunctionSpace;


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

      template <class DFSpace>
      struct LagrangePointSetExporter
      {
        LagrangePointSetExporter( const DFSpace& spc ) {}
      };

      template <class FunctionSpaceImp, class GridPartImp, int polOrder, template <class> class BaseFunctionStorageImp >
      struct LagrangePointSetExporter<
          LagrangeDiscreteFunctionSpace< FunctionSpaceImp, GridPartImp, polOrder, BaseFunctionStorageImp > >
      {
        typedef LagrangeDiscreteFunctionSpace< FunctionSpaceImp, GridPartImp, polOrder, BaseFunctionStorageImp > LagrangeSpaceType;
        typedef typename LagrangeSpaceType :: LagrangePointSetType LagrangePointSetType;

        const LagrangeSpaceType& lagrangeSpace_;

        LagrangePointSetExporter( const LagrangeSpaceType& spc ) : lagrangeSpace_( spc ) {}
        template <class Entity>
        const LagrangePointSetType& lagrangePointSet( const Entity& entity ) const
        {
          return lagrangeSpace_.lagrangePointSet( entity );
        }
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
      static const int codimension = 0;
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
      typedef typename ContainedDiscreteFunctionSpaceType :: BlockMapperType
        ContainedBlockMapperType ;
      typedef NonBlockMapper< ContainedBlockMapperType,
                              ContainedDiscreteFunctionSpaceType :: localBlockSize >  ContainedMapperType;

      typedef typename ContainedFunctionSpaceType::DomainFieldType 
      DomainFieldType;
      typedef typename ContainedFunctionSpaceType::RangeFieldType 
      RangeFieldType;
      typedef typename ContainedFunctionSpaceType::RangeType 
      ContainedRangeType;
      typedef typename ContainedFunctionSpaceType::JacobianRangeType
      ContainedJacobianRangeType;

      typedef typename ContainedSpaceTraits::GridType GridType;
      typedef typename ContainedSpaceTraits::GridPartType GridPartType;
      typedef typename ContainedSpaceTraits::IndexSetType IndexSetType;
      typedef typename ContainedSpaceTraits::IteratorType IteratorType;

      typedef CombinedSpace<
        DiscreteFunctionSpaceImp, N, policy> DiscreteFunctionSpaceType;

      typedef FunctionSpace<
        DomainFieldType, RangeFieldType, 
        ContainedDimDomain, ContainedDimRange*N > FunctionSpaceType;

      enum { dimLocal = GridType :: dimension };

      typedef typename ToLocalFunctionSpace< FunctionSpaceType, dimLocal > :: Type 
        BaseFunctionSpaceType; 

      // type of singleton factory 
      typedef VectorialBaseFunctionSet< BaseFunctionSpaceType, CachingStorage >
        BaseFunctionSetImp;
      typedef SimpleBaseFunctionProxy<BaseFunctionSetImp> BaseFunctionSetType;

      typedef CombinedMapper< ContainedMapperType, N, policy > MapperType;
      typedef CombinedSubMapper< ContainedMapperType, N, policy > SubMapperType;
      
      enum { localBlockSize = BlockTraits :: localBlockSize };
      typedef typename BlockTraits :: BlockMapperType BlockMapperType;
     
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

      typedef CombinedDofConversionUtility< ContainedMapperType, N, policy >
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
    template< class DiscreteFunctionSpaceImp, int N>
    class CombinedSpace<DiscreteFunctionSpaceImp, N, VariableBased>
    : public DiscreteFunctionSpaceDefault
          < CombinedSpaceTraits< DiscreteFunctionSpaceImp, N, VariableBased > >,
      public CombinedSpaceHelper::LagrangePointSetExporter< DiscreteFunctionSpaceImp > 
    {
      static const DofStoragePolicy policy = VariableBased ;
    public:

      typedef CombinedSpaceTraits< DiscreteFunctionSpaceImp, N, policy > Traits;
      
    private:
      typedef DiscreteFunctionSpaceDefault< Traits > BaseType;

      typedef CombinedSpaceHelper::LagrangePointSetExporter< DiscreteFunctionSpaceImp > LagrangePointSetExporterType;

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
      typedef typename Traits::SubMapperType  SubMapperType;

      enum { spaceId_ = 13 };
     
      dune_static_assert( (Traits::ContainedDimRange == 1),
                          "Use CombinedSpace only with scalar spaces." );
    public:
      //! default communication interface 
      static const InterfaceType defaultInterface =
            ContainedDiscreteFunctionSpaceType :: defaultInterface;

      //! default communication direction 
      static const CommunicationDirection defaultDirection =
            ContainedDiscreteFunctionSpaceType :: defaultDirection;

      //- Public methods
      //! constructor
      explicit CombinedSpace( GridPartType &gridpart,
          const InterfaceType commInterface = defaultInterface ,
          const CommunicationDirection commDirection = defaultDirection )
      : BaseType( gridpart, commInterface, commDirection  ),
        LagrangePointSetExporterType( containedSpace_ ),
        containedSpace_( gridpart ),
        containedMapper_( containedSpace_.blockMapper() ),
        mapper_( containedMapper_ ),
        blockMapper_( Traits :: BlockTraits :: containedBlockMapper( containedSpace_ ) ),
        baseSetMap_(),
        dm_( DofManagerType :: instance( containedSpace_.gridPart().grid() ) )
      {
        const std::vector<GeometryType>& geomTypes = containedSpace_.geomTypes(0);
        int maxNumDofs = -1;
        // create mappers and base sets for all existing geom types
        for(size_t i=0; i<geomTypes.size(); ++i)
        {
          if(baseSetMap_.find(geomTypes[i]) == baseSetMap_.end())
          {
            BaseFunctionSetImp* baseSet =
              & SingletonProviderType::getObject(geomTypes[i]);
            // store in map 
            baseSetMap_[ geomTypes[i] ] = baseSet;
            // calc max dofs 
            maxNumDofs = std::max(maxNumDofs,int( baseSet->size() ));
          }
        }
      }

    private:
      // prohibit copying
      CombinedSpace ( const ThisType & );

    public:
      //! destructor
      ~CombinedSpace()
      {
        typedef typename BaseFunctionMapType :: iterator iterator;
        iterator end = baseSetMap_.end();
        for (iterator it = baseSetMap_.begin(); it != end; ++it)
        {
          BaseFunctionSetImp * set = (BaseFunctionSetImp *) (*it).second;
          SingletonProviderType::removeObject(*set);
        }
      }


      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::contains(const int codim) const */
      inline bool contains ( const int codim ) const
      {
        return containedSpace().contains( codim );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous() const */
      inline bool continuous () const
      {
        return containedSpace().continuous();
      }

#if 0
      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::polynomOrder() const */
      inline int polynomOrder () const
      {
        return containedSpace().polynomOrder();
      }
#endif

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order() const */
      inline int order () const
      {
        return containedSpace().order();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order(const EntityType &entity) const */
      template< class EntityType >
      inline int order ( const EntityType &entity) const
      {
        return containedSpace().order( entity );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::begin() const */
      inline IteratorType begin () const
      {
        return containedSpace().begin();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::end() const */
      inline IteratorType end () const
      {
        return containedSpace().end();
      }

      //! Return the identifier
      inline DFSpaceIdentifier type () const
      {
        return CombinedSpace_id;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::baseFunctionSet(const EntityType &entity) const */
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
      inline 
      DUNE_VERSION_DEPRECATED(1,4,remove)
      MapperType &mapper () const
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
   
      //! return a reference to the contained space's mapper
      inline ContainedMapperType &containedMapper () const
      { 
        return containedMapper_;
      }

      const DiscreteFunctionSpaceImp& containedSpace() const { return containedSpace_; }

    protected:
      DiscreteFunctionSpaceImp containedSpace_;
      ContainedMapperType containedMapper_;
      mutable MapperType mapper_;
      mutable BlockMapperType blockMapper_;

      typedef std::map< const GeometryType, BaseFunctionSetImp* > BaseFunctionMapType; 
      mutable BaseFunctionMapType baseSetMap_; 
      const DofManagerType & dm_;
    }; // end class CombinedSpace  


    // new implementation (needed by AdaptiveDiscreteFunction to extract sub functions)
    template< class DiscreteFunctionSpaceImp, int N >
    class CombinedSpace< DiscreteFunctionSpaceImp, N, PointBased >
    : public DiscreteFunctionSpaceImp :: 
        template ToNewDimRange< DiscreteFunctionSpaceImp:: dimRange * N > :: Type 
    {
      static const DofStoragePolicy policy = PointBased;

      typedef CombinedSpace< DiscreteFunctionSpaceImp, N, policy > ThisType;
    public:  
      typedef typename DiscreteFunctionSpaceImp :: 
        template ToNewDimRange< DiscreteFunctionSpaceImp:: dimRange * N > :: Type   BaseType;

      typedef typename BaseType :: GridPartType GridPartType;

      typedef DiscreteFunctionSpaceImp ContainedDiscreteFunctionSpaceType;
      typedef CombinedSubMapper< typename ContainedDiscreteFunctionSpaceType ::
        MapperType, N, policy > SubMapperType;

      typedef CombinedDofConversionUtility< 
        typename ContainedDiscreteFunctionSpaceType :: MapperType, N, policy >   DofConversionType;

      explicit CombinedSpace( GridPartType &gridPart,
          const InterfaceType commInterface = BaseType :: defaultInterface ,
          const CommunicationDirection commDirection = BaseType :: defaultDirection )
       : BaseType( gridPart, commInterface, commDirection ),
         containedSpace_( gridPart, commInterface, commDirection )
      {}
    
      //- Additional methods
      //! number of components
      int numComponents() const { return N; }

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

    protected:  
      //- Member data  
      ContainedDiscreteFunctionSpaceType containedSpace_;
    };

    /** @} **/  
    
  } // namespace Fem

} // namespace Dune

// include implementation
#include "combinedadaptmanager.hh"

#endif // #ifndef DUNE_FEM_COMBINEDSPACE_HH
