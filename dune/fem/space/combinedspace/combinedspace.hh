#ifndef DUNE_FEM_COMBINEDSPACE_HH
#define DUNE_FEM_COMBINEDSPACE_HH

//- System includes
#include <vector>
#include <map>

//- Dune includes
#include <dune/fem/version.hh>

//- Local includes
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/combineddiscretefunctionspace.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>

#include <dune/fem/space/basisfunctionset/vectorial.hh>

#include "combineddofstorage.hh"
#include "mapper.hh"
#include "lagrangepointsetexporter.hh"

namespace Dune
{

  namespace Fem
  {


    namespace CombinedSpaceHelper
    {

      template< class ContainedSpace, int N, DofStoragePolicy policy >
      struct BlockTraits;

      template< class ContainedSpace, int N >
      struct BlockTraits< ContainedSpace, N, PointBased >
      {
        enum { containedLocalBlockSize = ContainedSpace :: localBlockSize };
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
        enum { localBlockSize = ContainedSpace :: localBlockSize };

        typedef typename ContainedSpace :: Traits :: BlockMapperType
          ContainedBlockMapperType;

        typedef CombinedMapper< ContainedBlockMapperType, N, VariableBased >
          BlockMapperType;

        inline static ContainedBlockMapperType &
        containedBlockMapper( const ContainedSpace &space )
        {
          return space.blockMapper();
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
    template< class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy >
    class CombinedSpace;

    template< class ContainedDiscreteFunctionSpaceImp,
              class NewFunctionSpace>
    struct DifferentDiscreteFunctionSpace;

    // specialization of DifferentDiscreteFunctionSpace for this CombinedSapce
    template <class ContainedSpace, int N, DofStoragePolicy policy, class NewFunctionSpace>
    struct DifferentDiscreteFunctionSpace< CombinedSpace< ContainedSpace, N, policy >, NewFunctionSpace >
    {
      typedef CombinedSpace< ContainedSpace, NewFunctionSpace::dimRange, policy > Type;
    };


    //! Traits class for CombinedSpace
    template< class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy >
    struct CombinedSpaceTraits
    {
      static const int codimension = DiscreteFunctionSpaceImp::Traits::codimension;

    private:
      typedef DiscreteFunctionSpaceImp ContainedDiscreteFunctionSpaceType;

      typedef typename ContainedDiscreteFunctionSpaceType::Traits
      ContainedSpaceTraits;
      typedef typename ContainedSpaceTraits::FunctionSpaceType
      ContainedFunctionSpaceType;
      typedef typename ContainedSpaceTraits::BasisFunctionSetType
      ContainedBasisFunctionSetType;

      enum { ContainedDimRange  = ContainedFunctionSpaceType::dimRange,
             ContainedDimDomain = ContainedFunctionSpaceType::dimDomain };

      typedef CombinedSpaceHelper :: BlockTraits< DiscreteFunctionSpaceImp, N, policy >
        BlockTraits;

    public:
      typedef typename ContainedDiscreteFunctionSpaceType::BlockMapperType
        ContainedBlockMapperType;

      typedef typename ContainedFunctionSpaceType::DomainFieldType
      DomainFieldType;
      typedef typename ContainedFunctionSpaceType::RangeFieldType
      RangeFieldType;
      typedef typename ContainedFunctionSpaceType::RangeType
      ContainedRangeType;
      typedef typename ContainedFunctionSpaceType::JacobianRangeType
      ContainedJacobianRangeType;

      typedef typename ContainedDiscreteFunctionSpaceType::GridType GridType;
      typedef typename ContainedDiscreteFunctionSpaceType::GridPartType GridPartType;
      typedef typename ContainedDiscreteFunctionSpaceType::IndexSetType IndexSetType;
      typedef typename ContainedDiscreteFunctionSpaceType::IteratorType IteratorType;

      typedef CombinedSpace<
        DiscreteFunctionSpaceImp, N, policy> DiscreteFunctionSpaceType;

      typedef CombinedSubMapper< typename ContainedDiscreteFunctionSpaceType ::
        BlockMapperType, N, policy > SubMapperType;

      typedef FunctionSpace<
        DomainFieldType, RangeFieldType,
        ContainedDimDomain, ContainedDimRange*N > FunctionSpaceType;

      // coordinates tpyes of this space
      typedef typename FunctionSpaceType::RangeType          RangeType;
      typedef typename FunctionSpaceType::DomainType         DomainType;
      typedef typename FunctionSpaceType::JacobianRangeType  JacobianRangeType;

      enum { dimRange  = FunctionSpaceType :: dimRange,
             dimDomain = FunctionSpaceType :: dimDomain };

      // type of Vectorial BasisFunctionSet
      // NOTE: local dof ordering is always point based !!!
      // VerticalDofAlignment   == PointBased 
      // HorizontalDofAlignment == VariableBased
      typedef VectorialBasisFunctionSet< ContainedBasisFunctionSetType, RangeType, VerticalDofAlignment >  BasisFunctionSetType;

      enum { localBlockSize = BlockTraits :: localBlockSize };
      typedef typename BlockTraits :: BlockMapperType BlockMapperType;

      // we will need further the SubBlockMapper and other stuff
      typedef CombinedDofConversionUtility< ContainedBlockMapperType, N, policy >
        DofConversionType;

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
      friend class CombinedMapper< ContainedBlockMapperType, N, policy >;
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

      enum { localBlockSize = Traits :: localBlockSize };

      typedef DofManager<typename Traits::GridType> DofManagerType;

      typedef typename Traits::IteratorType IteratorType;

      typedef typename Traits::RangeType RangeType;
      typedef typename Traits::DomainType DomainType;
      typedef typename Traits::RangeFieldType RangeFieldType;
      typedef typename Traits::DomainFieldType DomainFieldType;

      typedef typename Traits::ContainedRangeType ContainedRangeType;
      typedef typename Traits::ContainedJacobianRangeType ContainedJacobianRangeType;

      typedef typename Traits::BasisFunctionSetType  BasisFunctionSetType;

      typedef typename Traits::ContainedBlockMapperType ContainedBlockMapperType;
      typedef typename Traits::BlockMapperType BlockMapperType;

      typedef typename Traits::GridType GridType;
      typedef typename Traits::GridPartType GridPartType;
      typedef typename Traits::IndexSetType IndexSetType;

      typedef typename Traits::DofConversionType  DofConversionType;
      typedef typename Traits::SubMapperType      SubMapperType;

      enum { spaceId_ = 13 };

      static_assert( (Traits::ContainedDimRange == 1),
                      "Use CombinedSpace only with scalar spaces." );
    public:
      //- Public methods
      //! constructor
      explicit CombinedSpace( GridPartType &gridPart,
                              const InterfaceType commInterface = InteriorBorder_All_Interface,
                              const CommunicationDirection commDirection = ForwardCommunication )
      : BaseType( gridPart, commInterface, commDirection  ),
        LagrangePointSetExporterType( containedSpace_ ),
        containedSpace_( gridPart ),
        blockMapper_( Traits :: BlockTraits :: containedBlockMapper( containedSpace_ ) )
      {
        DofManagerType::instance( gridPart.grid() ).addIndexSet( blockMapper_ );
      }


      //! destructor removing mapper from DofManagers index set list
      ~CombinedSpace() 
      { 
        DofManagerType::instance( gridPart().grid() ).removeIndexSet( blockMapper_ );
      } 

    private:
      // prohibit copying
      CombinedSpace ( const ThisType & );

    public:
      using BaseType :: gridPart;

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::contains(const int codim) const */
      bool contains ( const int codim ) const
      {
        return containedSpace().contains( codim );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous() const */
      bool continuous () const
      {
        return containedSpace().continuous();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order() const */
      int order () const
      {
        return containedSpace().order();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order(const EntityType &entity) const */
      template< class EntityType >
      int order ( const EntityType &entity) const
      {
        return containedSpace().order( entity );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::begin() const */
      IteratorType begin () const
      {
        return containedSpace().begin();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::end() const */
      IteratorType end () const
      {
        return containedSpace().end();
      }

      //! Return the identifier
      DFSpaceIdentifier type () const
      {
        return CombinedSpace_id;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::baseFunctionSet(const EntityType &entity) const */
      template< class EntityType >
      const BasisFunctionSetType basisFunctionSet ( const EntityType &entity ) const
      {
        return BasisFunctionSetType( containedSpace().basisFunctionSet( entity ) );
      }

      //! access to mapper
      BlockMapperType &blockMapper () const
      {
        return blockMapper_;
      }

      //- Additional methods
      //! number of components
      int numComponents() const { return N; }

      //! policy of this space
      DofStoragePolicy myPolicy() const
      {
        return DofConversionType :: policy();
      }

      //! return a reference to the contained space's mapper
      ContainedBlockMapperType &containedBlockMapper () const
      {
        return containedSpace().blockMapper();
      }

      const ContainedDiscreteFunctionSpaceType& containedSpace() const
      {
        return containedSpace_;
      }

    protected:
      ContainedDiscreteFunctionSpaceType containedSpace_;
      mutable BlockMapperType blockMapper_;
    }; // end class CombinedSpace


    template < class DiscreteFunctionSpaceImp,
               int N>
    class DefaultLocalRestrictProlong< CombinedSpace<DiscreteFunctionSpaceImp, N, VariableBased> >
    : public DiscontinuousGalerkinLocalRestrictProlong< CombinedSpace<DiscreteFunctionSpaceImp, N, VariableBased>, false >
    {
      typedef DiscontinuousGalerkinLocalRestrictProlong< CombinedSpace<DiscreteFunctionSpaceImp, N, VariableBased>, false >
        BaseType;
    public:

      DefaultLocalRestrictProlong( const CombinedSpace<DiscreteFunctionSpaceImp,N, VariableBased> & space)
        : BaseType( space )
      {
        assert( !space.continuous() );
      }
    };


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
        BlockMapperType, N, policy > SubMapperType;

      typedef CombinedDofConversionUtility<
        typename ContainedDiscreteFunctionSpaceType :: BlockMapperType, N, policy >   DofConversionType;

      explicit CombinedSpace( GridPartType &gridPart,
          const InterfaceType commInterface = InteriorBorder_All_Interface,
          const CommunicationDirection commDirection = ForwardCommunication )
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

    // DefaultLocalRestrictProlong also simply derives from the base class for the PointBased approach
    template< class ContainedSpace, int N >
    class DefaultLocalRestrictProlong< CombinedSpace< ContainedSpace, N, PointBased > >
    : public DefaultLocalRestrictProlong< typename CombinedSpace< ContainedSpace, N, PointBased >::BaseType > 
    {
      typedef CombinedSpace< ContainedSpace, N, PointBased > SpaceType;
      typedef DefaultLocalRestrictProlong< typename SpaceType::BaseType > BaseType;
    public:
      DefaultLocalRestrictProlong ( const SpaceType& space )
      : BaseType( space )
      {}
    };

    /** @} **/

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_COMBINEDSPACE_HH
