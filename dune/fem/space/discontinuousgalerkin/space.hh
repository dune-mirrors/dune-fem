#ifndef DUNE_FEM_SPACE_DGSPACE_SPACE_HH
#define DUNE_FEM_SPACE_DGSPACE_SPACE_HH

// dune-common includes
#include <dune/common/static_assert.hh>

// dune-geometry includes
#include <dune/geometry/type.hh>

// dune-fem includes
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/mapper/codimensionmapper.hh>
#include <dune/fem/storage/singletonlist.hh>
#include <dune/fem/version.hh>

// local includes
#include "declaration.hh"


namespace Dune
{

  namespace Fem
  {

    // DiscontinuousGalerkinSpaceTraits
    // --------------------------------

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
    struct DiscontinuousGalerkinSpaceTraits
    {
      typedef FunctionSpace FunctionSpaceType;
      typedef GridPart GridPartType;
      static const int polynomialOrder = polOrder;

      // static const int localBlockSize = ;
      // typedef ... BlockMapperType;
      // typedef ... MapperType;

      static const int codimension = 0;

      // typedef ... BasisFunctionSetType;

      template< class DiscreteFunction, class Operation = DFCommunicationOperation::Copy >
      struct CommDataHandle
      {
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
        typedef Operation OperationType;
      };

      // typedef ... DiscreteFunctionSpaceType;
    };



    // DiscontinuousGalerkinSpace
    // --------------------------

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
    class DiscontinuousGalerkinSpace
    : public DiscreteFunctionSpaceDefault< DiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > >
    {
      typedef DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > ThisType;
      typedef DiscreteFunctionSpaceDefault< 
          DiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > 
        > BaseType;

    public:
      typedef typename BaseType::Traits Traits;

      static const int codimension = Traits::codimension;
      static const int polynomialOrder DUNE_DEPRECATED = 0;

      typedef typename BaseType::FunctionSpaceType FunctionSpaceType;

      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::GridType GridType;
      typedef typename BaseType::IndexSetType IndexSetType;
      typedef typename BaseType::IteratorType IteratorType;
      typedef typename IteratorType::Entity EntityType;
      typedef typename BaseType::IntersectionType IntersectionType;

      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;

      typedef typename BaseType::MapperType MapperType;
      typedef typename BaseType::BlockMapperType BlockMapperType;

    private:
      typedef CodimensionMapperSingletonFactory< GridPartType, codimension > BlockMapperSingletonFactoryType;
      typedef SingletonList< typename BlockMapperSingletonFactoryType::Key,
                             BlockMapperType, BlockMapperSingletonFactoryType 
                           > BlockMapperProviderType;

      static const InterfaceType defaultInterface = InteriorBorder_All_Interface;
      static const CommunicationDirection defaultDirection =  ForwardCommunication;

    public:
      using BaseType::order;

      DiscontinuousGalerkinSpace ( GridPartType &gridPart, 
                                   const std::vector< GeometryType > &types,
                                   const InterfaceType commInterface = defaultInterface,
                                   const CommunicationDirection commDirection = defaultDirection )
      : BaseType( gridPart, commInterface, commDirection ),
        blockMapper_( BlockMapperProviderType::getObject( gridPart ) ),
        mapper_( blockMapper_ )
      {}

      ~DiscontinuousGalerkinSpace ()
      {
        BlockMapperProviderType::removeObject( blockMapper_ );
      }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::type */
      DFSpaceIdentifier type () const
      {
        return DGSpace_id;
      }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::basisFunctionSet */
      BasisFunctionSetType basisFunctionSet ( const EntityType &entity ) const
      {
        // return BasisFunctionSetType( entity, &shapeFunctionSets_[ entity.type() ] );
      }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::contains */
      bool contains ( const int codimension ) const
      {
        return blockMapper_.contains( codimension );
      }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      bool continuous () const
      {
        return false;
      }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      bool continuous ( const IntersectionType &intersection ) const
      {
        return false;
      }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order */
      int order () const
      {
        return polynomialOrder;
      }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::mapper */
      DUNE_VERSION_DEPRECATED(1,4,remove)
      MapperType &mapper () const
      {
        return mapper_;
      }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::blockMapper */
      BlockMapperType &blockMapper () const
      {
        return blockMapper_;
      }

    private:
      BlockMapperType &blockMapper_;
      mutable MapperType mapper_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DGSPACE_SPACE_HH
