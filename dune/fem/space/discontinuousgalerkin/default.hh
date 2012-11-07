#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_DEFAULT_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_DEFAULT_HH

// dune-common includes
#include <dune/common/static_assert.hh>

// dune-geometry includes
#include <dune/geometry/type.hh>

// dune-fem includes
#include <dune/fem/misc/bartonnackmaninterface.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/mapper/codimensionmapper.hh>
#include <dune/fem/storage/singletonlist.hh>
#include <dune/fem/version.hh>


namespace Dune
{

  namespace Fem
  {

    // DiscontinuousGalerkinSpaceTraitsBase
    // ------------------------------------

    /* 
     * common base traits class for all Discontinuous Galerkin spaces
     */
    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
    struct DiscontinuousGalerkinSpaceTraitsBase
    {
      typedef FunctionSpace FunctionSpaceType;
      typedef GridPart GridPartType;

      static const int dimRange = FunctionSpaceType::dimRange;
      static const int dimLocal = GridPartType::dimension;
      dune_static_assert( (GridPartType::dimensionworld <= 3), "Use Legendre spaces for higher spatial dimensions." );
      
      static const int codimension = 0;
      static const int polynomialOrder = polOrder;
      dune_static_assert( (polOrder >= 0), "Negative polynomial order." );

      typedef CodimensionMapper< GridPartType, codimension > BlockMapperType;

      template <class DiscreteFunction, class Operation = DFCommunicationOperation::Copy >
      struct CommDataHandle
      {
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
        typedef Operation OperationType;
      };
    };



    // DiscontinuousGalerkinSpaceDefault
    // ---------------------------------

    /* 
     * Default implementation for discrete Discontinuous Galerkin spaces.
     *
     * Derived classes are expected to implement the following methods:
\code
  // return shape function set for given entity
  ShapeFunctionSetType shapeFunctionSet ( const EntityType &entity ) const;
  // return mapper
  MapperType &mapper () const;
\endcode
     */
    template< class Traits >
    class DiscontinuousGalerkinSpaceDefault
    : public DiscreteFunctionSpaceDefault< Traits >
    {
      typedef DiscontinuousGalerkinSpaceDefault< Traits > ThisType;
      typedef DiscreteFunctionSpaceDefault< Traits > BaseType;

    public:
      static const int codimension = Traits::codimension;
      static const int polynomialOrder = Traits::polynomialOrder;

      typedef typename BaseType::FunctionSpaceType FunctionSpaceType;

      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::GridType GridType;
      typedef typename BaseType::IndexSetType IndexSetType;
      typedef typename BaseType::IteratorType IteratorType;
      typedef typename IteratorType::Entity EntityType;
      typedef typename BaseType::IntersectionType IntersectionType;

      typedef typename Traits::ShapeFunctionSetType ShapeFunctionSetType;
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;

      typedef typename BaseType::MapperType MapperType;
      typedef typename BaseType::BlockMapperType BlockMapperType;

    protected:
      using BaseType::asImp;

      typedef CodimensionMapperSingletonFactory< GridPartType, codimension > BlockMapperSingletonFactoryType;
      typedef SingletonList< typename BlockMapperSingletonFactoryType::Key,
                             BlockMapperType, BlockMapperSingletonFactoryType 
                           > BlockMapperProviderType;

      static const InterfaceType defaultInterface = InteriorBorder_All_Interface;
      static const CommunicationDirection defaultDirection =  ForwardCommunication;

    public:
      using BaseType::order;

      DiscontinuousGalerkinSpaceDefault ( GridPartType &gridPart, 
                                          const InterfaceType commInterface,
                                          const CommunicationDirection commDirection )
      : BaseType( gridPart, commInterface, commDirection ),
        blockMapper_( BlockMapperProviderType::getObject( gridPart ) )
      {}

      ~DiscontinuousGalerkinSpaceDefault ()
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
        return BasisFunctionSetType( entity, shapeFunctionSet( entity ) );
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
        CHECK_INTERFACE_IMPLEMENTATION( asImp().mapper() );
        return asImp().mapper();
      }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::blockMapper */
      BlockMapperType &blockMapper () const
      {
        return blockMapper_;
      }

      /** \brief return shape function set for given entity
       *
       * \param[in]  entity  entity (of codim 0) for which shape function set 
       *                     is requested
       *
       * \returns  ShapeFunctionSetType  shape function set                     
       */
      ShapeFunctionSetType shapeFunctionSet ( const EntityType &entity ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().shapeFunctionSet( entity ) );
        return asImp().shapeFunctionSet( entity );
      }

    private:
      mutable BlockMapperType &blockMapper_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_DEFAULT_HH
