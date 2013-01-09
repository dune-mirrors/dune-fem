#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_DEFAULT_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_DEFAULT_HH

// dune-common includes
#include <dune/common/static_assert.hh>

// dune-fem includes
#include <dune/fem/misc/bartonnackmaninterface.hh>
#include <dune/fem/space/basisfunctionset/default.hh>
#include <dune/fem/space/common/basesetlocalkeystorage.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/mapper/codimensionmapper.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>
#include <dune/fem/space/shapefunctionset/proxy.hh>
#include <dune/fem/space/shapefunctionset/vectorial.hh>
#include <dune/fem/storage/singletonlist.hh>
#include <dune/fem/version.hh>

// local includes
#include "declaration.hh"
#include "localinterpolation.hh"
#include "localrestrictprolong.hh"


namespace Dune
{

  namespace Fem
  {

    // DiscontinuousGalerkinSpaceDefault
    // ---------------------------------

    /* 
     * Default implementation for discrete Discontinuous Galerkin spaces.
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

      typedef typename Traits::ScalarShapeFunctionSetType ScalarShapeFunctionSetType;
      typedef SingletonList< const GeometryType, ScalarShapeFunctionSetType, typename Traits::ScalarShapeFunctionSetFactoryType > SingletonProviderType;
      typedef BaseSetLocalKeyStorage< ScalarShapeFunctionSetType > ScalarShapeFunctionSetStorageType;
      
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
        blockMapper_( &BlockMapperProviderType::getObject( gridPart ) ),
        mapper_( blockMapper() )
      {
        // get geometry types
        std::vector< GeometryType > geomTypes = AllGeomTypes< IndexSetType, GridType >( gridPart.indexSet()).geomTypes( Traits::codimension );
        
        // store shape function sets per type
        const typename std::vector< GeometryType >::const_iterator end = geomTypes.end();
        for( typename std::vector< GeometryType >::const_iterator it = geomTypes.begin(); it != end; ++it )    
        {
          const GeometryType &type = *it;
          scalarShapeFunctionSets_.template insert< SingletonProviderType >( type );
        }
      }

      ~DiscontinuousGalerkinSpaceDefault ()
      {
        BlockMapperProviderType::removeObject( blockMapper() );
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
        return *blockMapper_;
      }

      ///////////////////////////
      // Non-interface methods //
      ///////////////////////////

      /** \brief local interpolation using discontinuous L2-projection
       *
       *  \param[in]  localFunction  local function to interpolate
       *  \param[in]  dofs           local degrees of freedom of the interpolation
       */
      template< class LocalFunction, class LocalDofVector >
      void interpolate ( const LocalFunction &localFunction, LocalDofVector &dofs ) const
      {
        typedef DiscontinuousGalerkinLocalInterpolation< typename BaseType::DiscreteFunctionSpaceType > LocalInterpolationType;
        LocalInterpolationType interpolation( asImp() );
        interpolation( localFunction, dofs );
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
        return shapeFunctionSet( entity.type() );
      }

      /** \brief return shape unique function set for geometry type 
       *
       * \param[in]  type  geometry type (must be a cube) for which 
       *                   shape function set is requested
       *
       * \returns  ShapeFunctionSetType  shape function set                     
       */
      ShapeFunctionSetType shapeFunctionSet ( const GeometryType &type ) const
      {
        return ShapeFunctionSetType( &scalarShapeFunctionSets_[ type ] );
      }

    private:
      BlockMapperType *blockMapper_;
      mutable MapperType mapper_;
      ScalarShapeFunctionSetStorageType scalarShapeFunctionSets_;
    };



    // DiscontinuousGalerkinSpaceDefaultTraits
    // ---------------------------------------

    /* 
     * Traits class for default Discontinuous Galerkin space implementation
     */
    template< class Traits >
    struct DiscontinuousGalerkinSpaceDefaultTraits
    {
      typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

      typedef typename Traits::FunctionSpaceType FunctionSpaceType;
      typedef typename Traits::GridPartType GridPartType;

      static const int codimension = Traits::codimension;
      static const int polynomialOrder = Traits::polynomialOrder;
      dune_static_assert( (polynomialOrder >= 0), "Negative polynomial order." );

      typedef typename GridPartType::template Codim< codimension >::EntityType EntityType;

      typedef typename Traits::ScalarShapeFunctionSetType ScalarShapeFunctionSetType;
      typedef typename Traits::ScalarShapeFunctionSetFactoryType ScalarShapeFunctionSetFactoryType;

      typedef ShapeFunctionSetProxy< ScalarShapeFunctionSetType > ScalarShapeFunctionSetProxyType;
      typedef VectorialShapeFunctionSet< ScalarShapeFunctionSetProxyType, typename FunctionSpaceType::RangeType > ShapeFunctionSetType;

      typedef Dune::Fem::DefaultBasisFunctionSet< EntityType, ShapeFunctionSetType > BasisFunctionSetType;

      typedef CodimensionMapper< GridPartType, codimension > BlockMapperType;
      static const int localBlockSize = Traits::localBlockSize;
      typedef NonBlockMapper< BlockMapperType, localBlockSize > MapperType;

      template <class DiscreteFunction, class Operation = DFCommunicationOperation::Copy >
      struct CommDataHandle
      {
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
        typedef Operation OperationType;
      };
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_DEFAULT_HH
