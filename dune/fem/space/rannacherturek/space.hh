#ifndef DUNE_FEM_SPACE_RANNACHERTUREK_HH
#define DUNE_FEM_SPACE_RANNACHERTUREK_HH

// C++ includes
#include <cassert>
#include <vector>

// dune-common includes
#include <dune/common/nullptr.hh>
#include <dune/common/static_assert.hh>

// dune-geometry types
#include <dune/geometry/type.hh>

// dune-fem includes
#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/space/basisfunctionset/default.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>
#include <dune/fem/space/shapefunctionset/localfunctions.hh>
#include <dune/fem/space/shapefunctionset/proxy.hh>
#include <dune/fem/space/shapefunctionset/selectcaching.hh>
#include <dune/fem/space/shapefunctionset/vectorial.hh>

// local includes
#include "declaration.hh"
#include "dofmappercode.hh"
#include "localfiniteelement.hh"
#include "localinterpolation.hh"
#include "localrestrictprolong.hh"

/**
  @file
  @author Christoph Gersbacher
  @brief  Provides space based on Rannacher-Turek finite element.
*/


namespace Dune
{

  namespace Fem
  {

    // RannacherTurekDiscreteFunctionSpaceTraits
    // -----------------------------------------

    template< class FunctionSpace, class GridPart, template< class > class Storage >
    struct RannacherTurekDiscreteFunctionSpaceTraits
    {
      dune_static_assert( Dune::Fem::GridPartCapabilities::hasSingleGeometryType< GridPart >::v,
                          "GridPart has more than one geometry type." );

      typedef RannacherTurekDiscreteFunctionSpace< FunctionSpace, GridPart, Storage > DiscreteFunctionSpaceType;

      typedef FunctionSpace FunctionSpaceType;
      typedef GridPart GridPartType;

      static const int codimension = 0;

    private:
      typedef typename GridPartType::template Codim< codimension >::EntityType EntityType;

      typedef typename FunctionSpaceType::ScalarFunctionSpaceType ScalarFunctionSpaceType;
      static const int dimLocal = GridPartType::dimension;
      typedef typename ToLocalFunctionSpace< ScalarFunctionSpaceType, dimLocal >::Type ScalarShapeFunctionSpaceType;

    public:
      typedef RannacherTurekLocalFiniteElement< ScalarShapeFunctionSpaceType > LocalFiniteElementType;
      typedef typename LocalFiniteElementType::LocalBasisType LocalBasisType;
      typedef typename LocalFiniteElementType::LocalCoefficientsType LocalCoefficientsType;
      typedef typename LocalFiniteElementType::LocalInterpolationType LocalInterpolationType;

      typedef RannacherTurekBlockMapperFactory< GridPartType, LocalCoefficientsType > BlockMapperFactoryType;
      typedef typename BlockMapperFactoryType::BlockMapperType BlockMapperType;

      static const int localBlockSize = FunctionSpaceType::dimRange;
      typedef NonBlockMapper< BlockMapperType, localBlockSize > MapperType;

      typedef LocalFunctionsShapeFunctionSet< LocalBasisType > LocalFunctionsShapeFunctionSetType;
      typedef SelectCachingShapeFunctionSet< LocalFunctionsShapeFunctionSetType, Storage > ScalarShapeFunctionSetType;

      struct ScalarShapeFunctionSetFactory
      {
        static ScalarShapeFunctionSetType *createObject ( const GeometryType &type )
        {
          assert( type.isCube() );
          return new ScalarShapeFunctionSetType( type, LocalFunctionsShapeFunctionSetType( LocalBasisType() ) );
        }
        
        static void deleteObject ( ScalarShapeFunctionSetType *object ) { delete object; }
      };

      typedef ScalarShapeFunctionSetFactory ScalarShapeFunctionSetFactoryType;

      typedef ShapeFunctionSetProxy< ScalarShapeFunctionSetType > ScalarShapeFunctionSetProxyType;
      typedef VectorialShapeFunctionSet< ScalarShapeFunctionSetProxyType, typename FunctionSpaceType::RangeType > ShapeFunctionSetType;

      typedef DefaultBasisFunctionSet< EntityType, ShapeFunctionSetType > BasisFunctionSetType;

      template< class DiscreteFunction, class Operation = DFCommunicationOperation::Add >
      struct CommDataHandle
      {
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
        typedef Operation OperationType;
      };
    };



    // RannacherTurekDiscreteFunctionSpace
    // -----------------------------------

    template< class FunctionSpace, class GridPart, template< class > class Storage = CachingStorage >
    struct RannacherTurekDiscreteFunctionSpace
    : public DiscreteFunctionSpaceDefault< RannacherTurekDiscreteFunctionSpaceTraits< FunctionSpace, GridPart, Storage > >
    {
      typedef RannacherTurekDiscreteFunctionSpace< FunctionSpace, GridPart, Storage > ThisType;
      typedef DiscreteFunctionSpaceDefault< RannacherTurekDiscreteFunctionSpaceTraits< FunctionSpace, GridPart, Storage > > BaseType;

      static const InterfaceType defaultInterface = InteriorBorder_All_Interface;
      static const CommunicationDirection defaultDirection =  ForwardCommunication;

    public:
      static const int polynomialOrder = 1;

      typedef typename BaseType::Traits Traits;

      typedef typename BaseType::FunctionSpaceType FunctionSpaceType;

      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::EntityType EntityType;
      typedef typename BaseType::IntersectionType IntersectionType;
     
      typedef typename BaseType::Traits::ShapeFunctionSetType ShapeFunctionSetType;
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;
      
      typedef typename BaseType::MapperType MapperType;
      typedef typename BaseType::BlockMapperType BlockMapperType;

    private:
      typedef typename Traits::ScalarShapeFunctionSetType ScalarShapeFunctionSetType;
      typedef SingletonList< GeometryType, ScalarShapeFunctionSetType, typename Traits::ScalarShapeFunctionSetFactoryType > ScalarShapeFunctionSetProviderType;

      typedef RannacherTurekBlockMapperSingletonKey< GridPartType > BlockMapperSingletonKeyType;
      typedef typename Traits::BlockMapperFactoryType BlockMapperFactoryType;
      typedef SingletonList< BlockMapperSingletonKeyType, BlockMapperType, BlockMapperFactoryType > BlockMapperProviderType;

    public:
      using BaseType::order;

      explicit RannacherTurekDiscreteFunctionSpace ( GridPartType &gridPart,
                                     const InterfaceType commInterface = defaultInterface,
                                     const CommunicationDirection commDirection = defaultDirection )
      : BaseType( gridPart, commInterface, commDirection ),
        scalarShapeFunctionSet_( nullptr ),
        blockMapper_( nullptr ),
        mapper_( nullptr )
      {
        // create scalar shape function set
        GeometryType type = GeometryType( Dune::Fem::GridPartCapabilities::hasSingleGeometryType< GridPartType >::topologyId, GridPartType::dimension );
        scalarShapeFunctionSet_ = &( ScalarShapeFunctionSetProviderType::getObject( type ) );

        // create block mapper
        BlockMapperSingletonKeyType key( gridPart );
        blockMapper_ = &( BlockMapperProviderType::getObject( key ) );

        // create mapper
        mapper_ = new MapperType( blockMapper() );
      }

      ~RannacherTurekDiscreteFunctionSpace ()
      {
        if( mapper_ )
          delete mapper_;
        mapper_ = nullptr;

        if( blockMapper_ )
          BlockMapperProviderType::removeObject( *blockMapper_ );
        blockMapper_ = nullptr;

        if( scalarShapeFunctionSet_ )
          ScalarShapeFunctionSetProviderType::removeObject( *scalarShapeFunctionSet_ );
        scalarShapeFunctionSet_ = nullptr;
      }
      
      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::type */
      DFSpaceIdentifier type () const { return GenericSpace_id; }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::basisFunctionSet */
      BasisFunctionSetType basisFunctionSet ( const EntityType &entity ) const
      {
        return BasisFunctionSetType( entity, shapeFunctionSet( entity ) );
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
        assert( scalarShapeFunctionSet_);
        return ShapeFunctionSetType( scalarShapeFunctionSet_ );
      }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      bool continuous () const { return false; }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      bool continuous ( const IntersectionType &intersection ) const { return false; }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order */
      int order () const { return polynomialOrder; }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::multipleGeometryTypes */
      inline bool multipleGeometryTypes () const { return false; }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::blockMapper */
      BlockMapperType &blockMapper () const
      {
        assert( blockMapper_ );
        return *blockMapper_;
      }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::mapper */
      DUNE_VERSION_DEPRECATED(1,4,remove)
      MapperType &mapper () const
      {
        assert( mapper_ );
        return *mapper_;
      }      

      /** \brief interpolate a function locally
       *
       *  \param[in]  f     local function to interpolate
       *  \param[out] dofs  local degrees of freedom of the interpolion
       */
      template< class LocalFunction, class LocalDofVector >
      static void interpolate ( const LocalFunction &f, LocalDofVector &dofs )
      {
        typedef typename Traits::LocalInterpolationType LocalInterpolationType;
        VectorialLocalInterpolation< LocalInterpolationType, typename FunctionSpaceType::RangeType > interpolation;
        interpolation( f, dofs );
      }

    private:
      RannacherTurekDiscreteFunctionSpace ( const ThisType & );
      ThisType &operator= ( const ThisType & );

      ScalarShapeFunctionSetType *scalarShapeFunctionSet_;
      BlockMapperType *blockMapper_;
      MapperType *mapper_;
    };

  } // namespace Fem

} // end namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_RANNACHERTUREK_HH
