#ifndef DUNE_FEM_SPACE_RANNACHERTUREK_SPACE_HH
#define DUNE_FEM_SPACE_RANNACHERTUREK_SPACE_HH

// ------------------------------------------------------------------------
// !!! RannacherTurekDiscreteFunctionSpace requires dune-localfunctions !!!
// ------------------------------------------------------------------------
#if HAVE_DUNE_LOCALFUNCTIONS

// C++ includes
#include <cassert>
#include <vector>

// dune-geometry types
#include <dune/geometry/type.hh>

// dune-localfunctions includes
#include <dune/localfunctions/rannacherturek.hh>

// dune-fem includes
#include <dune/fem/common/hybrid.hh>
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
#include <dune/fem/space/rannacherturek/capabilities.hh>
#include <dune/fem/space/rannacherturek/declaration.hh>
#include <dune/fem/space/rannacherturek/dofmappercode.hh>
#include <dune/fem/space/rannacherturek/localinterpolation.hh>

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

    template< class FunctionSpace, class GridPart, class Storage >
    struct RannacherTurekDiscreteFunctionSpaceTraits
    {
      static_assert( Dune::Fem::GridPartCapabilities::hasSingleGeometryType< GridPart >::v,
                          "GridPart has more than one geometry type." );

      typedef RannacherTurekDiscreteFunctionSpace< FunctionSpace, GridPart, Storage > DiscreteFunctionSpaceType;

      typedef GridPart GridPartType;
      typedef GridFunctionSpace< GridPartType, FunctionSpace > FunctionSpaceType;

      static const int codimension = 0;

    private:
      typedef typename GridPartType::template Codim< codimension >::EntityType EntityType;

      typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
      static const int dimLocal = GridPartType::dimension;

    public:
      typedef Dune::RannacherTurekLocalFiniteElement< DomainFieldType, RangeFieldType, dimLocal > LocalFiniteElementType;
      typedef typename LocalFiniteElementType::Traits::LocalBasisType LocalBasisType;
      typedef typename LocalFiniteElementType::Traits::LocalCoefficientsType LocalCoefficientsType;
      typedef typename LocalFiniteElementType::Traits::LocalInterpolationType LocalInterpolationType;

      typedef RannacherTurekBlockMapperFactory< GridPartType, LocalCoefficientsType > BlockMapperFactoryType;
      typedef typename BlockMapperFactoryType::BlockMapperType BlockMapperType;

      typedef Hybrid::IndexRange< int, FunctionSpaceType::dimRange > LocalBlockIndices;

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

    /** \brief Rannacher-Turek Space
     *  \class RannacherTurekDiscreteFunctionSpace
     *  \ingroup DiscreteFunctionSpace
     *
     *  \note The RannacherTurekDiscreteFunctionSpace depends on
     *        dune-localfunctions (see http://www.dune-project.org).
     *
     *  \todo please doc me
     */
    template< class FunctionSpace, class GridPart, class Storage = CachingStorage >
    struct RannacherTurekDiscreteFunctionSpace
    : public DiscreteFunctionSpaceDefault< RannacherTurekDiscreteFunctionSpaceTraits< FunctionSpace, GridPart, Storage > >
    {
      typedef RannacherTurekDiscreteFunctionSpace< FunctionSpace, GridPart, Storage > ThisType;
      typedef DiscreteFunctionSpaceDefault< RannacherTurekDiscreteFunctionSpaceTraits< FunctionSpace, GridPart, Storage > > BaseType;

      static const InterfaceType defaultInterface = InteriorBorder_All_Interface;
      static const CommunicationDirection defaultDirection =  ForwardCommunication;

    public:
      // must be 2 here since it contains x_i^2
      static const int polynomialOrder = 2;

      typedef typename BaseType::Traits Traits;

      typedef typename BaseType::FunctionSpaceType FunctionSpaceType;

      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::EntityType EntityType;
      typedef typename BaseType::IntersectionType IntersectionType;

      typedef typename BaseType::Traits::ShapeFunctionSetType ShapeFunctionSetType;
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;

      typedef typename BaseType::BlockMapperType BlockMapperType;

      typedef RannacherTurekLocalInterpolation< BasisFunctionSetType, typename Traits::LocalInterpolationType > InterpolationType;

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
        blockMapper_( nullptr )
      {
        // create scalar shape function set
        GeometryType type = GeometryType( Dune::Fem::GridPartCapabilities::hasSingleGeometryType< GridPartType >::topologyId, GridPartType::dimension );
        scalarShapeFunctionSet_ = &( ScalarShapeFunctionSetProviderType::getObject( type ) );

        // create block mapper
        BlockMapperSingletonKeyType key( gridPart );
        blockMapper_ = &( BlockMapperProviderType::getObject( key ) );
      }

      ~RannacherTurekDiscreteFunctionSpace ()
      {
        if( blockMapper_ )
          BlockMapperProviderType::removeObject( *blockMapper_ );
        blockMapper_ = nullptr;

        if( scalarShapeFunctionSet_ )
          ScalarShapeFunctionSetProviderType::removeObject( *scalarShapeFunctionSet_ );
        scalarShapeFunctionSet_ = nullptr;
      }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::type */
      DFSpaceIdentifier type () const { return RannacherTurekSpace_id; }

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

      /** \brief return local interpolation
       *
       *  \param[in]  entity  grid part entity
       */
      InterpolationType interpolation ( const EntityType &entity ) const
      {
        return InterpolationType( basisFunctionSet( entity ) );
      }

      RannacherTurekDiscreteFunctionSpace ( const ThisType & ) = delete;
      ThisType &operator= ( const ThisType & ) = delete;

    private:
      ScalarShapeFunctionSetType *scalarShapeFunctionSet_;
      BlockMapperType *blockMapper_;
    };


    // DifferentDiscreteFunctionSpace
    // ------------------------------

    template< class FunctionSpace, class GridPart, class Storage, class NewFunctionSpace >
    struct DifferentDiscreteFunctionSpace< RannacherTurekDiscreteFunctionSpace< FunctionSpace, GridPart, Storage >, NewFunctionSpace >
    {
      typedef RannacherTurekDiscreteFunctionSpace< NewFunctionSpace, GridPart, Storage > Type;
    };

  } // namespace Fem

} // end namespace Dune

#endif // #if HAVE_DUNE_LOCALFUNCTIONS

#endif // #ifndef DUNE_FEM_SPACE_RANNACHERTUREK_SPACE_HH
