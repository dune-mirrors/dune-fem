#ifndef DUNE_FEM_SPACE_BREZZIDOUGLASMARINI_SPACE_HH
#define DUNE_FEM_SPACE_BREZZIDOUGLASMARINI_SPACE_HH

// -------------------------------------------------------------
// !!! BDMDiscreteFunctionSpace requires dune-localfunctions !!!
// -------------------------------------------------------------
#if HAVE_DUNE_LOCALFUNCTIONS

// C++ includes
#include <cassert>
#include <vector>

// dune-geometry types
#include <dune/geometry/type.hh>

// dune-localfunctions includes
#include <dune/localfunctions/rannacherturek.hh>

// dune-fem includes
#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/space/basisfunctionset/default.hh>
#include <dune/fem/space/basisfunctionset/piolatransformation.hh>
#include <dune/fem/space/basisfunctionset/transformed.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>
#include <dune/fem/space/shapefunctionset/localfunctions.hh>
#include <dune/fem/space/shapefunctionset/proxy.hh>
#include <dune/fem/space/shapefunctionset/selectcaching.hh>
#include <dune/fem/space/shapefunctionset/tuple.hh>

// local includes
#include <dune/fem/common/tuplehelper.hh>
#include <dune/fem/space/brezzidouglasmarini/capabilities.hh>
#include <dune/fem/space/brezzidouglasmarini/declaration.hh>
#include <dune/fem/space/brezzidouglasmarini/dofmappercode.hh>
#include <dune/fem/space/brezzidouglasmarini/localfiniteelement.hh>
#include <dune/fem/space/brezzidouglasmarini/localinterpolation.hh>

/**
   @file
   @author Tobias Malkmus
   @brief  Provides space based on Brezzi-Douglas-Marini finite element.
 */


namespace Dune
{

  namespace Fem
  {

    // Forward decleration
    template< class ShapeFunctionSet, int blockSize >
    struct MakeTupledShapeFunctionSet;

    // BDMDiscreteFunctionSpaceTraits
    // ------------------------------

    template< class FunctionSpace, class GridPart, int order, template< class > class Storage >
    struct BDMDiscreteFunctionSpaceTraits
    {
      static_assert( Dune::Fem::GridPartCapabilities::hasSingleGeometryType< GridPart >::v,
                     "GridPart has more than one geometry type." );

      static const unsigned int topologyId = Dune::Fem::GridPartCapabilities::hasSingleGeometryType< GridPart >::topologyId;

      typedef BDMDiscreteFunctionSpace< FunctionSpace, GridPart, order, Storage > DiscreteFunctionSpaceType;

      typedef FunctionSpace FunctionSpaceType;
      typedef GridPart GridPartType;

      static const int codimension = 0;

    private:
      typedef typename GridPartType::template Codim< codimension >::EntityType EntityType;

      typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;

      static const int dimLocal = GridPartType::dimension;

      static_assert( RangeType::dimension == dimLocal, "Range of requested BDM Space is not supported." );

    public:
      typedef BDMLocalFiniteElement< topologyId, DomainFieldType, RangeFieldType, dimLocal, order > LocalFiniteElementType;
      typedef typename LocalFiniteElementType::Traits::LocalBasisType LocalBasisType;
      typedef typename LocalFiniteElementType::Traits::LocalCoefficientsType LocalCoefficientsType;
      typedef typename LocalFiniteElementType::Traits::LocalInterpolationType LocalInterpolationType;

      typedef BDMBlockMapperFactory< GridPartType, LocalCoefficientsType > BlockMapperFactoryType;
      typedef typename BlockMapperFactoryType::BlockMapperType BlockMapperType;

      static const int localBlockSize = 1;

      typedef LocalFunctionsShapeFunctionSet< LocalBasisType > LocalFunctionsShapeFunctionSetType;
      typedef SelectCachingShapeFunctionSet< LocalFunctionsShapeFunctionSetType, Storage > CachedShapeFunctionSetType;

      typedef std::pair< GeometryType, unsigned char > ShapeFunctionSetSingletonKey;
      struct ShapeFunctionSetFactory
      {
        static CachedShapeFunctionSetType *createObject ( const ShapeFunctionSetSingletonKey &key )
        {
          return new CachedShapeFunctionSetType( key.first, LocalFunctionsShapeFunctionSetType( LocalBasisType( key.second ) ) );
        }

        static void deleteObject ( CachedShapeFunctionSetType *object ) { delete object; }
      };

      typedef ShapeFunctionSetFactory ShapeFunctionSetFactoryType;
      typedef ShapeFunctionSetProxy< CachedShapeFunctionSetType > ShapeFunctionSetProxyType;

      typedef typename SameTypeTuple< ShapeFunctionSetProxyType, localBlockSize, TupleShapeFunctionSet >::Type ShapeFunctionSetType;
      typedef PiolaTransformation< typename EntityType::Geometry, RangeType::dimension > TransformationType;
      typedef TransformedBasisFunctionSet< EntityType, ShapeFunctionSetType, TransformationType > BasisFunctionSetType;

      template< class DiscreteFunction, class Operation = DFCommunicationOperation::Add >
      struct CommDataHandle
      {
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
        typedef Operation OperationType;
      };

      static ShapeFunctionSetType getShapeFunctionSet ( CachedShapeFunctionSetType *cache )
      {
        return SameTypeTuple< ShapeFunctionSetProxyType, localBlockSize, TupleShapeFunctionSet >::construct( cache );
      }
    };



    // BDMDiscreteFunctionSpace
    // ------------------------

    /** \brief Brezzi-Douglas-Marini Space
     *  \class BDMDiscreteFunctionSpace
     *  \ingroup DiscreteFunctionSpace
     *
     *  \note The BDMDiscreteFunctionSpace depends on
     *        dune-localfunctions (see http://www.dune-project.org).
     *
     *  \todo please doc me
     */
    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage = CachingStorage >
    struct BDMDiscreteFunctionSpace
      : public DiscreteFunctionSpaceDefault< BDMDiscreteFunctionSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > >
    {
      typedef BDMDiscreteFunctionSpace< FunctionSpace, GridPart, polOrder, Storage > ThisType;
      typedef DiscreteFunctionSpaceDefault< BDMDiscreteFunctionSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > > BaseType;

      static const InterfaceType defaultInterface = InteriorBorder_All_Interface;
      static const CommunicationDirection defaultDirection =  ForwardCommunication;

    public:
      static const int polynomialOrder = polOrder;

      typedef typename BaseType::Traits Traits;

      typedef typename BaseType::FunctionSpaceType FunctionSpaceType;

      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::EntityType EntityType;
      typedef typename BaseType::IntersectionType IntersectionType;

      typedef typename BaseType::Traits::ShapeFunctionSetType ShapeFunctionSetType;
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;

      typedef typename BaseType::BlockMapperType BlockMapperType;

      typedef BDMLocalInterpolation< BasisFunctionSetType, typename Traits::LocalInterpolationType > InterpolationType;

    private:
      typedef typename Traits::CachedShapeFunctionSetType CachedShapeFunctionSetType;
      typedef typename Traits::ShapeFunctionSetSingletonKey ShapeFunctionSetSingletonKey;
      typedef SingletonList< ShapeFunctionSetSingletonKey, CachedShapeFunctionSetType, typename Traits::ShapeFunctionSetFactoryType > ShapeFunctionSetProviderType;

      typedef BDMBlockMapperSingletonKey< GridPartType > BlockMapperSingletonKeyType;
      typedef typename Traits::BlockMapperFactoryType BlockMapperFactoryType;
      typedef SingletonList< BlockMapperSingletonKeyType, BlockMapperType, BlockMapperFactoryType > BlockMapperProviderType;

      typedef typename Traits::LocalInterpolationType LocalInterpolationType;
      static const int numOrientations = Traits::LocalFiniteElementType::numOrientations;

    public:
      using BaseType::order;

      explicit BDMDiscreteFunctionSpace ( GridPartType &gridPart,
                                          const InterfaceType commInterface = defaultInterface,
                                          const CommunicationDirection commDirection = defaultDirection )
        : BaseType( gridPart, commInterface, commDirection ),
        cachedShapeFunctionSet_( numOrientations, nullptr ),
        blockMapper_( nullptr )
      {
        // create scalar shape function set
        GeometryType type = GeometryType( Dune::Fem::GridPartCapabilities::hasSingleGeometryType< GridPartType >::topologyId, GridPartType::dimension );

        for( unsigned char i = 0; i < numOrientations; ++i )
        {
          ShapeFunctionSetSingletonKey key( type, i );
          cachedShapeFunctionSet_[ i ] = &( ShapeFunctionSetProviderType::getObject( key ) );
        }

        {
          // create block mapper
          BlockMapperSingletonKeyType key( gridPart );
          blockMapper_ = &( BlockMapperProviderType::getObject( key ) );
        }
      }

      ~BDMDiscreteFunctionSpace ()
      {
        if( blockMapper_ )
          BlockMapperProviderType::removeObject( *blockMapper_ );
        blockMapper_ = nullptr;

        for( unsigned char i = 0; i < numOrientations; ++i )
        {
          if( cachedShapeFunctionSet_[ i ] )
            ShapeFunctionSetProviderType::removeObject( *cachedShapeFunctionSet_[ i ] );
          cachedShapeFunctionSet_[ i ] = nullptr;
        }
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
        unsigned char orient = orientation( entity );
        assert( cachedShapeFunctionSet_[ orient ] );
        return ShapeFunctionSetType( Traits::getShapeFunctionSet( cachedShapeFunctionSet_[ orient ] ) );
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
        unsigned char orient = orientation( entity );
        assert( cachedShapeFunctionSet_[ orient ] );
        return InterpolationType( BasisFunctionSetType( entity, ShapeFunctionSetType( Traits::getShapeFunctionSet( cachedShapeFunctionSet_[ orient ] ) ) ), LocalInterpolationType( orient ) );
      }

      BDMDiscreteFunctionSpace ( const ThisType & ) = delete;
      ThisType &operator= ( const ThisType & ) = delete;

      using BaseType::gridPart;

    protected:
      template< class Entity >
      unsigned char orientation ( const Entity &entity ) const
      {
        unsigned char ret = 0;
        auto &idxSet = gridPart().indexSet();
        for( auto intersection : intersections( gridPart(), entity ) )
          if( intersection.neighbor() &&
              ( idxSet.index( entity ) < idxSet.index( intersection.outside() ) ) )
            ret |= 1 << intersection.indexInInside();
        return ret;
      }

    private:
      std::vector< CachedShapeFunctionSetType * > cachedShapeFunctionSet_;
      BlockMapperType *blockMapper_;
    };

  } // namespace Fem

} // end namespace Dune

#endif // #if HAVE_DUNE_LOCALFUNCTIONS

#endif // #ifndef DUNE_FEM_SPACE_BREZZIDOUGLASMARINI_SPACE_HH
