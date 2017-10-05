#ifndef DUNE_FEM_SPACE_LOCALFINITEELEMENT_DGSPACE_HH
#define DUNE_FEM_SPACE_LOCALFINITEELEMENT_DGSPACE_HH

#include <cassert>

#include <memory>
#include <utility>
#include <vector>

#include <dune/geometry/type.hh>

#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/space/basisfunctionset/default.hh>
#include <dune/fem/space/basisfunctionset/transformed.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/mapper/compile.hh>
#include <dune/fem/space/mapper/indexsetdofmapper.hh>
#include <dune/fem/space/shapefunctionset/proxy.hh>
#include <dune/fem/space/shapefunctionset/selectcaching.hh>
#include <dune/fem/space/shapefunctionset/vectorial.hh>

#include <dune/fem/space/localfiniteelement/shapefunctionset.hh>
#include <dune/fem/space/localfiniteelement/capabilities.hh>
#include <dune/fem/space/localfiniteelement/interpolation.hh>

namespace Dune
{

  namespace Fem
  {

    // DiscontinuousLocalFiniteElementSpaceTraits
    // ------------------------------------------

    template< class LFEMap, class FunctionSpace, template< class > class Storage >
    struct DiscontinuousLocalFiniteElementSpaceTraits
    {
      typedef DiscontinuousLocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > DiscreteFunctionSpaceType;

      typedef LFEMap LFEMapType;

      typedef typename LFEMapType::GridPartType GridPartType;
      typedef typename LFEMapType::LocalFiniteElementType LocalFiniteElementType;

      typedef GridFunctionSpace< GridPartType, FunctionSpace > FunctionSpaceType;

      static constexpr int codimension = 0;

      typedef Hybrid::IndexRange< int, 1 > LocalBlockIndices;

    private:
      typedef typename GridPartType::template Codim< codimension >::EntityType EntityType;

    public:
      typedef Dune::Fem::IndexSetDofMapper< GridPartType > BlockMapperType;

      typedef LocalFunctionsShapeFunctionSet< typename LocalFiniteElementType::Traits::LocalBasisType > LocalFunctionsShapeFunctionSetType;
      typedef SelectCachingShapeFunctionSet< LocalFunctionsShapeFunctionSetType, Storage > StoredShapeFunctionSetType;

      typedef ShapeFunctionSetProxy< StoredShapeFunctionSetType > ShapeFunctionSetType;
//      typedef VectorialShapeFunctionSet< ShapeFunctionSetProxy< StoredShapeFunctionSetType >, typename FunctionSpaceType::RangeType > ShapeFunctionSetType;

    private:
      template< class LFEM >
      static TransformedBasisFunctionSet< EntityType, ShapeFunctionSetType, typename LFEM::TransformationType > basisFunctionSet ( const LFEM & );

      static DefaultBasisFunctionSet< EntityType, ShapeFunctionSetType > basisFunctionSet ( ... );

    public:
      typedef decltype( basisFunctionSet( std::declval< const LFEMapType & >() ) ) BasisFunctionSetType;

      template< class DiscreteFunction, class Operation = DFCommunicationOperation::Copy >
      struct CommDataHandle
      {
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
        typedef Operation OperationType;
      };
    };



    // DiscontinuousLocalFiniteElementSpace
    // ------------------------------------

    /** \brief Rannacher-Turek Space
     *  \class DiscontinuousLocalFiniteElementSpace
     *  \ingroup DiscreteFunctionSpace
     *
     *  \note The DiscontinuousLocalFiniteElementSpace depends on
     *        dune-localfunctions (see http://www.dune-project.org).
     *
     *  \todo please doc me
     **/
    template< class LFEMap, class FunctionSpace, template< class > class Storage >
    class DiscontinuousLocalFiniteElementSpace
      : public DiscreteFunctionSpaceDefault< DiscontinuousLocalFiniteElementSpaceTraits< LFEMap, FunctionSpace, Storage > >
    {
      typedef DiscontinuousLocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > ThisType;
      typedef DiscreteFunctionSpaceDefault< DiscontinuousLocalFiniteElementSpaceTraits< LFEMap, FunctionSpace, Storage > > BaseType;

      typedef typename BaseType::Traits Traits;

    public:
      typedef typename BaseType::FunctionSpaceType FunctionSpaceType;

      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::EntityType EntityType;
      typedef typename BaseType::IntersectionType IntersectionType;

      typedef typename BaseType::Traits::ShapeFunctionSetType ShapeFunctionSetType;
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;

      typedef typename BaseType::BlockMapperType BlockMapperType;

      typedef typename Traits::LocalFiniteElementType LocalFiniteElementType;

      typedef typename Traits::LFEMapType LFEMapType;

      typedef typename LFEMapType::KeyType KeyType;

    private:
      typedef typename LocalFiniteElementType::Traits::LocalBasisType LocalBasisType;
      typedef typename LocalFiniteElementType::Traits::LocalInterpolationType LocalInterpolationType;
      typedef typename LocalFiniteElementType::Traits::LocalCoefficientsType LocalCoefficientsType;

      typedef typename Traits::LocalFunctionsShapeFunctionSetType LocalFunctionsShapeFunctionSetType;

      struct LFEMapFactory
      {
        static LFEMapType *createObject ( std::pair< GridPartType *, KeyType > key ) { return new LFEMapType( *key.first, key.second ); }
        static void deleteObject ( LFEMapType *object ) { delete object; }
      };

      typedef SingletonList< std::pair< GridPartType *, KeyType >, LFEMapType, LFEMapFactory > LFEMapProviderType;

      typedef typename Traits::StoredShapeFunctionSetType StoredShapeFunctionSetType;
      typedef std::vector< std::unique_ptr< StoredShapeFunctionSetType > > StoredShapeFunctionSetVectorType;

      struct StoredShapeFunctionSetVectorFactory
      {
        static StoredShapeFunctionSetVectorType *createObject ( LFEMapType *lfeMap ) { return new StoredShapeFunctionSetVectorType( lfeMap->size() ); }
        static void deleteObject ( StoredShapeFunctionSetVectorType *object ) { delete object; }
      };

      typedef SingletonList< LFEMapType *, StoredShapeFunctionSetVectorType, StoredShapeFunctionSetVectorFactory > StoredShapeFunctionSetVectorProviderType;

      struct BlockMapperSingletonFactory
      {
        static BlockMapperType *createObject ( LFEMapType *lfeMap )
        {
          return new BlockMapperType( lfeMap->gridPart(), [ lfeMap ] ( const auto &refElement ) {
            if( lfeMap->hasCoefficients( refElement.type() ) )
              return Dune::Fem::generateCodimensionCode( refElement, 0, lfeMap->localCoefficients( refElement.type() ).size() );
            else
              return Dune::Fem::DofMapperCode();
          } );
        }

        static void deleteObject ( BlockMapperType *object ) { delete object; }
      };

      typedef SingletonList< LFEMapType *, BlockMapperType, BlockMapperSingletonFactory > BlockMapperProviderType;

    public:
      typedef LocalFiniteElementInterpolation< BasisFunctionSetType, LocalInterpolationType > InterpolationType;

      using BaseType::order;

      template< class GridPart, std::enable_if_t< std::is_same< GridPart, GridPartType >::value &&std::is_same< KeyType, std::tuple<> >::value, int > = 0 >
      explicit DiscontinuousLocalFiniteElementSpace ( GridPart &gridPart,
                                                      const InterfaceType commInterface = InteriorBorder_All_Interface,
                                                      const CommunicationDirection commDirection = ForwardCommunication )
        : BaseType( gridPart, commInterface, commDirection ),
          lfeMap_( &LFEMapProviderType::getObject( std::make_pair( &gridPart, KeyType() ) ) ),
          storedShapeFunctionSetVector_( &StoredShapeFunctionSetVectorProviderType::getObject( lfeMap_.get() ) ),
          blockMapper_( &BlockMapperProviderType::getObject( lfeMap_.get() ) )
      {}

      template< class GridPart, std::enable_if_t< std::is_same< GridPart, GridPartType >::value && !std::is_same< KeyType, std::tuple<> >::value, int > = 0 >
      explicit DiscontinuousLocalFiniteElementSpace ( GridPart &gridPart, const KeyType &key,
                                                      const InterfaceType commInterface = InteriorBorder_All_Interface,
                                                      const CommunicationDirection commDirection = ForwardCommunication )
        : BaseType( gridPart, commInterface, commDirection ),
          lfeMap_( &LFEMapProviderType::getObject( std::make_pair( &gridPart, key ) ) ),
          storedShapeFunctionSetVector_( &StoredShapeFunctionSetVectorProviderType::getObject( lfeMap_.get() ) ),
          blockMapper_( &BlockMapperProviderType::getObject( lfeMap_.get() ) )
      {}

      DiscontinuousLocalFiniteElementSpace ( const ThisType & ) = delete;
      DiscontinuousLocalFiniteElementSpace ( ThisType && ) = delete;

      ThisType &operator= ( const ThisType & ) = delete;
      ThisType &operator= ( ThisType && ) = delete;

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::type **/
      DFSpaceIdentifier type () const { return DGSpace_id; }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::basisFunctionSet **/
      BasisFunctionSetType basisFunctionSet ( const EntityType &entity ) const
      {
        return BasisFunctionSetType( entity, shapeFunctionSet( entity ) );
      }

      /**
       * \brief return shape function set for given entity
       *
       * \param[in]  entity  entity for which shape function set is requested
       *
       * \returns  ShapeFunctionSetType  shape function set
       **/
      ShapeFunctionSetType shapeFunctionSet ( const EntityType &entity ) const
      {
        return getShapeFunctionSet( (*lfeMap_)( entity ), entity.type() );
      }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      bool continuous () const { return false; }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      bool continuous ( const IntersectionType &intersection ) const { return false; }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order */
      int order () const { return lfeMap_->order(); }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::multipleGeometryTypes */
      bool multipleGeometryTypes () const { return true; }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::blockMapper */
      BlockMapperType &blockMapper () const { assert( blockMapper_ ); return *blockMapper_; }

      /**
       * \brief return local interpolation
       *
       *  \param[in]  entity  grid part entity
       **/
      InterpolationType interpolation ( const EntityType &entity ) const
      {
        auto lfe = (*lfeMap_)( entity );
        return InterpolationType( BasisFunctionSetType( entity, getShapeFunctionSet( lfe, entity.type() ) ), std::get< 2 >( lfe ) );
      }

    private:
      ShapeFunctionSetType getShapeFunctionSet ( std::tuple< std::size_t, const LocalBasisType &, const LocalInterpolationType & > lfe, const GeometryType &type ) const
      {
        auto &storedShapeFunctionSet = (*storedShapeFunctionSetVector_)[ std::get< 0 >( lfe ) ];
        if( !storedShapeFunctionSet )
          storedShapeFunctionSet.reset( new StoredShapeFunctionSetType( type, LocalFunctionsShapeFunctionSetType( std::get< 1 >( lfe ) ) ) );
        return ShapeFunctionSetType( storedShapeFunctionSet.get() );
      }

      std::unique_ptr< LFEMapType, typename LFEMapProviderType::Deleter > lfeMap_;
      std::unique_ptr< StoredShapeFunctionSetVectorType, typename StoredShapeFunctionSetVectorProviderType::Deleter > storedShapeFunctionSetVector_;
      std::unique_ptr< BlockMapperType, typename BlockMapperProviderType::Deleter > blockMapper_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_LOCALFINITEELEMENT_DGSPACE_HH
