#ifndef DUNE_FEM_SPACE_LOCALFINITEELEMENT_SPACE_HH
#define DUNE_FEM_SPACE_LOCALFINITEELEMENT_SPACE_HH

#include <cassert>

#include <memory>
#include <utility>
#include <vector>

#include <dune/geometry/type.hh>

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

#include <dune/fem/space/localfiniteelement/capabilities.hh>
#include <dune/fem/space/rannacherturek/dofmappercode.hh>
#include <dune/fem/space/localfiniteelement/interpolation.hh>

namespace Dune
{

  namespace Fem
  {

    // LocalFiniteElementDiscreteFunctionSpaceTraits
    // ---------------------------------------------

    template< class LFEMap, class FunctionSpace, template< class > class Storage >
    struct LocalFiniteElementDiscreteFunctionSpaceTraits
    {
      typedef LocalFiniteElementDiscreteFunctionSpace< FunctionSpace, GridPart, Storage > DiscreteFunctionSpaceType;

      typedef typename LFEMap::GridPartType GridPartType;
      typedef typename LFEMap::LocalFiniteElementType LocalFiniteElementType;

      typedef FunctionSpace FunctionSpaceType;

      static constexpr int codimension = 0;
      static constexpr int localBlockSize = 1;

    private:
      typedef typename GridPartType::template Codim< codimension >::EntityType EntityType;

    public:
      typedef Dune::Fem::IndexSetDofMapper< GridPartType > BlockMapperType;

      typedef LocalFunctionsShapeFunctionSet< typename LocalFiniteElementType::Traits::LocalBasisType > LocalFunctionsShapeFunctionSetType;
      typedef SelectCachingShapeFunctionSet< LocalFunctionsShapeFunctionSetType, Storage > StoredShapeFunctionSetType;

      typedef VectorialShapeFunctionSet< ShapeFunctionSetProxy< StoredShapeFunctionSetType >, typename FunctionSpaceType::RangeType > ShapeFunctionSetType;

      typedef DefaultBasisFunctionSet< EntityType, ShapeFunctionSetType > BasisFunctionSetType;

      template< class DiscreteFunction, class Operation = DFCommunicationOperation::Add >
      struct CommDataHandle
      {
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
        typedef Operation OperationType;
      };
    };



    // LocalFiniteElementDiscreteFunctionSpace
    // ---------------------------------------

    /** \brief Rannacher-Turek Space
     *  \class LocalFiniteElementDiscreteFunctionSpace
     *  \ingroup DiscreteFunctionSpace
     *
     *  \note The LocalFiniteElementDiscreteFunctionSpace depends on
     *        dune-localfunctions (see http://www.dune-project.org).
     *
     *  \todo please doc me
     **/
    template< class LFEMap, class FunctionSpace, template< class > class Storage = CachingStorage >
    struct LocalFiniteElementDiscreteFunctionSpace
    : public DiscreteFunctionSpaceDefault< RannacherTurekDiscreteFunctionSpaceTraits< LFEMap, FunctionSpace, Storage > >
    {
      typedef LocalFiniteElementDiscreteFunctionSpace< LFEMap, FunctionSpace, Storage > ThisType;
      typedef DiscreteFunctionSpaceDefault< LocalFiniteElementDiscreteFunctionSpaceTraits< FunctionSpace, GridPart, Storage > > BaseType;

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

      typedef typename BaseType::BlockMapperType BlockMapperType;

      typedef typename Traits::LocalFiniteElementType LocalFiniteElementType;

    private:
      typedef typename LocalFiniteElementType::Traits::LocalBasisType LocalBasisType;
      typedef typename LocalFiniteElementType::Traits::LocalInterpolationType LocalInterpolationType;
      typedef typename LocalFiniteElementType::Traits::LocalCoefficientsType LocalCoefficientsType;

      typedef typename LFEMap::KeyType KeyType;

      struct LFEMapFactory
      {
        static LFEMap *createObject ( std::pair< GridPartType *, KeyType > key ) { return new LFEMap( *key.first, key.second ); }
        static void deleteObject ( LFEMap *object ) { delete object; }
      };

      typedef SingletonList< std::pair< GridPartType *, KeyType >, LFEMap, LFEMapFactory > LFEMapProviderType;

      typedef std::vector< std::unique_ptr< StoredShapeFunctionSetType > > StoredShapeFunctionSetVectorType;

      struct StoredShapeFunctionSetVectorFactory
      {
        static StoredShapeFunctionSetVectorType *createObject ( LFEMap *lfeMap ) { return new StoredShapeFunctionSetVectorType( lfeMap_->size() ); }
        static void deleteObject ( StoredShapeFunctionSetVectorType *object ) { delete object; }
      };

      typedef SingletonList< LFEMap *, StoredShapeFunctionSetVectorType, StoredShapeFunctionSetVectorFactory > StoredShapeFunctionSetVectorProviderType;

      struct BlockMapperSingletonFactory
      {
        static BlockMapperType *createObject ( LFEMap *lfeMap )
        {
          return new BlockMapperType( lfeMap->gridPart(), [ lfeMap ] ( const auto &refElement ) {
              return Dune::Fem::compile( refElement, lfeMap->localCoefficients( refElement.type() ) );
            } );
        }

        static void deleteObject ( BlockMapperType *object ) { delete object; }
      };

      typedef SingletonList< LFEMap *, BlockMapperType *, BlockMapperSingletonFactory > BlockMapperProviderType;

    public:
      typedef LocalFiniteElementInterpolation< BasisFunctionSetType, LocalInterpolationType > InterpolationType;

      using BaseType::order;

      template< class GridPart, std::enable_if_t< std::is_same< GridPart, GridPartType >::value && std::is_same< KeyType, std::tuple<> >::value, int > = 0 >
      explicit LocalFiniteElementDiscreteFunctionSpace ( GridPart &gridPart,
                                                         const InterfaceType commInterface = defaultInterface,
                                                         const CommunicationDirection commDirection = defaultDirection )
        : BaseType( gridPart, commInterface, commDirection ),
          lfeMap_( &LFEMapProviderType::getObject( std::make_pair( &gridPart, KeyType() ) ) ),
          storedShapeFunctionSetVector_( &StoredShapeFunctionSetVectorProviderType::getObject( lfeMap_.get() ) ),
          blockMapper_( &BlockMapperProviderType::getObject( lfeMap_.get() )
      {}

      template< class GridPart, std::enable_if_t< std::is_same< GridPart, GridPartType >::value && !std::is_same< KeyType, std::tuple<> >::value, int > = 0 >
      explicit LocalFiniteElementDiscreteFunctionSpace ( GridPart &gridPart, const KeyType &key,
                                                         const InterfaceType commInterface = defaultInterface,
                                                         const CommunicationDirection commDirection = defaultDirection )
        : BaseType( gridPart, commInterface, commDirection ),
          lfeMap_( &LFEMapProviderType::getObject( std::make_pair( &gridPart, key ) ) ),
          storedShapeFunctionSetVector_( &StoredShapeFunctionSetVectorProviderType::getObject( lfeMap_.get() ) ),
          blockMapper_( &BlockMapperProviderType::getObject( lfeMap_.get() )
      {}

      LocalFiniteElementDiscreteFunctionSpace ( const ThisType & ) = delete;
      LocalFiniteElementDiscreteFunctionSpace ( ThisType && ) = delete;

      ThisType &operator= ( const ThisType & ) = delete;
      ThisType &operator= ( ThisType && ) = delete;

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::type **/
      DFSpaceIdentifier type () const { return LocalFiniteElementSpace_id; }

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
        return shapeFunctionSet( lfeMap_( entity ) );
      }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      bool continuous () const { return false; }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      bool continuous ( const IntersectionType &intersection ) const { return false; }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order */
      int order () const { return polynomialOrder; }

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
        auto lfe = lfeMap_( entity );
        return InterpolationType( BasisFunctionSetType( entity, shapeFunctionSet( lfe ) ), std::get< 2 >( lfe ) );
      }

    private:
      ShapeFunctionSetType getShapeFunctionSet ( std::tuple< std::size_t, const LocalBasisType &, const LocalInterpolationType & > lfe )
      {
        auto &storedShapeFunctionSet = storedShapeFunctionSetVector_[ std::get< 0 >( lfe ) ];
        if( !storedShapeFunctionSet )
          storedShapeFunctionSet.reset( new StoredShapeFunctionSetType( std::get< 1 >( lfe ) ) );
        return ShapeFunctionSetType( storedShapeFunctionSet.get() );
      }

      std::unique_ptr< LFEMap, typename LFEMapProviderType::Deleter > lfeMap_;
      std::unique_ptr< StoredShapeFunctionSetVectorType, typename StoredShapeFunctionSetVectorProviderType::Delete > storedShapeFunctionSetVector_;
      BlockMapperType *blockMapper_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_LOCALFINITEELEMENT_SPACE_HH
