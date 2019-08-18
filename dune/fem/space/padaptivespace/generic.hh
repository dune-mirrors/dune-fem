#ifndef DUNE_FEM_SPACE_PADAPTIVE_GENERIC_HH
#define DUNE_FEM_SPACE_PADAPTIVE_GENERIC_HH

#include <cassert>
#include <list>
#include <vector>

#include <dune/common/math.hh>

#include <dune/fem/common/forloop.hh>
#include <dune/fem/function/vectorfunction.hh>
#include <dune/fem/space/basisfunctionset/default.hh>
#include <dune/fem/space/common/basesetlocalkeystorage.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/space/shapefunctionset/proxy.hh>
#include <dune/fem/space/shapefunctionset/selectcaching.hh>
#include <dune/fem/space/shapefunctionset/simple.hh>
#include <dune/fem/space/shapefunctionset/vectorial.hh>

#include <dune/fem/space/common/dataprojection/dataprojection.hh>
#include <dune/fem/space/common/dataprojection/tuple.hh>

namespace Dune
{

  namespace Fem
  {

    // GenericDiscreteFunctionSpace
    // ----------------------------

    /** \class   GenericDiscreteFunctionSpace
     *
     *  \ingroup PAdaptiveLagrangeSpace
     *
     *  \brief   Please doc me.
     */
    template< class Traits >
    class GenericDiscreteFunctionSpace
    : public DiscreteFunctionSpaceDefault< Traits >
    {
      typedef GenericDiscreteFunctionSpace< Traits > ThisType;
      typedef DiscreteFunctionSpaceDefault< Traits > BaseType;

    public:
      typedef ThisType GenericDiscreteFunctionSpaceType;
      typedef typename BaseType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

      typedef typename BaseType::FunctionSpaceType FunctionSpaceType;

      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::GridType GridType;
      typedef typename BaseType::IndexSetType IndexSetType;
      typedef typename BaseType::IteratorType IteratorType;
      typedef typename IteratorType::Entity EntityType;
      typedef typename BaseType::IntersectionType IntersectionType;

      typedef typename Traits::ShapeFunctionSetType ShapeFunctionSetType;
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;

      typedef typename BaseType::BlockMapperType BlockMapperType;

      /** \brief maximal available polynomial order */
      static const int polynomialOrder = Traits::polynomialOrder;

    protected:
      // single type for shape function sets of all polynomial orders
      typedef typename Traits::ScalarShapeFunctionSetType ScalarShapeFunctionSetType;
      // storage for scalar shape function set per polynomial order
      typedef BaseSetLocalKeyStorage< ScalarShapeFunctionSetType > ScalarShapeFunctionSetStorageType;
      // factory for shape function set of static order
      template< int pOrd >
      struct ScalarShapeFunctionSetFactory
      {
        typedef typename Traits::template ScalarShapeFunctionSetFactory< pOrd >::Type Type;
      };

    protected:
      template< int pOrd >
      struct Initialize;

      // DoF manager
      using BaseType :: dofManager_;

    public:
      typedef typename Traits::CompiledLocalKeyType CompiledLocalKeyType;
      typedef BaseSetLocalKeyStorage< CompiledLocalKeyType > LocalKeyStorageType;

      // key that identifies the basis function set, here the polynomial order
      typedef int KeyType;

      //! type of identifier for this discrete function space
      typedef int IdentifierType;
      //! identifier of this discrete function space
      static const IdentifierType id = 665;

      using BaseType::asImp;
      using BaseType::gridPart;

      /** \brief constructor
       *
       *  \param[in]  gridPart       grid part
       *  \param[in]  commInterface  communication interface to use
       *  \param[in]  commDirection  communication direction to use
       */
      GenericDiscreteFunctionSpace ( GridPartType &gridPart,
                                     const int order,
                                     const InterfaceType commInterface,
                                     const CommunicationDirection commDirection )
      : BaseType( gridPart, commInterface, commDirection ),
        order_( order ),
        scalarShapeFunctionSets_( order_+1 ),
        compiledLocalKeys_( order_+1 ),
        blockMapper_( initialize() )
      {
        assert( Capabilities::isAdaptiveDofMapper< BlockMapperType>::v );
      }

    protected:
      // copy constructor needed for p-adaptation
      GenericDiscreteFunctionSpace ( const GenericDiscreteFunctionSpace &other )
      : BaseType( other.gridPart_, other.commInterface_, other.commDirection_ ),
        order_( other.order_ ),
        scalarShapeFunctionSets_( order_+1 ),
        compiledLocalKeys_( order_+1 ),
        blockMapper_( initialize( &other.blockMapper() ) )
      {}

    public:

      ///////////////////////
      // Interface methods //
      ///////////////////////

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::type */
      inline DFSpaceIdentifier type () const { return GenericSpace_id; }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::basisFunctionSet */
      BasisFunctionSetType basisFunctionSet ( const EntityType &entity ) const
      {
        return BasisFunctionSetType( entity, shapeFunctionSet( entity ) );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      inline bool continuous () const { return Traits::continuousSpace; }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order */
      inline int order () const { return order_; }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order */
      inline int order (const typename BaseType::EntityType &entity) const
      {
        return blockMapper().polynomOrder( entity );
      }

      /** \brief this space has more than one base function set */
      inline bool multipleBaseFunctionSets () const { return (polynomialOrder > 1); }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::blockMapper */
      BlockMapperType &blockMapper () const
      {
        assert( blockMapper_ );
        return *blockMapper_;
      }


      ///////////////////////////
      // Non-interface methods //
      ///////////////////////////

      /** \brief return shape function set for given entity
       *
       *  \param[in]  entity  entity (of codim 0) for which shape function set
       *                      is requested
       *
       * \returns  ShapeFunctionSetType  shape function set
       */
      ShapeFunctionSetType shapeFunctionSet ( const EntityType &entity ) const
      {
        return shapeFunctionSet( entity.type(), order( entity ) );
      }

      /** \brief return shape unique function set for geometry type
       *
       *  \param[in]  type   geometry type (must be a cube) for which
       *                     shape function set is requested
       *  \param[in]  order  polynomial order
       *
       * \returns  ShapeFunctionSetType  shape function set
       */
      ShapeFunctionSetType shapeFunctionSet ( const GeometryType &type, const int order = polynomialOrder ) const
      {
        return ShapeFunctionSetType( &scalarShapeFunctionSets_[ order ][ type ] );
      }

      /** \brief provide access to the compiled local keys for an entity
       *
       *  \note This method is not part of the DiscreteFunctionSpaceInterface. It
       *        is unique to the GenericDiscreteFunctionSpace.
       *
       *  \param[in]  entity  entity the Lagrange point set is requested for
       *
       *  \returns CompiledLocalKey
       */
      template< class EntityType >
      inline const CompiledLocalKeyType &compiledLocalKey ( const EntityType &entity ) const
      {
        return compiledLocalKey( entity.type(), order( entity ) );
      }

      /** \brief provide access to the compiled local keys for a geometry type and polynomial order
       *
       *  \note This method is not part of the DiscreteFunctionSpaceInterface. It
       *        is unique to the GenericDiscreteFunctionSpace.
       *
       *  \param[in]  type  type of geometry the compiled local key is requested for
       *  \param[in]  order polynomial order for given geometry type
       *
       *  \returns CompiledLocalKey
       */
      inline const CompiledLocalKeyType &compiledLocalKey ( const GeometryType type, const int order = polynomialOrder ) const
      {
        return compiledLocalKeys_[ order ][ type ];
      }


      ////////////////////////////////
      // Adaptive interface methods //
      ////////////////////////////////
      /** \name Adaptation
       *  \{
       */

      /** \brief get identifiying basis function set key assigned to given entity
       *
       *  \param[in]  entity  grid part entity
       *
       *  \returns key
       */
      KeyType key ( const EntityType &entity ) const
      {
        return blockMapper().order( entity );
      }

      /** \brief assign new key to given entity
       *
       *  \param[in]  key  key identifying basis function set
       *  \param[in]  entity  grid part entity
       */
      void mark ( const KeyType &key, const EntityType &entity )
      {
        return blockMapper().suggestPolynomOrder( entity, key );
      }

      /** \brief get key to be assigned to an entity after next call to adapt()
       *
       *  \param[in]  entity  grid part entity
       *
       *  \returns key
       */
      KeyType getMark ( const EntityType &entity ) const
      {
        return blockMapper().suggestPolynomOrder( entity );
      }

      /** \brief p adaptation
       *
       *  \param[in] function  oialOrders  vector containing polynomial orders for each cell
       *  \param[in]  polOrderShift     possible shift of polynomial order (i.e. in case of
       *                                Taylor-Hood put -1 for the pressure) (default = 0)
       */
      void adapt ()
      {
        // adjust mapper by using previously set new polynomial orders
        blockMapper().adapt();

        // resize discrete functions (only functions belonging
        // to this space will be affected ), for convenience
        dofManager_.resize();
        dofManager_.compress();
      }

      template< class DiscreteFunctionSpace, class Implementation >
      void adapt ( DataProjection< DiscreteFunctionSpace, Implementation > &projection )
      {
        // create a copy of this space (to be improved, avoid DofManager involvement)
        DiscreteFunctionSpaceType oldSpace( asImp() );

        // adjust mapper by using previously set new polynomial orders
        blockMapper().adapt();

        // possibly enlarge memory attached to this space
        dofManager_.enlargeMemory();

        // create temporary storage for projection of discrete functions
        typedef std::vector< typename BaseType::RangeFieldType > TmpDofVectorType;
        TmpDofVectorType tmpVector( oldSpace.size() );

        // type of intermediate storage
        typedef VectorDiscreteFunction< DiscreteFunctionSpaceType, TmpDofVectorType > IntermediateStorageFunctionType;

        // Adapt space and then discrete functions
        IntermediateStorageFunctionType tmp( "padapt-temp", oldSpace, tmpVector );

        // go through list and adjust discrete functions
        // see DefaultDataProjectionTuple and DefaultDataProjection for possible implementation
        projection( tmp );

        // resize discrete functions (only functions belonging
        // to this space will be affected ), for convenience
        dofManager_.resize();
        dofManager_.compress();
      }

      /** \} */

    protected:
      // initialize space and create block mapper
      BlockMapperType *initialize ( const BlockMapperType *otherMapper = 0 )
      {
        const IndexSetType &indexSet = gridPart().indexSet();

        AllGeomTypes< IndexSetType, GridType > allGeometryTypes( indexSet );
        const std::vector< GeometryType > &geometryTypes
          = allGeometryTypes.geomTypes( 0 );

        for( unsigned int i = 0; i < geometryTypes.size(); ++i )
        {
          Fem::ForLoop< Initialize, 1, polynomialOrder >::
            apply( order_, scalarShapeFunctionSets_, compiledLocalKeys_, geometryTypes[ i ] );
        }

        if( otherMapper )
        {
          // make a copy of the other block mapper
          return new BlockMapperType( *otherMapper, order_, compiledLocalKeys_ );
        }
        else
        {
          // create new block mapper, this mapper is unique for each space since
          // the polynomial degrees might be different for each element
          return new BlockMapperType( gridPart(), order_, compiledLocalKeys_ );
        }
      }

    protected:
      // dynamically set maximal polynomial order
      const int order_;

      // storage for base function sets
      std::vector< ScalarShapeFunctionSetStorageType > scalarShapeFunctionSets_;
      // storage for compiled local keys
      std::vector< LocalKeyStorageType > compiledLocalKeys_;

      // corresponding mapper
      std::unique_ptr< BlockMapperType > blockMapper_;
    };



    // Implementation of GenericDiscreteFunctionSpace::Initialize
    // ----------------------------------------------------------

    template< class Traits >
    template <int pOrd>
    struct GenericDiscreteFunctionSpace< Traits >::Initialize
    {
      struct CompiledLocalKeyFactory
      {
        static CompiledLocalKeyType *createObject ( const GeometryType &type )
        {
          return new CompiledLocalKeyType( type, pOrd );
        }
        static void deleteObject ( CompiledLocalKeyType *obj )
        {
          delete obj;
        }
      };

      static void apply ( const int maxOrder,
                          std::vector< ScalarShapeFunctionSetStorageType > &scalarShapeFunctionSets,
                          std::vector< LocalKeyStorageType > &compiledLocalKeys,
                          const GeometryType &type )
      {
        // avoid creating shape function sets for polynomial orders that are not used
        if( pOrd > maxOrder ) return ;

        typedef typename ScalarShapeFunctionSetFactory< pOrd >::Type ScalarShapeFunctionSetFactoryType;
        typedef SingletonList< const GeometryType, ScalarShapeFunctionSetType, ScalarShapeFunctionSetFactoryType > SingletonProviderType;
        scalarShapeFunctionSets[ pOrd ].template insert< SingletonProviderType >( type );

        typedef SingletonList< GeometryType, CompiledLocalKeyType, CompiledLocalKeyFactory > CompiledLocalKeySingletonProviderType;
        compiledLocalKeys[ pOrd ].template insert< CompiledLocalKeySingletonProviderType >( type );
      }
    };

  } // namespace Fem

} // Dune namespace

#endif // #ifndef DUNE_FEM_SPACE_PADAPTIVE_GENERIC_HH
