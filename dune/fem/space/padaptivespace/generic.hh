#ifndef DUNE_FEM_SPACE_PADAPTIVE_GENERIC_HH
#define DUNE_FEM_SPACE_PADAPTIVE_GENERIC_HH

#include <cassert>
#include <list>
#include <vector>

#include <dune/common/math.hh>

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/space/basisfunctionset/default.hh>
#include <dune/fem/space/common/basesetlocalkeystorage.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/shapefunctionset/proxy.hh>
#include <dune/fem/space/shapefunctionset/selectcaching.hh>
#include <dune/fem/space/shapefunctionset/simple.hh>
#include <dune/fem/space/shapefunctionset/vectorial.hh>


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

    public:

      typedef typename Traits::CompiledLocalKeyType CompiledLocalKeyType;
      typedef BaseSetLocalKeyStorage< CompiledLocalKeyType > LocalKeyStorageType;

    protected:
      template< int pOrd > 
      struct Initialize;

      // type of intermediate storage 
      typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > IntermediateStorageFunctionType;
      
      /// interface for list of p-adaptive functions 
      struct PAdaptiveDiscreteFunctionEntryInterface
      {
        virtual ~PAdaptiveDiscreteFunctionEntryInterface() {}
        virtual bool equals( void * ) const = 0 ;
        virtual void adaptFunction( IntermediateStorageFunctionType &tmp ) = 0;

      protected:
        PAdaptiveDiscreteFunctionEntryInterface () {}
      };

      template < class DF, class LocalInterpolation > 
      class PAdaptiveDiscreteFunctionEntry;

      typedef std::list< PAdaptiveDiscreteFunctionEntryInterface * > PAdaptiveDiscreteFunctionListType; 
      typedef typename PAdaptiveDiscreteFunctionListType::iterator DFListIteratorType;
     
      // type of DoF manager
      typedef DofManager< GridType > DofManagerType;

    public:
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
                                     const InterfaceType commInterface,
                                     const CommunicationDirection commDirection )
      : BaseType( gridPart, commInterface, commDirection ),
        scalarShapeFunctionSets_( polynomialOrder+1 ),
        compiledLocalKeys_( polynomialOrder+1 ),
        blockMapper_( initialize() )
      {}

    protected:
      // copy constructor needed for p-adaptation
      GenericDiscreteFunctionSpace ( const GenericDiscreteFunctionSpace &other ) 
      : BaseType( other.gridPart_, other.commInterface_, other.commDirection_ ),
        scalarShapeFunctionSets_( polynomialOrder+1 ),
        compiledLocalKeys_( polynomialOrder+1 ),
        blockMapper_( initialize( &other.blockMapper() ) )
      {}

    public:
      // Destructor (freeing base functions pointers and block mapper)
      ~GenericDiscreteFunctionSpace ()
      {
        assert( dfList_.empty() );
        delete blockMapper_;
      }


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
      inline int order () const { return polynomialOrder; }

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

      /** \brief p adaptation 
       *
       *  \param[in]  polynomialOrders  vector containing polynomial orders for each cell 
       *  \param[in]  polOrderShift     possible shift of polynomial order (i.e. in case of
       *                                Taylor-Hood put -1 for the pressure) (default = 0)
       */
      template< class Vector > 
      void adapt ( const Vector &polynomialOrders, const int polOrderShift = 0 ) const
      {
        typedef typename IteratorType::Entity EntityType ;
        const IteratorType endit = this->end();

        // create a copy of this space (to be improved)
        DiscreteFunctionSpaceType oldSpace( asImp() );

        //std::cout << "Old space size = " << oldSpace.size() << std::endl;

        // set new polynomial order for space  
        for( IteratorType it = this->begin(); it != endit; ++it )
        {
          const EntityType &entity = *it;
          const int polOrder = polynomialOrders[ this->indexSet().index( entity ) ] + polOrderShift ;
          blockMapper().setPolynomOrder( entity, polOrder );
        }

        // adjust mapper 
        blockMapper().adapt();

        //std::cout << "New space size = " << blockMapper().size() << std::endl;
        //for(size_t i=0; i<polynomialOrders.size(); ++i)
        //  std::cout << "   " << polynomialOrders[ i ] << std::endl;

        // Adapt space and then discrete functions 
        IntermediateStorageFunctionType tmp( "padapt-temp", oldSpace );
        //std::cout << "created tmp with size = " << oldSpace.size() << " " << tmp.space().size() << std::endl;

        const DFListIteratorType endDF = dfList_.end();
        for( DFListIteratorType it = dfList_.begin(); it != endDF; ++it ) 
        {
          (*it)->adaptFunction( tmp );
        }

        DofManagerType &dm = DofManagerType::instance( this->grid() );
        // resize discrete functions (only functions belonging 
        // to this space will be affected ), for convenience 
        dm.resize();
        dm.compress();

        //std::cout <<"This spaces size = " << this->size() << std::endl;
        //std::cout << std::endl;
      }


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
          ForLoop< Initialize, 1, polynomialOrder >::
            apply( scalarShapeFunctionSets_, compiledLocalKeys_, geometryTypes[ i ] );
        }

        if( otherMapper ) 
        {
          // make a copy of the other block mapper 
          return new BlockMapperType( *otherMapper, compiledLocalKeys_ );
        }
        else 
        {
          // create new block mapper, this mapper is unique for each space since 
          // the polynomial degrees might be different for each element 
          return new BlockMapperType( gridPart(), compiledLocalKeys_ );
        }
      }

      template< class DiscreteFunction >
      DFListIteratorType searchFunction ( const DiscreteFunction &df ) const
      {
        assert( &df.space() == this );
        const DFListIteratorType endDF = dfList_.end();
        for( DFListIteratorType it = dfList_.begin(); it != endDF; ++it ) 
        {
          if( (*it)->equals( (void *) &df ) ) 
            return it;
        }
        return endDF;
      }

    public:
      template< class DiscreteFunction >
      void removeFunction ( const DiscreteFunction &df ) const
      {
        DFListIteratorType it = searchFunction( df );
        if( it != dfList_.end() )
        {
          delete (*it);
          dfList_.erase( it );
        }
      }

    protected:
      // storage for base function sets 
      std::vector< ScalarShapeFunctionSetStorageType > scalarShapeFunctionSets_;
      // storage for compiled local keys  
      std::vector< LocalKeyStorageType > compiledLocalKeys_; 

      // corresponding mapper
      BlockMapperType *blockMapper_;
      // list of registered discrete functions 
      mutable PAdaptiveDiscreteFunctionListType dfList_;
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

      static void apply ( std::vector< ScalarShapeFunctionSetStorageType > &scalarShapeFunctionSets, 
                          std::vector< LocalKeyStorageType > &compiledLocalKeys, 
                          const GeometryType &type ) 
      {
        typedef typename ScalarShapeFunctionSetFactory< pOrd >::Type ScalarShapeFunctionSetFactoryType;
        typedef SingletonList< const GeometryType, ScalarShapeFunctionSetType, ScalarShapeFunctionSetFactoryType > SingletonProviderType;
        scalarShapeFunctionSets[ pOrd ].template insert< SingletonProviderType >( type );

        typedef SingletonList< GeometryType, CompiledLocalKeyType, CompiledLocalKeyFactory > CompiledLocalKeySingletonProviderType;
        compiledLocalKeys[ pOrd ].template insert< CompiledLocalKeySingletonProviderType >( type );
      }
    };



    // Implementation of GenericDiscreteFunctionSpace::PAdaptiveDiscreteFunctionEntry
    // ------------------------------------------------------------------------------

    template< class Traits >
    template < class DF, class LocalInterpolation > 
    class GenericDiscreteFunctionSpace< Traits >::PAdaptiveDiscreteFunctionEntry 
    : public PAdaptiveDiscreteFunctionEntryInterface
    {
      DF &df_;
      DofManagerType &dm_;

    public:
      PAdaptiveDiscreteFunctionEntry ( DF &df ) 
      : df_( df ),
      dm_( DofManagerType::instance( df.space().grid() ) )
      {}

      virtual bool equals ( void *ptr ) const { return (((void *) &df_) == ptr ); }

      virtual void adaptFunction ( IntermediateStorageFunctionType &tmp ) 
      {
        //const int oldSize = tmp.space().size() ;
        //const int newSize = df_.space().size() ;

        //std::cout << " Start adaptFct: old size " << tmp.space().size() 
        //          << "                 new size " << df_.space().size() << std::endl;

        typedef typename IntermediateStorageFunctionType::DofIteratorType TmpIteratorType;
        typedef typename DF::DofIteratorType  DFIteratorType;

        // copy dof to temporary storage 
        DFIteratorType dfit = df_.dbegin();
        const TmpIteratorType endtmp = tmp.dend();
        for( TmpIteratorType it = tmp.dbegin(); it != endtmp; ++it, ++dfit )
        {
          assert( dfit != df_.dend() );
          *it = *dfit;
        }
        assert( dfit == df_.dend() );

        // adjust size of discrete function 
        df_.resize(); 

        //std::cout << " End adaptFct: old size " << tmp.space().size() 
        //          << "               new size " << df_.space().size() << std::endl;

        // interpolate to new space, this can be a 
        // Lagrange interpolation or a L2 projection
        LocalInterpolation::apply( tmp, df_ );
      }
    };

  } // namespace Fem
    
} // Dune namespace  

#endif // #ifndef DUNE_FEM_SPACE_PADAPTIVE_GENERIC_HH
