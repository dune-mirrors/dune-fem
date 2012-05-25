#ifndef DUNE_FEM_PADATPTIVELAGRANGESPACE_HH
#define DUNE_FEM_PADATPTIVELAGRANGESPACE_HH

//- system includes 
#include <algorithm>

//- Dune includes 
#include <dune/common/misc.hh>
#include <dune/grid/common/grid.hh>

//- Dune-Fem includes 
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/padaptivespace/restrictprolong.hh>
#include <dune/fem/space/basefunctions/basefunctionstorage.hh>
#include <dune/fem/space/basefunctions/basefunctionsets.hh>
#include <dune/fem/space/basefunctions/basefunctionproxy.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/lagrangeinterpolation.hh>
#include <dune/fem/operator/projection/dgl2projection.hh>

#include <dune/fem/space/common/basesetlocalkeystorage.hh>

//- local includes 
#include <dune/fem/space/lagrangespace/basefunctions.hh>

#include "mapper.hh"

namespace Dune
{

  namespace Fem 
  {

    /** \class   GenericDiscreteFunctionSpace 
     *  \ingroup PAdaptiveLagrangeSpace
     *  \brief   Lagrange discrete function space
     */
    template< class SpaceImpTraits >
    class GenericDiscreteFunctionSpace 
      : public DiscreteFunctionSpaceDefault< SpaceImpTraits > 
    {
    public:
      typedef SpaceImpTraits  Traits;

    protected:  
      //! traits for the discrete function space
      typedef GenericDiscreteFunctionSpace< Traits > ThisType;
      typedef DiscreteFunctionSpaceDefault< Traits > BaseType;

    public:  
      typedef ThisType GenericDiscreteFunctionSpaceType;

      typedef typename Traits :: GridPartType GridPartType;
      typedef typename Traits :: GridType GridType;
      typedef typename Traits :: IndexSetType IndexSetType;
      typedef typename Traits :: IteratorType IteratorType;

      //! dimension of the grid (not the world)
      enum { dimension = GridType :: dimension };

      typedef typename Traits :: FunctionSpaceType FunctionSpaceType;
      //! field type for function space's domain
      typedef typename Traits :: DomainFieldType DomainFieldType;
      //! type for function space's domain
      typedef typename Traits :: DomainType DomainType;
      //! field type for function space's range
      typedef typename Traits :: RangeFieldType RangeFieldType;
      //! type for function space's range
      typedef typename Traits :: RangeType RangeType;
      //! dimension of function space's range
      enum { dimRange = FunctionSpaceType :: dimRange };
      //! type of scalar function space
      typedef typename Traits :: BaseFunctionSpaceType BaseFunctionSpaceType;
     
      //! maximum polynomial order of functions in this space
      enum { polynomialOrder = Traits :: polynomialOrder };
      
      //! type of the shape function set (on reference element)
      typedef typename Traits :: ShapeFunctionSetType ShapeFunctionSetType;

      // deprecated name 
      typedef ShapeFunctionSetType BaseFunctionSetImp ;
      
      //! type of BaseFunctionSet (entity dependent)
      typedef typename Traits :: BaseFunctionSetType BaseFunctionSetType;

      //! type of compiled local key 
      typedef typename Traits :: CompiledLocalKeyType  CompiledLocalKeyType;

      //! mapper used to implement mapToGlobal
      typedef typename Traits :: MapperType MapperType;

      //! mapper used to for block vector function 
      typedef typename Traits :: BlockMapperType BlockMapperType;

      //! size of local blocks
      enum { localBlockSize = Traits :: localBlockSize };

      //! type for DoF
      typedef RangeFieldType DofType;
      //! dimension of a value
      enum { dimVal = 1 };

      //! type of storage class for base function sets 
      typedef BaseSetLocalKeyStorage< ShapeFunctionSetType > ShapeSetStorageType;

      //! type of storage class for compiled local keys 
      typedef BaseSetLocalKeyStorage< CompiledLocalKeyType > LocalKeyStorageType;

      // vector containing storages for each polynomial order 
      typedef std::vector< ShapeSetStorageType > BaseSetVectorType;

      // vector containing storages for each polynomial order 
      typedef std::vector< LocalKeyStorageType > LocalKeyVectorType;

      template <int pOrd> 
      struct ConstructBaseFunctionSets
      {
       /** HelperClasses 
           \brief 
           CompiledLocalKeyFactory method createObject and
           deleteObject for the SingletonList  
        */
        class CompiledLocalKeyFactory
        {
        public:
          //! create new BaseFunctionSet 
          static CompiledLocalKeyType* createObject( const GeometryType& type )
          {
            return new CompiledLocalKeyType( type, pOrd );
          }

          //! delete BaseFunctionSet 
          static void deleteObject( CompiledLocalKeyType* obj )
          {
            delete obj;
          }
        };

        static void apply( BaseSetVectorType& baseFunctionSets, 
                           LocalKeyVectorType& compiledLocalKeys,
                           const GeometryType& geometryType ) 
        {
          typedef LagrangeBaseFunctionFactory
            < typename BaseFunctionSpaceType :: ScalarFunctionSpaceType, dimension, pOrd >
            ScalarFactoryType;

          //! type of singleton base function factory
          typedef BaseFunctionSetSingletonFactory
            < GeometryType, ShapeFunctionSetType, ScalarFactoryType >
            BaseFunctionSetSingletonFactoryType;

          //! type of singleton list (singleton provider) for base functions
          typedef SingletonList
            < GeometryType, ShapeFunctionSetType, BaseFunctionSetSingletonFactoryType >
          BaseFunctionSetSingletonProviderType;

          const size_t k = pOrd ;

          // insert base function set into list 
          baseFunctionSets[ k ].template insert< BaseFunctionSetSingletonProviderType> ( geometryType );

          //! type of singleton list (singleton provider) for compiled local keys 
          typedef SingletonList
            < GeometryType, CompiledLocalKeyType, CompiledLocalKeyFactory >
          CompiledLocalKeySingletonProviderType;

          // insert compiled local key 
          compiledLocalKeys[ k ].template insert< CompiledLocalKeySingletonProviderType > ( geometryType );
        }
      };

    public:
      //! type of identifier for this discrete function space
      typedef int IdentifierType;
      //! identifier of this discrete function space
      static const IdentifierType id = 665;
      
      //! type of DoF manager
      typedef DofManager< GridType > DofManagerType;

    protected:
      typedef typename Traits :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      //! type of intermediate storage 
      typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > 
          IntermediateStorageFunctionType;

      /// interface for list of p-adaptive functions 
      class PAdaptiveDiscreteFunctionEntryInterface   
      {
      protected:
        PAdaptiveDiscreteFunctionEntryInterface () {}
      public:
        virtual ~PAdaptiveDiscreteFunctionEntryInterface() {}
        virtual bool equals( void * ) const = 0 ;
        virtual void adaptFunction( IntermediateStorageFunctionType& tmp ) = 0;
      };

      template < class DF, class LocalInterpolation > 
      class PAdaptiveDiscreteFunctionEntry 
        : public PAdaptiveDiscreteFunctionEntryInterface
      {
        DF& df_;
        DofManagerType& dm_;

      public:
        PAdaptiveDiscreteFunctionEntry ( DF& df ) 
          : df_( df ),
            dm_( DofManagerType :: instance( df.space().grid() ) )
        {
        }

        virtual bool equals( void* ptr ) const  
        {
          return (((void *) &df_) == ptr );
        }

        virtual void adaptFunction( IntermediateStorageFunctionType& tmp ) 
        {
          //const int oldSize = tmp.space().size() ;
          //const int newSize = df_.space().size() ;

          //std::cout << " Start adaptFct: old size " << tmp.space().size() 
          //          << "                 new size " << df_.space().size() << std::endl;

          typedef typename IntermediateStorageFunctionType :: DofIteratorType
            TmpIteratorType;
          typedef typename DF :: DofIteratorType  DFIteratorType;

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
          LocalInterpolation :: apply( tmp, df_ );
        }
      };

      typedef std::list< PAdaptiveDiscreteFunctionEntryInterface* >
        PAdaptiveDiscreteFunctionListType; 
      typedef typename PAdaptiveDiscreteFunctionListType :: iterator DFListIteratorType;

    public:
      using BaseType :: gridPart;
      using BaseType :: asImp;

    public:
      /** \brief constructor
       *
       *  \param[in]  gridPart       grid part for the Lagrange space
       *  \param[in]  commInterface  communication interface to use 
       *  \param[in]  commDirection  communication direction to use
       */
      explicit GenericDiscreteFunctionSpace
        ( GridPartType &gridPart,
          const InterfaceType commInterface,
          const CommunicationDirection commDirection )
      : BaseType( gridPart, commInterface, commDirection ),
        baseFunctionSets_( polynomialOrder+1 ),
        compiledLocalKeys_( polynomialOrder+1 ),
        blockMapper_( initialize() ),
        mapper_( blockMapper() )
      {
      }

    protected:
      //! copy constructor needed for p-adaptation 
      GenericDiscreteFunctionSpace( const GenericDiscreteFunctionSpace& other ) 
      : BaseType( other.gridPart_, other.commInterface_, other.commDirection_ ),
        baseFunctionSets_( polynomialOrder+1 ),
        compiledLocalKeys_( polynomialOrder+1 ),
        blockMapper_( initialize( &other.blockMapper() ) ),
        mapper_( blockMapper() )
      {
      }

      //! initialize space and create block mapper 
      BlockMapperType* initialize( const BlockMapperType* otherMapper = 0 ) 
      {
        const IndexSetType &indexSet = gridPart().indexSet();

        AllGeomTypes< IndexSetType, GridType > allGeometryTypes( indexSet );
        const std :: vector< GeometryType >& geometryTypes
          = allGeometryTypes.geomTypes( 0 );

        for( unsigned int i = 0; i < geometryTypes.size(); ++i )
        {
          ForLoop< ConstructBaseFunctionSets, 1, polynomialOrder > :: 
            apply( baseFunctionSets_, compiledLocalKeys_, geometryTypes[ i ] );
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

    public:
      /** \brief Destructor (freeing base functions pointers and block mapper)
          \return 
      **/
      ~GenericDiscreteFunctionSpace ()
      {
        assert( dfList_.empty() );
        delete blockMapper_;
      }

      template <class DiscreteFunction> 
      DFListIteratorType searchFunction( const DiscreteFunction& df ) const
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

      template <class DiscreteFunction> 
      void removeFunction( const DiscreteFunction& df ) const
      {
        DFListIteratorType it = searchFunction( df );
        if( it != dfList_.end() )
        {
          delete (*it);
          dfList_.erase( it );
        }
      }

      /** \brief pAdaptation 
          \param polynomialOrders  vector containing polynomial orders for each cell 
          \param polOrderShift possible shift of polynomial order (i.e. in case of
                               Taylor-Hood put -1 for the pressure) (default = 0)
      */
      //-  --adapt 
      template <class Vector> 
      void adapt( const Vector& polynomialOrders, const int polOrderShift = 0 ) const
      {
        typedef typename IteratorType :: Entity EntityType ;
        const IteratorType endit = this->end();

        // create a copy of this space (to be improved)
        DiscreteFunctionSpaceType oldSpace( asImp() );

        //std::cout << "Old space size = " << oldSpace.size() << std::endl;

        // set new polynomial order for space  
        for( IteratorType it = this->begin(); it != endit; ++it )
        {
          const EntityType& entity = *it;
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

        DofManagerType& dm = DofManagerType :: instance( this->grid() );
        // resize discrete functions (only functions belonging 
        // to this space will be affected ), for convenience 
        dm.resize();
        dm.compress();

        //std::cout <<"This spaces size = " << this->size() << std::endl;
        //std::cout << std::endl;
      }

      /** \copydoc Dune::DiscreteFunctionSpaceInterface::contains */
      inline bool contains ( const int codim ) const
      {
        // forward to mapper since this information is held there 
        return blockMapper().contains( codim );
      }

      /** \copydoc Dune::DiscreteFunctionSpaceInterface::continuous */
      inline bool continuous () const
      {
        return Traits :: continuousSpace;
      }

      /** \brief this space has more than one base function set */
      inline bool multipleBaseFunctionSets () const
      {
        return (polynomialOrder > 1);
      }

      /** \brief get the type of this discrete function space 
          \return DFSpaceIdentifier
      **/
      inline DFSpaceIdentifier type () const
      {
        return LagrangeSpace_id;
      }

      /** \copydoc Dune::DiscreteFunctionSpaceInterface::order */
      inline int order () const 
      {
        return polynomialOrder;
      }
      /** \copydoc Dune::DiscreteFunctionSpaceInterface::order */
      inline int order (const typename BaseType::EntityType &entity) const 
      {
        return blockMapper().polynomOrder( entity );
      }

      /** \copydoc Dune::DiscreteFunctionSpaceInterface::baseFunctionSet(const EntityType &entity) const */
      template< class EntityType >
      inline const BaseFunctionSetType baseFunctionSet ( const EntityType &entity ) const
      {
        return baseFunctionSet( entity.type(), 
                                blockMapper().polynomOrder( entity ) );
      }

      /** \brief provide access to the base function set for a geometry type
       *
       *  \param[in]  type  type of geometry the base function set is requested for
       *
       *  \returns base function set for the specified geometry
       */
      inline const BaseFunctionSetType baseFunctionSet ( const GeometryType type ) const
      {
        return baseFunctionSet( type, polynomialOrder );
      }

      inline const BaseFunctionSetType baseFunctionSet ( const GeometryType type, const int k ) const
      {
        assert( k <= polynomialOrder );
        assert( k > 0 );
        return BaseFunctionSetType( &baseFunctionSets_[ k ][ type ] );
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
      inline const CompiledLocalKeyType &compiledLocalKey( const EntityType &entity ) const
      {
        return compiledLocalKey( entity.type(),
                                 blockMapper().polynomOrder( entity ) );
      }

      /** \brief provide access to the compiled local keys for a geometry type 
       *
       *  \note This method is not part of the DiscreteFunctionSpaceInterface. It
       *        is unique to the GenericDiscreteFunctionSpace.
       *
       *  \param[in]  type  type of geometry the compiled local key is requested for
       *
       *  \returns CompiledLocalKey 
       */
      inline const CompiledLocalKeyType &compiledLocalKey( const GeometryType type ) const
      {
        return compiledLocalKey( type, polynomialOrder );
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
      inline const CompiledLocalKeyType &compiledLocalKey( const GeometryType type, const int order ) const
      {
        return compiledLocalKeys_[ order ][ type ];
      }

      /** \brief get dimension of value
          \return int
      **/
      inline int dimensionOfValue () const
      {
        return dimVal;
      }

      /** \copydoc Dune::DiscreteFunctionSpaceInterface::mapper */
      MapperType &mapper () const
      {
        return mapper_;
      }

      /** \brief obtain the DoF block mapper of this space
          \return BlockMapperType
      **/
      BlockMapperType &blockMapper () const
      {
        assert( blockMapper_ != 0 );
        return *blockMapper_;
      }

    protected:
      //! storage for base function sets 
      mutable BaseSetVectorType baseFunctionSets_;
      
      //! storage for compiled local keys  
      mutable LocalKeyVectorType compiledLocalKeys_;
      
      //! corresponding mapper
      BlockMapperType *blockMapper_;

      //! corresponding mapper
      mutable MapperType mapper_;

      //! list of registered discrete functions 
      mutable PAdaptiveDiscreteFunctionListType dfList_;
    };



    /////////////////////////////////////////////////////////////////////
    //
    // --padaptive Lagrange space 
    //
    /////////////////////////////////////////////////////////////////////
    template< class FunctionSpaceImp,
              class GridPartImp,
              int polOrder,
              template< class > class BaseFunctionStorageImp = CachingStorage >
    class PAdaptiveLagrangeSpace;

    /** \addtogroup PAdaptiveLagrangeSpace
     *
     *  Provides access to base function sets for different element types in
     *  one grid and size of function space and maps from local to global dof
     *  number.
     *
     *  \note This space can only be used with special index sets. If you want
     *  to use the PAdaptiveLagrangeSpace with an index set only
     *  supporting the index set interface you will have to use the
     *  IndexSetWrapper class to provide the required functionality.
     *
     *  \note For adaptive calculations one has to use index sets that are
     *  capable of adaption (i.e. the method adaptive returns true). See also
     *  AdaptiveLeafIndexSet.
     */
    //- --padaptivetraits 
    template< class FunctionSpace, class GridPart, unsigned int polOrder,
              template< class > class BaseFunctionStorage = CachingStorage >
    struct PAdaptiveLagrangeSpaceTraits
    {
      dune_static_assert((polOrder > 0), "LagrangeSpace only defined for polOrder > 0" );
      
      typedef FunctionSpace FunctionSpaceType;
      typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
      typedef typename FunctionSpaceType :: DomainType DomainType;
      typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;
      typedef typename FunctionSpaceType :: RangeType RangeType;
      typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

      enum { dimRange = FunctionSpaceType :: dimRange };
      
      typedef GridPart GridPartType;
      typedef typename GridPartType :: GridType GridType;
      typedef typename GridPartType :: IndexSetType IndexSetType;
      typedef typename GridPartType :: template Codim< 0 > :: IteratorType
        IteratorType;

      // get dimension of local coordinate 
      enum { dimLocal = GridType :: dimension };

      typedef typename ToLocalFunctionSpace< FunctionSpaceType, dimLocal > :: Type 
        BaseFunctionSpaceType;

      enum { polynomialOrder = polOrder };
      
      typedef PAdaptiveLagrangeSpace
        < FunctionSpaceType, GridPartType, polynomialOrder, BaseFunctionStorage >
        DiscreteFunctionSpaceType;

      enum { localBlockSize = dimRange };

      //! this is a continuous space 
      static const bool continuousSpace = true ;

      // mapper for block
      typedef PAdaptiveLagrangeMapper< GridPartType, polynomialOrder > BlockMapperType;
      typedef NonBlockMapper< BlockMapperType, localBlockSize > MapperType;
      
      // implementation of shapefunction set 
      typedef VectorialBaseFunctionSet< BaseFunctionSpaceType, BaseFunctionStorage >
        ShapeFunctionSetType ;

      // exported type of base function set 
      typedef SimpleBaseFunctionProxy< ShapeFunctionSetType > BaseFunctionSetType;

      //! type of a compiled local key 
      typedef LagrangePointSet< GridPartType, polynomialOrder >
        CompiledLocalKeyType;

      /** \brief defines type of communication data handle for this type of space
       */
      template< class DiscreteFunction,
                class Operation = DFCommunicationOperation :: Add >
      struct CommDataHandle
      {
        //! type of data handle 
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
        //! type of operatation to perform on scatter 
        typedef Operation OperationType;
      };
    };
    

    /** \class   PAdaptiveLagrangeSpace
     *  \ingroup PAdaptiveLagrangeSpace
     *  \brief   Lagrange discrete function space
     */
    template< class FunctionSpaceImp,
              class GridPartImp,
              int polOrder,
              template< class > class BaseFunctionStorageImp >
    class PAdaptiveLagrangeSpace
    : public GenericDiscreteFunctionSpace
             < PAdaptiveLagrangeSpaceTraits< FunctionSpaceImp,
                                             GridPartImp,
                                             polOrder,
                                             BaseFunctionStorageImp > >
    {
    public:
      //! traits for the discrete function space
      typedef PAdaptiveLagrangeSpaceTraits< FunctionSpaceImp,
                                            GridPartImp,
                                            polOrder,
                                            BaseFunctionStorageImp >
        Traits;

      typedef GenericDiscreteFunctionSpace< Traits > BaseType ;

      //! type of the discrete function space
      typedef PAdaptiveLagrangeSpace< FunctionSpaceImp,
                                      GridPartImp,
                                      polOrder,
                                      BaseFunctionStorageImp >
              PAdaptiveLagrangeSpaceType;

      typedef typename Traits :: GridPartType GridPartType;
      typedef typename Traits :: GridType GridType;
      typedef typename Traits :: IndexSetType IndexSetType;
      typedef typename Traits :: IteratorType IteratorType;

      //! maximum polynomial order of functions in this space
      enum { polynomialOrder = Traits :: polynomialOrder };
      
      //! type of compiled local key 
      typedef typename Traits :: CompiledLocalKeyType  CompiledLocalKeyType;

      // deprecated name 
      typedef CompiledLocalKeyType LagrangePointSetType;

      //! mapper used to implement mapToGlobal
      typedef typename Traits :: MapperType MapperType;

      //! mapper used to for block vector function 
      typedef typename Traits :: BlockMapperType BlockMapperType;

      //! size of local blocks
      enum { localBlockSize = Traits :: localBlockSize };

      //! dimension of a value
      enum { dimVal = 1 };

      //! type of DoF manager
      typedef DofManager< GridType > DofManagerType;

      //! type of intersections
      typedef typename BaseType ::IntersectionType IntersectionType;

    public:
      using BaseType :: gridPart;
      using BaseType :: blockMapper;
      using BaseType :: compiledLocalKey;
      using BaseType :: order;

      //! default communication interface 
      static const InterfaceType defaultInterface = InteriorBorder_InteriorBorder_Interface;

      //! default communication direction 
      static const CommunicationDirection defaultDirection = ForwardCommunication;

      /** \brief constructor
       *
       *  \param[in]  gridPart       grid part for the Lagrange space
       *  \param[in]  commInterface  communication interface to use (optional)
       *  \param[in]  commDirection  communication direction to use (optional)
       */
      explicit PAdaptiveLagrangeSpace
        ( GridPartType &gridPart,
          const InterfaceType commInterface = defaultInterface,
          const CommunicationDirection commDirection = defaultDirection )
      : BaseType( gridPart, commInterface, commDirection )
      {
      }

      //! copy constructor needed for p-adaption 
      PAdaptiveLagrangeSpace( const PAdaptiveLagrangeSpace& other ) 
      : BaseType( other )
      {
      }

    protected:
      using BaseType :: dfList_ ;
      using BaseType :: searchFunction ;

    public:
      /*! \brief add function to discrete function space for p-adaptation 
          (currently only supported by AdaptiveDiscreteFunction )
       */
      template <class DiscreteFunction> 
      void addFunction( DiscreteFunction& df ) const
      {
        assert( searchFunction( df ) == dfList_.end() );
        // select LagrangeInterpolation to be the LocalInterpolation 
        typedef typename BaseType :: template PAdaptiveDiscreteFunctionEntry< 
            DiscreteFunction, LagrangeInterpolation< DiscreteFunction > > RealEntryType ;
        typedef typename BaseType :: PAdaptiveDiscreteFunctionEntryInterface
          EntryInterface;
        EntryInterface* entry = new RealEntryType( df );

        assert( entry );
        dfList_.push_front( entry );
      }

      //! deprecated method 
      template< class EntityType >
      inline const CompiledLocalKeyType &lagrangePointSet( const EntityType &entity ) const
      {
        return compiledLocalKey( entity.type(),
                                 blockMapper().polynomOrder( entity ) );
      }

      //! deprecated method 
      inline const CompiledLocalKeyType &lagrangePointSet( const GeometryType type ) const
      {
        return compiledLocalKey( type, polynomialOrder );
      }

      //! deprecated method 
      inline const CompiledLocalKeyType &lagrangePointSet( const GeometryType type, const int order ) const
      {
        return compiledLocalKey( type, order );
      }

      using BaseType::continuous;
      /** @copydoc Dune::DiscreteFunctionSpaceInterface::continuous */
      inline bool continuous (const IntersectionType &intersection) const
      { 
        if ( order() > 0 && intersection.conforming())
        {
          return true;
          if (intersection.neighbor())
            return (order(*(intersection.inside())) == order(*(intersection.outside())));
          else
            return true;
        }
        return false;
      }
    };


    /////////////////////////////////////////////////////////////////////
    //
    // --padaptive DG space 
    //
    /////////////////////////////////////////////////////////////////////

    template< class FunctionSpaceImp,
              class GridPartImp,
              int polOrder,
              template< class > class BaseFunctionStorageImp = CachingStorage >
    class PAdaptiveDGSpace;

    /** \addtogroup PAdaptiveDGSpace
     *
     *  Provides access to base function sets for different element types in
     *  one grid and size of function space and maps from local to global dof
     *  number.
     *
     *  \note This space can only be used with special index sets. If you want
     *  to use the PAdaptiveDGSpace with an index set only
     *  supporting the index set interface you will have to use the
     *  IndexSetWrapper class to provide the required functionality.
     *
     *  \note For adaptive calculations one has to use index sets that are
     *  capable of adaption (i.e. the method adaptive returns true). See also
     *  AdaptiveLeafIndexSet.
     */
    //- --padaptivetraits 
    template< class FunctionSpace, class GridPart, unsigned int polOrder,
              template< class > class BaseFunctionStorage = CachingStorage >
    struct PAdaptiveDGSpaceTraits 
      : public PAdaptiveLagrangeSpaceTraits
          < FunctionSpace, GridPart, polOrder, BaseFunctionStorage >
    {
      typedef PAdaptiveDGSpace
        < FunctionSpace, GridPart, polOrder, BaseFunctionStorage >
        DiscreteFunctionSpaceType;

      enum { localBlockSize = FunctionSpace :: dimRange };

      //! this is a continuous space 
      static const bool continuousSpace = false ;

      // mapper for block
      typedef PAdaptiveDGMapper< GridPart, polOrder > BlockMapperType;
      typedef NonBlockMapper< BlockMapperType, localBlockSize > MapperType;
      
      /** \brief defines type of communication data handle for this type of space
       */
      template< class DiscreteFunction,
                class Operation = DFCommunicationOperation :: Copy >
      struct CommDataHandle
      {
        //! type of data handle 
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
        //! type of operatation to perform on scatter 
        typedef Operation OperationType;
      };
    };
    

    /** \class   PAdaptiveDGSpace
     *  \ingroup PAdaptiveDGSpace
     *  \brief   adaptive DG discrete function space
     */
    template< class FunctionSpaceImp,
              class GridPartImp,
              int polOrder,
              template< class > class BaseFunctionStorageImp >
    class PAdaptiveDGSpace
    : public GenericDiscreteFunctionSpace
             < PAdaptiveDGSpaceTraits< FunctionSpaceImp,
                                       GridPartImp,
                                       polOrder,
                                       BaseFunctionStorageImp > >
    {
    public:
      //! traits for the discrete function space
      typedef PAdaptiveDGSpaceTraits< FunctionSpaceImp,
                                      GridPartImp,
                                      polOrder,
                                      BaseFunctionStorageImp >
        Traits;
      //! type of the discrete function space
      typedef PAdaptiveDGSpace< FunctionSpaceImp,
                                GridPartImp,
                                polOrder,
                                BaseFunctionStorageImp >
              PAdaptiveDGSpaceType;

    private:
      typedef GenericDiscreteFunctionSpace< Traits > BaseType;
      typedef PAdaptiveDGSpaceType ThisType;
    public:

      typedef typename Traits :: GridPartType GridPartType;
      typedef typename Traits :: GridType GridType;
      typedef typename Traits :: IndexSetType IndexSetType;
      typedef typename Traits :: IteratorType IteratorType;
      typedef typename BaseType ::IntersectionType IntersectionType;

      //! maximum polynomial order of functions in this space
      enum { polynomialOrder = Traits :: polynomialOrder };
      
      //! type of compiled local key 
      typedef typename Traits :: CompiledLocalKeyType  CompiledLocalKeyType;

      // deprecated name 
      typedef CompiledLocalKeyType LagrangePointSetType;

      //! mapper used to implement mapToGlobal
      typedef typename Traits :: MapperType MapperType;

      //! mapper used to for block vector function 
      typedef typename Traits :: BlockMapperType BlockMapperType;

      //! size of local blocks
      enum { localBlockSize = Traits :: localBlockSize };

      //! dimension of a value
      enum { dimVal = 1 };

      //! type of DoF manager
      typedef DofManager< GridType > DofManagerType;

      using BaseType :: gridPart;
      using BaseType :: blockMapper;
      using BaseType :: compiledLocalKey;

      //! default communication interface 
      static const InterfaceType defaultInterface = InteriorBorder_All_Interface;

      //! default communication direction 
      static const CommunicationDirection defaultDirection = ForwardCommunication;

      /** \brief constructor
       *
       *  \param[in]  gridPart       grid part for the Lagrange space
       *  \param[in]  commInterface  communication interface to use (optional)
       *  \param[in]  commDirection  communication direction to use (optional)
       */
      explicit PAdaptiveDGSpace
        ( GridPartType &gridPart,
          const InterfaceType commInterface = defaultInterface,
          const CommunicationDirection commDirection = defaultDirection )
      : BaseType( gridPart, commInterface, commDirection )
      {
      }

      //! copy constructor needed for p-adaption 
      PAdaptiveDGSpace( const PAdaptiveDGSpace& other ) 
      : BaseType( other )
      {
      }

    protected:
      using BaseType :: dfList_;
      using BaseType :: searchFunction;

    public:
      /*! \brief add function to discrete function space for p-adaptation 
          (currently only supported by AdaptiveDiscreteFunction )
       */
      template <class DiscreteFunction> 
      void addFunction( DiscreteFunction& df ) const
      {
        assert( searchFunction( df ) == dfList_.end() );
        // select L2Porjection to be the LocalInterpolation 
        typedef typename BaseType :: template PAdaptiveDiscreteFunctionEntry< 
            DiscreteFunction, DGL2ProjectionImpl > RealEntryType ;
        typedef typename BaseType :: PAdaptiveDiscreteFunctionEntryInterface
          EntryInterface;

        EntryInterface* entry = new RealEntryType( df );
        assert( entry );
        dfList_.push_front( entry );
      }

      //! deprecated method 
      template< class EntityType >
      inline const CompiledLocalKeyType &lagrangePointSet( const EntityType &entity ) const
      {
        return compiledLocalKey( entity.type(),
                                 blockMapper().polynomOrder( entity ) );
      }

      //! deprecated method 
      inline const CompiledLocalKeyType &lagrangePointSet( const GeometryType type ) const
      {
        return compiledLocalKey( type, polynomialOrder );
      }

      //! deprecated method 
      inline const CompiledLocalKeyType &lagrangePointSet( const GeometryType type, const int order ) const
      {
        return compiledLocalKey( type, order );
      }

      using BaseType::continuous;
      /** @copydoc Dune::DiscreteFunctionSpaceInterface::continuous */
      inline bool continuous (const IntersectionType &intersection) const
      { 
        return false;
      }

    };
  } // end namespace Fem
    
} // end Dune namespace  

// include definition of RestrictProlongDefault for Lagrange Space.
#include "adaptmanager.hh"

#endif // #ifndef DUNE_LAGRANGESPACE_LAGRANGESPACE_HH
