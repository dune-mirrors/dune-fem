#ifndef DUNE_LAGRANGESPACE_LAGRANGESPACE_HH
#define DUNE_LAGRANGESPACE_LAGRANGESPACE_HH

//- system includes 
#include <algorithm>

//- Dune includes 
#include <dune/common/misc.hh>
#include <dune/grid/common/grid.hh>

//- Dune-Fem includes 
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/basefunctionfactory.hh>
#include <dune/fem/space/basefunctions/basefunctionstorage.hh>
#include <dune/fem/space/basefunctions/basefunctionsets.hh>
#include <dune/fem/space/basefunctions/basefunctionproxy.hh>

//- local includes 
#include "basefunctions.hh"
#include "mapper.hh"

namespace Dune
{

  template< class FunctionSpaceImp,
            class GridPartImp,
            int polOrder,
            template< class > class BaseFunctionStorageImp = CachingStorage >
  class LagrangeDiscreteFunctionSpace;



  template< class FunctionSpaceImp,
            class GridPartImp,
            unsigned int polOrder,
            template< class > class BaseFunctionStorageImp = CachingStorage >
  struct LagrangeDiscreteFunctionSpaceTraits
  {
    CompileTimeChecker< (polOrder > 0) > __LagrangeSpace_only_defined_for_polOrder_greater_zero__;
    
    typedef FunctionSpaceImp FunctionSpaceType;
    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;
    typedef typename FunctionSpaceType :: RangeType RangeType;
    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;
    typedef typename FunctionSpaceType :: ScalarFunctionSpaceType
      ScalarFunctionSpaceType;
    enum { DimRange = FunctionSpaceType :: DimRange };
    
    typedef GridPartImp GridPartType;
    typedef typename GridPartType :: GridType GridType;
    typedef typename GridPartType :: IndexSetType IndexSetType;
    typedef typename GridPartType :: template Codim< 0 > :: IteratorType
      IteratorType;

    enum { polynomialOrder = polOrder };
    
    typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType,
                                           GridPartType,
                                           polynomialOrder,
                                           BaseFunctionStorageImp >
      DiscreteFunctionSpaceType;
    typedef LagrangeMapper< GridPartType, polynomialOrder, DimRange >
      MapperType;
    
    // mapper for block 
    typedef LagrangeMapper< GridPartType, polynomialOrder, 1 >
      BlockMapperType;
    
    // implementation of basefunction set 
    typedef VectorialBaseFunctionSet< FunctionSpaceType,
                                      BaseFunctionStorageImp >
        BaseFunctionSetImp;

    // exported type 
    typedef SimpleBaseFunctionProxy<BaseFunctionSetImp>  BaseFunctionSetType;

    enum { localBlockSize = DimRange };
  };



  /** @addtogroup LagrangeDiscreteFunctionSpace 
   *
   *  Provides access to bse function sets for different element types in
   *  one grid and size of function space and maps from local to global dof
   *  number.
   *
   *  \note This space can only be used with special index sets. If you want
   *  to use the LagrangeDiscreteFunctionSpace with an index set only
   *  supporting the index set interface you will have to use the
   *  IndexSetWrapper class to provide the required functionality.
   *
   *  \note For adaptive calculations one has to use index sets that are
   *  capable of adaption (i.e. the method adaptive returns true). See also
   *  AdaptiveLeafIndexSet.
   *
   *  @{
   **/

  /** @brief
   *  Lagrange Discrete Function Space
   **/
  template< class FunctionSpaceImp,
            class GridPartImp,
            int polOrder,
            template< class > class BaseFunctionStorageImp >
  class LagrangeDiscreteFunctionSpace
  : public DiscreteFunctionSpaceDefault
           < LagrangeDiscreteFunctionSpaceTraits< FunctionSpaceImp,
                                                  GridPartImp,
                                                  polOrder,
                                                  BaseFunctionStorageImp > >
  {
  public:
    //! traits for the discrete function space
    typedef LagrangeDiscreteFunctionSpaceTraits< FunctionSpaceImp,
                                                 GridPartImp,
                                                 polOrder,
                                                 BaseFunctionStorageImp >
      Traits;

    //! type of the discrete function space
    typedef LagrangeDiscreteFunctionSpace< FunctionSpaceImp,
                                           GridPartImp,
                                           polOrder,
                                           BaseFunctionStorageImp >
     LagrangeDiscreteFunctionSpaceType;

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
    enum { DimRange = FunctionSpaceType :: DimRange };
    //! type of scalar function space
    typedef typename Traits :: ScalarFunctionSpaceType ScalarFunctionSpaceType;
   
    //! maximum polynomial order of functions in this space
    enum { polynomialOrder = Traits :: polynomialOrder };
    
    //! type of the base function set(s)
    typedef typename Traits :: BaseFunctionSetImp BaseFunctionSetImp;
    
    typedef typename Traits :: BaseFunctionSetType BaseFunctionSetType;
    //! type of the base function set map
    typedef std :: map< const GeometryType, const BaseFunctionSetImp* >
      BaseFunctionMapType;
    //! type of base function factory
    typedef LagrangeBaseFunctionFactory< ScalarFunctionSpaceType,
                                         dimension,
                                         polynomialOrder >
      ScalarFactoryType;
    //! type of singleton base function factory
    typedef BaseFunctionSetSingletonFactory< GeometryType,
                                             BaseFunctionSetImp,
                                             ScalarFactoryType >
      BaseFunctionSetSingletonFactoryType;
    //! type of singleton list (singleton provider) for base functions
    typedef SingletonList< GeometryType,
                           BaseFunctionSetImp,
                           BaseFunctionSetSingletonFactoryType >
      BaseFunctionSetSingletonProviderType;

    //! type of a Lagrange point set
    typedef LagrangePointSet< GridPartType, polynomialOrder >
      LagrangePointSetType;
    //! type of Lagrange point set map
    typedef std :: map< const GeometryType, const LagrangePointSetType* >
      LagrangePointSetMapType;

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
    //! type of DoF manager
    typedef DofManager< GridType > DofManagerType;
    //! type of DoF manager factory
    typedef DofManagerFactory< DofManagerType > DofManagerFactoryType;

  public:
    //! type of identifier for this discrete function space
    typedef int IdentifierType;
    //! identifier of this discrete function space
    static const IdentifierType id = 665;
    
  private:
    typedef LagrangeDiscreteFunctionSpaceType ThisType;
    typedef DiscreteFunctionSpaceDefault< Traits > BaseType;

  private:
    //! map for the base function sets
    mutable BaseFunctionMapType baseFunctionSet_;

    //! map for the langrage point sets
    mutable LagrangePointSetMapType lagrangePointSet_;
    
    //! corresponding mapper
    MapperType *mapper_;

    //! corresponding mapper
    mutable BlockMapperType* blockMapper_;

    //! reference to the DoF manager
    DofManagerType &dofManager_;

  public:
    /** \brief Constructor generating a LagrangeBaseFunctionSet of the requested polynomial order for each element type of the grid 
        \param[in] gridPart 
        \return 
    **/
    inline explicit LagrangeDiscreteFunctionSpace ( GridPartType &gridPart )
    : BaseType( gridPart )
    , baseFunctionSet_()
    , lagrangePointSet_()
    , mapper_(0)
    , blockMapper_(0)
    , dofManager_( DofManagerFactoryType :: getDofManager( gridPart.grid() ) )
    {
      const IndexSetType &indexSet = gridPart.indexSet();
      GridType &grid = gridPart.grid();
      
      dofManager_.addIndexSet( grid, const_cast< IndexSetType& >( indexSet ) );

      AllGeomTypes< IndexSetType, GridType > allGeometryTypes( indexSet );
      const std :: vector< GeometryType >& geometryTypes
        = allGeometryTypes.geomTypes( 0 );
      for( unsigned int i = 0; i < geometryTypes.size(); ++i )
      {
        const GeometryType &geometryType = geometryTypes[ i ];
        
        if( baseFunctionSet_.find( geometryType ) == baseFunctionSet_.end() ) {
          const BaseFunctionSetImp *baseFunctionSet
            = &(BaseFunctionSetSingletonProviderType
                :: getObject( geometryType ));
          assert( baseFunctionSet != NULL );
          baseFunctionSet_[ geometryType ] = baseFunctionSet;
        }

        if( lagrangePointSet_.find( geometryType )
            == lagrangePointSet_.end() ) {
          const LagrangePointSetType *lagrangePointSet
            = new LagrangePointSetType( geometryType );
          assert( lagrangePointSet != NULL );
          lagrangePointSet_[ geometryType ] = lagrangePointSet;
        }
      }

      mapper_ = new LagrangeMapper< GridPartType, polynomialOrder, DimRange >
                  ( this->gridPart_, lagrangePointSet_ );
      assert( mapper_ != NULL );
    }

  private:
    // forbid the copy constructor
    LagrangeDiscreteFunctionSpace ( const ThisType& );

  public:
    /** \brief Destructor (freeing base functions and mapper)
        \return 
    **/
    inline ~LagrangeDiscreteFunctionSpace ()
    {
      delete blockMapper_;
      delete mapper_;

      typedef typename BaseFunctionMapType :: iterator BFIteratorType;
      BFIteratorType bfend = baseFunctionSet_.end();
      for( BFIteratorType it = baseFunctionSet_.begin(); it != bfend; ++it ) 
      {
        const BaseFunctionSetImp *baseFunctionSet = (*it).second;
        if( baseFunctionSet != NULL )
          BaseFunctionSetSingletonProviderType
          :: removeObject( *baseFunctionSet );
      }

      typedef typename LagrangePointSetMapType :: iterator LPIteratorType;
      LPIteratorType lpend = lagrangePointSet_.end();
      for( LPIteratorType it = lagrangePointSet_.begin(); it != lpend; ++it ) 
      {
        const LagrangePointSetType *lagrangePointSet = (*it).second;
        if( lagrangePointSet != NULL )
          delete lagrangePointSet;
      }
      
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::continuous */
    inline bool continuous () const
    {
      return (polynomialOrder > 0);
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

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::baseFunctionSet */
    template< class EntityType >
    inline const BaseFunctionSetType 
      baseFunctionSet ( const EntityType &entity ) const
    {
      return this->baseFunctionSet( entity.geometry().type() );
    }

    /** \brief provide access to the base function set for a geometry type
        \param[in] type
        \return BaseFunctionSetType
    **/
    inline const BaseFunctionSetType
      baseFunctionSet ( const GeometryType type ) const
    {
      assert( baseFunctionSet_.find( type ) != baseFunctionSet_.end() );
      assert( baseFunctionSet_[ type ] != NULL );
      return BaseFunctionSetType(baseFunctionSet_[ type ]);
    }

    /** \brief provide access to the Lagrange point set for an entity 
        \param[in] entity Entity the Lagrange point set is requested
        \return LagrangePointSet
    **/
    template< class EntityType >
    inline const LagrangePointSetType&
      lagrangePointSet ( const EntityType &entity ) const
    {
      return this->lagrangePointSet( entity.geometry().type() );
    }

    /** \brief provide access to the Lagrange point set for a geometry type
        \param[in] type
        \return LagrangePointSetType
    **/
    inline const LagrangePointSetType&
      lagrangePointSet ( const GeometryType type ) const
    {
      assert( lagrangePointSet_.find( type ) != lagrangePointSet_.end() );
      assert( lagrangePointSet_[ type ] != NULL );
      return *lagrangePointSet_[ type ];
    }

    /** \brief get dimension of value
        \return int
    **/
    inline int dimensionOfValue () const
    {
      return dimVal;
    }

    /** \brief obtain the DoF mapper of this space
        \return MapperType
    **/
    inline MapperType& mapper () const
    {
      assert( mapper_ != 0 );
      return *mapper_;
    }

    /** \brief obtain the DoF mapper of this space
        \return MapperType
    **/
    inline BlockMapperType& blockMapper () const
    {
      if( ! blockMapper_ )
      {
        blockMapper_ = new BlockMapperType( this->gridPart_, lagrangePointSet_ ); 
      }
      return *blockMapper_;
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::mapToGlobal */
    template< class EntityType >
    inline int mapToGlobal( const EntityType &entity, const int localDof ) const
    {
      return mapper().mapToGlobal( entity, localDof );
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::sequence */
    inline int sequence () const
    {
      return dofManager_.sequence();
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::size */
    inline int size () const
    {
      return mapper().size();
    }
  };
  /** @}
   **/
} // end Dune namespace  
#endif
