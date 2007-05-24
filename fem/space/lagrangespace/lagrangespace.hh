#ifndef DUNE_LAGRANGESPACE_LAGRANGESPACE_HH
#define DUNE_LAGRANGESPACE_LAGRANGESPACE_HH

#include <algorithm>

#include <dune/grid/common/grid.hh>

#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/allgeomtypes.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/common/basefunctionfactory.hh>
#include <dune/fem/space/basefunctions/basefunctionstorage.hh>
#include <dune/fem/space/basefunctions/basefunctionsets.hh>

#include "basefunctions.hh"
#include "mapper.hh"

namespace Dune
{

  template< class FunctionSpaceImp,
            class GridPartImp,
            int polOrder,
            template< class > class BaseFunctionStorageImp = SimpleStorage >
  class LagrangeDiscreteFunctionSpace;



  template< class FunctionSpaceImp,
            class GridPartImp,
            unsigned int polOrder,
            template< class > class BaseFunctionStorageImp = SimpleStorage >
  struct LagrangeDiscreteFunctionSpaceTraits
  {
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
    
    typedef VectorialBaseFunctionSet< FunctionSpaceType,
                                      BaseFunctionStorageImp >
      BaseFunctionSetType;

    enum { localBlockSize = 1 };
  };



  /** @defgroup LagrangeDiscreteFunctionSpace Lagrange Discrete Function Space
   *  @ingroup DiscreteFunctionSpace
   *
   *  Provides access to bse function sets for different element types in
   *  one grid and size of function space and maps from local to global dof
   *  number.
   *
   *  NOTE: This space can only be used with special index sets. If you want
   *  to use the LagrangeDiscreteFunctionSpace with an index set only
   *  supporting the index set interface you will have to use the
   *  IndexSetWrapper class to provide the required functionality.
   *
   *  NOTE: For adaptive calculations one has to use index sets that are
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

    //! type of the grid partition
    typedef typename Traits :: GridPartType GridPartType;
    //! type of the grid
    typedef typename Traits :: GridType GridType;
    //! type of the index set
    typedef typename Traits :: IndexSetType IndexSetType;
    //! type of iterator for entites of codim 0
    typedef typename Traits :: IteratorType IteratorType;
    //! dimension of the grid (not the world)
    enum { dimension = GridType :: dimension };

    //! type of function space
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
    typedef typename Traits :: BaseFunctionSetType BaseFunctionSetType;
    //! type of the base function set map
    typedef std :: map< const GeometryType, const BaseFunctionSetType* >
      BaseFunctionMapType;
    //! type of base function factory
    typedef LagrangeBaseFunctionFactory< ScalarFunctionSpaceType,
                                         dimension,
                                         polynomialOrder >
      ScalarFactoryType;
    //! type of singleton base function factory
    typedef BaseFunctionSetSingletonFactory< GeometryType,
                                             BaseFunctionSetType,
                                             ScalarFactoryType >
      SingletonFactoryType;
    //! type of singleton list (singleton provider)
    typedef SingletonList< GeometryType,
                           BaseFunctionSetType,
                           SingletonFactoryType >
      SingletonProviderType;

    //! type of the Lagrange point set factory
    typedef LagrangePointSetFactory< DomainFieldType,
                                     dimension,
                                     polynomialOrder >
      LagrangePointSetFactoryType;
    //! type of a Lagrange point set
    typedef typename LagrangePointSetFactoryType :: LagrangePointSetType
      LagrangePointSetType;
    typedef std :: map< const GeometryType, const LagrangePointSetType* >
      LagrangePointSetMapType;

    //! mapper used to implement mapToGlobal
    typedef typename Traits :: MapperType MapperType;

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
    typedef int IdentifierType;
    //! identifier of this discrete function space
    static const IdentifierType id = 665;
    
  private:
    typedef LagrangeDiscreteFunctionSpaceType ThisType;
    typedef DiscreteFunctionSpaceDefault< Traits > BaseType;

  private:
    //! corresponding grid partition
    GridPartType &gridPart_;

    //! map for the base function sets
    mutable BaseFunctionMapType baseFunctionSet_;
    //std :: vector< BaseFunctionSetType* > baseFunctionSet_;

    //! map for the langrage point sets
    mutable LagrangePointSetMapType lagrangePointSet_;
    //std :: vector< LagrangePointSetType* > lagrangePointSet_;
    
    //! corresponding mapper
    MapperType *mapper_;

    //! reference to the DoF manager
    DofManagerType &dofManager_;

  public:
    //! Constructor generating a LagrangeBaseFunctionSet of the requested
    //! polynomial order for each element type of the grid
    inline LagrangeDiscreteFunctionSpace ( GridPartType &gridPart )
    : BaseType( gridPart ),
      gridPart_( gridPart ),
      baseFunctionSet_(),
      lagrangePointSet_(),
      dofManager_( DofManagerFactoryType :: getDofManager( gridPart_.grid() ) )
    {
      const IndexSetType &indexSet = gridPart_.indexSet();
      GridType &grid = gridPart_.grid();
      
      dofManager_.addIndexSet( grid, const_cast< IndexSetType& >( indexSet ) );

      AllGeomTypes< IndexSetType, GridType > allGeometryTypes( indexSet );
      const std :: vector< GeometryType >& geometryTypes
        = allGeometryTypes.geomTypes( 0 );
      for( unsigned int i = 0; i < geometryTypes.size(); ++i )
      {
        const GeometryType &geometryType = geometryTypes[ i ];
        /*
        GeometryIdentifier :: IdentifierType id
          = GeometryIdentifier :: fromGeo( geometryType );
        */
        
        /*
        BaseFunctionSetType* &baseFunctionSet = baseFunctionSet_[ id ];
        if( baseFunctionSet == NULL )
          baseFunctionSet = &(SingletonProviderType
                              :: getObject( geometryType ));
        */
        if( baseFunctionSet_.find( geometryType ) == baseFunctionSet_.end() ) {
          const BaseFunctionSetType *baseFunctionSet
            = &(SingletonProviderType :: getObject( geometryType ));
          assert( baseFunctionSet != NULL );
          baseFunctionSet_[ geometryType ] = baseFunctionSet;
        }

        /*
        LagrangePointSetType* &lagrangePointSet = lagrangePointSet_[ id ];
        if( lagrangePointSet == NULL )
          lagrangePointSet = LagrangePointSetFactoryType
                             :: pointSet( geometryType.basicType() );
        */
        if( lagrangePointSet_.find( geometryType )
            == lagrangePointSet_.end() ) {
          const LagrangePointSetType *lagrangePointSet
            = LagrangePointSetFactoryType :: pointSet( geometryType );
          assert( lagrangePointSet != NULL );
          lagrangePointSet_[ geometryType ] = lagrangePointSet;
        }
      }

      mapper_ = new LagrangeMapper< GridPartType, polynomialOrder, DimRange >
                  ( gridPart_, lagrangePointSet_ );
      assert( mapper_ != NULL );
    }

    //! Destructor (freeing base functions and mapper)
    inline ~LagrangeDiscreteFunctionSpace ()
    {
      if( mapper_ != NULL )
        delete mapper_;

      typedef typename BaseFunctionMapType :: iterator BFIteratorType;
      BFIteratorType bfend = baseFunctionSet_.end();
      for( BFIteratorType it = baseFunctionSet_.begin(); it != bfend; ++it ) {
        const BaseFunctionSetType *baseFunctionSet = (*it).second;
        if( baseFunctionSet != NULL )
          SingletonProviderType :: removeObject( *baseFunctionSet );
      }

      typedef typename LagrangePointSetMapType :: iterator LPIteratorType;
      LPIteratorType lpend = lagrangePointSet_.end();
      for( LPIteratorType it = lagrangePointSet_.begin(); it != lpend; ++it ) {
        const LagrangePointSetType *lagrangePointSet = (*it).second;
        if( lagrangePointSet != NULL )
          delete lagrangePointSet;
      }

      /*
      for( unsigned int i = 0; i < baseFunctionSet_.size(); ++i ) {
        BaseFunctionSetType* &baseFunctionSet = baseFunctionSet_[ i ];
        if( baseFunctionSet != NULL )
          SingletonProviderType :: removeObject( *baseFunctionSet );
      }
      */
    }

    //! are the functions continuous?
    inline bool continuous () const
    {
      return (polynomialOrder > 0);
    }

    //! get the type of this discrete function space
    inline DFSpaceIdentifier type () const
    {
      return LagrangeSpace_id;
    }

    //! get the polynomial order of this discrete function space
    inline int order () const
    {
      return polynomialOrder;
    }

    inline int polynomOrder () const DUNE_DEPRECATED
    {
      return order();
    }

    //! begin iterator
    inline IteratorType begin () const
    {
      return gridPart_.template begin< 0 >();
    }

    //! end iterator
    inline IteratorType end () const
    {
      return gridPart_.template end< 0 >();
    }

    //! provide access to the base function set for an entity
    template< class EntityType >
    inline const BaseFunctionSetType& 
      baseFunctionSet ( const EntityType &entity ) const
    {
      return this->baseFunctionSet( entity.geometry().type() );
    }

    //! provide access to the base function set for a geometry type
    inline const BaseFunctionSetType&
      baseFunctionSet ( const GeometryType type ) const
    {
      assert( baseFunctionSet_.find( type ) != baseFunctionSet_.end() );
      assert( baseFunctionSet_[ type ] != NULL );
      return *baseFunctionSet_[ type ];
    }
    
    //! provide access to the Lagrange point set for an entity
    template< class EntityType >
    inline const LagrangePointSetType&
      lagrangePointSet ( const EntityType &entity ) const
    {
      return this->lagrangePointSet( entity.geometry().type() );
    }

    //! provide access to the Lagrange point set for a geometry type
    inline const LagrangePointSetType&
      lagrangePointSet ( const GeometryType type ) const
    {
      assert( lagrangePointSet_.find( type ) != lagrangePointSet_.end() );
      assert( lagrangePointSet_[ type ] != NULL );
      return *lagrangePointSet_[ type ];
    }

    //! get dimension of balue
    inline int dimensionOfValue () const
    {
      return dimVal;
    }

    //! obtain the associated grid
    inline const GridType& grid () const
    {
      return gridPart_.grid();
    }
    
    //! obtain the associated grid partition
    inline const GridPartType& gridPart () const
    {
      return gridPart_;
    }

    //! obtain the associated grid partition
    inline GridPartType& gridPart ()
    {
      return gridPart_;
    }

    //! obtain the associated index set
    inline const IndexSetType& indexSet () const
    {
      return gridPart_.indexSet();
    }

    //! obtain the DoF mapper of this space
    inline const MapperType& mapper () const
    {
      return *mapper_;
    }

    //! map local DoF number to global DoF number
    template< class EntityType >
    inline int mapToGlobal( EntityType &entity, int localDof ) const
    {
      return mapper_->mapToGlobal( entity, localDof );
    }

    //! get index in grid sequences
    inline int sequence () const
    {
      return dofManager_.sequence();
    }

    //! number of DoFs in the function space
    inline int size () const
    {
      return mapper_->size();
    }
  };
  /** @}
   **/

}

#endif
