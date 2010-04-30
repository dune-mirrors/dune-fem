#ifndef  DUNE_FEM_PK2DLAGRANGESPACE_HH
#define  DUNE_FEM_PK2DLAGRANGESPACE_HH

#if HAVE_DUNE_LOCALFUNCTIONS

#include <dune/common/misc.hh>
#include <dune/grid/common/grid.hh>

#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/mapper/genericdofmapper.hh>
#include <dune/fem/space/basefunctions/genericbasefunctionsets.hh>
#include <dune/fem/space/basefunctions/basefunctionproxy.hh>
#include <dune/fem/space/lagrangespace/lagrangedatahandle.hh>

#include  <dune/localfunctions/lagrange/pk2d.hh>

namespace Dune
{

  template< class FunctionSpace, class GridPart, int polOrder = 1 >
  class Pk2DLagrangeSpace;



  // Pk2DLocalCoefficientsMap
  // ------------------------

  template< class IndexSet, class LocalFiniteElement >
  struct Pk2DLocalCoefficientsMap
  {
    typedef typename LocalFiniteElement::Traits::LocalCoefficientsType
      LocalCoefficientsType;

    typedef GenericBaseFunctionSet< typename LocalFiniteElement::Traits::LocalBasisType >
      BaseFunctionSetType;

    Pk2DLocalCoefficientsMap ( const IndexSet &indexSet )
    : indexSet_( indexSet )
    {
      for( int i = 0; i < 8; ++i )
      {
        fem_[ i ] = LocalFiniteElement( i );
        baseFunctionSet_[ i ] = new BaseFunctionSetType( fem_[ i ].localBasis(), fem_[ i ].type() );
      }
    }

    template< class Entity >
    unsigned int operator () ( const Entity &entity ) const
    {
      assert( entity.type().isSimplex() );
      const unsigned int n0 = indexSet_.subIndex( entity, 0, 2 );
      const unsigned int n1 = indexSet_.subIndex( entity, 1, 2 );
      const unsigned int n2 = indexSet_.subIndex( entity, 2, 2 );
      return int( n0 > n1 ) + 2*int( n0 > n2 ) + 4*int( n1 > n2 );
    }

    template< class Entity >
    const BaseFunctionSetType &baseFunctionSet ( const Entity &entity ) const
    {
      return *baseFunctionSet_[ (*this)( entity ) ];
    }

    template< class Topology >
    unsigned int size () const
    {
      return GenericGeometry::IsSimplex< Topology >::value ? 8 : 0;
    }

    template< class Topology >
    const LocalCoefficientsType &localCoefficients ( const unsigned int i ) const
    {
      assert( i < size< Topology >() );
      return fem_[ i ].localCoefficients();
    }

  private:
    const IndexSet &indexSet_;
    LocalFiniteElement fem_[ 8 ];
    const BaseFunctionSetType *baseFunctionSet_[ 8 ];
  };



  // Pk2DLagrangeSpaceTraits
  // -----------------------

  template< class FunctionSpace, class GridPart, int polOrder >
  struct Pk2DLagrangeSpaceTraits
  {
    typedef FunctionSpace FunctionSpaceType;
    typedef GridPart GridPartType;

    typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
    typedef typename FunctionSpaceType::ScalarFunctionSpaceType ScalarFunctionSpaceType;

    static const int dimRange  = FunctionSpaceType::dimRange;
    static const int dimDomain = FunctionSpaceType::dimDomain;
    dune_static_assert( (dimRange == 1), "Pk2DLagrangeSpace expects dimRange == 1" );
    dune_static_assert( (dimDomain == 2), "Pk2DLagrangeSpace expects dimDomain == 2" );
    
    typedef typename GridPartType::GridType GridType;
    typedef typename GridPartType::IndexSetType IndexSetType;
    typedef typename GridPartType::template Codim< 0 >::IteratorType IteratorType;

    static const int polynomialOrder = polOrder;

    typedef Pk2DLocalFiniteElement< DomainFieldType, RangeFieldType, polynomialOrder > LocalFiniteElementType;

    typedef Pk2DLocalCoefficientsMap< IndexSetType, LocalFiniteElementType > LocalCoefficientsMapType;
    typedef GenericDofMapper< GridPartType, LocalCoefficientsMapType > MapperType;
    typedef MapperType BlockMapperType;

    typedef GenericBaseFunctionSet< typename LocalFiniteElementType::Traits::LocalBasisType > BaseFunctionSetImp;
    typedef SimpleBaseFunctionProxy< BaseFunctionSetImp > BaseFunctionSetType;

    static const int localBlockSize = dimRange;

    template< class DiscreteFunction, class Operation = DFCommunicationOperation::Add >
    struct CommDataHandle
    {
      typedef LagrangeCommunicationHandler< DiscreteFunction, Operation > Type;
      typedef Operation OperationType;
    };

    typedef Pk2DLagrangeSpace< FunctionSpace, GridPart, polOrder > DiscreteFunctionSpaceType;
  };



  // Pk2DLagrangeSpace
  // -----------------

  template< class FunctionSpace, class GridPart, int polOrder >
  class Pk2DLagrangeSpace
  : public DiscreteFunctionSpaceDefault< Pk2DLagrangeSpaceTraits< FunctionSpace, GridPart, polOrder > >,
    public GenericDiscreteFunctionSpace
  {
    typedef Pk2DLagrangeSpace< FunctionSpace, GridPart, polOrder > ThisType;
    typedef DiscreteFunctionSpaceDefault< Pk2DLagrangeSpaceTraits< FunctionSpace, GridPart, polOrder > > BaseType;

  public:
    typedef typename BaseType::Traits Traits;

    typedef typename Traits :: GridPartType                          GridPartType;
    typedef typename Traits :: GridType                              GridType;
    typedef typename Traits :: IndexSetType                          IndexSetType;
    typedef typename Traits :: IteratorType                          IteratorType;
    //! dimension of the grid (not the world)
    enum { dimension = GridType :: dimension };

    typedef typename Traits :: FunctionSpaceType                     FunctionSpaceType;
    //! field type for function space's domain
    typedef typename Traits :: DomainFieldType                       DomainFieldType;
    //! type for function space's domain
    typedef typename Traits :: DomainType                            DomainType;
    //! field type for function space's range
    typedef typename Traits :: RangeFieldType                        RangeFieldType;
    //! type for function space's range
    typedef typename Traits :: RangeType                             RangeType;
    //! dimension of function space's range
    enum { dimRange = FunctionSpaceType :: dimRange };
    //! type of scalar function space
    typedef typename Traits :: ScalarFunctionSpaceType               ScalarFunctionSpaceType;
   
    //! maximum polynomial order of functions in this space
    enum { polynomialOrder = Traits :: polynomialOrder };
    
    //! type of the base function set(s)
    typedef typename Traits :: BaseFunctionSetImp                    BaseFunctionSetImp;
    
    typedef typename Traits :: BaseFunctionSetType                   BaseFunctionSetType;

    typedef typename Traits::LocalCoefficientsMapType LocalCoefficientsMapType;

    //! mapper used to implement mapToGlobal
    typedef typename Traits :: MapperType                            MapperType;

    //! local finite element type
    typedef typename Traits :: LocalFiniteElementType                LocalFiniteElementType;

    //! type for DoF
    typedef RangeFieldType DofType;

    //! dimension of a value
    enum { dimVal = 1 };

    //! type of DoF manager
    typedef DofManager< GridType > DofManagerType;

  public:
    //! type of identifier for this discrete function space
    typedef int IdentifierType;
    //! identifier of this discrete function space
    static const IdentifierType id = 664;

  public:
    using BaseType::gridPart;

  public:
    /** \brief constructor
     *
     *  \param[in]  gridPart  grid part for the Lagrange space
     */
    explicit Pk2DLagrangeSpace ( GridPartType &gridPart )
    : BaseType( gridPart ),
      localCoefficientsMap_( gridPart.indexSet() ),
      mapper_( gridPart, localCoefficientsMap_ )
    {}

  private:
    // forbid the copy constructor
    Pk2DLagrangeSpace ( const ThisType & );

  public:
    /** \copydoc Dune::DiscreteFunctionSpaceInterface::contains */
    bool contains ( const int codim ) const
    {
      return mapper().contains( codim );
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::continuous */
    bool continuous () const
    {
      return true;
    }

    /** \brief get the type of this discrete function space 
        \return DFSpaceIdentifier
    **/
    DFSpaceIdentifier type () const
    {
      return GenericSpace_id;
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::order */
    int order () const
    {
      return polynomialOrder;
    }

    template< class Entity >
    BaseFunctionSetType baseFunctionSet ( const Entity &entity ) const
    {
      return BaseFunctionSetType( &localCoefficientsMap_.baseFunctionSet( entity ) );
    }

    int dimensionOfValue () const
    {
      return dimVal;
    }

    MapperType &mapper () const
    {
      return mapper_;
    }

    MapperType &blockMapper () const
    {
      return mapper_;
    }

  private:
    LocalCoefficientsMapType localCoefficientsMap_;
    mutable MapperType mapper_;
  };

}

#endif // #if HAVE_DUNE_LOCALFUNCTIONS

#endif // #ifndef DUNE_FEM_PK2DLAGRANGESPACE_HH
