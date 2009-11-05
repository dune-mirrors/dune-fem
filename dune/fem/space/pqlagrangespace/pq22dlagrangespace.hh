#ifndef  DUNE_FEM_PQ22DSPACE_HH
#define  DUNE_FEM_PQ22DSPACE_HH

#if HAVE_DUNE_LOCALFUNCTIONS

//- Dune includes 
#include <dune/common/misc.hh>
#include <dune/grid/common/grid.hh>

//- Dune-Fem includes 
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/genericdofmapper.hh>
#include <dune/fem/space/basefunctions/genericbasefunctionsets.hh>
#include <dune/fem/space/basefunctions/basefunctionproxy.hh>
#include <dune/fem/space/lagrangespace/lagrangedatahandle.hh>

//- Dune-Localfunctions include
#include  <dune/finiteelements/pq22d.hh>

namespace Dune
{

  template< class FunctionSpaceImp, class GridPartImp >
  class PQ22DLagrangeSpace;

  // PQLocalCoefficientsMap
  // ----------------------

  template <class LocalFiniteElement>
  struct PQ22DLocalCoefficientsMap
  {
    typedef typename LocalFiniteElement::Traits::LocalCoefficientsType
      LocalCoefficientsType;
    typedef GenericBaseFunctionSet< typename LocalFiniteElement::Traits::LocalBasisType >
      BaseFunctionSetType;

    PQ22DLocalCoefficientsMap ( )
    : femSimplex_( GeometryType(GeometryType::simplex,2) ),
      femCube_( GeometryType(GeometryType::cube,2) ),
      basisSimplex_( femSimplex_.localBasis(), femSimplex_.type() ),
      basisCube_( femCube_.localBasis(), femCube_.type() )
    {}

    template< class Entity >
    int operator () ( const Entity &entity ) const
    {
      return 0;
    }

    template< class Topology >
    unsigned int size () const
    {
      return GenericGeometry::IsSimplex< Topology >::value || 
             GenericGeometry::IsCube< Topology >::value ? 1 : 0;
    }

    template< class Entity >
    const BaseFunctionSetType &baseFunctionSet ( const Entity &entity ) const
    {
      if (entity.type().isTriangle())
        return basisSimplex_;
      if (entity.type().isQuadrilateral())
        return basisCube_;
      abort();
    }

    template< class Topology >
    const LocalCoefficientsType &localCoefficients ( const unsigned int i ) const
    {
      if (GenericGeometry::IsSimplex< Topology >::value) 
        return femSimplex_.localCoefficients();
      if (GenericGeometry::IsCube< Topology >::value) 
        return femCube_.localCoefficients();
    }

  private:
    const LocalFiniteElement femSimplex_,femCube_;
    const BaseFunctionSetType basisSimplex_, basisCube_;
  };



  template< class FunctionSpace, class GridPart >
  struct PQ22DLagrangeSpaceTraits
  {
    typedef FunctionSpace                                            FunctionSpaceType;
    typedef typename FunctionSpaceType :: DomainFieldType            DomainFieldType;
    typedef typename FunctionSpaceType :: DomainType                 DomainType;
    typedef typename FunctionSpaceType :: RangeFieldType             RangeFieldType;
    typedef typename FunctionSpaceType :: RangeType                  RangeType;
    typedef typename FunctionSpaceType :: JacobianRangeType          JacobianRangeType;
    typedef typename FunctionSpaceType :: ScalarFunctionSpaceType    ScalarFunctionSpaceType;
    static const unsigned int dimRange  = FunctionSpaceType :: dimRange;
    static const unsigned int dimDomain = FunctionSpaceType :: dimDomain;
    dune_static_assert(( dimRange == 1 ), "P12DSpace expects range dimension == 1");
    
    typedef GridPart                                                 GridPartType;
    typedef typename GridPartType :: GridType                        GridType;
    typedef typename GridPartType :: IndexSetType                    IndexSetType;
    typedef typename GridPartType :: template Codim< 0 >
              :: IteratorType                                        IteratorType;

    static const int polynomialOrder = 2;
    typedef Dune::PQ22DLocalFiniteElement< DomainFieldType, RangeFieldType > LocalFiniteElementType;

  public:
   
    typedef PQ22DLocalCoefficientsMap< LocalFiniteElementType >        LocalCoefficientsMapType;
    typedef GenericDofMapper< GridPartType, LocalCoefficientsMapType > MapperType;
    typedef MapperType BlockMapperType;

    // implementation of basefunction set 
    typedef GenericBaseFunctionSet< typename LocalFiniteElementType::Traits::LocalBasisType > BaseFunctionSetImp;

    // exported type 
    typedef SimpleBaseFunctionProxy< BaseFunctionSetImp >            BaseFunctionSetType;

    enum { localBlockSize = dimRange };

    /** \brief defines type of communication data handle for this type of space
     */
    template< class DiscreteFunction,
	      class Operation = DFCommunicationOperation :: Add >
    struct CommDataHandle
    {
      //! type of data handle 
      typedef LagrangeCommunicationHandler< DiscreteFunction,
                                            Operation >              Type;
      //! type of operatation to perform on scatter 
      typedef Operation                                              OperationType;
    };

    typedef PQ22DLagrangeSpace< FunctionSpace, GridPart > DiscreteFunctionSpaceType;
  };



  // PQ22DLagrangeSpace
  // ---------------

  template< class FunctionSpaceImp, class GridPartImp >
  class PQ22DLagrangeSpace
  : public DiscreteFunctionSpaceDefault< PQ22DLagrangeSpaceTraits< FunctionSpaceImp, GridPartImp > >,
    public GenericDiscreteFunctionSpace
  {
    typedef PQ22DLagrangeSpace< FunctionSpaceImp, GridPartImp> ThisType;
    typedef DiscreteFunctionSpaceDefault
              < PQ22DLagrangeSpaceTraits< FunctionSpaceImp, GridPartImp > >        BaseType;

  public:
    //! traits for the discrete function space
    typedef PQ22DLagrangeSpaceTraits< FunctionSpaceImp, GridPartImp >   Traits;

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

    //! size of local blocks
    enum { localBlockSize = Traits :: localBlockSize };

    //! type for DoF
    typedef RangeFieldType                                           DofType;
    //! dimension of a value
    enum { dimVal = 1 };
    //! type of DoF manager
    typedef DofManager< GridType >                                   DofManagerType;

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
    explicit PQ22DLagrangeSpace ( GridPartType &gridPart )
    : BaseType( gridPart ),
      localCoefficientsMap_( ),
      mapper_( gridPart, localCoefficientsMap_ )
    {}

  private:
    // forbid the copy constructor
    PQ22DLagrangeSpace ( const ThisType & );

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

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::baseFunctionSet(const EntityType &entity) const */
    template< class Entity >
    const BaseFunctionSetType baseFunctionSet ( const Entity &entity ) const
    {
      return BaseFunctionSetType( &localCoefficientsMap_.baseFunctionSet( entity ) );
    }

    /** \brief get dimension of value
        \return int
    **/
    int dimensionOfValue () const
    {
      return dimVal;
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::mapper */
    MapperType &mapper () const
    {
      return mapper_;
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::blockMapper */
    MapperType &blockMapper () const
    {
      return mapper_;
    }

  private:
    LocalCoefficientsMapType localCoefficientsMap_;
    mutable MapperType mapper_;
  };

} // end Dune namespace  

#endif // #if HAVE_DUNE_LOCALFUNCTIONS

#endif // #ifndef  DUNE_FEM_P12DSPACE_HH

/* vim: set sw=2 et: */
