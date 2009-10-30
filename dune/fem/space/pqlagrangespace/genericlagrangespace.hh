#ifndef  DUNE_FEM_GENERICLAGRANGESPACE_HH
#define  DUNE_FEM_GENERICLAGRANGESPACE_HH

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
#include  <dune/finiteelements/generic/math/field.hh>

#include  <dune/finiteelements/generic/lagrangefiniteelement.hh>
#include  <dune/finiteelements/generic/lagrangebasis/equidistantpoints.hh>
#if HAVE_ALGLIB
#include  <dune/finiteelements/generic/lagrangebasis/lobattopoints.hh>
#endif
#include "pqlagrangespace.hh"

namespace Dune
{

  template< class FunctionSpaceImp, class GridPartImp, int polOrder = 1 >
  class LagrangeSpace;



  template< class FunctionSpace, class GridPart >
  struct GenericLagrangeSpaceTraits
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

    /*
    static const int polynomialOrder = polOrder;
    */

    /*
    // TODO: extract this information from the grid(part) and merge
    // PLagrangeSpace and QLagrangeSpace implementations.
    static const bool isSimplex = isSimplexP;
    static const bool isCube    = isCubeP;
    */
  private:
    template< int dimDomain >
    struct LocalFiniteElementFactoryTraits
    {
      typedef LagrangeLocalFiniteElement< EquidistantPointSet, dimDomain, double, double > LocalFiniteElementType;
      // typedef LagrangeLocalFiniteElement< LobattoPointSet, dimDomain,double,amp::ampf<256>> LocalFiniteElementType;
      // typedef LagrangeLocalFiniteElement< EquidistantPointSet, dimDomain,double,double, GMPField<128>,GMPField<512> > LocalFiniteElementType;
      // typedef LagrangeLocalFiniteElement< EquidistantPointSet, dimDomain,double,double, amp::ampf<128>,amp::ampf<512> > LocalFiniteElementType;
      // typedef LagrangeLocalFiniteElement< EquidistantPointSet, dimDomain,double,double, double,GMPField<512> > LocalFiniteElementType;
      static const bool geometryTypeIsFixed = false;
    };

    template<class Traits, bool geometryTypeIsFixed = Traits :: geometryTypeIsFixed>
    class LocalFiniteElementFactory
    {
    private:
      typedef typename Traits::LocalFiniteElementType LocalFiniteElementType;
      // typedef LagrangeLocalFiniteElement< dimDomain, DomainFieldType, RangeFieldType > LocalFiniteElementType;

    public:
      LocalFiniteElementFactory(unsigned int order)
      : localFiniteElement_(0,order) { }  // simplex grid

      const LocalFiniteElementType & getObject() const
      {
        return localFiniteElement_;
      }

    private:
      LocalFiniteElementType             localFiniteElement_;
    };

  public:

    typedef LocalFiniteElementFactoryTraits < dimDomain >            LocalFEFactoryTraitsType;

    typedef typename LocalFEFactoryTraitsType
              :: LocalFiniteElementType                              LocalFiniteElementType;

    typedef LocalFiniteElementFactory< LocalFEFactoryTraitsType >    LocalFEFactoryType;

    typedef PQLocalCoefficientsMap< LocalFiniteElementType, true, false >
      LocalCoefficientsMapType;
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
  };

  template< class FunctionSpace, class GridPart, int polOrder >
  struct LagrangeSpaceTraits
  : public GenericLagrangeSpaceTraits< FunctionSpace, GridPart >
  {
    typedef LagrangeSpace< FunctionSpace, GridPart, polOrder >      DiscreteFunctionSpaceType;
  };



  // GenericLagrangeSpace
  // ---------

  /** \class   GenericLagrangeSpace
   *  \ingroup LocalFunctionSpaces
   *  \brief   langrange discrete function space for simplex and cube grids
   */
  template< class FunctionSpaceImp, class GridPartImp, class SpaceTraits >
  class GenericLagrangeSpace
  : public DiscreteFunctionSpaceDefault< SpaceTraits >,
    public GenericDiscreteFunctionSpace
  {
    typedef GenericLagrangeSpace< FunctionSpaceImp, GridPartImp,
                                  SpaceTraits >                      ThisType;
    typedef DiscreteFunctionSpaceDefault
              < SpaceTraits >       BaseType;

  public:
    //! traits for the discrete function space
    typedef SpaceTraits            Traits;

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
   
#if 0
    //! maximum polynomial order of functions in this space
    enum { polynomialOrder = -1 };
#endif    

    //! type of the base function set(s)
    typedef typename Traits :: BaseFunctionSetImp                    BaseFunctionSetImp;
    
    typedef typename Traits :: BaseFunctionSetType                   BaseFunctionSetType;

    typedef typename Traits::LocalCoefficientsMapType LocalCoefficientsMapType;

    //! mapper used to implement mapToGlobal
    typedef typename Traits :: MapperType                            MapperType;

    //! local finite element type
    typedef typename Traits :: LocalFiniteElementType                LocalFiniteElementType;

    typedef typename Traits :: LocalFEFactoryType                    LocalFEFactoryType;

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
    GenericLagrangeSpace ( GridPartType &gridPart, unsigned int polOrder )
    : BaseType( gridPart ),
      polynomialOrder(polOrder),
      finiteElementFactory_(polynomialOrder),
      finiteElement_( finiteElementFactory_.getObject() ),
      localCoefficientsMap_( finiteElement_ ),
      mapper_( gridPart, localCoefficientsMap_ ),
      baseFunctionSet_( finiteElement_.localBasis(), finiteElement_.type() )
    {}

  private:
    // forbid the copy constructor
    GenericLagrangeSpace ( const ThisType & );

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

    /** returns the underlying LocalFiniteElement class from dune-localfunctions */
    const LocalFiniteElementType &localFiniteElement () const
    {
      return finiteElement_;
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::baseFunctionSet(const EntityType &entity) const */
    template< class Entity >
    const BaseFunctionSetType baseFunctionSet ( const Entity &entity ) const
    {
      return BaseFunctionSetType( &baseFunctionSet_ );
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
    //! corresponding mapper
    unsigned int polynomialOrder;

    LocalFEFactoryType finiteElementFactory_;
    const LocalFiniteElementType &finiteElement_;
    LocalCoefficientsMapType localCoefficientsMap_;
    mutable MapperType mapper_;
    BaseFunctionSetImp baseFunctionSet_;
  };

  template<class FunctionSpaceImp, class GridPartImp, int polOrder>
  class LagrangeSpace
    : public GenericLagrangeSpace<FunctionSpaceImp, GridPartImp, LagrangeSpaceTraits< FunctionSpaceImp, GridPartImp, polOrder > >
  {
  public:
    typedef LagrangeSpaceTraits< FunctionSpaceImp, GridPartImp, polOrder >  Traits;
    typedef GenericLagrangeSpace< FunctionSpaceImp, GridPartImp,
                                  Traits >                  BaseType;
    //! maximum polynomial order of functions in this space
    enum { polynomialOrder = polOrder };

    explicit LagrangeSpace ( GridPartImp &gridPart )
      : BaseType(gridPart,polOrder) {}
  };
  
} // end Dune namespace  

#endif // #if HAVE_DUNE_GENERICLOCALFUNCTIONS

#endif // #ifndef  DUNE_FEM_GENERICLAGRANGESPACE_HH

/* vim: set sw=2 et: */
