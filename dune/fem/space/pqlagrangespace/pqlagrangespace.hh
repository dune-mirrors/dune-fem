#ifndef  DUNE_FEM_P12DSPACE_HH
#define  DUNE_FEM_P12DSPACE_HH

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
#include  <dune/localfunctions/lagrange/p0.hh>
#include  <dune/localfunctions/lagrange/p1.hh>
#include  <dune/localfunctions/lagrange/p2.hh>
#include  <dune/localfunctions/lagrange/pk2d.hh>
#include  <dune/localfunctions/lagrange/pk3d.hh>
#include  <dune/localfunctions/lagrange/q1.hh>
#include  <dune/localfunctions/lagrange/q22d.hh>

namespace Dune
{

  template< class FunctionSpaceImp, class GridPartImp, int polOrder = 1 >
  class PLagrangeSpace;

  template< class FunctionSpaceImp, class GridPartImp, int polOrder = 1 >
  class QLagrangeSpace;



  // PQLocalCoefficientsMap
  // ----------------------

  template< class LocalFiniteElement, bool isSimplex, bool isCube >
  struct PQLocalCoefficientsMap
  {
    typedef typename LocalFiniteElement::Traits::LocalCoefficientsType
      LocalCoefficientsType;

    PQLocalCoefficientsMap ( const LocalFiniteElement &fem )
    : fem_( fem )
    {}

    template< class Entity >
    int operator () ( const Entity &entity ) const
    {
      return 0;
    }

    template< class Topology >
    unsigned int size () const
    {
      int size = 0;
      if( isSimplex )
        size += GenericGeometry::IsSimplex< Topology >::value ? 1 : 0;
      if( isCube )
        size += GenericGeometry::IsCube< Topology >::value ? 1 : 0;
      return size;
    }

    template< class Topology >
    const LocalCoefficientsType &localCoefficients ( const unsigned int i ) const
    {
      return fem_.localCoefficients();
    }

  private:
    const LocalFiniteElement &fem_;
  };



  template< class FunctionSpace, class GridPart, bool isSimplexP, bool isCubeP, int polOrder = 1 >
  struct PQLagrangeSpaceTraits
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

    static const int polynomialOrder = polOrder;

    // TODO: extract this information from the grid(part) and merge
    // PLagrangeSpace and QLagrangeSpace implementations.
    static const bool isSimplex = isSimplexP;
    static const bool isCube    = isCubeP;
  private:
    template<int polOrd, 
             bool isSimplex, bool isCube,
             int dimDomain, bool polOrderGreaterThanOne = (polOrd > 1) >
    struct LocalFiniteElementFactoryTraits
    {
      struct NotImplemented
      {
        NotImplemented() 
        { 
          DUNE_THROW(NotImplemented, "Shape function implementation missing for" << \
                     "selected GridType at dimension " << dimDomain << " and for" << \
                     "polynomial order " << polOrd << ".");
        }
      };
      typedef NotImplemented                                         LocalFiniteElementType;
      static const bool geometryTypeIsFixed = true;
    };

    template<bool isCube, int dimDomain, bool dummy>
    struct LocalFiniteElementFactoryTraits<0, true, isCube, dimDomain, dummy>
    {
      typedef P0LocalFiniteElement< DomainFieldType, RangeFieldType,
                                    dimDomain >                      LocalFiniteElementType;
      static const bool geometryTypeIsFixed = false;
    };

    template<bool isCube, int dimDomain, bool dummy>
    struct LocalFiniteElementFactoryTraits<1, true, isCube, dimDomain, dummy>
    {
      typedef P1LocalFiniteElement< DomainFieldType, RangeFieldType,
                                    dimDomain >                      LocalFiniteElementType;
      static const bool geometryTypeIsFixed = true;
    };

    template<bool isCube, int polOrd>
    struct LocalFiniteElementFactoryTraits<polOrd, true, isCube, 2, true>
    {
      typedef Pk2DLocalFiniteElement< DomainFieldType,
                                      RangeFieldType, polOrd >     LocalFiniteElementType;
      static const bool geometryTypeIsFixed = true;
    };

    template<bool isCube, int polOrd>
    struct LocalFiniteElementFactoryTraits<polOrd, true, isCube, 3, true>
    {
      typedef Pk3DLocalFiniteElement< DomainFieldType,
                                      RangeFieldType, polOrd >     LocalFiniteElementType;
      static const bool geometryTypeIsFixed = true;
    };

    template<int dimDomain, bool dummy>
    struct LocalFiniteElementFactoryTraits<1, false, true, dimDomain, dummy>
    {
      typedef Q1LocalFiniteElement< DomainFieldType, RangeFieldType,
                                    dimDomain >                      LocalFiniteElementType;
      static const bool geometryTypeIsFixed = true;
    };

    template<bool isSimplex, bool dummy>
    struct LocalFiniteElementFactoryTraits<2, isSimplex, true, 2, dummy>
    {
      typedef Q22DLocalFiniteElement< DomainFieldType,
                                      RangeFieldType >               LocalFiniteElementType;
      static const bool geometryTypeIsFixed = true;
    };

    template<class Traits, bool geometryTypeIsFixed = Traits :: geometryTypeIsFixed>
    class LocalFiniteElementFactory;

    template<class Traits>
    class LocalFiniteElementFactory<Traits, true>
    {
    private:
      typedef typename Traits :: LocalFiniteElementType              LocalFiniteElementType;

    public:
      LocalFiniteElementFactory()
        : localFiniteElement_() { }

      const LocalFiniteElementType & getObject() const
      {
        return localFiniteElement_;
      }

    private:
      LocalFiniteElementType                localFiniteElement_;
    };

    template<class Traits>
    class LocalFiniteElementFactory<Traits, false>
    {
    private:
      typedef P0LocalFiniteElement< DomainFieldType, RangeFieldType,
              dimDomain >                        LocalFiniteElementType;

    public:
      LocalFiniteElementFactory()
        : basicType_(isCube ? GeometryType :: cube : GeometryType :: simplex),
        localFiniteElement_(basicType_) { }

      const LocalFiniteElementType & getObject() const
      {
        return localFiniteElement_;
      }

    private:
      typename GeometryType :: BasicType basicType_;
      LocalFiniteElementType             localFiniteElement_;
    };

  public:

    typedef LocalFiniteElementFactoryTraits< polynomialOrder,
                                             isSimplex, isCube,
                                             dimDomain >             LocalFEFactoryTraitsType;

    typedef typename LocalFEFactoryTraitsType
              :: LocalFiniteElementType                              LocalFiniteElementType;

    typedef LocalFiniteElementFactory< LocalFEFactoryTraitsType >    LocalFEFactoryType;

    typedef PQLocalCoefficientsMap< LocalFiniteElementType, isSimplex, isCube >
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



  template< class FunctionSpace, class GridPart, int polOrder = 1 >
  struct PLagrangeSpaceTraits
  : public PQLagrangeSpaceTraits< FunctionSpace, GridPart, true, false, polOrder>
  {
    typedef PLagrangeSpace< FunctionSpace, GridPart, polOrder >      DiscreteFunctionSpaceType;
  };



  template< class FunctionSpace, class GridPart, int polOrder = 1 >
  struct QLagrangeSpaceTraits
  : public PQLagrangeSpaceTraits< FunctionSpace, GridPart, false, true, polOrder>
  {
    typedef QLagrangeSpace< FunctionSpace, GridPart, polOrder >      DiscreteFunctionSpaceType;
  };



  // PQLagrangeSpace
  // ---------------

  /** \class   PQLagrangeSpace
   *  \ingroup LocalFunctionSpaces
   *  \brief   langrange discrete function space for simplex and cube grids
   */
  template< class FunctionSpaceImp, class GridPartImp, int polOrder, template<class,class,int> class SpaceTraits >
  class PQLagrangeSpace
  : public DiscreteFunctionSpaceDefault< SpaceTraits< FunctionSpaceImp, GridPartImp, polOrder > >,
    public GenericDiscreteFunctionSpace
  {
    typedef PQLagrangeSpace< FunctionSpaceImp, GridPartImp, polOrder,
                             SpaceTraits >                           ThisType;
    typedef DiscreteFunctionSpaceDefault
              < SpaceTraits< FunctionSpaceImp, GridPartImp,
                             polOrder > >                            BaseType;

  public:
    //! traits for the discrete function space
    typedef SpaceTraits< FunctionSpaceImp, GridPartImp, polOrder >   Traits;

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

    typedef typename Traits :: LocalFEFactoryType                    LocalFEFactoryType;

    //! size of local blocks
    enum { localBlockSize = Traits :: localBlockSize };

    static const bool isSimplex = Traits :: isSimplex;
    static const bool isCube    = Traits :: isCube;

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
    explicit PQLagrangeSpace ( GridPartType &gridPart )
    : BaseType( gridPart ),
      finiteElementFactory_(),
      finiteElement_( finiteElementFactory_.getObject() ),
      localCoefficientsMap_( finiteElement_ ),
      mapper_( gridPart, localCoefficientsMap_ ),
      baseFunctionSet_( finiteElement_.localBasis(), finiteElement_.type() )
    {}

  private:
    // forbid the copy constructor
    PQLagrangeSpace ( const ThisType & );

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
    LocalFEFactoryType finiteElementFactory_;
    const LocalFiniteElementType &finiteElement_;
    LocalCoefficientsMapType localCoefficientsMap_;
    mutable MapperType mapper_;
    BaseFunctionSetImp baseFunctionSet_;
  };



  template<class FunctionSpaceImp, class GridPartImp, int polOrder>
  class PLagrangeSpace
    : public PQLagrangeSpace<FunctionSpaceImp, GridPartImp, polOrder, PLagrangeSpaceTraits>
  {
  private:
    typedef PQLagrangeSpace< FunctionSpaceImp, GridPartImp, polOrder,
                             PLagrangeSpaceTraits >                  BaseType;
  public:
    explicit PLagrangeSpace ( GridPartImp &gridPart )
      : BaseType(gridPart) {}
  };



  template<class FunctionSpaceImp, class GridPartImp, int polOrder>
  class QLagrangeSpace
    : public PQLagrangeSpace<FunctionSpaceImp, GridPartImp, polOrder, QLagrangeSpaceTraits>
  {
  private:
    typedef PQLagrangeSpace< FunctionSpaceImp, GridPartImp, polOrder,
                             QLagrangeSpaceTraits >                  BaseType;
  public:
    explicit QLagrangeSpace ( GridPartImp &gridPart )
      : BaseType(gridPart) {}
  };

} // end Dune namespace  

#endif // #if HAVE_DUNE_LOCALFUNCTIONS

#endif // #ifndef  DUNE_FEM_P12DSPACE_HH

/* vim: set sw=2 et: */
