#ifndef  DUNE_FEM_GENERICLAGRANGESPACE_HH
#define  DUNE_FEM_GENERICLAGRNAGESPACE_HH

#if HAVE_DUNE_GENERICLOCALFUNCTIONS

//- Dune includes 
#include <dune/common/misc.hh>
#include <dune/grid/common/grid.hh>

//- Dune-Fem includes 
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/basefunctionfactory.hh>
#include <dune/fem/space/common/genericdofmapper.hh>
#include <dune/fem/space/basefunctions/basefunctionstorage.hh>
#include <dune/fem/space/basefunctions/basefunctionsets.hh>
#include <dune/fem/space/basefunctions/genericbasefunctionsets.hh>
#include <dune/fem/space/basefunctions/basefunctionproxy.hh>
#include <dune/fem/space/lagrangespace/lagrangedatahandle.hh>

//- Dune-Localfunctions include
#include  <dune/finiteelements/common/field.hh>

#include  <dune/finiteelements/lagrangefiniteelement.hh>
#include  <dune/finiteelements/lagrangebasis/equidistantpoints.hh>
#if HAVE_ALGLIB
#include  <dune/finiteelements/lagrangebasis/lobattopoints.hh>
#endif
#include "pqlagrangespace.hh"

namespace Dune
{

  template< class FunctionSpaceImp,
            class GridPartImp, int polOrder = 1 >
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
      typedef LagrangeLocalFiniteElement<EquidistantCoefficients,dimDomain,double,double> LocalFiniteElementType;
      // typedef LobattoLocalFiniteElement<dimDomain,double,double> LocalFiniteElementType;
      // typedef LagrangeLocalFiniteElement<dimDomain,double,double, GMPField<128>,GMPField<512> > LocalFiniteElementType;
      // typedef LagrangeLocalFiniteElement<dimDomain,double,double, amp::ampf<128>,amp::ampf<512> > LocalFiniteElementType;
      // typedef LagrangeLocalFiniteElement<dimDomain,double,double, double,GMPField<512> > LocalFiniteElementType;
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

    typedef GenericDofMapper< GridPartType >                         MapperType;
    typedef MapperType                                               BlockMapperType;

    // implementation of basefunction set 
    typedef GenericBaseFunctionSet< LocalFiniteElementType,
                                    FunctionSpaceType >              BaseFunctionSetImp;

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
  struct LagrangeSpaceTraits
  : public GenericLagrangeSpaceTraits< FunctionSpace, GridPart >
  {
    typedef LagrangeSpace< FunctionSpace, GridPart, polOrder >      DiscreteFunctionSpaceType;
  };

#if 0
  // GenericMapperSingletonKey
  // -------------------------

  //! Key for Mapper singleton list 
  template< class GridPartImp, class FiniteElementImp >
  struct GenericMapperSingletonKey 
  {
    //! constructor taking index set and numDofs 
    GenericMapperSingletonKey(const GridPartImp & gridPart,
                              const FiniteElementImp & finiteElement,
                              const int polOrd)
      : gridPart_(gridPart),  finiteElement_(finiteElement) , polOrd_(polOrd) 
    {}
    //! copy constructor 
    GenericMapperSingletonKey(const GenericMapperSingletonKey &org) 
      : gridPart_(org.gridPart_), finiteElement_(org.finiteElement_), polOrd_(org.polOrd_)
    {}
    //! returns true if indexSet pointer and numDofs are equal 
    bool operator == (const GenericMapperSingletonKey & otherKey) const 
    {
      return ((&(gridPart_.indexSet()) == &(otherKey.gridPart().indexSet())) 
	      && (polOrd_ == otherKey.polOrd_));
    }

    //! return reference to index set 
    const GridPartImp & gridPart() const { return gridPart_; }
    //! return local finite element object
    const FiniteElementImp& finiteElement() const { return finiteElement_; }

  private:
    const GridPartImp      &gridPart_;
    const FiniteElementImp &finiteElement_;
    const int               polOrd_;
  };



  // GenericMapperSingletonFactory
  // -----------------------------

  // Factory class for SingletonList to tell how objects are created and
  // how compared.
  template< class Key, class Object, bool isSimplex=true, bool isCube=false >
  class GenericMapperSingletonFactory;

  template< class GridPart, class LocalFiniteElement, class Object, bool isSimplex, bool isCube >
  class GenericMapperSingletonFactory< GenericMapperSingletonKey< GridPart, LocalFiniteElement >, Object, isSimplex, isCube >
  {
    struct LocalCoefficientsProvider
    {
      typedef typename LocalFiniteElement::Traits::LocalCoefficientsType
        LocalCoefficientsType;

      LocalCoefficientsProvider ( const LocalFiniteElement &fem )
      : fem_( fem )
      {}

      template< class Topology >
      unsigned int size () const
      {
        int size = 0;
        if( isSimplex ) {
          size+=GenericGeometry::IsSimplex< Topology >::value ? 1 : 0;
        }
        if( isCube ) {
          size+=GenericGeometry::IsCube< Topology >::value ? 1 : 0;
        }
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

  public:
    typedef GenericMapperSingletonKey< GridPart, LocalFiniteElement > KeyType;

    static Object *createObject( const KeyType &key )
    {
      LocalCoefficientsProvider localCoefficientsProvider( key.finiteElement() );
      return new Object( key.gridPart(), localCoefficientsProvider );
    }
    
    static void deleteObject( Object *obj )
    {
      delete obj;
    }
  };
#endif


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

    //! mapper singleton key
    typedef GenericMapperSingletonKey< GridPartType,
                                       LocalFiniteElementType >      MapperSingletonKeyType;

    //! mapper factory
    typedef GenericMapperSingletonFactory< MapperSingletonKeyType,
                                           MapperType >              MapperSingletonFactoryType;

    //! singleton list of mappers
    typedef SingletonList< MapperSingletonKeyType, MapperType,
                           MapperSingletonFactoryType >              MapperProviderType;

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
     *  \param[in]  gridPart       grid part for the Lagrange space
     *  \param[in]  commInterface  communication interface to use (optional)
     *  \param[in]  commDirection  communication direction to use (optional)
     */
    explicit GenericLagrangeSpace ( GridPartType &gridPart, unsigned int polOrder )
    : BaseType( gridPart ),
      polynomialOrder(polOrder),
      mapper_( 0 ),
      finiteElementFactory_(polynomialOrder),
      finiteElement_( finiteElementFactory_.getObject() ),
      genericBaseFunctionSet_( finiteElement_ ),
      baseFunctionSet_( &genericBaseFunctionSet_ )
    {
      MapperSingletonKeyType key( gridPart, finiteElement_, polynomialOrder );
      mapper_ = &MapperProviderType::getObject( key );
      assert( mapper_ != 0 );
    }

  private:
    // forbid the copy constructor
    GenericLagrangeSpace ( const ThisType & );

  public:
    /** \brief Destructor (freeing mapper)
        \return 
    **/
    ~GenericLagrangeSpace ()
    {
      MapperProviderType::removeObject( *mapper_ );
    }

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
/*      assert( entity.type().isSimplex() && (Entity::dimension == 2) );*/
      return baseFunctionSet_;
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
      assert( mapper_ != 0 );
      return *mapper_;
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::blockMapper */
    MapperType &blockMapper () const
    {
      assert( mapper_ != 0 );
      return *mapper_;
    }

  private:
    //! corresponding mapper
    unsigned int polynomialOrder;
    MapperType *mapper_;

    LocalFEFactoryType            finiteElementFactory_;
    const LocalFiniteElementType &finiteElement_;
    BaseFunctionSetImp            genericBaseFunctionSet_;
    BaseFunctionSetType           baseFunctionSet_;
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

#endif // #ifndef  _HH

/* vim: set sw=2 et: */
