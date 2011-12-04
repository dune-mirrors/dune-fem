#ifndef DUNE_DISCONTINUOUSGALERKINSPACE_HH
#define DUNE_DISCONTINUOUSGALERKINSPACE_HH

/*************************************************************
1) static method to create a space<dimRange> from a given Grid
 *************************************************************/

#include <vector>

#include <dune/common/static_assert.hh>

#include <dune/grid/common/grid.hh>

#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/common/allgeomtypes.hh>

#include <dune/fem/space/basefunctions/basefunctionsets.hh>
#include <dune/fem/space/basefunctions/basefunctionstorage.hh>
#include <dune/fem/space/basefunctions/basefunctionproxy.hh>

#include <dune/fem/space/mapper/codimensionmapper.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>

#include <dune/fem/space/dgspace/dgbasefunctions.hh>
#include <dune/fem/space/dgspace/legendredgbasefunctions.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/lagrangespace/basefunctions.hh>
#include <dune/fem/space/common/basesetlocalkeystorage.hh>

namespace Dune
{

  template< class FunctionSpace, unsigned int topologyId,
            unsigned int dim, unsigned int pOrder >
  class LagrangeBaseFunction;

  /** @addtogroup DGDSpace 
   DiscreteFunctionSpace for discontinuous functions 
   NOTE: To use this space for adaptive calcuations one has to
   use an index set that is capable for adaptive calculations, e.g
   DGAdaptiveLeafIndexSet and DGAdaptiveLeafGridPart.
   
   @{

  */
  
  //! A discontinuous Galerkin space base for DGSpaces 
  template <class SpaceImpTraits>
  class DiscontinuousGalerkinSpaceBase : 
    public DiscreteFunctionSpaceDefault <SpaceImpTraits>
  {
    // - Local enums
    enum { DGFSpaceId = 0 };
  
    enum { MaxNumElType = 9 };

    // - Local typedefs
    // the type of this class
    typedef DiscontinuousGalerkinSpaceBase<SpaceImpTraits> 
    ThisType;

    //! the actual space implementation 
    typedef typename SpaceImpTraits::DiscreteFunctionSpaceType 
            DiscreteFunctionSpaceImp;
  public:
    //- Public typedefs

    //! The traits class
    typedef SpaceImpTraits Traits;

    typedef typename Traits::GridPartType GridPartType;

    enum { polOrd = Traits :: polynomialOrder };

    //! Exporting the interface type
    typedef DiscreteFunctionSpaceDefault<Traits> BaseType;
    //! Base function set type
    typedef typename Traits::ShapeFunctionSetType ShapeFunctionSetType;
    //! Base function set type
    typedef typename Traits::BaseFunctionSetType  BaseFunctionSetType;

    // deprecated name 
    typedef ShapeFunctionSetType  BaseFunctionSetImp;

    //! Domain vector type
    typedef typename Traits::DomainType DomainType;
    //! Range vector type
    typedef typename Traits::RangeType RangeType;
    //! Space iterator
    typedef typename Traits::IteratorType IteratorType;
    //! Grid type
    typedef typename Traits::GridType GridType;

    //! Index set of space
    typedef typename Traits::IndexSetType IndexSetType;

    //! factory type for base functions 
    typedef DiscontinuousGalerkinBaseFunctionFactory<
      typename Traits::BaseFunctionSpaceType, polOrd> FactoryType;

    //! Dimension of the range vector field
    enum { dimRange = Traits::FunctionSpaceType::dimRange };

    //! Dimension of the domain vector field
    enum { dimDomain = Traits::FunctionSpaceType::dimDomain };

    //! codimension of space 
    enum { codimension = Traits :: codimension };

    //! The polynom order of the base functions
    enum { polynomialOrder = polOrd };
    dune_static_assert( (polOrd >= 0), "Only use DGSpace with positive polOrd." );

    //! size of local blocks 
    enum { localBlockSize = Traits::localBlockSize };

    //! dimensionworld of grid 
    enum { dimensionworld = GridType :: dimensionworld };
    dune_static_assert( (dimensionworld <= 3), "Use Legendre spaces for higher dimensions." );

    //! type of entity of this space 
    typedef typename GridPartType :: template Codim< codimension > :: EntityType  EntityType;

    //! mapper used to implement mapToGlobal 
    typedef typename Traits::MapperType MapperType; 

    //! mapper used to implement mapToGlobal (for each block of DoFs)
    typedef typename Traits::BlockMapperType BlockMapperType; 

    //! mapper factory 
    typedef CodimensionMapperSingletonFactory< GridPartType, codimension > BlockMapperSingletonFactoryType;

    //! singleton list of mappers 
    typedef SingletonList
      < typename BlockMapperSingletonFactoryType::Key, BlockMapperType, BlockMapperSingletonFactoryType >
      BlockMapperProviderType;

    // type of base function factory 
    typedef typename Traits :: ScalarFactoryType ScalarFactoryType ;

    // type of singleton factory 
    typedef BaseFunctionSetSingletonFactory<const GeometryType, ShapeFunctionSetType,
                ScalarFactoryType> SingletonFactoryType; 

    // type of singleton list  
    typedef SingletonList<const GeometryType, ShapeFunctionSetType,
            SingletonFactoryType > SingletonProviderType;

  public:
    using BaseType :: order ;

    //! default communication interface 
    static const InterfaceType defaultInterface = InteriorBorder_All_Interface;

    //! default communication direction 
    static const CommunicationDirection defaultDirection =  ForwardCommunication;

    //- Constructors and destructors
    /** \brief Constructor taking grid part */
    explicit DiscontinuousGalerkinSpaceBase( GridPartType &gridPart, 
                                             const std::vector<GeometryType>& geomTypes,
                                             const InterfaceType commInterface,
                                             const CommunicationDirection commDirection)
    : BaseType( gridPart , commInterface, commDirection ),
      blockMapper_( BlockMapperProviderType::getObject( gridPart ) ),
      mapper_( blockMapper_ ),
      shapeFunctionSets_()
    {
      // create mappers and base sets for all existing geom types
      for(size_t i=0; i<geomTypes.size(); ++i)
      {
        // insert shape function set for given geometry type 
        shapeFunctionSets_.template insert< SingletonProviderType >( geomTypes[ i ] );
      }

#ifndef NDEBUG
      int maxNumDofs = shapeFunctionSets_.maxSize();

      // check maxNumDofs 
      assert( maxNumDofs > 0 );
      // this should be the same 
      assert( maxNumDofs == mapper().maxNumDofs() );
#endif
    }

    /** \brief Destructor */
    virtual ~DiscontinuousGalerkinSpaceBase () 
    {
      BlockMapperProviderType::removeObject( blockMapper_ );
    }
  
    /** @copydoc Dune::DiscreteFunctionSpaceInterface::type */
    DFSpaceIdentifier type () const 
    {
      return DGSpace_id;
    }
  
    /** \brief obtain the base functions corresponding to an entity
     *
     *  \param[in]  entity  entity for which the base function set is requested
     *
     *  \return a reference to the base function set
     */
    const BaseFunctionSetType
    baseFunctionSet (const EntityType & entity) const 
    {
      return this->baseFunctionSet(entity.type());
    }


    /**\brief return reference to base functions set according to the  geometry type 
       \param[in] geomType geometry type 
       \return BasefunctipnSetType
    */
    const BaseFunctionSetType
    baseFunctionSet (const GeometryType geomType) const 
    {
      return BaseFunctionSetType( &shapeFunctionSets_[ geomType ] );
    }

    /** @copydoc Dune::DiscreteFunctionSpaceInterface::contains */
    bool contains (const int codim) const
    {
      return blockMapper_.contains( codim );
    }
  
    /** @copydoc Dune::DiscreteFunctionSpaceInterface::continuous */
    bool continuous () const
    {
      return false;
    }
  
    /** @copydoc Dune::DiscreteFunctionSpaceInterface::order */
    int order () const
    {
      return polynomialOrder;
    }
    
    //! Return dof mapper of the space
    /** \brief Return dof mapper of the space
        \return MapperType 
    */
    MapperType& mapper() const 
    {
      return mapper_;
    }

    /** \brief Return dof mapper for block located one elements 
    */
    BlockMapperType& blockMapper() const 
    {
      return blockMapper_;
    }

  protected:  
    //! \brief prohibited empty constructor  
    DiscontinuousGalerkinSpaceBase();
    //! \brief prohibited empty copy constructor  
    DiscontinuousGalerkinSpaceBase(const DiscontinuousGalerkinSpaceBase&);
    //! \brief prohibited empty assignment operator   
    DiscontinuousGalerkinSpaceBase& operator=(const DiscontinuousGalerkinSpaceBase&);

  protected:
    //! mapper for blocks 
    BlockMapperType& blockMapper_;

    //! mapper for function space 
    mutable MapperType mapper_; 

    // storage of pointers to shape function sets 
    typedef Fem :: BaseSetLocalKeyStorage< ShapeFunctionSetType > ShapeSetStorageType;

    // storage of pointers to shape function sets 
    mutable ShapeSetStorageType shapeFunctionSets_;
  };


  /////////////////////////////////////////////////
  //! Common traits class for DG spaces 
  /////////////////////////////////////////////////
  template <class FunctionSpaceImp, class GridPartImp, 
            int polOrd, template <class> class BaseFunctionStorageImp>
  struct DiscontinuousGalerkinSpaceTraitsBase  
  {
    enum { polynomialOrder = polOrd };
    typedef FunctionSpaceImp FunctionSpaceType;
    typedef GridPartImp GridPartType;
    typedef typename GridPartType::GridType GridType;

    enum { dimRange  = FunctionSpaceType::dimRange };
    enum { dimDomain = FunctionSpaceType::dimDomain };
    enum { dimLocal  = GridType :: dimension };
    enum { codimension = GridType :: dimensionworld - dimDomain }; 

    typedef typename GridPartType::IndexSetType IndexSetType;
    typedef typename GridPartType::template Codim<codimension>::IteratorType IteratorType;
    typedef typename IteratorType :: Entity  EntityType ;

    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
    typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

    typedef typename ToLocalFunctionSpace< FunctionSpaceType, dimLocal> :: Type 
      BaseFunctionSpaceType;
 
    typedef VectorialBaseFunctionSet<BaseFunctionSpaceType, BaseFunctionStorageImp > ShapeFunctionSetType;
    typedef SimpleBaseFunctionProxy< ShapeFunctionSetType > BaseFunctionSetType;
    
    //! mapper for block vector function 
    typedef CodimensionMapper< GridPartType, codimension > BlockMapperType;

    /** \brief defines type of data handle for communication 
        for this type of space.
    */
    template <class DiscreteFunctionImp,
              class OperationImp = DFCommunicationOperation :: Copy>
    struct CommDataHandle
    {
      //! type of data handle 
      typedef DefaultCommunicationHandler<DiscreteFunctionImp, OperationImp> Type;
      //! type of operation to perform on scatter 
      typedef OperationImp OperationType;
    };
  };


  //********************************************************
  // DG Space with orthonormal basis functions 
  //********************************************************
  template <class FunctionSpaceImp, class GridPartImp, int polOrd,
            template<class> class BaseFunctionStorageImp = CachingStorage >
  class DiscontinuousGalerkinSpace;

  //! Traits class for DiscontinuousGalerkinSpace
  template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp>
  struct DiscontinuousGalerkinSpaceTraits 
  : public DiscontinuousGalerkinSpaceTraitsBase< FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp >
  {
    typedef DiscontinuousGalerkinSpaceTraitsBase< FunctionSpaceImp, GridPartImp, polOrd,
                                                  BaseFunctionStorageImp > BaseType ;

    //! number of base functions * dimRange (use dimLocal here)
    enum { localBlockSize = BaseType :: dimRange * 
        DGNumberOfBaseFunctions<polOrd, BaseType :: dimLocal>::numBaseFunctions }; 

    //! type of DG mapper (based on BlockMapper)
    typedef NonBlockMapper< typename BaseType :: BlockMapperType, localBlockSize > MapperType;

    // type of base function factory 
    typedef DiscontinuousGalerkinBaseFunctionFactory< 
      typename BaseType :: BaseFunctionSpaceType :: ScalarFunctionSpaceType, polOrd> ScalarFactoryType;   

    // type of DG space 
    typedef DiscontinuousGalerkinSpace<
      FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp > DiscreteFunctionSpaceType;
  };

  //! \brief A discontinuous Galerkin space
  template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp >
  class DiscontinuousGalerkinSpace : 
    public DiscontinuousGalerkinSpaceBase 
  <DiscontinuousGalerkinSpaceTraits<FunctionSpaceImp, GridPartImp,polOrd,BaseFunctionStorageImp> >
  {
    // - Local enums
    enum { DGFSpaceId = 1 };

    // - Local typedefs
    // the type of this class
    typedef DiscontinuousGalerkinSpace<FunctionSpaceImp, GridPartImp, polOrd,BaseFunctionStorageImp> 
    ThisType;

  public:
    //- Public typedefs

    //! The traits class
    typedef DiscontinuousGalerkinSpaceTraits<
      FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp> Traits;

    //! Exporting the interface type
    typedef DiscreteFunctionSpaceDefault<Traits> BaseType;

    // for dimensions higher then 3 use LegendreDiscontinuousGalerkinSpace 
    dune_static_assert( (GridPartImp :: GridType :: dimensionworld <= 3), 
        "Use Legendre spaces for higher dimensions." );

  protected:
    typedef DiscontinuousGalerkinSpaceBase <Traits> BaseImpType;

    // set of geometry types 
    typedef AllGeomTypes<typename Traits :: IndexSetType,
                         typename Traits :: GridType>   GeometryTypes;
  public:
    //- Constructors and destructors
    /** Constructor */
    DiscontinuousGalerkinSpace(GridPartImp& gridPart,
                               const InterfaceType commInterface = BaseImpType :: defaultInterface,
                               const CommunicationDirection commDirection = BaseImpType :: defaultDirection ) :
      BaseImpType (gridPart, GeometryTypes(gridPart.indexSet()).geomTypes(BaseImpType :: codimension),
                   commInterface, commDirection)
    {}
  };



  //********************************************************
  // DG Space with Legendre basis functions 
  //********************************************************
  template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp = CachingStorage >
  class LegendreDiscontinuousGalerkinSpace; 
    
  //! Traits class for DiscontinuousGalerkinSpace
  template <class FunctionSpaceImp, class GridPartImp, int polOrd, template
    <class> class BaseFunctionStorageImp >
  struct LegendreDiscontinuousGalerkinSpaceTraits 
  : public DiscontinuousGalerkinSpaceTraitsBase< FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp >
  {
    typedef DiscontinuousGalerkinSpaceTraitsBase< FunctionSpaceImp, GridPartImp, polOrd,
                                                  BaseFunctionStorageImp > BaseType ;

    //! number of base functions * dimRange 
    enum { localBlockSize = BaseType :: dimRange * 
        Fem :: NumLegendreBaseFunctions<polOrd, BaseType :: dimLocal>::numBaseFct }; 

    //! type of DG mapper (based on BlockMapper)
    typedef NonBlockMapper< typename BaseType :: BlockMapperType, localBlockSize > MapperType;

    typedef Fem :: LegendreDGBaseFunctionFactory<
      typename BaseType :: BaseFunctionSpaceType :: ScalarFunctionSpaceType, polOrd> ScalarFactoryType;

    // type of DG space 
    typedef LegendreDiscontinuousGalerkinSpace<
      FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp> DiscreteFunctionSpaceType;
  };

  /** \brief A discontinuous Galerkin space using tensor product base
    functions 
  */
  template <class FunctionSpaceImp, class GridPartImp, int polOrd,
           template<class> class BaseFunctionStorageImp >
  class LegendreDiscontinuousGalerkinSpace : 
    public DiscontinuousGalerkinSpaceBase
      <LegendreDiscontinuousGalerkinSpaceTraits<FunctionSpaceImp,GridPartImp,polOrd,BaseFunctionStorageImp> > 
  {
    // - Local enums
    enum { DGFSpaceId = 2 };
  
    // - Local typedefs
    // the type of this class
    typedef LegendreDiscontinuousGalerkinSpace<FunctionSpaceImp, GridPartImp, polOrd,BaseFunctionStorageImp> 
    ThisType;

  public:
    //- Public typedefs

    //! The traits class
    typedef LegendreDiscontinuousGalerkinSpaceTraits<
      FunctionSpaceImp, GridPartImp, polOrd,BaseFunctionStorageImp> Traits;

    //! Exporting the interface type
    typedef DiscreteFunctionSpaceDefault<Traits> BaseType;

  protected:
    typedef DiscontinuousGalerkinSpaceBase <Traits> BaseImpType;

  public:
    //- Constructors and destructors
    /** Constructor */
    LegendreDiscontinuousGalerkinSpace(GridPartImp& gridPart,
                               const InterfaceType commInterface = BaseImpType :: defaultInterface ,
                               const CommunicationDirection commDirection = BaseImpType :: defaultDirection ) :
      BaseImpType (gridPart, gridPart.indexSet().geomTypes( BaseImpType :: codimension ) ,  
                   commInterface, commDirection ) 
    {
    }

    // return cube for all codims 
    const std::vector<GeometryType>& geomTypes(int codim) const 
    {
      // make sure grid only contains cube elements 
      assert( this->indexSet().geomTypes(codim).size() == 1 );
      assert( this->indexSet().geomTypes(codim)[0].isCube() );
      return this->indexSet().geomTypes(codim); 
    }
  };

  //********************************************************
  // DG Space with Lagrange basis functions 
  //********************************************************
  template <class FunctionSpaceImp, class GridPartImp, int polOrd,
            template<class> class BaseFunctionStorageImp = CachingStorage >
  class LagrangeDiscontinuousGalerkinSpace;

  //! Traits class for DiscontinuousGalerkinSpace
  template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp>
  struct LagrangeDiscontinuousGalerkinSpaceTraits 
  : public DiscontinuousGalerkinSpaceTraitsBase< FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp >
  {
    typedef DiscontinuousGalerkinSpaceTraitsBase< FunctionSpaceImp, GridPartImp, polOrd,
                                                  BaseFunctionStorageImp > BaseType ;

    typedef typename BaseType :: GridType  GridType;

    // define Lagrange generic base function to get numBaseFunctions 
    typedef typename GeometryWrapper< 
          Capabilities :: hasSingleGeometryType < GridType > :: topologyId, 
          GridType :: dimension >::GenericGeometryType   GenericGeometryType;

    typedef GenericLagrangeBaseFunction< 
          typename BaseType :: BaseFunctionSpaceType :: ScalarFunctionSpaceType, 
          GenericGeometryType, polOrd >  GenericBaseFunctionType;

    //! number of base functions * dimRange (use dimLocal here)
    enum { localBlockSize = BaseType :: dimRange * GenericBaseFunctionType :: numBaseFunctions };

    //! type of DG mapper (based on BlockMapper)
    typedef NonBlockMapper< typename BaseType :: BlockMapperType, localBlockSize > MapperType;

    // type of base function factory 
    typedef LagrangeBaseFunctionFactory< 
      typename BaseType :: BaseFunctionSpaceType :: ScalarFunctionSpaceType, GridType :: dimension, polOrd> ScalarFactoryType;   

    // type of DG space 
    typedef LagrangeDiscontinuousGalerkinSpace<
      FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp > DiscreteFunctionSpaceType;
  };

  //! \brief A discontinuous Galerkin space
  template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp >
  class LagrangeDiscontinuousGalerkinSpace : 
    public DiscontinuousGalerkinSpaceBase 
  < LagrangeDiscontinuousGalerkinSpaceTraits<FunctionSpaceImp, GridPartImp,polOrd,BaseFunctionStorageImp> >
  {
    // - Local enums
    enum { DGFSpaceId = 3 };
  
    // - Local typedefs
    // the type of this class
    typedef LagrangeDiscontinuousGalerkinSpace<FunctionSpaceImp, GridPartImp, polOrd,BaseFunctionStorageImp> 
    ThisType;

  public:
    //- Public typedefs

    //! The traits class
    typedef LagrangeDiscontinuousGalerkinSpaceTraits<
      FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp> Traits;

    //! Exporting the interface type
    typedef DiscreteFunctionSpaceDefault<Traits> BaseType;

    //! Index set of space
    typedef typename Traits::IndexSetType IndexSetType;

  protected:
    typedef DiscontinuousGalerkinSpaceBase <Traits> BaseImpType;

    // set of geometry types 
    typedef AllGeomTypes<typename Traits :: IndexSetType,
                         typename Traits :: GridType>   GeometryTypes;
    
  public:
    //- Constructors and destructors
    /** Constructor */
    LagrangeDiscontinuousGalerkinSpace(GridPartImp& gridPart,
                                       const InterfaceType commInterface = BaseImpType :: defaultInterface,
                                       const CommunicationDirection commDirection = BaseImpType :: defaultDirection ) :
      BaseImpType (gridPart, GeometryTypes(gridPart.indexSet()).geomTypes( BaseImpType :: codimension ),
                   commInterface, commDirection)
    {}

    /** @copydoc Dune::DiscreteFunctionSpaceInterface::continuous */
    bool continuous () const
    {
      return false;
    }
  };

  //@}
} // end namespace Dune 

#include <dune/fem/space/dgspace/localdgmassmatrix.hh>
#include <dune/fem/space/dgspace/localrestrictprolong.hh>

#endif
