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

#include <dune/fem/space/dgspace/dgmapper.hh>
#include <dune/fem/space/dgspace/dgbasefunctions.hh>
#include <dune/fem/space/dgspace/legendredgbasefunctions.hh>
#include <dune/fem/space/dgspace/dgdatahandle.hh>

namespace Dune
{

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
    enum { DGFSpaceId = 1 };
  
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
    typedef typename Traits::BaseFunctionSetImp  BaseFunctionSetImp;
    //! Base function set type
    typedef typename Traits::BaseFunctionSetType BaseFunctionSetType;
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
    enum { dimRange = Traits::FunctionSpaceType :: dimRange };

    //! codimension of space 
    enum { codimension = Traits :: codimension };

    //! type of entity of this space 
    typedef typename GridType :: template Codim< codimension > :: Entity  EntityType;

    //! The polynom order of the base functions (deprecated)
    enum { PolOrd = polOrd };
    
    //! The polynom order of the base functions
    enum { polynomialOrder = polOrd };
    dune_static_assert( (polOrd >= 0), "Only use DGSpace with positive polOrd." );

    //! mapper used to implement mapToGlobal 
    typedef typename Traits::MapperType MapperType; 

    //! mapper used to implement mapToGlobal 
    typedef typename Traits::BlockMapperType BlockMapperType; 

    //! mapper singleton key 
    typedef MapperSingletonKey< GridPartType > MapperSingletonKeyType;
    //! mapper factory 
    typedef MapperSingletonFactory< MapperSingletonKeyType , 
              MapperType > MapperSingletonFactoryType;

    //! singleton list of mappers 
    typedef SingletonList< MapperSingletonKeyType , MapperType ,
            MapperSingletonFactoryType > MapperProviderType;

    //! mapper factory 
    typedef MapperSingletonFactory< MapperSingletonKeyType , 
              BlockMapperType > BlockMapperSingletonFactoryType;

    //! singleton list of mappers 
    typedef SingletonList< MapperSingletonKeyType , BlockMapperType ,
            BlockMapperSingletonFactoryType > BlockMapperProviderType;

  public:
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
      mapper_( 0 ),
      blockMapper_( BlockMapperProviderType::getObject(
            MapperSingletonKeyType (this->gridPart(),1) )),
      baseFuncSet_()
    {
      int maxNumDofs = -1;

      // create mappers and base sets for all existing geom types
      for(size_t i=0; i<geomTypes.size(); ++i)
      {
        // get geometry type 
        const GeometryType type = geomTypes[i];

        if(baseFuncSet_.find( type ) == baseFuncSet_.end())
        {
          const BaseFunctionSetImp* set = & setBaseFuncSetPointer( type );
          assert( set );
          baseFuncSet_[ type ] = set;
          maxNumDofs = std::max(maxNumDofs,set->numBaseFunctions());
        }
      }

      // create mapper 
      assert( maxNumDofs > 0 );
      {
        MapperSingletonKeyType key( gridPart, maxNumDofs );
        mapper_ = & MapperProviderType::getObject(key);
      }
      assert( mapper_ );
      assert( mapper_->maxNumDofs() == maxNumDofs );
    }

    /** \brief Destructor */
    virtual ~DiscontinuousGalerkinSpaceBase () 
    {
      typedef typename BaseFunctionMapType :: iterator iterator;
      iterator end = baseFuncSet_.end();
      for (iterator it = baseFuncSet_.begin(); it != end; ++it)
      {
        BaseFunctionSetImp * set = (BaseFunctionSetImp *) (*it).second; 
        if( set ) removeBaseFuncSetPointer( *set );
      }

      MapperProviderType::removeObject( *mapper_ );
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
      assert(baseFuncSet_.find( geomType ) != baseFuncSet_.end());
      return BaseFunctionSetType(baseFuncSet_[geomType]);
    }

    /** @copydoc Dune::DiscreteFunctionSpaceInterface::contains */
    bool contains (const int codim) const
    {
      return (codim == codimension);
    }
  
    /** @copydoc Dune::DiscreteFunctionSpaceInterface::continuous */
    bool continuous () const
    {
      return false;
    }
  
    /** @copydoc Dune::DiscreteFunctionSpaceInterface::order */
    int order () const
    {
      return polOrd;
    }
    
    //! Return dof mapper of the space
    /** \brief Return dof mapper of the space
        \return MapperType 
    */
    MapperType& mapper() const 
    {
      assert( mapper_ );
      return *mapper_;
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

    /** \brief return BaseFunctionSetPointer 
        \param en entity for which base function set is collected 
        \return BaseFunctionSetImp pointer 
    */
    template <class EntityType>
    BaseFunctionSetImp& setBaseFuncSetPointer(EntityType& en) 
    {
      // calls static method of actual implementation to create set
      return DiscreteFunctionSpaceImp::setBaseFuncSetPointer(en);
    }

    /** \brief remove BaseFunctionSetPointer in singleton list (if no
         other references exist, pointer is deleted  
        \param set pointer to base function set 
    */
    void removeBaseFuncSetPointer(BaseFunctionSetImp& set) 
    {
      // calls static method of actual implementation to remove set
      DiscreteFunctionSpaceImp::removeBaseFuncSetPointer(set);
    }

  protected:
    //! mapper for function space 
    MapperType* mapper_; 
    // mapper for blocks 
    BlockMapperType& blockMapper_;

    //! map holding base function sets
    typedef std::map < const GeometryType, const BaseFunctionSetImp* > BaseFunctionMapType;
    //! base function set map 
    mutable BaseFunctionMapType baseFuncSet_;
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
  {
    enum { polynomialOrder = polOrd };
    typedef FunctionSpaceImp FunctionSpaceType;
    typedef GridPartImp GridPartType;
    typedef typename GridPartType::GridType GridType;

    enum { dimRange  = FunctionSpaceType::dimRange };
    enum { dimDomain = FunctionSpaceType::dimDomain };
    enum { codimension = GridType :: dimensionworld - dimDomain }; 

    typedef typename GridPartType::IndexSetType IndexSetType;
    typedef typename GridPartType::template Codim<codimension>::IteratorType IteratorType;

    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
    typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

    typedef DiscontinuousGalerkinSpace<
      FunctionSpaceType, GridPartType, polOrd, BaseFunctionStorageImp > DiscreteFunctionSpaceType;

    enum { dimLocal = GridType :: dimension };

    typedef typename ToLocalFunctionSpace< FunctionSpaceType, dimLocal> :: Type 
      BaseFunctionSpaceType;
 
    typedef VectorialBaseFunctionSet<BaseFunctionSpaceType, BaseFunctionStorageImp > BaseFunctionSetImp;
    typedef SimpleBaseFunctionProxy<BaseFunctionSetImp> BaseFunctionSetType;
    
    //! type of DG mapper 
    typedef DGMapper< GridPartType, polOrd, dimRange > MapperType;

    //! mapper for block vector function 
    typedef DGMapper< GridPartType, polOrd, 1 > BlockMapperType;
    
    //! number of base functions * dimRange 
    enum { localBlockSize = dimRange * 
        DGNumberOfBaseFunctions<polOrd,dimDomain>::numBaseFunctions }; 
    
    /** \brief defines type of data handle for communication 
        for this type of space.
    */
    template <class DiscreteFunctionImp,
              class OperationImp = DFCommunicationOperation :: Copy>
    struct CommDataHandle
    {
      //! type of data handle 
      typedef DGCommunicationHandler<DiscreteFunctionImp, OperationImp> Type;
      //! type of operation to perform on scatter 
      typedef OperationImp OperationType;
    };
  };

  //! \brief A discontinuous Galerkin space
  template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp >
  class DiscontinuousGalerkinSpace : 
    public DiscontinuousGalerkinSpaceBase 
  <DiscontinuousGalerkinSpaceTraits<FunctionSpaceImp, GridPartImp,polOrd,BaseFunctionStorageImp> >
  {
    // - Local enums
    enum { DGFSpaceId = 1 };
  
    enum { MaxNumElType = 9 };

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
    typedef typename Traits::BaseFunctionSetImp BaseFunctionSetImp;
    //! Base function set type
    typedef typename Traits::BaseFunctionSetType BaseFunctionSetType;
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

    //! Dimension of the range vector field
    enum { dimRange = Traits::FunctionSpaceType::dimRange };

    //! Dimension of the domain vector field
    enum { dimDomain = Traits::FunctionSpaceType::dimDomain };

    //! codimension of space 
    enum { codimension = Traits :: codimension };

    //! The polynom order of the base functions (deprecated)
    enum { PolOrd = polOrd };
    
    //! The polynom order of the base functions
    enum { polynomialOrder = polOrd };

    //! size of local blocks 
    enum { localBlockSize = Traits::localBlockSize };

    //! dimensionworld of grid 
    enum { dimensionworld = GridType :: dimensionworld };
    dune_static_assert( (dimensionworld <= 3), "Use Legendre spaces for higher dimensions." );

    //! scalar space type 
    typedef typename Traits::BaseFunctionSpaceType  BaseFunctionSpaceType;

    // type of base function factory 
    typedef DiscontinuousGalerkinBaseFunctionFactory< 
      typename BaseFunctionSpaceType :: ScalarFunctionSpaceType, polOrd> ScalarFactoryType;   

    // type of singleton factory 
    typedef BaseFunctionSetSingletonFactory<const GeometryType,BaseFunctionSetImp,
                ScalarFactoryType> SingletonFactoryType; 

    // type of singleton list  
    typedef SingletonList<const GeometryType, BaseFunctionSetImp,
            SingletonFactoryType > SingletonProviderType;

  protected:
    typedef DiscontinuousGalerkinSpaceBase <Traits> BaseImpType;

    // set of geometry types 
    typedef AllGeomTypes<IndexSetType,typename Traits::GridType> GeometryTypes;
    
  public:
    //- Constructors and destructors
    /** Constructor */
    DiscontinuousGalerkinSpace(GridPartImp& gridPart,
                               const InterfaceType commInterface = BaseImpType :: defaultInterface,
                               const CommunicationDirection commDirection = BaseImpType :: defaultDirection ) :
      BaseImpType (gridPart, GeometryTypes(gridPart.indexSet()).geomTypes(codimension),
                   commInterface, commDirection)
    {}

    /** \brief ! get object from singleton list 
      \param[in] type 
      \return BasefunctionSetImp
    */
    static BaseFunctionSetImp& setBaseFuncSetPointer(const GeometryType type) 
    {
      return SingletonProviderType::getObject(type);
    }
    //! remove object from singleton list 
    static void removeBaseFuncSetPointer(BaseFunctionSetImp& set) 
    {
      SingletonProviderType::removeObject(set);
    }
  };

  template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp = CachingStorage >
  class LegendreDiscontinuousGalerkinSpace; 
    
  //! Traits class for DiscontinuousGalerkinSpace
  template <class FunctionSpaceImp, class GridPartImp, int polOrd, template
    <class> class BaseFunctionStorageImp >
  struct LegendreDiscontinuousGalerkinSpaceTraits {
    typedef FunctionSpaceImp FunctionSpaceType;
    typedef GridPartImp GridPartType;
    enum { polynomialOrder = polOrd };

    typedef typename GridPartType::GridType GridType;
    enum { dimRange  = FunctionSpaceType::dimRange };
    enum { dimDomain = FunctionSpaceType::dimDomain };
    enum { codimension = GridType :: dimensionworld - dimDomain }; 

    typedef typename GridPartType::IndexSetType IndexSetType;
    typedef typename GridPartType::template Codim<codimension>::IteratorType IteratorType;

    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
    typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

    typedef LegendreDiscontinuousGalerkinSpace<
      FunctionSpaceType, GridPartType, polOrd, BaseFunctionStorageImp> DiscreteFunctionSpaceType;
 
    enum { dimLocal = GridType :: dimension };

    typedef typename ToLocalFunctionSpace< FunctionSpaceType, dimLocal> :: Type 
      BaseFunctionSpaceType;
 
    typedef VectorialBaseFunctionSet<BaseFunctionSpaceType, BaseFunctionStorageImp > BaseFunctionSetImp;
    typedef SimpleBaseFunctionProxy< BaseFunctionSetImp > BaseFunctionSetType;

    typedef DGMapper< GridPartType, polOrd, dimRange > MapperType;

    //! mapper with only one dof 
    typedef DGMapper< GridPartType, polOrd, 1 > BlockMapperType;

    //! number of base functions * dimRange 
    enum { localBlockSize = dimRange * 
        NumLegendreBaseFunctions<polOrd,dimDomain>::numBaseFct }; 

    /** \brief defines type of data handle for communication 
        for this type of space.
    */
    template <class DiscreteFunctionImp, 
              class OperationImp = DFCommunicationOperation :: Copy>
    struct CommDataHandle
    {
      //! type of data handle 
      typedef DGCommunicationHandler<DiscreteFunctionImp,OperationImp> Type;
      //! type of operation to perform on scatter 
      typedef OperationImp OperationType;
    };
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
  
    enum { MaxNumElType = 9 };

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
    //! Base function set type
    typedef typename Traits::BaseFunctionSetImp BaseFunctionSetImp;
    //! Base function set type
    typedef typename Traits::BaseFunctionSetType BaseFunctionSetType;
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

    typedef LegendreDGBaseFunctionFactory<
      typename Traits::FunctionSpaceType, polOrd> FactoryType;

    //! Dimension of the range vector field
    enum { dimRange  = Traits::FunctionSpaceType::dimRange };
    enum { dimDomain = Traits::FunctionSpaceType::dimDomain };

    //! codimension of space 
    enum { codimension = Traits :: codimension };

    //! The polynom order of the base functions (deprecated)
    enum { PolOrd = polOrd };
    
    //! The polynom order of the base functions
    enum { polynomialOrder = polOrd };

    //! type of scalar space 
    typedef typename Traits::BaseFunctionSpaceType BaseFunctionSpaceType;
    
    // type of base function factory 
    typedef LegendreDGBaseFunctionFactory< 
        typename BaseFunctionSpaceType :: ScalarFunctionSpaceType, 
        polOrd
          > ScalarFactoryType;   

    // type of singleton factory  
    typedef BaseFunctionSetSingletonFactory<GeometryType,BaseFunctionSetImp,
                ScalarFactoryType> SingletonFactoryType; 

    // type of singleton list  
    typedef SingletonList< GeometryType, BaseFunctionSetImp,
            SingletonFactoryType > SingletonProviderType;

    //! size of local blocks 
    enum { localBlockSize = Traits::localBlockSize };

  protected:
    typedef DiscontinuousGalerkinSpaceBase <Traits> BaseImpType;
  public:
    //- Constructors and destructors
    /** Constructor */
    LegendreDiscontinuousGalerkinSpace(GridPartImp& gridPart,
                               const InterfaceType commInterface = BaseImpType :: defaultInterface ,
                               const CommunicationDirection commDirection = BaseImpType :: defaultDirection ) :
      BaseImpType (gridPart, gridPart.indexSet().geomTypes(codimension) ,  
                   commInterface, commDirection ) 
    {
    }

    //! get base set from singleton list 
    static BaseFunctionSetImp& setBaseFuncSetPointer(GeometryType type) 
    {
      return SingletonProviderType::getObject(type);
    }
    //! remove base set 
    static void removeBaseFuncSetPointer(BaseFunctionSetImp& set) 
    {
      SingletonProviderType::removeObject(set);
    }

    // return cueb for all codims 
    const std::vector<GeometryType>& geomTypes(int codim) const 
    {
      // make sure grid only contains cube elements 
      assert( this->indexSet().geomTypes(codim).size() == 1 );
      assert( this->indexSet().geomTypes(codim)[0].isCube() );
      return this->indexSet().geomTypes(codim); 
    }
  };
  //@}
} // end namespace Dune 

#include "dgadaptmanager.hh"

#endif
