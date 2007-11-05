#ifndef DUNE_DISCONTINUOUSGALERKINSPACE_HH
#define DUNE_DISCONTINUOUSGALERKINSPACE_HH

/*************************************************************
1) static method to create a space<dimRange> from a given Grid
 *************************************************************/

//- system includes 
#include <vector>

// -Dune includes 
#include <dune/grid/common/grid.hh>


//- local includes 
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/common/allgeomtypes.hh>
#include <dune/fem/space/basefunctions/basefunctionsets.hh>
#include <dune/fem/space/basefunctions/basefunctionstorage.hh>
#include <dune/fem/space/basefunctions/basefunctionproxy.hh>

#include "dgmapper.hh"
#include "dgbasefunctions.hh"
#include "legendredgbasefunctions.hh"

namespace Dune {

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
      typename Traits::FunctionSpaceType, polOrd> FactoryType;

    //! type of dof manager 
    typedef DofManager<typename Traits::GridType> DofManagerType;
    //! type of dof manager factory 
    typedef DofManagerFactory<DofManagerType> DofManagerFactoryType;

    //! Dimension of the range vector field
    enum { dimRange = Traits::FunctionSpaceType::DimRange };

    //! The polynom order of the base functions
    enum { PolOrd = polOrd };
    
    //! The polynom order of the base functions
    enum { polynomialOrder = polOrd };
    CompileTimeChecker<(polOrd>=0)> only_use_DGSpace_with_positive_polOrd;

    //! mapper used to implement mapToGlobal 
    typedef typename Traits::MapperType MapperType; 

    //! mapper used to implement mapToGlobal 
    typedef typename Traits::BlockMapperType BlockMapperType; 

    //! mapper singleton key 
    typedef MapperSingletonKey< IndexSetType > MapperSingletonKeyType;
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
    //- Constructors and destructors
    /** \brief Constructor taking grid part */
    explicit DiscontinuousGalerkinSpaceBase(GridPartType& gridPart) :
      BaseType (gridPart),
      mapper_(0),
      blockMapper_(0),
      baseFuncSet_(),
      dm_(DofManagerFactoryType::getDofManager(gridPart.grid()))
    {
      // add index set to list of indexset of dofmanager 
      typedef typename Traits::GridPartType::IndexSetType IndexSetType;
      DofManagerType & dm = 
        DofManagerFactoryType::getDofManager(gridPart.grid());

      dm.addIndexSet(gridPart.grid(), 
                     const_cast<IndexSetType&>(gridPart.indexSet()));

      // create info for all geom types 
      AllGeomTypes<IndexSetType,typename Traits::GridType>
        allGeomTypes(gridPart.indexSet());
      
      int maxNumDofs = -1;
      for(int cd=0; cd<2; ++cd)
      {
        // get types for codim 0 and 1   
        const std::vector<GeometryType>& geomTypes = allGeomTypes.geomTypes(cd);

        // create mappers and base sets for all existing geom types
        for(size_t i=0; i<geomTypes.size(); ++i)
        {
	        if(baseFuncSet_.find( geomTypes[i] ) == baseFuncSet_.end() )
          {
      	    const BaseFunctionSetImp* set = & setBaseFuncSetPointer(geomTypes[i]);
            assert( set );
            baseFuncSet_[ geomTypes[i] ] = set;
            maxNumDofs = std::max(maxNumDofs,set->numBaseFunctions());
          }
        }
      }

      assert( maxNumDofs > 0 );
      MapperSingletonKeyType key(gridPart.indexSet(),maxNumDofs);
      mapper_ = & MapperProviderType::getObject(key);

      assert( mapper_ );
      assert( mapper_->numDofs() == maxNumDofs );
    }

    /** @copydoc DiscreteFunctionSpaceInterface::sequence const */
    int sequence () const { return dm_.sequence(); }
    
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
      if( blockMapper_ ) 
      {
        BlockMapperProviderType::removeObject( *blockMapper_ );
      }
    }
  
    /** @copydoc DiscreteFunctionSpaceInterface::type */
    DFSpaceIdentifier type () const 
    {
      return DGSpace_id;
    }
  
    //! return reference to base functions set according to the geometry's geometry type 
    /**\brief return reference to base functions set according to the geometry's geometry type 
       \param[in] geo
       \return BasefunctipnSetType
    */
    template<class Geometry>
    const BaseFunctionSetType
    subBaseFunctionSet (const Geometry & geo) const 
    {
      return this->baseFunctionSet(geo);
    }
    
    //! return reference to base functions set according to the geometry's geometry type 
    /**\brief return reference to base functions set according to the geometry's geometry type 
       \param[in] type geometry type 
       \param[in] dummy for function overloading 
       \return BasefunctipnSetType
    */
    const BaseFunctionSetType
    subBaseFunctionSet (const GeometryType & type, bool dummy ) const 
    {
      return this->baseFunctionSet(type);
    }
    
    /**\brief return reference to base functions set according to the entity's geometry type 
       \param[in] en entity type 
       \return BaseFunctionSetType
    */
    template <class Entity>
    const BaseFunctionSetType
    baseFunctionSet (const Entity& en) const 
    {
      return this->baseFunctionSet(en.geometry().type());
    }


    /**\brief return reference to base functions set according to the  geometry type 
       \param[in] geomType geometry type 
       \return BasefunctipnSetType
    */
    const BaseFunctionSetType
    baseFunctionSet (const GeometryType geomType) const 
    {
      //assert(baseFuncSet_.find( geomType ) != baseFuncSet_.end());
      //return *baseFuncSet_[geomType];
      assert(baseFuncSet_.find( geomType ) != baseFuncSet_.end());
      return BaseFunctionSetType(baseFuncSet_[geomType]);
    }

    /** @copydoc DiscreteFunctionSpaceInterface::continuous */
    bool continuous () const
    {
      return false;
    }
  
    /** @copydoc DiscreteFunctionSpaceInterface::order */
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
      // only access mapper if really needed 
      if( !blockMapper_ )
      {
        // create mapper with 1 dof for locating the block per element 
        MapperSingletonKeyType key(this->gridPart().indexSet(),1);
        blockMapper_ = & BlockMapperProviderType::getObject(key);
      }
      return *blockMapper_;
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
    mutable BlockMapperType* blockMapper_;

    //! map holding base function sets
    typedef std::map < const GeometryType, const BaseFunctionSetImp* > BaseFunctionMapType;
    //! base function set map 
    mutable BaseFunctionMapType baseFuncSet_;

    //! dof manager 
    const DofManagerType & dm_;
  };



  //********************************************************
  // DG Space with orthonormal basis functions 
  //********************************************************
  template <class FunctionSpaceImp, class GridPartImp, int polOrd,
            template<class> class BaseFunctionStorageImp = CachingStorage >
  class DiscontinuousGalerkinSpace;

  //! Traits class for DiscontinuousGalerkinSpace
  template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp>
  struct DiscontinuousGalerkinSpaceTraits {
    enum { polynomialOrder = polOrd };
    typedef FunctionSpaceImp FunctionSpaceType;
    typedef GridPartImp GridPartType;

    typedef typename GridPartType::GridType GridType;
    typedef typename GridPartType::IndexSetType IndexSetType;
    typedef typename GridPartType::template Codim<0>::IteratorType IteratorType;

    enum { DimRange  = FunctionSpaceType::DimRange };
    enum { DimDomain = FunctionSpaceType::DimDomain };
    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
    typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
    typedef typename FunctionSpaceType::ScalarFunctionSpaceType ScalarFunctionSpaceType;

    typedef DiscontinuousGalerkinSpace<
      FunctionSpaceType, GridPartType, polOrd, BaseFunctionStorageImp > DiscreteFunctionSpaceType;
 
    typedef VectorialBaseFunctionSet<FunctionSpaceType, BaseFunctionStorageImp > BaseFunctionSetImp;
    //typedef BaseFunctionSetImp BaseFunctionSetType;
    typedef SimpleBaseFunctionProxy<BaseFunctionSetImp> BaseFunctionSetType;
    //
    typedef DGMapper<IndexSetType, polOrd, DimRange> MapperType;

    //! mapper for block vector function 
    typedef DGMapper<IndexSetType, polOrd, 1> BlockMapperType;
    
    //! number of base functions * dimRange 
    enum { localBlockSize = DimRange * 
        DGNumberOfBaseFunctions<polOrd,DimDomain>::numBaseFunctions }; 
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
    enum { dimRange = Traits::FunctionSpaceType::DimRange };

    //! Dimension of the domain vector field
    enum { dimDomain = Traits::FunctionSpaceType::DimDomain };

    //! The polynom order of the base functions
    enum { PolOrd = polOrd };
    
    //! The polynom order of the base functions
    enum { polynomialOrder = polOrd };

    //! size of local blocks 
    enum { localBlockSize = Traits::localBlockSize };

    //! dimensionworld of grid 
    enum { dimensionworld = GridType :: dimensionworld };
    //! only use with dimension <= 3 
    CompileTimeChecker<(dimensionworld<=3)> use_Legendre_Spaces_for_higher_worlddims;

    //! scalar space type 
    typedef typename Traits::ScalarFunctionSpaceType ScalarFunctionSpaceType;

    // type of base function factory 
    typedef DiscontinuousGalerkinBaseFunctionFactory<
      ScalarFunctionSpaceType, polOrd> ScalarFactoryType;   

    // type of singleton factory 
    typedef BaseFunctionSetSingletonFactory<const GeometryType,BaseFunctionSetImp,
                ScalarFactoryType> SingletonFactoryType; 

    // type of singleton list  
    typedef SingletonList<const GeometryType, BaseFunctionSetImp,
            SingletonFactoryType > SingletonProviderType;

  public:
    //- Constructors and destructors
    /** Constructor */
    DiscontinuousGalerkinSpace(GridPartImp& gridPart) :
    DiscontinuousGalerkinSpaceBase <Traits> (gridPart) {}

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
    typedef typename GridPartType::IndexSetType IndexSetType;
    typedef typename GridPartType::template Codim<0>::IteratorType IteratorType;

    enum { DimRange  = FunctionSpaceType::DimRange };
    enum { DimDomain = FunctionSpaceType::DimDomain };
    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
    typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
    typedef typename FunctionSpaceType::ScalarFunctionSpaceType ScalarFunctionSpaceType;

    typedef LegendreDiscontinuousGalerkinSpace<
      FunctionSpaceType, GridPartType, polOrd, BaseFunctionStorageImp> DiscreteFunctionSpaceType;
 
    typedef VectorialBaseFunctionSet<FunctionSpaceType, BaseFunctionStorageImp > BaseFunctionSetImp;
    typedef SimpleBaseFunctionProxy< BaseFunctionSetImp > BaseFunctionSetType;

    typedef DGMapper<IndexSetType, polOrd, DimRange> MapperType;

    //! mapper with only one dof 
    typedef DGMapper<IndexSetType, polOrd, 1> BlockMapperType;

    //! number of base functions * dimRange 
    enum { localBlockSize = DimRange * 
        NumLegendreBaseFunctions<polOrd,DimDomain>::numBaseFct }; 

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
    typedef typename GridType :: template Codim<0> :: Entity EntityCodim0Type; 
    //! Index set of space
    typedef typename Traits::IndexSetType IndexSetType;

    typedef LegendreDGBaseFunctionFactory<
      typename Traits::FunctionSpaceType, polOrd> FactoryType;

    //! Dimension of the range vector field
    enum { dimRange  = Traits::FunctionSpaceType::DimRange };
    enum { DimRange  = Traits::FunctionSpaceType::DimRange };
    enum { DimDomain = Traits::FunctionSpaceType::DimDomain };

    //! The polynom order of the base functions
    enum { PolOrd = polOrd };
    
    //! The polynom order of the base functions
    enum { polynomialOrder = polOrd };

    //! type of scalar space 
    typedef typename Traits::ScalarFunctionSpaceType ScalarFunctionSpaceType;
    
    // type of base function factory 
    typedef LegendreDGBaseFunctionFactory<
      ScalarFunctionSpaceType, polOrd> ScalarFactoryType;   

    // type of singleton factory  
    typedef BaseFunctionSetSingletonFactory<GeometryType,BaseFunctionSetImp,
                ScalarFactoryType> SingletonFactoryType; 

    // type of singleton list  
    typedef SingletonList< GeometryType, BaseFunctionSetImp,
            SingletonFactoryType > SingletonProviderType;

    //! size of local blocks 
    enum { localBlockSize = Traits::localBlockSize };

  public:
    //- Constructors and destructors
    /** Constructor */
    LegendreDiscontinuousGalerkinSpace(GridPartImp& gridPart) :
      DiscontinuousGalerkinSpaceBase<Traits> (gridPart) {}

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
  };
  //@}
} // end namespace Dune 

#endif
