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

#include "dgmapper.hh"
#include "dgbasefunctions.hh"
#include "legendredgbasefunctions.hh"

namespace Dune {

  /** @defgroup DGDSpace Discontinuous Galerkin Space
      @ingroup DiscreteFunctionSpace
  
   DiscreteFunctionSpace for discontinuous functions 
   NOTE: To use this space for adaptive calcuations one has to
   use an index set that is capable for adaptive calculations, e.g
  DGAdaptiveLeafIndexSet and DGAdaptiveLeafGridPart.
  **/
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

    typedef DiscontinuousGalerkinBaseFunctionFactory<
      typename Traits::FunctionSpaceType, polOrd> FactoryType;

    typedef DofManager<typename Traits::GridType> DofManagerType;
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
    //! mapper singleton key 
    typedef MapperSingletonKey< IndexSetType > MapperSingletonKeyType;
    //! mapper factory 
    typedef MapperSingletonFactory< MapperSingletonKeyType , 
              MapperType > MapperSingletonFactoryType;

    //! singleton list of mappers 
    typedef SingletonList< MapperSingletonKeyType , MapperType ,
            MapperSingletonFactoryType > MapperProviderType;

  public:
    //- Constructors and destructors
    /** Constructor */
    DiscontinuousGalerkinSpaceBase(GridPartType& gridPart) :
      BaseType (gridPart),
      gridPart_(gridPart),
      mapper_(0),
      // baseFuncSet_(MaxNumElType, 0),
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
          // GeometryIdentifier::IdentifierType id =
//                       GeometryIdentifier::fromGeo(geomTypes[i]);
        
          //if(baseFuncSet_[id] == 0 )
	  if(baseFuncSet_.find( geomTypes[i] ) == baseFuncSet_.end() )
          {
	    const BaseFunctionSetType* set = & setBaseFuncSetPointer(geomTypes[i]);
            baseFuncSet_[ geomTypes[i] ] = set;
            //baseFuncSet_[id] = & setBaseFuncSetPointer(geomTypes[i]);
            maxNumDofs = std::max(maxNumDofs,set->numBaseFunctions());
          }
        }
      }

      assert( maxNumDofs > 0 );
      MapperSingletonKeyType key(gridPart_.indexSet(),maxNumDofs);
      mapper_ = & MapperProviderType::getObject(key);

      assert( mapper_ );
      assert( mapper_->numDofs() == maxNumDofs );
    }

    //! returns index of sequence in grid sequences 
    int sequence () const { return dm_.sequence(); }
    
    /** Destructor */
    virtual ~DiscontinuousGalerkinSpaceBase () 
    {
      for (unsigned int i = 0; i < baseFuncSet_.size(); ++i) 
      {
	/*
        BaseFunctionSetType * set = (BaseFunctionSetType *) baseFuncSet_[i]; 
        if( set ) removeBaseFuncSetPointer( *set );
        baseFuncSet_[i] = 0;
	*/
      }

      MapperProviderType::removeObject( *mapper_ );
    }
  
    //- Methods
    IteratorType begin() const { return gridPart_.template begin<0>(); }

    IteratorType end() const { return gridPart_.template end<0>(); }

    const GridType& grid() const { return gridPart_.grid(); }

    const IndexSetType& indexSet() const { return gridPart_.indexSet(); }

    //! Return the identifier
    DFSpaceIdentifier type () const 
    {
      return DGSpace_id;
    }
  
    //! Get base function set for a given entity
    // template<class Geometry>
//     const BaseFunctionSetType&
//     subBaseFunctionSet (const Geometry & geo) const 
//     {
//       GeometryIdentifier::IdentifierType id = 
//         GeometryIdentifier::fromGeometry(geo);
//       return this->baseFunctionSet(id);
//     }
  
    template<class Geometry>
    const BaseFunctionSetType&
    subBaseFunctionSet (const Geometry & geo) const 
    {
      
      return this->baseFunctionSet(geo);
    }






    //! Get base function set for a given entity
    // const BaseFunctionSetType&
//     subBaseFunctionSet (const GeometryType & type, bool ) const 
//     {
//       GeometryIdentifier::IdentifierType id = 
//         GeometryIdentifier::fromGeo(type.dim(),type);
//       return this->baseFunctionSet(id);
//     }
  const BaseFunctionSetType&
    subBaseFunctionSet (const GeometryType & type, bool ) const 
    {
      // GeometryIdentifier::IdentifierType id = 
//         GeometryIdentifier::fromGeo(type.dim(),type);
      return this->baseFunctionSet(type);
    }
    //! Get base function set for a given entity
    // template <class Entity>
//     const BaseFunctionSetType&
//     baseFunctionSet (const Entity& en) const 
//     {
//       GeometryIdentifier::IdentifierType id = 
//         GeometryIdentifier::fromGeometry(en.geometry());
//       return this->baseFunctionSet(id);
//     }
    template <class Entity>
    const BaseFunctionSetType&
    baseFunctionSet (const Entity& en) const 
    {
      // GeometryIdentifier::IdentifierType id = 
//         GeometryIdentifier::fromGeometry(en.geometry());
      return this->baseFunctionSet(en.geometry().type());
    }
    //! Get base function set for a given entity
    template <class Entity>
    const BaseFunctionSetType&
    getBaseFunctionSet (const Entity& en) const
    {
      return this->baseFunctionSet( en );
    }
  
    //! Get base function set for a given id of geom type (mainly used by
    //! CombinedSpace) 
   //  const BaseFunctionSetType&
//     baseFunctionSet (const GeometryIdentifier::IdentifierType id) const 
//     {
//       assert(id < (int) baseFuncSet_.size());
//       assert(id >= 0);
//       assert(baseFuncSet_[id]);
//       return *baseFuncSet_[id];
//     }
   const BaseFunctionSetType&
    baseFunctionSet (const GeometryType geo) const 
    {
      //assert(id < (int) baseFuncSet_.size());
      //assert(id >= 0);
      assert(baseFuncSet_.find( geo ) != baseFuncSet_.end());
      return *baseFuncSet_[geo];
    }
    // const BaseFunctionSetType&
//     getBaseFunctionSet (const GeometryIdentifier::IdentifierType id) const
//     {
//       return this->baseFunctionSet(id);
//     }
  
    //! return true if we have continuous discrete functions 
    bool continuous () const
    {
      return false;
    }
  
    //! get global order of space  
    int order () const
    {
      return polOrd;
    }
    
    //! get global order of space  
    int polynomOrder () const DUNE_DEPRECATED
    {
      return order();
    }
  
    //! length of the dof vector  
    //! size knows the correct way to calculate the size of the functionspace
    int size () const 
    {
      return mapper().size();
    }

    //! for given entity map local dof number to global dof number 
    template <class EntityType>
    int mapToGlobal ( EntityType &en, int localNum ) const
    {
      return mapper().mapToGlobal ( en , localNum );
    }

    //! Return dof mapper of the space
    const MapperType& mapper() const 
    {
      assert( mapper_ );
      return *mapper_;
    }

    //! return reference to the spaces grid part
    GridPartType & gridPart () { return gridPart_; }
    //! return reference to the spaces grid part
    const GridPartType & gridPart () const { return gridPart_; }

  protected:
    DiscontinuousGalerkinSpaceBase();
    DiscontinuousGalerkinSpaceBase(const DiscontinuousGalerkinSpaceBase&);
    DiscontinuousGalerkinSpaceBase& operator=(const DiscontinuousGalerkinSpaceBase&);

    template <class EntityType>
    BaseFunctionSetType& setBaseFuncSetPointer(EntityType& en) 
    {
      // calls static method of actual implementation to create set
      return DiscreteFunctionSpaceImp::setBaseFuncSetPointer(en);
    }

    void removeBaseFuncSetPointer(BaseFunctionSetType& set) 
    {
      // calls static method of actual implementation to remove set
      DiscreteFunctionSpaceImp::removeBaseFuncSetPointer(set);
    }

  protected:
    // grid part
    GridPartType& gridPart_;

    // mapper for function space 
    mutable MapperType* mapper_; 

    // vector of base function sets
    typedef std::map < const GeometryType ,const BaseFunctionSetType* > BaseFunctionMapType;
    mutable BaseFunctionMapType baseFuncSet_;
    //std::vector<const BaseFunctionSetType*> baseFuncSet_;

    const DofManagerType & dm_;
  };

  //********************************************************
  // DG Space with orthonormal basis functions 
  //********************************************************
  template <class FunctionSpaceImp, class GridPartImp, int polOrd,
            template<class> class BaseFunctionStorageImp = SimpleStorage >
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
 
    typedef VectorialBaseFunctionSet<FunctionSpaceType, BaseFunctionStorageImp > BaseFunctionSetType;
    typedef DGMapper<IndexSetType, polOrd, DimRange> MapperType;
    
    //! number of base functions 
    enum { localBlockSize = DimRange * 
        DGNumberOfBaseFunctions<polOrd,DimDomain>::numBaseFunctions }; 

  };

  //! A discontinuous Galerkin space
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

    //! scalar space type 
    typedef typename Traits::ScalarFunctionSpaceType ScalarFunctionSpaceType;

    // type of base function factory 
    typedef DiscontinuousGalerkinBaseFunctionFactory<
      ScalarFunctionSpaceType, polOrd> ScalarFactoryType;   

    // type of singleton factory 
    typedef BaseFunctionSetSingletonFactory<GeometryType,BaseFunctionSetType,
                ScalarFactoryType> SingletonFactoryType; 

    // type of singleton list  
    typedef SingletonList< GeometryType, BaseFunctionSetType,
            SingletonFactoryType > SingletonProviderType;

  public:
    //- Constructors and destructors
    /** Constructor */
    DiscontinuousGalerkinSpace(GridPartImp& gridPart) :
    DiscontinuousGalerkinSpaceBase <Traits> (gridPart) {}

    //! get object from singleton list 
    static BaseFunctionSetType& setBaseFuncSetPointer(GeometryType type) 
    {
      return SingletonProviderType::getObject(type);
    }
    //! remove object from singleton list 
    static void removeBaseFuncSetPointer(BaseFunctionSetType& set) 
    {
      SingletonProviderType::removeObject(set);
    }
  };

  template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp = SimpleStorage >
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

    enum { DimRange = FunctionSpaceType::DimRange };
    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
    typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
    typedef typename FunctionSpaceType::ScalarFunctionSpaceType ScalarFunctionSpaceType;

    typedef LegendreDiscontinuousGalerkinSpace<
      FunctionSpaceType, GridPartType, polOrd, BaseFunctionStorageImp> DiscreteFunctionSpaceType;
 
    typedef VectorialBaseFunctionSet<FunctionSpaceType, BaseFunctionStorageImp > BaseFunctionSetType;
    typedef DGMapper<IndexSetType, polOrd, DimRange> MapperType;
  };

  //! A discontinuous Galerkin space
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
    enum { dimRange = Traits::FunctionSpaceType::DimRange };

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
    typedef BaseFunctionSetSingletonFactory<GeometryType,BaseFunctionSetType,
                ScalarFactoryType> SingletonFactoryType; 

    // type of singleton list  
    typedef SingletonList< GeometryType, BaseFunctionSetType,
            SingletonFactoryType > SingletonProviderType;

  public:
    //- Constructors and destructors
    /** Constructor */
    LegendreDiscontinuousGalerkinSpace(GridPartImp& gridPart) :
      DiscontinuousGalerkinSpaceBase<Traits> (gridPart) {}

    //! get base set from singleton list 
    static BaseFunctionSetType& setBaseFuncSetPointer(GeometryType type) 
    {
      return SingletonProviderType::getObject(type);
    }
    //! remove base set 
    static void removeBaseFuncSetPointer(BaseFunctionSetType& set) 
    {
      SingletonProviderType::removeObject(set);
    }
  };
  
} // end namespace Dune 

#endif
