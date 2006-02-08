#ifndef DUNE_DISCONTINUOUSGALERKINSPACE_HH
#define DUNE_DISCONTINUOUSGALERKINSPACE_HH

#include <vector>

#include <dune/grid/common/grid.hh>

#include <dune/fem/dofmanager.hh>
#include <dune/fem/common/discretefunctionspace.hh>
//#include <dune/fem/common/basefunctionsets.hh>
#include "../basefunctions/basefunctionsets.hh"
#include "../basefunctions/basefunctionstorage.hh"

#include "dgmapper.hh"
#include "dgbasefunctions.hh"
#include "legendredgbasefunctions.hh"

namespace Dune {

  //**********************************************************************
  //
  //!  DiscreteFunctionSpace for discontinuous functions 
  //  
  //**********************************************************************
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

    //! Dimension of the range vector field
    enum { dimRange = Traits::FunctionSpaceType::DimRange };

    //! The polynom order of the base functions
    enum { PolOrd = polOrd };
    
    //! The polynom order of the base functions
    enum { polynomialOrder = polOrd };

  public:
    //- Constructors and destructors
    /** Constructor */
    DiscontinuousGalerkinSpaceBase(GridPartType& gridPart) :
      BaseType (DGFSpaceId),
      gridPart_(gridPart),
      mapper_(0),
      baseFuncSet_(MaxNumElType, 0)
    {
      // add index set to list of indexset of dofmanager 
      typedef DofManager<typename Traits::GridType> DofManagerType;
      typedef DofManagerFactory<DofManagerType> DofManagerFactoryType;
      typedef typename Traits::GridPartType::IndexSetType IndexSetType;
      DofManagerType & dm = 
        DofManagerFactoryType::getDofManager(gridPart.grid());
      dm.addIndexSet(gridPart.grid(), 
                     const_cast<IndexSetType&>(gridPart.indexSet()));
  
      typedef typename GridType :: template Codim<0>::Entity EntityType;
      // search the macro grid for different element types 
      IteratorType endit  = gridPart.template end<0>();
      for(IteratorType it = gridPart.template begin<0>(); it != endit; ++it) {
        GeometryType geo = (*it).geometry().type(); // Hack
        int dimension = static_cast<int>( EntityType::mydimension);
        GeometryIdentifier::IdentifierType id = 
          GeometryIdentifier::fromGeo(dimension, geo);
        if(baseFuncSet_[id] == 0 ) {
          baseFuncSet_[id] = setBaseFuncSetPointer(*it);
          mapper_ = 
            new typename Traits::MapperType(const_cast<IndexSetType&>(gridPart_.indexSet()),
                                   baseFuncSet_[id]->numBaseFunctions());
        }
      }
    }

    
    /** Destructor */
    virtual ~DiscontinuousGalerkinSpaceBase () {
      for (unsigned int i = 0; i < baseFuncSet_.size(); ++i) {
        delete baseFuncSet_[i];
        baseFuncSet_[i] = 0;
      }
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
    template <class Entity>
    BaseFunctionSetType&
    getBaseFunctionSet (const Entity& en) const {
      GeometryType geom = en.geometry().type();
      int dimension = static_cast<int>(Entity::mydimension);
      assert(GeometryIdentifier::fromGeo(dimension,geom)
             <(int) baseFuncSet_.size());
      assert(GeometryIdentifier::fromGeo(dimension, geom) >= 0);
      assert(baseFuncSet_[GeometryIdentifier::fromGeo(dimension, geom)]);
      
      return *baseFuncSet_[GeometryIdentifier::fromGeo(dimension, geom)];
    }
  
    //! return true if we have continuous discrete functions 
    bool continuous () const
    {
      return false;
    }
  
    //! get maximal global polynom order 
    int polynomOrder () const
    {
      return polOrd;
    }
  
    //! length of the dof vector  
    //! size knows the correct way to calculate the size of the functionspace
    int size () const 
    {
      return mapper_->size();
    }

    //! for given entity map local dof number to global dof number 
    template <class EntityType>
    int mapToGlobal ( EntityType &en, int localNum ) const
    {
      return mapper_->mapToGlobal ( en , localNum );
    }

    //! Return dof mapper of the space
    const typename Traits::MapperType& mapper() const 
    {
      return *mapper_;
    }

    //! default for polOrd 0
    template <class EntityType> 
    bool evaluateLocal(int baseFunc, const EntityType &en, 
                       const DomainType &local, RangeType & ret) const 
    {
      enum { dim = EntityType::dimension };
      const BaseFunctionSetType & baseSet = getBaseFunctionSet(en);  
      baseSet.eval(baseFunc, local, ret);
      return true;
    }
    
    //! default for polOrd 0
    template <class EntityType, class QuadratureType> 
    bool evaluateLocal ( int baseFunc, const EntityType &en, 
                         const QuadratureType &quad, 
                         int quadPoint, RangeType & ret) const 
    {
      enum { dim = EntityType::dimension };
      const BaseFunctionSetType & baseSet = getBaseFunctionSet(en);  
      
      baseSet.eval(baseFunc, quad, quadPoint, ret);
      return true;
    }
 
  protected:
    DiscontinuousGalerkinSpaceBase();
    DiscontinuousGalerkinSpaceBase(const DiscontinuousGalerkinSpaceBase&);
    DiscontinuousGalerkinSpaceBase& operator=(const DiscontinuousGalerkinSpaceBase&);

    template <class EntityType>
    BaseFunctionSetType* setBaseFuncSetPointer(EntityType& en) 
    {
      // calls static method of actual implementation to create set
      return DiscreteFunctionSpaceImp::setBaseFuncSetPointer(en);
    }

  protected:
    // grid part
    GridPartType& gridPart_;

    // mapper for function space 
    mutable typename Traits::MapperType* mapper_; 

    // vector of base function sets
    std::vector<BaseFunctionSetType*> baseFuncSet_;
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

    enum { DimRange = FunctionSpaceType::DimRange };
    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
    typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

    typedef DiscontinuousGalerkinSpace<
      FunctionSpaceType, GridPartType, polOrd, BaseFunctionStorageImp > DiscreteFunctionSpaceType;
 
    typedef VecBaseFunctionSet<FunctionSpaceType, BaseFunctionStorageImp > BaseFunctionSetType;
    typedef DGMapper<IndexSetType, polOrd, DimRange> MapperType;
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

   typedef DiscontinuousGalerkinBaseFunctionFactory<
      typename Traits::FunctionSpaceType, polOrd> FactoryType;

    //! Dimension of the range vector field
    enum { dimRange = Traits::FunctionSpaceType::DimRange };

    //! The polynom order of the base functions
    enum { PolOrd = polOrd };
    
    //! The polynom order of the base functions
    enum { polynomialOrder = polOrd };

  public:
    //- Constructors and destructors
    /** Constructor */
    DiscontinuousGalerkinSpace(GridPartImp& gridPart) :
    DiscontinuousGalerkinSpaceBase <Traits> (gridPart) {}

    template <class EntityType>
    static BaseFunctionSetType* setBaseFuncSetPointer(EntityType& en) 
    {
      typedef typename ToScalarFunctionSpace<
        typename Traits::FunctionSpaceType>::Type ScalarFunctionSpaceType;
      
      DiscontinuousGalerkinBaseFunctionFactory<
        ScalarFunctionSpaceType, polOrd> fac(en.geometry().type());
      return new BaseFunctionSetType(fac);
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

    typedef LegendreDiscontinuousGalerkinSpace<
      FunctionSpaceType, GridPartType, polOrd, BaseFunctionStorageImp> DiscreteFunctionSpaceType;
 
    typedef VecBaseFunctionSet<FunctionSpaceType, BaseFunctionStorageImp > BaseFunctionSetType;
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

  public:
    //- Constructors and destructors
    /** Constructor */
    LegendreDiscontinuousGalerkinSpace(GridPartImp& gridPart) :
      DiscontinuousGalerkinSpaceBase<Traits> (gridPart) {}

    template <class EntityType>
    static BaseFunctionSetType* setBaseFuncSetPointer(EntityType& en) 
    {
      typedef typename ToScalarFunctionSpace<
        typename Traits::FunctionSpaceType>::Type ScalarFunctionSpaceType;
      
      LegendreDGBaseFunctionFactory<
        ScalarFunctionSpaceType, polOrd> fac(en.geometry().type());
      return new BaseFunctionSetType(fac);
    }
  };
  
} // end namespace Dune 

#endif
