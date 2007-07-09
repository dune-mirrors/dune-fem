#ifndef DUNE_FVSPACE_HH
#define DUNE_FVSPACE_HH

//- system includes 
#include <map>

//- Dune includes 
#include <dune/grid/common/grid.hh>

//- Dune-Fem includes 
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/common/basefunctionfactory.hh>
#include <dune/fem/space/basefunctions/basefunctionstorage.hh>
#include <dune/fem/space/basefunctions/basefunctionsets.hh>
#include <dune/fem/space/basefunctions/basefunctionproxy.hh>

//- local includes 
#include "fvspacebasefunctions.hh"
#include "fvspacemapper.hh"

// * Note: the dofmanager could be removed from the space altogether now.
// (Maybe this wouldn't be a clever move, though. In my view of a perfect Dune,
// there would be one DofManager per space and the DiscreteFunctions wouldn't
// need to fiddle with the DofMapper anymore...

namespace Dune {

  // Forfard declarations
  template <class FunctionSpaceImp, class GridPartImp, int polOrd,
            template<class> class BaseFunctionStorageImp = SimpleStorage >
  class FiniteVolumeSpace;

  template <class FunctionSpaceImp,class GridPartImp, int polOrd,
            template <class> class BaseFunctionStorageImp > 
  struct FiniteVolumeSpaceTraits 
  {   
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
    typedef typename FunctionSpaceType::ScalarFunctionSpaceType ScalarFunctionSpaceType;

    typedef FiniteVolumeSpace<
      FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp > DiscreteFunctionSpaceType;
    
    typedef VectorialBaseFunctionSet<FunctionSpaceImp, BaseFunctionStorageImp > BaseFunctionSetImp;
    typedef SimpleBaseFunctionProxy<BaseFunctionSetImp> BaseFunctionSetType;
    typedef FiniteVolumeMapper<IndexSetType,polOrd,DimRange> MapperType;

    // always 1 for FVSpace Base , to be revised 
    enum { localBlockSize = 1 };
  };
  //
  //  --FiniteVolumeSpace
  //

  /** @defgroup FVDFSpace Finie Volume Function Space
      @ingroup DiscreteFunctionSpace

   Provides access to base function set for different element 
   type in one grid and size of functionspace 
   and map from local to global dof number 
   NOTE: This space can only be used with a special set of index sets.
   If you want to use the FiniteVolumeSpace with an index set only
   supportting the index set interface, then use the IndexSetWrapper
   class which will add the needed functionalty.
   NOTE: For adaptive calculations one have to use Index Sets that are
   capable for adaptation, i.e. the method adaptive returns true, see 
   AdaptiveLeafIndexSet. 
   @{
  **/

  /** @brief 
      Finite Volume Function Space 
      **/
  template<class FunctionSpaceImp, class GridPartImp, int polOrd, 
           template <class> class BaseFunctionStorageImp >
  class FiniteVolumeSpace : 
    public DiscreteFunctionSpaceDefault
  <
    FiniteVolumeSpaceTraits<FunctionSpaceImp, GridPartImp, 
                            polOrd, BaseFunctionStorageImp > 
  >
  {
 public:
    //! my Grid's type 
    typedef typename GridPartImp::GridType GridType;

    /** \todo Please doc me! */
    typedef FiniteVolumeSpace< 
          FunctionSpaceImp, GridPartImp, polOrd , BaseFunctionStorageImp
      > FiniteVolumeSpaceType;

    //! type of this pointer 
    typedef FiniteVolumeSpaceType ThisType; 
 
    //! my Traits 
    typedef FiniteVolumeSpaceTraits<
      FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp
      > Traits;

    /** \todo Please doc me! */
    typedef DiscreteFunctionSpaceDefault<Traits> DefaultType;
  
    /** type of base function set implementation  */
    typedef typename Traits::BaseFunctionSetImp BaseFunctionSetImp;
    /** exported type of base function set */
    typedef typename Traits::BaseFunctionSetType  BaseFunctionSetType;
    /** \todo Please doc me! */
    typedef typename Traits::IndexSetType IndexSetType;

    /** \todo Please doc me! */
    typedef typename Traits::GridPartType GridPartType;
    
    /** \todo Please doc me! */
    typedef typename Traits::IteratorType IteratorType;

    /** \todo Please doc me! */
    typedef typename Traits::FunctionSpaceType FunctionSpaceType;

    //! type of id 
    typedef int IdentifierType;
    
    //! id is neighbor of the beast
    static const IdentifierType id = 665;

    /** \copydoc FunctionSpace::DimRange */
    enum { DimRange = FunctionSpaceType::DimRange };

    //! size of local blocks, here always 1 
    enum { localBlockSize = Traits :: localBlockSize };
  
    //! mapper used to implement mapToGlobal */
    typedef typename Traits :: MapperType MapperType; 

    //! mapper singleton key  
    typedef MapperSingletonKey< IndexSetType > MapperSingletonKeyType;
    //! mapper factory 
    typedef MapperSingletonFactory< MapperSingletonKeyType ,    
              MapperType > MapperSingletonFactoryType;

    //! mapper singleton list 
    typedef SingletonList< MapperSingletonKeyType , MapperType ,
            MapperSingletonFactoryType > MapperProviderType;

    /** \todo Please doc me! */
    typedef typename FunctionSpaceType::DomainType DomainType;
    /** \todo Please doc me! */
    typedef typename FunctionSpaceType::RangeType RangeType;
    /** \todo Please doc me! */
    typedef typename FunctionSpaceType::RangeFieldType DofType;
    /** \todo Please doc me! */
    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
    /** \todo Please doc me! */
    typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;

    //! scalar space type 
    typedef typename Traits::ScalarFunctionSpaceType ScalarFunctionSpaceType;

    // type of base function factory 
    typedef FVBaseFunctionFactory<
      ScalarFunctionSpaceType, polOrd> ScalarFactoryType;

    // type of singleton factory 
    typedef BaseFunctionSetSingletonFactory<GeometryType,BaseFunctionSetImp,
                ScalarFactoryType> SingletonFactoryType;

    // type of singleton list  
    typedef SingletonList< GeometryType, BaseFunctionSetImp,
            SingletonFactoryType > SingletonProviderType;

    typedef DofManager<GridType> DofManagerType;
    typedef DofManagerFactory<DofManagerType> DofManagerFactoryType;

  public:
    //! remember polynomial order 
    enum { polynomialOrder =  polOrd };

    //! Constructor generating Finite Volume Space 
    FiniteVolumeSpace(GridPartType & g);

    //! Desctructor 
    virtual ~FiniteVolumeSpace (); 

    //! continuous
    bool continuous() const { return (polOrd == 0) ? false : true; }
 
    //! return type of this fucntion space 
    DFSpaceIdentifier type () const;

    //! returns polynomial order
    int order() const { return polynomialOrder; }

    int polynomOrder () const DUNE_DEPRECATED { return order(); } 

    //! begin iterator
    IteratorType begin() const { return gridPart_.template begin<0>(); }

    //! end iterator
    IteratorType end() const { return gridPart_.template end<0>(); }

    //! provide the access to the base function set for a given entity
    template <class EntityType>
    const BaseFunctionSetType 
    baseFunctionSet ( const EntityType &en ) const;

    //! Get base function set for a given id of geom type (mainly used by
    //! CombinedSpace) 
    const BaseFunctionSetType
    baseFunctionSet (const GeometryType& geoType) const;

    //! get dimension of value 
    int dimensionOfValue () const;
  
    //! Return grid
    const GridType& grid() const { return gridPart_.grid(); }

    //! Return index set
    const IndexSetType& indexSet() const { return gridPart_.indexSet(); }

    //! return reference to the spaces grid part
    GridPartType & gridPart () { return gridPart_; }
    
    //! return reference to the spaces grid part
    const GridPartType & gridPart () const { return gridPart_; }

    //! number of unknows for this function space   
    int size () const;

    //! for given entity map local dof number to global dof number 
    template <class EntityType>
    int mapToGlobal ( EntityType &en, int localNum ) const;

    //! Return the dof mapper of the space
    const MapperType& mapper() const;

    //! \brief return index in grid sequences 
    int sequence () const { return dm_.sequence(); }

  protected:
    // create functions space with basefunction set for given level
    void makeFunctionSpace (GridPartType& gridPart); 
  
  protected:
    //! the corresponding map of base function sets
    typedef std::map< const GeometryType, const BaseFunctionSetImp *> BaseFunctionMapType;
    mutable BaseFunctionMapType baseFuncSet_;

    //! the index set, used by the mapper for mapping between grid and space 
    GridPartType& gridPart_;

  private:
    //! the corresponding FiniteVolumeMapper 
    MapperType* mapper_; 

    // reference to dof manager 
    const DofManagerType & dm_;
  }; // end class FiniteVolumeSpace

/** @} **/  
} // end namespace Dune

// contains the implementation of FiniteVolumeSpace
#include "fvspace_inline.hh"
#endif
