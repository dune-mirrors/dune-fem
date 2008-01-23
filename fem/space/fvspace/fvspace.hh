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

// use dg data handle 
#include <dune/fem/space/dgspace/dgdatahandle.hh>

// * Note: the dofmanager could be removed from the space altogether now.
// (Maybe this wouldn't be a clever move, though. In my view of a perfect Dune,
// there would be one DofManager per space and the DiscreteFunctions wouldn't
// need to fiddle with the DofMapper anymore...

namespace Dune {

  // Forward declarations
  template <class FunctionSpaceImp, class GridPartImp, int polOrd,
            template<class> class BaseFunctionStorageImp = CachingStorage >
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
    typedef FiniteVolumeMapper< GridPartType, polOrd, DimRange > MapperType;

    // type of mapper for block vector functions 
    typedef FiniteVolumeMapper< GridPartType, polOrd, 1 > BlockMapperType;

    enum { localBlockSize = DimRange };
    
    /** \brief defines type of data handle for communication 
        for this type of space.
    */
    template< class DiscreteFunction,
              class Operation = DFCommunicationOperation :: Copy >
    struct CommDataHandle
    {
      //! type of data handle 
      typedef DGCommunicationHandler< DiscreteFunction, Operation > Type;
      //! type of operation to perform on scatter 
      typedef Operation OperationType;
    };
  };
  //
  //  --FiniteVolumeSpace
  //

  /** @addtogroup FVDFSpace

   Provides access to base function set for different element 
   type in one grid and size of functionspace 
   and map from local to global dof number

   \note This space can only be used with a special set of index sets.
   If you want to use the FiniteVolumeSpace with an index set only
   supportting the index set interface, then use the IndexSetWrapper
   class which will add the needed functionalty.

   \note For adaptive calculations one have to use Index Sets that are
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
    typedef typename GridPartImp::GridType GridType;

    typedef FiniteVolumeSpace< 
          FunctionSpaceImp, GridPartImp, polOrd , BaseFunctionStorageImp
      > FiniteVolumeSpaceType;

    //! type of this pointer 
    typedef FiniteVolumeSpaceType ThisType; 
 
    //! my Traits 
    typedef FiniteVolumeSpaceTraits<
      FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp
      > Traits;

    typedef DiscreteFunctionSpaceDefault<Traits> DefaultType;
  
    /** type of base function set implementation  */
    typedef typename Traits::BaseFunctionSetImp BaseFunctionSetImp;

    typedef typename Traits::BaseFunctionSetType  BaseFunctionSetType;

    typedef typename Traits::IndexSetType IndexSetType;

    typedef typename Traits::GridPartType GridPartType;
    
    typedef typename Traits::IteratorType IteratorType;

    typedef typename Traits::FunctionSpaceType FunctionSpaceType;

    //! type of id 
    typedef int IdentifierType;
    
    //! id is neighbor of the beast
    static const IdentifierType id = 665;

    /** \copydoc FunctionSpace::DimRange */
    enum { DimRange = FunctionSpaceType::DimRange };

    //! size of local blocks, here always 1 
    enum { localBlockSize = Traits :: localBlockSize };
  
    typedef typename Traits :: MapperType MapperType; 

    //! block mapper 
    typedef typename Traits :: BlockMapperType BlockMapperType; 

    //! mapper singleton key  
    typedef MapperSingletonKey< IndexSetType > MapperSingletonKeyType;

    //! mapper factory 
    typedef MapperSingletonFactory< MapperSingletonKeyType ,    
              MapperType > MapperSingletonFactoryType;

    //! mapper singleton list 
    typedef SingletonList< MapperSingletonKeyType , MapperType ,
            MapperSingletonFactoryType > MapperProviderType;

    //! mapper factory 
    typedef MapperSingletonFactory< MapperSingletonKeyType ,    
              BlockMapperType > BlockMapperSingletonFactoryType;

    //! mapper singleton list 
    typedef SingletonList< MapperSingletonKeyType , BlockMapperType ,
            BlockMapperSingletonFactoryType > BlockMapperProviderType;

    /** \copydoc FunctionSpace::DomainType */
    typedef typename FunctionSpaceType::DomainType DomainType;
    /** \copydoc FunctionSpace::RangeType */
    typedef typename FunctionSpaceType::RangeType RangeType;

    typedef typename FunctionSpaceType::RangeFieldType DofType;
    /** \copydoc FunctionSpace::RangeFieldType */
    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
    /** \copydoc FunctionSpace::DomainFieldType */
    typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;

    //! scalar space type 
    typedef typename Traits::ScalarFunctionSpaceType ScalarFunctionSpaceType;

    //! type of base function factory 
    typedef FVBaseFunctionFactory<
      ScalarFunctionSpaceType, polOrd> ScalarFactoryType;

    //! type of singleton factory 
    typedef BaseFunctionSetSingletonFactory<GeometryType,BaseFunctionSetImp,
                ScalarFactoryType> SingletonFactoryType;

    //! type of singleton list  
    typedef SingletonList< GeometryType, BaseFunctionSetImp,
            SingletonFactoryType > SingletonProviderType;

    typedef DofManager<GridType> DofManagerType;
    typedef DofManagerFactory<DofManagerType> DofManagerFactoryType;

  public:
    //! remember polynomial order 
    enum { polynomialOrder =  polOrd };

    //! Constructor generating Finite Volume Space 
    inline explicit FiniteVolumeSpace(GridPartType & g);

    //! Desctructor 
    ~FiniteVolumeSpace (); 

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::contains */
    inline bool contains(const int codim) const
    {   
      return (polynomialOrder == 0) ? (codim == 0) : true;
    }
    
    /** \copydoc Dune::DiscreteFunctionSpaceInterface::continuous */
    bool continuous() const { return (polynomialOrder == 0) ? false : true; }
 
    //! return type of this function space 
    DFSpaceIdentifier type () const;

    //! returns polynomial order
    int order() const { return polynomialOrder; }

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

    //! Return the dof mapper of the space
    MapperType& mapper() const;

    //! Return the dof mapper of the space
    BlockMapperType& blockMapper() const;

    //! \brief return index in grid sequences 
    int sequence () const { return dm_.sequence(); }

  protected:
    //! create functions space
    void makeFunctionSpace (GridPartType& gridPart); 
  
  protected:
    //! type of corresponding map of base function sets
    typedef std::map< const GeometryType, const BaseFunctionSetImp *> BaseFunctionMapType;

    //! the corresponding map of base function sets
    mutable BaseFunctionMapType baseFuncSet_;

  private:
    //! the corresponding FiniteVolumeMapper 
    MapperType* mapper_; 

    //! mapper for block vector functions 
    BlockMapperType& blockMapper_;

    //! reference to dof manager 
    const DofManagerType & dm_;
  }; // end class FiniteVolumeSpace

/** @} **/  
} // end namespace Dune

// contains the implementation of FiniteVolumeSpace
#include "fvspace_inline.hh"
#endif
