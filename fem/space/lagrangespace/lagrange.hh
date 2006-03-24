#ifndef DUNE_LAGRANGEBASE_HH
#define DUNE_LAGRANGEBASE_HH

//- Dune includes 
#include <dune/grid/common/grid.hh>

//- local includes 
#include "../common/discretefunctionspace.hh"
#include "lagrangebasefunctions.hh"
#include "lagrangemapper.hh"

#include "../common/dofmanager.hh"

#include "../../basefunctions/common/basefunctionsets.hh"
#include "../../basefunctions/common/basefunctionfactory.hh"

// * Note: the dofmanager could be removed from the space altogether now.
// (Maybe this wouldn't be a clever move, though. In my view of a perfect Dune,
// there would be one DofManager per space and the DiscreteFunctions wouldn't
// need to fiddle with the DofMapper anymore...

namespace Dune {

  // Forfard declarations
  template <class FunctionSpaceImp, class GridPartImp, int polOrd>
  class LagrangeDiscreteFunctionSpace;

  template <class FunctionSpaceImp,class GridPartImp, int polOrd> 
  struct LagrangeDiscreteFunctionSpaceTraits 
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

    typedef LagrangeDiscreteFunctionSpace<
      FunctionSpaceImp, GridPartImp, polOrd> DiscreteFunctionSpaceType;
    //typedef FastBaseFunctionSet<DiscreteFunctionSpaceType> BaseFunctionSetType;
    typedef VectorialBaseFunctionSet<FunctionSpaceImp> BaseFunctionSetType;
    typedef LagrangeMapper<IndexSetType,polOrd,DimRange> MapperType;
  };

  //****************************************************************
  //
  //  --LagrangeDiscreteFunctionSpace
  //
  //! Provides access to base function set for different element 
  //! type in one grid and size of functionspace 
  //! and map from local to global dof number 
  //! NOTE: This space can only be used with a special set of index sets.
  //! If you want to use the Lagrangespace with an index set only
  //! supportting the index set interface, then use the IndexSetWrapper
  //! class which will add the needed functionalty.
  //!
  //****************************************************************
  template<class FunctionSpaceImp, class GridPartImp, int polOrd>
  class LagrangeDiscreteFunctionSpace : 
    public DiscreteFunctionSpaceDefault
  <
    LagrangeDiscreteFunctionSpaceTraits<FunctionSpaceImp, GridPartImp, 
                                        polOrd> 
  >
  {
 public:
    //! my Grid's type 
    typedef typename GridPartImp::GridType GridType;

    /** \todo Please doc me! */
    typedef LagrangeDiscreteFunctionSpace< 
      FunctionSpaceImp, GridPartImp, polOrd
      > LagrangeDiscreteFunctionSpaceType;
 
    //! my Traits 
    typedef LagrangeDiscreteFunctionSpaceTraits<
      FunctionSpaceImp, GridPartImp, polOrd
      > Traits;

    /** \todo Please doc me! */
    typedef DiscreteFunctionSpaceDefault<Traits> DefaultType;
  
    /** \todo Please doc me! */
    typedef typename Traits::BaseFunctionSetType  BaseFunctionSetType;
    /** \todo Please doc me! */
    typedef typename Traits::IndexSetType IndexSetType;

    /** \todo Please doc me! */
    typedef typename Traits::GridPartType GridPartType;
    
    /** \todo Please doc me! */
    typedef typename Traits::IteratorType IteratorType;

    /** \todo Please doc me! */
    typedef typename Traits::FunctionSpaceType FunctionSpaceType;

    //! id is neighbor of the beast
    static const IdentifierType id = 665;

    // Lagrange 1 , to be revised in this matter 
    /** \todo Please doc me! */
    enum { DimRange = FunctionSpaceType::DimRange };
  
    /** \todo Please doc me! */
    typedef LagrangeMapper<IndexSetType,polOrd,DimRange> MapperType; 

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

    //! The base function factory
    typedef LagrangeBaseFunctionFactory<
      typename Traits::FunctionSpaceType, polOrd> FactoryType;

    //! dimension of value 
    enum { dimVal = 1 };
  
  public:
    //! remember polynomial order 
    enum { polynomialOrder =  polOrd };
  
    //! Constructor generating for each different element type of the grid a 
    //! LagrangeBaseSet with polOrd 
    LagrangeDiscreteFunctionSpace(GridPartType & g);

    //! Desctructor 
    virtual ~LagrangeDiscreteFunctionSpace (); 

    //! continuous
    bool continuous() const { return (polOrd == 0) ? false : true; }
 
    //! return type of this fucntion space 
    DFSpaceIdentifier type () const;

    //! returns polynomial order
    int polynomOrder() const { return polynomialOrder; }

    //! begin iterator
    IteratorType begin() const { return grid_.template begin<0>(); }

    //! end iterator
    IteratorType end() const { return grid_.template end<0>(); }

    //! provide the access to the base function set for a given entity
    template <class EntityType>
    const BaseFunctionSetType& getBaseFunctionSet ( EntityType &en ) const;

    //! default for polOrd 0
    template <class EntityType> 
    bool evaluateLocal (int baseFunc, EntityType &en, const DomainType &local, RangeType & ret) const; 

    //! default for polOrd 0
    template <class EntityType, class QuadratureType> 
    bool evaluateLocal ( int baseFunc, EntityType &en, QuadratureType &quad, 
                         int quadPoint, RangeType & ret) const;


    //! get dimension of value 
    int dimensionOfValue () const;
  
    //! Return grid
    const GridType& grid() const { return grid_.grid(); }

    //! Return index set
    const IndexSetType& indexSet() const { return grid_.indexSet(); }

    //! level
    int level() const { return grid_.level(); }

    //! number of unknows for this function space   
    int size () const;

    //! for given entity map local dof number to global dof number 
    template <class EntityType>
    int mapToGlobal ( EntityType &en, int localNum ) const;

    //! Return the dof mapper of the space
    const MapperType& mapper() const;

  protected:
    // create functions space with basefunction set for given level
    void makeFunctionSpace (GridPartType& gridPart); 
  
    //! get the right BaseFunctionSet for a given Entity 
    BaseFunctionSetType* setBaseFuncSetPointer(GeometryType type);

  protected:
    //! the corresponding vector of base function sets
    //! length is different element types 
    std::vector < BaseFunctionSetType * > baseFuncSet_;

    //! the index set, used by the mapper for mapping between grid and space 
    GridPartType& grid_;
  private:
    //! the corresponding LagrangeMapper 
    MapperType *mapper_; 

  }; // end class LagrangeDiscreteFunctionSpace
  
} // end namespace Dune

// contains the implementation of LagrangeSpace
#include "lagrangespace_imp.cc"
#endif
