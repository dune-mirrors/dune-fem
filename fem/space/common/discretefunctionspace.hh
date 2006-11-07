#ifndef DUNE_DISCRETEFUNCTIONSPACE_HH
#define DUNE_DISCRETEFUNCTIONSPACE_HH

//- system includes
#include <assert.h>

//- Dune includes 
#include <dune/fem/space/common/geometryconversion.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/basefunctions/common/basefunctions.hh>

//- local includes 

namespace Dune{

  /** @defgroup DiscreteFunctionSpace DiscreteFunctionSpace
      @ingroup DiscreteFunction
      This provides the interfaces for discrete function spaces. 
  
      @{
  */

  enum DFSpaceIdentifier {  LagrangeSpace_id , DGSpace_id , RaviartThomasSpace_id, PerLagrangeSpace_id };
    
 
  //**************************************************************************
  //
  //  --DiscreteFunctionSpaceInterface
  //
  /*! This is the interface for discrete function spaces. All methods
    declared here have to be implemented by the implementation class.
    The discrete function space always depends on a given grid. 
    For all diffrent element types of the grid the function space provides 
    a set of base functions for the different element types. 
    Because of the knowledge of on the one hand the grid an on the other
    hand the base functions sets, the discrete function space provides the size
    of the function space and a mapping from entity and local dof number
    to global dof number of the level of the entity.
    NOTE: A FunctionSpace is defined on a certain grid part.
  */
  template<class FunctionSpaceTraits>
  class DiscreteFunctionSpaceInterface : 
    public FunctionSpaceTraits::FunctionSpaceType  
  {
  public:
    //- Typedefs and enums
    typedef typename FunctionSpaceTraits::FunctionSpaceType FunctionSpaceType;
    //! type of DiscretefunctionSapce implementation (Barton-Nackman)
    typedef typename FunctionSpaceTraits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    //! type of BaseFunctionSet implementation (Barton-Nackman)
    typedef typename FunctionSpaceTraits::BaseFunctionSetType BaseFunctionSetType;
    //! type of underlying grid part 
    typedef typename FunctionSpaceTraits::GridPartType GridPartType;
    //! type of underlying grid  
    typedef typename GridPartType::GridType GridType;
    //! type of used index set 
    typedef typename GridPartType::IndexSetType IndexSetType;
    //! iterator type traversing the set of entities defining the discrete 
    //! function space (only codim 0 at the moment, to be revised)
    typedef typename GridPartType:: template Codim<0>::IteratorType IteratorType;
    
  public:
    //- Public methods
    //! Constructor 
    DiscreteFunctionSpaceInterface(int ident) DUNE_DEPRECATED :
      FunctionSpaceType(ident) {};

    DiscreteFunctionSpaceInterface() :
      FunctionSpaceType() {};

    //! Get base function set for given entity. 
    //! For a type of element the base function set is unique.
    template <class EntityType>
    const BaseFunctionSetType&  baseFunctionSet ( EntityType &en ) const 
    {
      return asImp().baseFunctionSet( en );
    }
  
    //! deprecated 
    template <class EntityType>
    const BaseFunctionSetType&  getBaseFunctionSet ( EntityType &en ) const {
      return asImp().getBaseFunctionSet( en );
    }
  
    //! return true if space is continuous 
    bool continuous() const { return asImp().continuous(); }

    //! Return the corresponding Grid (const version) 
    const GridType& grid() const { return asImp().grid(); }

    //! Return the corresponding Grid 
    GridType& grid() { return asImp().grid(); }

    //! Return the corresponding grid part (const version) 
    const GridPartType& gridPart() const { return asImp().gridPart(); }

    //! Return the corresponding grid part  
    GridPartType& gridPart() { return asImp().gridPart(); }

    //! Return the index set corresponding to the iterator
    const IndexSetType& indexSet() const { return asImp().indexSet(); }

    //! Return number of degrees of freedom for specified grid part 
    int size () const { return asImp().size(); }

    //! For given entity map local dof number to global dof number 
    //! at the level of the given entity. 
    template <class EntityType>
    int mapToGlobal ( EntityType &en, int localNum ) const
    {
      return asImp().mapToGlobal ( en , localNum );
    }
  
    //! Iterator over the entities of a space 
    //! The index set specifies the subset of grid entities of all available 
    //! codimensions; usually, the elements of the space normally belong only
    //! to one codimension, which is selected by the space
    IteratorType begin() const {
      return asImp().begin();
    }

    //! End iterator
    IteratorType end() const {
      return asImp().end();
    }

  protected:
    //! Barton-Nackman trick 
    DiscreteFunctionSpaceType& asImp() 
    { 
      return static_cast<DiscreteFunctionSpaceType&>(*this); 
    }

    //! Barton-Nackman trick 
    const DiscreteFunctionSpaceType& asImp() const  
    { 
      return static_cast<const DiscreteFunctionSpaceType&>(*this); 
    }
  
  }; // end class DiscreteFunctionSpaceInterface

  //**************************************************************************
  //
  // --DiscreteFunctionSpaceDefault
  //
  //! This is the class with default implementations for discrete function
  //! space. 
  //!
  //**************************************************************************
  template <class FunctionSpaceTraits>
  class DiscreteFunctionSpaceDefault :
    public DiscreteFunctionSpaceInterface<FunctionSpaceTraits>
  {
  public:
    //! The implementation type
    typedef typename FunctionSpaceTraits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename FunctionSpaceTraits::GridPartType  GridPartType;

  public:
    //! Constructor
    DiscreteFunctionSpaceDefault(GridPartType & gridPart, int id) 
      : DiscreteFunctionSpaceInterface<FunctionSpaceTraits>(id) 
      , multipleGeometryTypes_( gridPart.indexSet().geomTypes(0).size() > 1 ) 
    {
    }

    //! returns true if grid has more than one geometry type (hybrid grid)
    bool multipleGeometryTypes() const { return multipleGeometryTypes_; }

  protected: 
    //! true if grid has more than one geometry type (hybrid grids)
    const bool multipleGeometryTypes_;
 
  private:
    //! Barton-Nackman trick 
    DiscreteFunctionSpaceType& asImp() 
    { 
      return static_cast<DiscreteFunctionSpaceType&>(*this); 
    }

    const DiscreteFunctionSpaceType &asImp() const  
    { 
      return static_cast<const DiscreteFunctionSpaceType&>(*this); 
    }
  };
  
  /** @} end documentation group */

} // end namespace Dune 

#endif
