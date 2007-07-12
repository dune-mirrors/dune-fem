#ifndef DUNE_DISCRETEFUNCTIONSPACE_HH
#define DUNE_DISCRETEFUNCTIONSPACE_HH

//- system includes
#include <assert.h>

//- Dune includes 
#include <dune/fem/space/common/geometryconversion.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/common/basefunctioninterface.hh>


//- local includes 
#include "allgeomtypes.hh"
#include "singletonlist.hh"

namespace Dune{

  /** @defgroup DiscreteFunctionSpace DiscreteFunctionSpace
      @ingroup FunctionCommon
      This provides the interfaces for discrete function spaces. 
  
      @{
  */

  enum DFSpaceIdentifier {  LagrangeSpace_id , DGSpace_id , 
    CombinedSpace_id , FiniteVolumeSpace_id , DFAdapter_id };
 
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
    DiscreteFunctionSpaceInterface() :
      FunctionSpaceType() {};

    //! Get base function set for given entity. 
    //! For a type of element the base function set is unique.
    template< class EntityType >
    const BaseFunctionSetType baseFunctionSet ( const EntityType &entity ) const 
    {
      return asImp().baseFunctionSet( entity );
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

    //! returns index of sequence in grid sequences 
    int sequence () const { return asImp().sequence(); }

    //! get global order of space  
    int order () const { return asImp().order(); } 
  
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
    DiscreteFunctionSpaceDefault(const GridPartType & gridPart) 
      : DiscreteFunctionSpaceInterface<FunctionSpaceTraits>() 
      , multipleGeometryTypes_( 
          AllGeomTypes< typename GridPartType::IndexSetType ,
                        typename GridPartType::GridType > :: multipleGeomTypes() )
    {
    }

    //! returns true if grid has more than one geometry type (hybrid grid)
    bool multipleGeometryTypes() const { return multipleGeometryTypes_; }

    //! default returns false  
    bool multipleBaseFunctionSets() const { return false; }

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

  ////////////////////////////////////////////////////////////
  //
  //  DiscreteFunctionSpaceAdapter 
  //
  ////////////////////////////////////////////////////////////
  /** \brief Create Obejct that behaves like a discrete function space 
      without to provide functions with the iterator facilities. 
  */
  template <class FunctionSpaceImp, class GridPartImp>
  class DiscreteFunctionSpaceAdapter : public FunctionSpaceImp
  {
  public:  
    enum { polynomialOrder = 111 };
    
    //! type of function space 
    typedef FunctionSpaceImp FunctionSpaceType;
    //! grid part type 
    typedef GridPartImp GridPartType;
    //! grid type 
    typedef typename GridPartType :: GridType GridType;
    //! type of used entity
    typedef typename GridType :: template Codim<0> :: Entity EntityType;
    //! type of iterator 
    typedef typename GridPartType :: template Codim<0> :: IteratorType IteratorType; 
    //! type of IndexSet 
    typedef typename GridPartType :: IndexSetType IndexSetType; 
    
    //! constructor taking grid Part 
    DiscreteFunctionSpaceAdapter(const GridPartType& gridPart) 
      : gridPart_(gridPart) 
    {
    }

    //! copy constructor
    DiscreteFunctionSpaceAdapter(const DiscreteFunctionSpaceAdapter& org) 
      : gridPart_(org.gridPart_) 
    {
    }

    //! return begin iterator 
    IteratorType begin () const { return gridPart_.template begin<0> (); }
    //! return end iterator 
    IteratorType end () const { return gridPart_.template end<0> (); }

    //! return reference to grid part 
    const GridPartType& gridPart() const { return gridPart_; }
    //! return reference to index set 
    const IndexSetType& indexSet() const { return gridPart_.indexSet(); }
    //! return reference to grid 
    const GridType& grid () const { return gridPart_.grid(); }

    //! space for continuous functions 
    bool continuous () const { return true; }

    //! return order which is infinity 
    int order () const { return polynomialOrder; }

    //! return type of space 
    DFSpaceIdentifier type () const { return DFAdapter_id; }

  protected:
    //! grid part to select view of grid 
    const GridPartType& gridPart_;
  };


  //! BaseFunctionSetSingletonFactory provides method createObject and
  //! deleteObject for the SingletonList  
  template <class KeyImp, class ObjectImp, class ObjectFactoryImp>
  class BaseFunctionSetSingletonFactory
  { 
  public:
    //! create new BaseFunctionSet 
    static ObjectImp * createObject( const KeyImp & key )
    {
      ObjectFactoryImp fac(key); 
      return new ObjectImp(fac); 
    }
    
    //! delete BaseFunctionSet 
    static void deleteObject( ObjectImp * obj ) 
    {
      delete obj;
    }
  };

  /** @} end documentation group */

} // end namespace Dune 
#endif
