#ifndef DUNE_DISCRETEFUNCTIONSPACE_HH
#define DUNE_DISCRETEFUNCTIONSPACE_HH

//- system includes
#include <assert.h>

//- Dune includes 
#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/fem/space/common/geometryconversion.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/common/basefunctioninterface.hh>


//- local includes 
#include "allgeomtypes.hh"
#include "singletonlist.hh"

namespace Dune{

  /** @addtogroup DiscreteFunctionSpace 
      This provides the interfaces for discrete function spaces. 
      Discrete function spaces contain functions
      from a \ref FunctionSpaceInterface "function space"
      but the domain is defined by a grid
      or more precisly by a \ref GridPart "grid part".
 
      \remarks The interface for using a DiscreteFunctionSpace is
      defined by the class DiscreteFunctionSpaceInterface.

      @{
  */

  //! \brief enumerator for identification of spaces 
  enum DFSpaceIdentifier {  
    LagrangeSpace_id , //!< id for Lagrange Space 
    DGSpace_id ,       //!< id for Discontinuous Galerkin Space 
    CombinedSpace_id , //!< id for Combined Space 
    FiniteVolumeSpace_id , //!< id for Finite Volume Space 
    DFAdapter_id    //!< id for DiscreteFunctionSpace Adapter
  };
 
  //**************************************************************************
  //
  //  --DiscreteFunctionSpaceInterface
  //
  /**  
    \brief This is the interface for discrete function spaces. All methods
    declared here have to be implemented by the implementation class.

    The discrete function space always depends on a given grid. 
    For all diffrent element types of the grid the function space provides 
    a set of base functions for the different element types. 
    Because of the knowledge of on the one hand the grid an on the other
    hand the base functions sets, the discrete function space provides the size
    of the function space and a mapping from entity and local dof number
    to global dof number of the level of the entity.
    \note A DiscreteFunctionSpace is defined on a certain grid part.

  */
  template<class FunctionSpaceTraits>
  class DiscreteFunctionSpaceInterface : 
    public FunctionSpaceTraits::FunctionSpaceType  
  {
  public:
    //- Typedefs and enums
    //! type of traits class 
    typedef FunctionSpaceTraits Traits; 
    //! type of DiscretefunctionSapce implementation (Barton-Nackman)
    typedef typename FunctionSpaceTraits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    //! type of \ref FunctionSpaceInterface "function space" (define domain and range types)
    typedef typename FunctionSpaceTraits::FunctionSpaceType FunctionSpaceType;
    //! type of \ref BaseFunctionSetInterface "base function set" of this space 
    typedef typename FunctionSpaceTraits::BaseFunctionSetType BaseFunctionSetType;
    //! type of \ref DofMapperInterface "mapper" of this space 
    typedef typename FunctionSpaceTraits::MapperType MapperType;
    //! type of underlying \ref GridPart "grid part" 
    typedef typename FunctionSpaceTraits::GridPartType GridPartType;
    //! type of underlying dune grid  
    typedef typename GridPartType::GridType GridType;
    //! type of used dune index set 
    typedef typename GridPartType::IndexSetType IndexSetType;
    
    /** \brief iterator type traversing the set of 
        entities defining the discrete  function space 
        (only codim 0 at the moment, to be revised)
    */
    typedef typename GridPartType:: template Codim<0>::IteratorType IteratorType;
    
  public:
    //- Public methods
    //! Constructor 
    DiscreteFunctionSpaceInterface() : FunctionSpaceType() 
    {}

    //! Method provided by implementation

    /** \brief return type identifier of discrete function space 
        \return return type identifier of discrete function space
    */
    DFSpaceIdentifier type () const 
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().type());
      return asImp().type();
    }

    /** \brief get base function set for given entity 
        \param[in] entity Entity for which base function is requested 
        \return BaseFunctionSet for Entity 
    */
    template< class EntityType >
    const BaseFunctionSetType baseFunctionSet ( const EntityType &entity ) const 
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().baseFunctionSet( entity ));
      return asImp().baseFunctionSet( entity );
    }
  
    /** \brief return true if space contains global continuous functions 
       (i.e. for LagrangeSpace \b true is returned, for DiscontinuousGalerkinSpace \b false is returned. 
       \return \b true if space contians global continous functions, \b false> otherwise 
    */
    bool continuous() const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().continuous());
      return asImp().continuous(); 
    }

    /** \brief return reference to grid which belongs to discrete function space 
        \return reference to grid  
    */
    const GridType& grid() const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().grid());
      return asImp().grid(); 
    }

    /** \brief return reference to grid which belongs to discrete function space 
        \return reference to grid  
    */
    GridType& grid() 
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().grid());
      return asImp().grid(); 
    }

    /** \brief Return the corresponding grid part (const version) 
        \return reference to grid part 
    */ 
    const GridPartType& gridPart() const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().gridPart());
      return asImp().gridPart(); 
    }

    /** \brief Return the corresponding grid part (const version) 
        \return reference to grid part 
    */ 
    GridPartType& gridPart() 
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().gridPart());
      return asImp().gridPart(); 
    }

    /** \brief Return reference to the corresponding index set of the space 
        \return reference to index set  
    */ 
    const IndexSetType& indexSet() const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().indexSet());
      return asImp().indexSet(); 
    }

    /** \brief Return number of degrees of freedom for this space 
        \return number of degrees of freedom 
    */
    int size () const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().size());
      return asImp().size(); 
    }

    /** \brief returns index of sequence in grid sequences 
        \return number of current sequence 
    */
    int sequence () const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().sequence());
      return asImp().sequence(); 
    }

    /** \brief get global order of space  
        \return order of space, i.e. polynomial order of base functions 
    */
    int order () const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().order());
      return asImp().order(); 
    } 
  
    /** \brief return the instance of the mapper
        \return refernce to mapper
    */    
    MapperType& mapper () const  
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().mapper());
      return asImp().mapper();
    }

    /** \brief Iterator pointing to first entity associated 
               with this discrete function space 
        \return Iterator pointing to first Entity
    */
    IteratorType begin() const {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().begin());
      return asImp().begin();
    }

    /** \brief Iterator pointing behind last entity associated 
               with this discrete function space 
        \return Iterator pointing behind last Entity
    */
    IteratorType end() const {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().end());
      return asImp().end();
    }

    /** \brief returns true if grid has more than one geometry type (hybrid grid)
        \return \b true  if  grid has more than one geometry type
        (hybrid grid), \b false otherwise 
        \hasdefault
    */
    bool multipleGeometryTypes() const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().multipleGeometryTypes());
      return asImp().multipleGeometryTypes();
    }

    /** \brief returns true if base function sets depend on entity 
        \return \b true if base function set depend on entities, \b false otherwise 
        \hasdefault
    */
    bool multipleBaseFunctionSets() const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().multipleBaseFunctionSets());
      return asImp().multipleBaseFunctionSets(); 
    }
    
    /** \brief For given entity map local dof number to global dof number 
        \param[in] entity   Entity for which mapping is done 
        \param[in] localDof local dof number 
        \return global dof number, i.e. position in dof array 
        \hasdefault
    */    
    template <class EntityType>
    int mapToGlobal ( const EntityType &entity, 
                      const int localDof ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().mapToGlobal(entity, localDof ));
      return asImp().mapToGlobal ( entity , localDof );
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

  //---------------------------------------------------------------------------
  //-
  //-  --DiscreteFunctionSpaceDefault
  //-
  //-
  //---------------------------------------------------------------------------
  /** \brief This is the class with default implementations for discrete function */
  template <class FunctionSpaceTraits>
  class DiscreteFunctionSpaceDefault :
    public DiscreteFunctionSpaceInterface<FunctionSpaceTraits>
  {
  public:
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

    /** \brief @copydoc DiscreteFunctionSpaceInterface::multipleGeometryTypes */
    bool multipleGeometryTypes() const { return multipleGeometryTypes_; }

    /** \brief @copydoc DiscreteFunctionSpaceInterface::multipleBaseFunctionSets 
        \note The default implementation returns false.
    */
    bool multipleBaseFunctionSets() const { return false; }

    /** \brief @copydoc DiscreteFunctionSpaceInterface::mapToGlobal */
    template <class EntityType>
    int mapToGlobal ( const EntityType &entity, 
                      const int localDof ) const
    {
      return this->mapper().mapToGlobal ( entity , localDof );
    }
  protected: 
    //! true if grid has more than one geometry type (hybrid grids)
    const bool multipleGeometryTypes_;
 
  private:
    //! Barton-Nackman trick 
    DiscreteFunctionSpaceType& asImp() 
    { 
      return static_cast<DiscreteFunctionSpaceType&>(*this); 
    }

    //! Barton-Nackman trick 
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
    
    //- type of function space 
    typedef FunctionSpaceImp FunctionSpaceType;
    //- grid part type 
    typedef GridPartImp GridPartType;
    //- grid type 
    typedef typename GridPartType :: GridType GridType;
    //- type of used entity
    typedef typename GridType :: template Codim<0> :: Entity EntityType;
    //- type of iterator 
    typedef typename GridPartType :: template Codim<0> :: IteratorType IteratorType; 
    //- type of IndexSet 
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

    /** \brief @copydoc DiscreteFunctionSpaceInterface::begin */
    IteratorType begin () const { return gridPart_.template begin<0> (); }
    /** \brief @copydoc DiscreteFunctionSpaceInterface::end */
    IteratorType end () const { return gridPart_.template end<0> (); }

    /** \brief @copydoc DiscreteFunctionSpaceInterface::gridPart */
    const GridPartType& gridPart() const { return gridPart_; }
    /** \brief @copydoc DiscreteFunctionSpaceInterface::indexSet */
    const IndexSetType& indexSet() const { return gridPart_.indexSet(); }
    /** \brief @copydoc DiscreteFunctionSpaceInterface::grid */
    const GridType& grid () const { return gridPart_.grid(); }

    /** \brief @copydoc DiscreteFunctionSpaceInterface::continuous */
    bool continuous () const { return true; }

    /** \brief @copydoc DiscreteFunctionSpaceInterface::order */
    int order () const { return polynomialOrder; }

    /** \brief @copydoc DiscreteFunctionSpaceInterface::type */
    DFSpaceIdentifier type () const { return DFAdapter_id; }

  protected:
    //! grid part to select view of grid 
    const GridPartType& gridPart_;
  };
///@}  

  /**\ingroup HelperClasses 
     \brief 
     BaseFunctionSetSingletonFactory provides method createObject and
     deleteObject for the SingletonList  
  */
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
  
} // end namespace Dune 
#endif
