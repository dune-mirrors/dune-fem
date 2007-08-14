#ifndef DUNE_DISCRETEFUNCTIONSPACE_HH
#define DUNE_DISCRETEFUNCTIONSPACE_HH

//- system includes
#include <assert.h>

//- Dune includes 
#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/fem/space/common/geometryconversion.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/common/basefunctioninterface.hh>
#include <dune/fem/function/common/temporarylocalfunction.hh>
#include <dune/fem/function/common/localfunctionwrapper.hh>


//- local includes 
#include "allgeomtypes.hh"
#include "singletonlist.hh"

namespace Dune
{

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
   
    \interfaceclass 
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
    
    /** \brief return reference to grid which belongs to discrete function space 
        \hasdefault
        \return reference to grid  
    */
    const GridType& grid() const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().grid());
      return asImp().grid(); 
    }

    /** \brief return reference to grid which belongs to discrete function space 
        \hasdefault
        \return reference to grid  
    */
    GridType& grid() 
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().grid());
      return asImp().grid(); 
    }

    /** \brief Return the corresponding grid part (const version) 
        \hasdefault
        \return reference to grid part 
    */ 
    const GridPartType& gridPart() const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().gridPart());
      return asImp().gridPart(); 
    }

    /** \brief Return the corresponding grid part (const version) 
        \hasdefault
        \return reference to grid part 
    */ 
    GridPartType& gridPart() 
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().gridPart());
      return asImp().gridPart(); 
    }

    /** \brief Return reference to the corresponding index set of the space 
        \hasdefault
        \return reference to index set  
    */ 
    const IndexSetType& indexSet() const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().indexSet());
      return asImp().indexSet(); 
    }

    /** \brief Return number of degrees of freedom for this space 
        \hasdefault
        \return number of degrees of freedom 
    */
    int size () const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().size());
      return asImp().size(); 
    }

    /** \brief Iterator pointing to first entity associated 
               with this discrete function space 
        \hasdefault
        \return Iterator pointing to first Entity
    */
    IteratorType begin() const {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().begin());
      return asImp().begin();
    }

    /** \brief Iterator pointing behind last entity associated 
               with this discrete function space 
        \hasdefault
        \return Iterator pointing behind last Entity
    */
    IteratorType end() const {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().end());
      return asImp().end();
    }

    /** \brief returns true if grid has more than one geometry type (hybrid grid)
        \hasdefault
        \return \b true  if  grid has more than one geometry type
        (hybrid grid), \b false otherwise 
        
    */
    bool multipleGeometryTypes() const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().multipleGeometryTypes());
      return asImp().multipleGeometryTypes();
    }

    /** \brief returns true if base function sets depend on entity 
        \hasdefault
        \return \b true if base function set depend on entities, \b false
        otherwise
    */
    bool multipleBaseFunctionSets() const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().multipleBaseFunctionSets());
      return asImp().multipleBaseFunctionSets(); 
    }
    
    /** \brief For given entity map local dof number to global dof number 
        \hasdefault
        \param[in] entity   Entity for which mapping is done 
        \param[in] localDof local dof number 
        \return global dof number, i.e. position in dof array 
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
  /** \brief This is the class with default implementations for discrete
     function.
     The  methods not marked with having a default
     in the interface class must be provided by
     the implementation; all other methods
     have a default implementation here.
     
     \remark An reference to the GridPart is
             stored in the default implementation. */
  template< class FunctionSpaceTraits >
  class DiscreteFunctionSpaceDefault
  : public DiscreteFunctionSpaceInterface< FunctionSpaceTraits >
  {
  public:
    typedef FunctionSpaceTraits Traits;
    
  private:
    typedef DiscreteFunctionSpaceDefault< Traits > ThisType;
    typedef DiscreteFunctionSpaceInterface< Traits > BaseType;
    
  public:
    typedef typename FunctionSpaceTraits :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;
    typedef typename FunctionSpaceTraits::GridPartType  GridPartType;
    typedef typename GridPartType:: template Codim<0>::IteratorType IteratorType;
    typedef typename GridPartType :: GridType GridType;
    typedef typename GridPartType :: IndexSetType IndexSetType;

    typedef TemporaryLocalFunctionFactory< DiscreteFunctionSpaceType >
      LocalFunctionFactoryType;
    typedef LocalFunctionStack< LocalFunctionFactoryType > LocalFunctionStorageType;

    typedef typename LocalFunctionStorageType :: LocalFunctionType
      LocalFunctionType;

  protected:
    GridPartType &gridPart_;

    const LocalFunctionFactoryType lfFactory_;
    mutable LocalFunctionStorageType lfStorage_;
 
    // true if grid has more than one geometry type (hybrid grids)
    const bool multipleGeometryTypes_;
 
  public:
    //! constructor
    inline explicit DiscreteFunctionSpaceDefault( GridPartType &gridPart )
    : BaseType(),
      gridPart_( gridPart ),
      lfFactory_( asImp() ),
      lfStorage_( lfFactory_ ),
      multipleGeometryTypes_( AllGeomTypes< IndexSetType, GridType > 
                                :: multipleGeomTypes() )
    {
    }

    /** obtain a local function for an entity (to store intermediate values)
     *  
     *  \param[in]  entity  entity for which a local function is desired
     *
     *  \returns a local function backed by a small, fast array
     */
    template< class EntityType >
    LocalFunctionType localFunction ( const EntityType &entity ) const
    {
      return lfStorage_.localFunction( entity );
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

    /** \brief @copydoc DiscreteFunctionSpaceInterface::size */
    int size () const 
    {
      return this->mapper().size();
    }
   
    /** \brief @copydoc DiscreteFunctionSpaceInterface::begin const */
    IteratorType begin() const { return gridPart_.template begin<0>(); }

    /** \brief @copydoc DiscreteFunctionSpaceInterface::end const */
    IteratorType end() const { return gridPart_.template end<0>(); }

    /** \brief @copydoc DiscreteFunctionSpaceInterface::grid const */
    const GridType& grid() const { return gridPart_.grid(); }

    /** \brief @copydoc DiscreteFunctionSpaceInterface::indexSet const */ 
    const IndexSetType& indexSet() const { return gridPart_.indexSet(); }

    /** \brief @copydoc DiscreteFunctionSpaceInterface::gridPart */
    GridPartType & gridPart () { return gridPart_; }
    /** \brief @copydoc DiscreteFunctionSpaceInterface::gridPart const */
    const GridPartType & gridPart () const { return gridPart_; }

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

      \note DiscreteFunctionSpaceAdapter is itself derived from the template
            argument FunctionSpaceImp. Hence, the constructor will call the
            default constructor of FunctionSpaceImp when this class is
            instanciated. So do not use discrete function spaces for the
            first template argument.
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
    explicit DiscreteFunctionSpaceAdapter ( const GridPartType &gridPart )
    : gridPart_( gridPart ) 
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
