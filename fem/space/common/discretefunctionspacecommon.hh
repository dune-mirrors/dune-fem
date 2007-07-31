#ifndef DUNE_DISCRETEFUNCTIONSPACECOMMON_HH
#define DUNE_DISCRETEFUNCTIONSPACECOMMON_HH

//- system includes
#include <cassert>

//- local includes
#include <dune/fem/space/common/discretefunctionspace.hh>

namespace Dune{

  /** @ingroup DiscreteFunctionSpace 
      @{
  */

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
  
//@}  
} // end namespace Dune 
#endif
