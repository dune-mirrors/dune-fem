#ifndef DUNE_LOCALFUNCTIONWRAPPER_HH
#define DUNE_LOCALFUNCTIONWRAPPER_HH

//-s system includes 
#include <cassert>

//- Dune includes 
#include <dune/fem/space/common/objectstack.hh>

namespace Dune{

template < class DFTraits > class DiscreteFunctionDefault;
template < class DiscreteFunctionSpaceType, class LocalFunctionImp > class LocalFunctionDefault;
//**************************************************************************
//
//  --LocalFunctionWrapper 
//
//**************************************************************************
//! Manages the getting and deleting of local function pointers and 
//! acts like a local functions 
template < class DiscreteFunctionImp > 
class LocalFunctionWrapper : 
    public LocalFunctionDefault < 
  typename DiscreteFunctionImp::DiscreteFunctionSpaceType,
  LocalFunctionWrapper < DiscreteFunctionImp > > 
{
public:
  //! type of local function implementation 
  typedef typename DiscreteFunctionImp::LocalFunctionImp  LocalFunctionImp; 
  //! type of discrete functions space 
  typedef typename DiscreteFunctionImp::DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;
  //! Iterator over the space
  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  //! The codim 0 entity
  typedef typename IteratorType::Entity Entity;

  //! type of discrete function 
  typedef DiscreteFunctionImp  DiscreteFunctionType;  
  //! type of base function set 
  typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;

  enum { dimrange = DiscreteFunctionSpaceType::DimRange };
  enum { dimRange = DiscreteFunctionSpaceType::DimRange };
  enum { DimRange = DiscreteFunctionSpaceType::DimRange };
  enum { dimDomain = DiscreteFunctionSpaceType::DimDomain };
  enum { DimDomain = DiscreteFunctionSpaceType::DimDomain };

  // no docu here, then docu is copied from base class 
  typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
  typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
  typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef DiscreteFunctionDefault< typename DiscreteFunctionImp::Traits > 
  DiscreteFunctionDefaultType;
      
  typedef typename DiscreteFunctionDefaultType::LocalFunctionStorageType LFStorage;

private:
  // local function storage stack 
  LFStorage & storage_;

  // type of stack entry 
  typedef typename LFStorage :: StackStorageType StackStorageType;
  
  // pair storing pointer to local function and poiner to ref-counter 
  StackStorageType obj_;

  // reference to local function 
  LocalFunctionImp & lf_;

public:

  //! Constructor initializing the underlying local function 
  template < class EntityType > 
  LocalFunctionWrapper(const EntityType & en, const DiscreteFunctionImp & df) 
    : storage_( df.localFunctionStorage() )
    , obj_ ( storage_.getObject() )
    , lf_ ( *obj_.first )
  {
    // init real local function with entity
    localFunc().init( en );
  }

  //! Constructor creating empty local function 
  LocalFunctionWrapper (const DiscreteFunctionImp & df) 
    : storage_( df.localFunctionStorage() ) 
    , obj_( storage_.getObject() )
    , lf_ ( *obj_.first )
  {}

  //! Copy constructor
  LocalFunctionWrapper(const LocalFunctionWrapper& org) 
    : storage_(org.storage_)
    , obj_( org.obj_ )
    , lf_ ( *obj_.first )
  {
    assert(*obj_.second == 1);
    ++(*(obj_.second));
  }

  //! Destructor , push local function to stack if there are no other 
  //! to it references
  ~LocalFunctionWrapper () 
  { 
    removeObj();
  }

  /** \brief @copydoc LocalFunctionInterface::operator [] */
  RangeFieldType & operator [] (const int num) { return localFunc()[num]; }
  
  /** \brief @copydoc LocalFunctionInterface::operator [] const */
  const RangeFieldType & operator [] (const int num) const { return localFunc()[num]; }

  /** \brief @copydoc LocalFunctionInterface::numDofs */
  int numDofs () const { return localFunc().numDofs(); }
  
  /** \brief @copydoc LocalFunctionInterface::evaluate */
  void evaluate (const DomainType & x, RangeType & ret) const
  {
    localFunc().evaluate( x , ret );
  }
  
  /** \brief @copydoc LocalFunctionInterface::evaluate */
  template <class QuadratureType> 
  void evaluate (const QuadratureType &quad, int quadPoint , RangeType & ret) const
  {
    localFunc().evaluate( quad, quadPoint , ret );
  }
  
  /** \brief @copydoc LocalFunctionInterface::jacobian */
  template <class QuadratureType> 
  void jacobian (const QuadratureType &quad, 
		             const int quadPoint , 
                 JacobianRangeType & grad) const
  {
    localFunc().jacobian(quad, quadPoint, grad ); 
  }
 
  /** \brief @copydoc LocalFunctionInterface::jacobian */
  void jacobian(const DomainType& x, 
            		JacobianRangeType& grad) const
  {
    localFunc().jacobian( x , grad ); 
  }

  /** \brief update local function for given Entity */
  template <class EntityType > 
  void init ( const EntityType &en )
  { 
    localFunc().init(en);
  } 

  /** \brief @copydoc LocalFunctionInterface::axpy */
  template <class QuadratureType> 
  inline void axpy(const QuadratureType &quad, 
                   const int quadPoint , const RangeType & factor)
  {
    localFunc().axpy( quad, quadPoint , factor );
  }
  
  /** \brief @copydoc LocalFunctionInterface::axpy */
  template <class QuadratureType> 
  inline void axpy(const QuadratureType &quad, 
                   const int quadPoint , const JacobianRangeType & factor)
  {
    localFunc().axpy( quad, quadPoint , factor );
  }
  
  /** \brief @copydoc LocalFunctionInterface::axpy */
  template <class QuadratureType> 
  inline void axpy(const QuadratureType &quad, 
                   const int quadPoint ,
                   const RangeType& factor1,
                   const JacobianRangeType & factor2)
  {
    localFunc().axpy( quad, quadPoint, factor1, factor2 );
  }
  
  /** \brief @copydoc LocalFunctionInterface::baseFunctionSet */ 
  const BaseFunctionSetType& baseFunctionSet() const 
  {
    return localFunc().baseFunctionSet();
  }

private:
  // prohibit assignment 
  LocalFunctionWrapper& operator=(const LocalFunctionWrapper);

private:
  //! return reference to local function 
  LocalFunctionImp & localFunc() 
  { 
    assert( obj_.first );
    assert( &lf_ == obj_.first );
    return lf_;
  } 
  //! return reference to local function 
  const LocalFunctionImp & localFunc() const 
  { 
    assert( obj_.first );
    assert( &lf_ == obj_.first );
    return lf_;
  } 
  
  //! method remove the obj by using the storage 
  void removeObj () 
  {
    assert( obj_.first );
    assert( *(obj_.second) > 0);
    // the second counter is left at a value of 1, so we dont have to
    // initialize when getting the object again 
    if( (*(obj_.second)) == 1) {
      storage_.freeObject ( obj_ );
    } 
    else 
      --(*(obj_.second));
  }
}; // end LocalFunctionWrapper  

} // end namespace Dune
#endif
