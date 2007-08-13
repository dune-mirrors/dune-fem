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
class LocalFunctionWrapper
: public LocalFunctionDefault
  < typename DiscreteFunctionImp :: DiscreteFunctionSpaceType,
    LocalFunctionWrapper< DiscreteFunctionImp >
  >
{
public:
  //! type of discrete function
  typedef DiscreteFunctionImp DiscreteFunctionType;

  //! type of discrete function space 
  typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;
  
private:
  typedef LocalFunctionWrapper< DiscreteFunctionType > ThisType;
  typedef LocalFunctionDefault< DiscreteFunctionSpaceType, ThisType > BaseType;

public:
  //! type of local function implementation 
  typedef typename DiscreteFunctionImp :: LocalFunctionImp  LocalFunctionImp;
  
  //! Iterator over the space
  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  //! The codim 0 entity
  typedef typename IteratorType::Entity EntityType;

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

protected:
  // type of stack entry 
  typedef typename LFStorage :: ObjectPointerType LocalFunctionPtrType;
 
private:
  // pair storing pointer to local function and poiner to ref-counter 
  LocalFunctionPtrType lfptr_;

  // reference to local function 
  LocalFunctionImp &lf_;

public:
  //! Constructor initializing the underlying local function 
  template< class EntityType > 
  LocalFunctionWrapper( const EntityType &entity,
                        const DiscreteFunctionImp &discreteFunction )
  : lfptr_( discreteFunction.localFunctionStorage().getObject() ),
    lf_( *lfptr_ )
  {
    // init real local function with entity
    localFunction().init( entity );
  }

  //! Constructor creating empty local function 
  LocalFunctionWrapper( const DiscreteFunctionImp &discreteFunction ) 
  : lfptr_( discreteFunction.localFunctionStorage().getObject() ),
    lf_( *lfptr_ )
  {
  }

  //! Copy constructor
  LocalFunctionWrapper( const LocalFunctionWrapper &org )
  : lfptr_( org.lfptr_ ),
    lf_( *lfptr_ )
  {
  }

  //! destructor pushing local funciton back to the stack (just forget the pointer)
  ~LocalFunctionWrapper () 
  { 
  }
  
private:
  // prohibit assignment 
  inline ThisType &operator= ( const ThisType& );

public:
  /** \brief @copydoc LocalFunctionInterface::operator [] */
  RangeFieldType & operator [] (const int num) { return localFunction()[num]; }
  
  /** \brief @copydoc LocalFunctionInterface::operator [] const */
  const RangeFieldType & operator [] (const int num) const { return localFunction()[num]; }

  /** \brief @copydoc LocalFunctionInterface::numDofs */
  int numDofs () const { return localFunction().numDofs(); }
  
  /** \brief @copydoc LocalFunctionInterface::evaluate */
  void evaluate (const DomainType & x, RangeType & ret) const
  {
    localFunction().evaluate( x , ret );
  }
  
  /** \brief @copydoc LocalFunctionInterface::evaluate */
  template <class QuadratureType> 
  void evaluate (const QuadratureType &quad, int quadPoint , RangeType & ret) const
  {
    localFunction().evaluate( quad, quadPoint , ret );
  }
  
  /** \brief @copydoc LocalFunctionInterface::jacobian */
  template <class QuadratureType> 
  void jacobian (const QuadratureType &quad, 
		             const int quadPoint , 
                 JacobianRangeType & grad) const
  {
    localFunction().jacobian(quad, quadPoint, grad ); 
  }
 
  /** \brief @copydoc LocalFunctionInterface::jacobian */
  void jacobian(const DomainType& x, 
            		JacobianRangeType& grad) const
  {
    localFunction().jacobian( x , grad ); 
  }

  /** \brief update local function for given Entity */
  template <class EntityType > 
  void init ( const EntityType &en )
  { 
    localFunction().init(en);
  } 

  /** \brief @copydoc LocalFunctionInterface::axpy */
  template <class QuadratureType> 
  inline void axpy(const QuadratureType &quad, 
                   const int quadPoint , const RangeType & factor)
  {
    localFunction().axpy( quad, quadPoint , factor );
  }
  
  /** \brief @copydoc LocalFunctionInterface::axpy */
  template <class QuadratureType> 
  inline void axpy(const QuadratureType &quad, 
                   const int quadPoint , const JacobianRangeType & factor)
  {
    localFunction().axpy( quad, quadPoint , factor );
  }
  
  /** \brief @copydoc LocalFunctionInterface::axpy */
  template <class QuadratureType> 
  inline void axpy(const QuadratureType &quad, 
                   const int quadPoint ,
                   const RangeType& factor1,
                   const JacobianRangeType & factor2)
  {
    localFunction().axpy( quad, quadPoint, factor1, factor2 );
  }
  
  /** \brief @copydoc LocalFunctionInterface::baseFunctionSet */ 
  const BaseFunctionSetType& baseFunctionSet() const 
  {
    return localFunction().baseFunctionSet();
  }

private:
  //! return reference to local function 
  LocalFunctionImp &localFunction () 
  {
    return lf_;
  } 
  //! return reference to local function 
  const LocalFunctionImp &localFunction () const 
  { 
    return lf_;
  } 
}; // end LocalFunctionWrapper  

} // end namespace Dune
#endif
