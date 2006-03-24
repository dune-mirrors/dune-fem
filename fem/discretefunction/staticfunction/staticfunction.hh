#ifndef DUNE_STATICFUNCTION_HH
#define DUNE_STATICFUNCTION_HH

//- system includes
#include <fstream>
#include <rpc/xdr.h>

//- Dune inlcudes 
#include <dune/common/array.hh>

#include <dune/fem/common/discretefunction.hh>
#include <dune/fem/common/fastbase.hh>
#include <dune/fem/common/localfunction.hh>
#include <dune/fem/common/dofiterator.hh>



namespace Dune{

template <class DiscreteFunctionSpaceType,  class DofStorageImp = Array<typename DiscreteFunctionSpaceType::RangeFieldType> > class StaticDiscreteFunction;
template <class DiscreteFunctionSpaceType,  class DofStorageImp> class StaticDiscreteLocalFunction;
template <class DofStorageImp,class DofImp> class DofIteratorStaticDiscreteFunction;


template <class DiscreteFunctionSpaceImp, class DofStorageImp>
struct StaticDiscreteFunctionTraits {
  typedef DofStorageImp DofStorageType;
  typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;
  typedef StaticDiscreteFunction<DiscreteFunctionSpaceImp,DofStorageType> DiscreteFunctionType;
  typedef StaticDiscreteLocalFunction<DiscreteFunctionSpaceImp,DofStorageType> LocalFunctionImp;
  typedef LocalFunctionWrapper<DiscreteFunctionType> LocalFunctionType;
  typedef DofIteratorStaticDiscreteFunction<DofStorageType,typename DiscreteFunctionSpaceImp::RangeFieldType> DofIteratorType;
  typedef ConstDofIteratorDefault<DofIteratorType> ConstDofIteratorType;
};

namespace DofTypeWrapper { 
  
template <class DofType>
static DofType & convert(FieldVector<DofType,1> & val) { return val[0]; }

template <class DofType>
static const DofType & convert(DofType & val) { return val; }

}

//**********************************************************************
//
//  --StaticDiscreteFunction 
//
//! this is one special implementation of a discrete function using an
//! array for storing the dofs.  
//!
//**********************************************************************
template<class DiscreteFunctionSpaceType , class DofStorageImp> 
class StaticDiscreteFunction 
: public DiscreteFunctionDefault <StaticDiscreteFunctionTraits<DiscreteFunctionSpaceType,DofStorageImp> > 
{
  typedef DiscreteFunctionDefault<StaticDiscreteFunctionTraits <DiscreteFunctionSpaceType,DofStorageImp> >
  DiscreteFunctionDefaultType;

  friend class DiscreteFunctionDefault< StaticDiscreteFunctionTraits <DiscreteFunctionSpaceType,DofStorageImp > > ;

  enum { myId_ = 0};
public:
  //! type of underlying array
  typedef DofStorageImp DofStorageType;

  //! my type 
  typedef StaticDiscreteFunction < DiscreteFunctionSpaceType, DofStorageType > DiscreteFunctionType;

  //! Type of the range field
  typedef typename DiscreteFunctionSpaceType::Traits::RangeFieldType RangeFieldType;
  
  /** \brief For ISTL-compatibility */
  typedef FieldVector<RangeFieldType,1> block_type;

    //! Type of the grid
  typedef typename DiscreteFunctionSpaceType::Traits::GridType GridType;

  //! the local function implementation e 
  typedef StaticDiscreteLocalFunction<DiscreteFunctionSpaceType,DofStorageType> LocalFunctionImp;

  //! LocalFunctionType is the exported lf type 
  typedef LocalFunctionWrapper < DiscreteFunctionType > LocalFunctionType;

  // the storage of the local functions 
  typedef LocalFunctionStorage< DiscreteFunctionType > LocalFunctionStorageType;
  
  //! the dof iterator type of this function
  typedef DofIteratorStaticDiscreteFunction <DofStorageType,typename DiscreteFunctionSpaceType::RangeFieldType> DofIteratorType;
  typedef ConstDofIteratorDefault<DofIteratorType> ConstDofIteratorType;

  //! The associated discrete function space
  typedef DiscreteFunctionSpaceType FunctionSpaceType;

  //! our traits, like DofIterator etc. 
  typedef StaticDiscreteFunctionTraits 
    <DiscreteFunctionSpaceType,DofStorageType > Traits;

  //! the type of the unknowns 
  typedef RangeFieldType DofType;
  
  //! Constructor makes Discrete Function  
  StaticDiscreteFunction ( const DiscreteFunctionSpaceType & f ) ;
  
  //! Constructor makes Discrete Function  
  StaticDiscreteFunction ( const DiscreteFunctionSpaceType & f, const DofStorageType & org ) ;
  
  //! Constructor makes Discrete Function with name 
  StaticDiscreteFunction ( const char * name, const DiscreteFunctionSpaceType & f ) ;
  
  //! Constructor makes Discrete Function from copy 
  StaticDiscreteFunction (const StaticDiscreteFunction <DiscreteFunctionSpaceType,DofStorageType> & df); 

  //! delete stack of free local functions belonging to this discrete
  //! function 
  ~StaticDiscreteFunction ();

  // ***********  Interface  *************************
  //! return object of type LocalFunctionType 
  LocalFunctionType newLocalFunction () DUNE_DEPRECATED;

  //! return local function for given entity
  template <class EntityType>
  LocalFunctionType localFunction(const EntityType& en) const;

  //! update LocalFunction to given Entity en  
  template <class EntityType> 
  void localFunction ( const EntityType &en, LocalFunctionType & lf) DUNE_DEPRECATED; 

  //! return reference to this 
  //! this methods is only to fullfill the interface as parameter classes 
  DiscreteFunctionType & argument    () { return *this; }

  //! return reference to this 
  //! this methods is only to fullfill the interface as parameter classes 
  const DiscreteFunctionType & argument () const { return *this; }

  //! return reference to this 
  //! this methods is only to fullfill the interface as parameter classes 
  DiscreteFunctionType & destination () { return *this; }
 
  //! we use the default implementation 
  DofIteratorType dbegin ();
  
  //! points behind the last dof of type cc
  DofIteratorType dend   ();

  //! the const versions 
  //! we use the default implementation 
  ConstDofIteratorType dbegin () const;
  
  //! points behind the last dof of type cc
  ConstDofIteratorType dend   () const;

  //! Return the name of the discrete function
  const std::string& name() const {return name_;}

  //! return size of this discrete function
  int size() const { return dofVec_.size(); }

  //! set all dofs to zero  
  void clear( );

  //! set all dof to value x 
  void set( RangeFieldType x ); 

  //! add g * scalar to discrete function 
  void addScaled ( const DiscreteFunctionType & g,
      const RangeFieldType &scalar); 
  
  /** \todo Please to me! */
  template <class GridIteratorType>
  void addScaledLocal (GridIteratorType &it, 
      const DiscreteFunctionType & g,
      const RangeFieldType &scalar); 
 
  //! add g to this on local entity
  template <class GridIteratorType>
  void addLocal (GridIteratorType &it, 
      const DiscreteFunctionType & g); 
  
  //! add g to this on local entity
  template <class GridIteratorType>
  void substractLocal (GridIteratorType &it, 
      const DiscreteFunctionType & g); 
  
    /** \todo Please to me! */
  template <class GridIteratorType>
  void setLocal (GridIteratorType &it, const RangeFieldType &scalar);
  
  //! print all dofs 
  void print(std::ostream& s) const;

  //! write data of discrete function to file filename 
  //! with xdr methods 
  bool write_xdr( const char *filename );

  //! write data of discrete function to file filename 
  //! with xdr methods 
  bool read_xdr( const char *filename );

  //! write function data to file filename in ascii Format
  bool write_ascii(const char *filename);

  //! read function data from file filename in ascii Format
  bool read_ascii(const char *filename);

  //! write function data in pgm fromat file
  bool write_pgm(const char *filename) ;

  //! read function data from pgm fromat file
  bool read_pgm(const char *filename); 

  //! return pointer to internal array for use of BLAS routines 
  DofType * leakPointer () { return &(dofVec_[0]);  };
  //! return pointer to internal array for use of BLAS routines 
  const DofType * leakPointer () const { return &(dofVec_[0]); };

private:  
  //! return object pointer of type LocalFunctionImp 
  LocalFunctionImp * newLocalFunctionObject () const;

  // get memory for discrete function
  void getMemory(); 

  //! the name of the function
  std::string name_;

  //! true if memory was allocated
  bool built_;

  //! hold one object for addLocal and setLocal and so on 
  LocalFunctionImp localFunc_;

  //! the dofs stored in an array
  mutable DofStorageType dofVec_;
}; // end class StaticDiscreteFunction 


//**************************************************************************
//
//  --StaticDiscreteLocalFunction
//
//! Implementation of the local functions 
//
//**************************************************************************
template < class DiscreteFunctionSpaceType , class DofStorageImp > 
class StaticDiscreteLocalFunction 
: public LocalFunctionDefault <DiscreteFunctionSpaceType ,
  StaticDiscreteLocalFunction < DiscreteFunctionSpaceType, DofStorageImp >  >
{
  typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;
  typedef StaticDiscreteLocalFunction < DiscreteFunctionSpaceType, DofStorageImp > MyType;
  typedef StaticDiscreteFunction <DiscreteFunctionSpaceType,DofStorageImp> DiscFuncType;

  enum { dimrange = DiscreteFunctionSpaceType::DimRange };
  typedef typename DiscreteFunctionSpaceType::Traits::DomainType DomainType;
  typedef typename DiscreteFunctionSpaceType::Traits::RangeType RangeType;
  typedef typename DiscreteFunctionSpaceType::Traits::RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionSpaceType::Traits::JacobianRangeType JacobianRangeType;
  typedef DofStorageImp DofStorageType;

  friend class StaticDiscreteFunction <DiscreteFunctionSpaceType,DofStorageType>;
  friend class LocalFunctionWrapper < StaticDiscreteFunction <DiscreteFunctionSpaceType,DofStorageType> >;
public:
  //! Constructor 
  StaticDiscreteLocalFunction ( const DiscreteFunctionSpaceType &f , DofStorageType & dofVec );

  //! Destructor 
  ~StaticDiscreteLocalFunction ();

  //! access to dof number num, all dofs of the dof entity
  RangeFieldType & operator [] (int num);
  
  //! access to dof number num, all dofs of the dof entity
  const RangeFieldType & operator [] (int num) const;

  //! return number of degrees of freedom 
  int numberOfDofs () const DUNE_DEPRECATED;

  //! return number of degrees of freedom 
  int numDofs () const;

  //! sum over all local base functions 
  template <class EntityType> 
  void evaluate (EntityType &en, const DomainType & x, RangeType & ret) const ;
  
  template <class EntityType>
  void evaluateLocal(EntityType &en, const DomainType & x, RangeType & ret) const ;
  //! sum over all local base functions evaluated on given quadrature point
  template <class EntityType, class QuadratureType>
  void evaluate (EntityType &en, QuadratureType &quad, int quadPoint , RangeType & ret) const;

  //! sum over all local base functions evaluated on given quadrature point
  template <class EntityType, class QuadratureType>
  void jacobian (EntityType &en, QuadratureType &quad, int quadPoint , JacobianRangeType & ret) const;

  template <class EntityType>
  void jacobianLocal(EntityType& en, const DomainType& x, JacobianRangeType& ret) const ;

  template <class EntityType>
  void jacobian(EntityType& en, const DomainType& x, JacobianRangeType& ret) const;
  
protected:
  //! update local function for given Entity  
  template <class EntityType > 
  void init ( const EntityType &en ) const;

  //! the corresponding function space which provides the base function set
  const DiscreteFunctionSpaceType &fSpace_;
  
  //! dofVec from all levels of the discrete function 
  DofStorageType & dofVec_;

  //! Array holding pointers to the local dofs 
  mutable Array < RangeFieldType * > values_ ;

  //! needed once 
  mutable RangeType tmp_;
  mutable DomainType xtmp_;

  //! needed once 
  mutable JacobianRangeType tmpGrad_;

  //! diffVar for evaluate, is empty 
  const DiffVariable<0>::Type diffVar;

  //! number of all dofs 
  mutable int numOfDof_;

  //! for example number of corners for linear elements 
  mutable int numOfDifferentDofs_;
 
  //! do we have the same base function set for all elements
  bool uniform_;

  //! is it initialised?
  mutable bool init_;
  
}; // end StaticDiscreteLocalFunction 


//***********************************************************************
//
//  --DofIteratorStaticDiscreteFunction
//
//***********************************************************************
  /** \brief Iterator over an array of dofs 
      \todo Please doc me!
  */
template < class DofStorageImp, class DofImp >
class DofIteratorStaticDiscreteFunction : public
DofIteratorDefault < DofImp , DofIteratorStaticDiscreteFunction < DofStorageImp, DofImp > >
{
public:
  typedef DofImp DofType;
  typedef DofStorageImp DofStorageType;
  typedef DofIteratorStaticDiscreteFunction<DofStorageType,DofType> ThisType;
  
  //! Default constructor
  DofIteratorStaticDiscreteFunction() :
    dofArray_ (0) ,
    count_() {}

  //! Constructor (const)
  DofIteratorStaticDiscreteFunction ( const DofStorageType & dofArray , int count )
    :  dofArray_ (const_cast<DofStorageType*>(&dofArray)) ,
       count_ (count) {}
  
  //! Constructor
  DofIteratorStaticDiscreteFunction(DofStorageType& dofArray, int count)
    : dofArray_(&dofArray),
      count_(count) {}

  //! Copy Constructor
  DofIteratorStaticDiscreteFunction (const ThisType& other)
    : dofArray_(other.dofArray_)
    , count_(other.count_) 
  {}

  //! Assignment operator
  ThisType& operator=(const ThisType& other)
  {
    if (&other != this) 
    {
      dofArray_ = other.dofArray_;
      count_ = other.count_;
    }
    return *this;
  }

  //! return dof
  DofType & operator *()
  {
    assert((count_ >=0) && (count_ < dofArray_->size()));
    //return GetValue::get((*dofArray_) [ count_ ]);
    return DofTypeWrapper::template convert<DofType>((*dofArray_)[count_]);
    //return (*dofArray_) [ count_ ];
  }

  //! return dof read only 
  const DofType & operator * () const
  {
    assert((count_ >=0) && (count_ < dofArray_->size()));
    //return GetValue::get((*dofArray_) [ count_ ]);
    return DofTypeWrapper::template convert<DofType>((*dofArray_)[count_]);
    //return (*dofArray_) [ count_ ];
  }

  //! go next dof
  ThisType & operator++ ()
  {
    ++count_;
    return (*this);
  }
  
  //! random access 
  DofType& operator[] (int i)
  {
    assert((i >=0) && (i < dofArray_->size()));
    return (*dofArray_)[i];
  }

  //! random access read only 
  const DofType& operator[] (int i) const
  {
    assert((i >=0) && (i < dofArray_->size()));
    return (*dofArray_)[i];
  }

  //! compare
  bool operator == (const ThisType & I ) const
  {
    return count_ == I.count_;
  }

  //! compare 
  bool operator != (const ThisType & I ) const
  {
    return !((*this) == I);
  }

  //! return actual index 
  int index () const { return count_; }

  //! set dof iterator back to begin , for const and not const Iterators
  void reset () { count_ = 0; }
  
private:
  //! the array holding the dofs 
  DofStorageType * dofArray_;  
  
  //! index 
  mutable int count_;

}; // end DofIteratorStaticDiscreteFunction 


} // end namespace Dune

#include "staticfunction.cc"

#endif


