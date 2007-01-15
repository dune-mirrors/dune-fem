#ifndef DUNE_STATICFUNCTION_HH
#define DUNE_STATICFUNCTION_HH

//- system includes
#include <fstream>
#include <rpc/xdr.h>

//- Dune inlcudes 
#include <dune/common/array.hh>

//- local includes 
#include "../common/discretefunction.hh"
#include "../common/localfunction.hh"
#include "../common/dofiterator.hh"

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

template <class DofType>
struct DofTypeWrapper
{
  enum { blockSize = DofType :: dimension };
  static typename DofType::field_type & convert(DofType & val, int idx) 
  { 
    return val[idx]; 
  }
};

template <>
struct DofTypeWrapper<double>
{
  enum { blockSize = 1 };
  static double & convert(double  & val, int idx) { return val; }
};

template < class DofStorageImp, class DofImp >
class DoubleArrayWrapper
{
public:
  typedef DofImp DofType;
  typedef DofStorageImp DofStorageType;
  typedef typename DofStorageType :: block_type DofBlockType;
  typedef DoubleArrayWrapper<DofStorageType,DofType> ThisType;

  typedef DofTypeWrapper<DofBlockType> DofWrapperType;
  enum { blockSize = DofWrapperType :: blockSize };
  
  //! Constructor (const)
  DoubleArrayWrapper(DofStorageType & dofArray)
    :  dofArray_ (dofArray) {}
  
  DofType& operator [] (int idx) 
  {
    int newIdx = idx%blockSize; 
    int count  = idx/blockSize; 
    return dofArray_[count][newIdx];
  }

  const DofType& operator [] (int idx) const
  {
    int newIdx = idx%blockSize; 
    int count  = idx/blockSize; 
    return dofArray_[count][newIdx];
  }

private:
  //! the array holding the dofs 
  DofStorageType & dofArray_;  
}; // end DofIteratorStaticDiscreteFunction 


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
  //typedef FieldVector<RangeFieldType,1> block_type;
  typedef typename DofStorageImp :: block_type block_type; 

  //! Type of the grid
  typedef typename DiscreteFunctionSpaceType::Traits::GridType GridType;

  //! dof manager 
  typedef DofManager<GridType> DofManagerType;
  typedef DofManagerFactory<DofManagerType> DofManagerFactoryType;

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
  typedef block_type DofBlockType;
  typedef DoubleArrayWrapper<DofStorageType,DofType> LeakPointerType;
  
  //! Constructor makes Discrete Function  
  StaticDiscreteFunction ( const DiscreteFunctionSpaceType & f ) ;
  
  //! Constructor makes Discrete Function  
  StaticDiscreteFunction ( const DiscreteFunctionSpaceType & f, const DofStorageType & org ) ;
  
  //! Constructor makes Discrete Function with name 
  StaticDiscreteFunction ( const std::string name, const DiscreteFunctionSpaceType & f ) ;
  
  //! Constructor makes Discrete Function from copy 
  StaticDiscreteFunction (const StaticDiscreteFunction <DiscreteFunctionSpaceType,DofStorageType> & df); 

  //! delete stack of free local functions belonging to this discrete
  //! function 
  ~StaticDiscreteFunction ();

  //! return local function for given entity
  template <class EntityType>
  LocalFunctionType localFunction(const EntityType& en) const;

  //! update LocalFunction to given Entity en  
  template <class EntityType> 
  void localFunction ( const EntityType &en, LocalFunctionType & lf) DUNE_DEPRECATED; 

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

  //! add g * scalar to discrete function 
  void addScaled ( const DiscreteFunctionType & g,
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
  LeakPointerType& leakPointer () { return leakPointer_;  }
  //! return pointer to internal array for use of BLAS routines 
  const LeakPointerType& leakPointer () const { return leakPointer_; }

  //! return reference to internal block vector 
  DofStorageType& blockVector () const { return dofVec_; }

private:  
  //! return object pointer of type LocalFunctionImp 
  LocalFunctionImp * newLocalFunctionObject () const;

  //! the name of the function
  std::string name_;

  typedef DGMapper<typename DiscreteFunctionSpaceType :: IndexSetType ,0,1> MapperType;
  MapperType mapper_;

  DofManagerType & dm_;

  // MemObject that manages the memory for the dofs of this function
  std::pair<MemObjectInterface*, DofStorageType*> memPair_;

  //! true if memory was allocated
  bool built_;

  //! the dofs stored in an array
  mutable DofStorageType& dofVec_;

  //! hold one object for addLocal and setLocal and so on 
  LocalFunctionImp localFunc_;

  //! pretend to be double array
  mutable LeakPointerType leakPointer_;
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
  //CompileTimeChecker<dimrange == 1> check; 
  typedef typename DiscreteFunctionSpaceType::Traits::DomainType DomainType;
  typedef typename DiscreteFunctionSpaceType::Traits::RangeType RangeType;
  typedef typename DiscreteFunctionSpaceType::Traits::RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionSpaceType::Traits::JacobianRangeType JacobianRangeType;

  typedef typename DiscreteFunctionSpaceType::Traits::GridType :: template
    Codim<0> :: Entity EntityType;
  typedef DofStorageImp DofStorageType;

  typedef typename DofStorageType :: block_type DofBlockType;

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
  int numDofs () const;

  //! sum over all local base functions 
  void evaluate (EntityType &en, const DomainType & x, RangeType & ret) const ;
  
  //! sum over all local base functions 
  void evaluate (const DomainType & x, RangeType & ret) const ;
  
  void evaluateLocal(const DomainType & x, RangeType & ret) const ;
  void evaluateLocal(EntityType &en, const DomainType & x, RangeType & ret) const ;
  //! sum over all local base functions evaluated on given quadrature point
  template <class QuadratureType>
  void evaluate (EntityType &en, QuadratureType &quad, int quadPoint , RangeType & ret) const;

  template <class QuadratureType>
  void evaluate (QuadratureType &quad, int quadPoint , RangeType & ret) const;

  //! sum over all local base functions evaluated on given quadrature point
  template < class QuadratureType>
  void jacobian (EntityType &en, QuadratureType &quad, int quadPoint , JacobianRangeType & ret) const;

  //! evaluate jacobian of local function 
  void jacobianLocal(EntityType& en, const DomainType& x, JacobianRangeType& ret) const ;

  //! evaluate jacobian of local function 
  void jacobian(EntityType& en, const DomainType& x, JacobianRangeType& ret) const;
  
  //! return reference to base function set
  const BaseFunctionSetType& baseFunctionSet() const;
protected:
  //! update local function for given Entity  
  template <class EntityImp> 
  void init ( const EntityImp &en ) const;

  //! return reference to entity
  const EntityType& en() const 
  {
    assert( en_ ); 
    return *en_;
  }

  //! actual entity 
  mutable const EntityType* en_;

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

  //! do we have the same base function set for all elements
  bool uniform_;

  //! is it initialised?
  mutable bool init_;

  //! corresponding base function set 
  mutable const BaseFunctionSetType* baseSet_; 
  
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

  typedef typename DofStorageType :: block_type DofBlockType;
  typedef DofTypeWrapper<DofBlockType> DofWrapperType;
  enum { blockSize = DofWrapperType :: blockSize };
  
  //! Default constructor
  DofIteratorStaticDiscreteFunction() :
    dofArray_ (0) ,
    count_(0),
    idx_(0)
  {}

  //! Constructor (const)
  DofIteratorStaticDiscreteFunction ( const DofStorageType & dofArray , int count )
    :  dofArray_ (const_cast<DofStorageType*>(&dofArray)) ,
       count_(count),
       idx_(0) {}
  
  //! Constructor
  DofIteratorStaticDiscreteFunction(DofStorageType& dofArray, int count)
    : dofArray_(&dofArray),
      count_(count),
      idx_(0) {}

  //! Copy Constructor
  DofIteratorStaticDiscreteFunction (const ThisType& other)
    : dofArray_(other.dofArray_)
    , count_(other.count_) 
    , idx_(other.idx_)
  {}

  //! Assignment operator
  ThisType& operator=(const ThisType& other)
  {
    if (&other != this) 
    {
      dofArray_ = other.dofArray_;
      count_ = other.count_;
      idx_   = other.idx_;
    }
    return *this;
  }

  //! return dof
  DofType & operator *()
  {
    assert((count_ >=0) && (count_ < dofArray_->size()));
    //return DofWrapperType::convert((*dofArray_)[count_],idx_);
    return ((*dofArray_)[count_][idx_]);
  }

  //! return dof read only 
  const DofType & operator * () const
  {
    assert((count_ >=0) && (count_ < dofArray_->size()));
    //return DofWrapperType::convert((*dofArray_)[count_],idx_);
    return ((*dofArray_)[count_][idx_]);
  }

  //! go next dof
  ThisType & operator++ ()
  {
    ++idx_;
    if(idx_ >= blockSize) 
    {
      idx_ = 0;
      ++count_;
    }
    return (*this);
  }
  
  //! compare
  bool operator == (const ThisType & I ) const
  {
    return (count_ == I.count_) && (idx_ == I.idx_);
  }

  //! compare 
  bool operator != (const ThisType & I ) const
  {
    return !((*this) == I);
  }

  //! return actual index 
  int index () const { return count_; }

  //! set dof iterator back to begin , for const and not const Iterators
  void reset () { count_ = 0; idx_ = 0; }
  
private:
  //! the array holding the dofs 
  DofStorageType * dofArray_;  
  
  //! index 
  mutable int count_;
  mutable int idx_;

}; // end DofIteratorStaticDiscreteFunction 

} // end namespace Dune

#include "staticfunction.cc"

#endif


