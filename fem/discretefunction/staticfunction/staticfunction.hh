#ifndef DUNE_STATICFUNCTION_HH
#define DUNE_STATICFUNCTION_HH

//- system includes
#include <fstream>
#include <rpc/xdr.h>

//- Dune inlcudes 
#include <dune/common/array.hh>
#include <dune/istl/bvector.hh>

#include <dune/fem/space/dgspace/dgmapper.hh>

//- local includes 
#include "../common/discretefunction.hh"
#include "../common/localfunction.hh"
#include "../common/dofiterator.hh"

namespace Dune{

template <class DiscreteFunctionSpaceType,  class DofStorageImp = Array<typename DiscreteFunctionSpaceType::RangeFieldType> > class StaticDiscreteFunction;
template <class DiscreteFunctionType> class StaticDiscreteLocalFunction;
template <class DofStorageImp,class DofImp> class DofIteratorStaticDiscreteFunction;


template <class DiscreteFunctionSpaceImp, class DofStorageImp>
struct StaticDiscreteFunctionTraits 
{
  typedef DofStorageImp DofStorageType;
  typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;
  
  typedef StaticDiscreteFunction<DiscreteFunctionSpaceType,DofStorageType> DiscreteFunctionType;
  typedef StaticDiscreteLocalFunction<DiscreteFunctionType> LocalFunctionImp;
  
  typedef LocalFunctionWrapper<DiscreteFunctionType> LocalFunctionType;
  typedef DofIteratorStaticDiscreteFunction<DofStorageType,
            typename DofStorageType::field_type> DofIteratorType;
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
  //! traits of this type 
  typedef StaticDiscreteFunctionTraits<DiscreteFunctionSpaceType,DofStorageImp> Traits;
  
  //! type of underlying array
  typedef DofStorageImp DofStorageType;

  //! my type 
  typedef StaticDiscreteFunction < DiscreteFunctionSpaceType, DofStorageType > DiscreteFunctionType;

  //! Type of the range field
  typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;

  /** \brief For ISTL-compatibility */
  typedef typename DofStorageType :: block_type block_type; 

  //! Type of the grid
  typedef typename DiscreteFunctionSpaceType::GridType GridType;

  //! dof manager 
  typedef DofManager<GridType> DofManagerType;
  typedef DofManagerFactory<DofManagerType> DofManagerFactoryType;

  //! the local function implementation e 
  typedef typename Traits :: LocalFunctionImp LocalFunctionImp;

  //! LocalFunctionType is the exported lf type 
  typedef typename Traits :: LocalFunctionType LocalFunctionType;

  // the storage of the local functions 
  typedef LocalFunctionStorage< DiscreteFunctionType > LocalFunctionStorageType;
  
  //! the dof iterator type of this function
  typedef typename Traits :: DofIteratorType DofIteratorType;
  typedef typename Traits :: ConstDofIteratorType ConstDofIteratorType;

  //! The associated discrete function space
  typedef DiscreteFunctionSpaceType FunctionSpaceType;

  //! the type of the unknowns 
  typedef RangeFieldType DofType;
  typedef block_type DofBlockType;
  
  typedef typename DiscreteFunctionSpaceType :: IndexSetType IndexSetType; 
  //! needs additional mapper 
  typedef DGMapper<IndexSetType,0,1> MapperType;

  //! Constructor makes Discrete Function  
  StaticDiscreteFunction ( const DiscreteFunctionSpaceType & f ) ;
  
  //! Constructor makes Discrete Function with name 
  StaticDiscreteFunction ( const std::string name, const DiscreteFunctionSpaceType & f ) ;
  
  //! Constructor makes Discrete Function  
  StaticDiscreteFunction ( const std::string name, const DiscreteFunctionSpaceType & f, const DofStorageType & data ) ;
  
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
  bool write_xdr( std::string filename ) const;

  //! write data of discrete function to file filename 
  //! with xdr methods 
  bool read_xdr( std::string filename );

  //! write function data to file filename in ascii Format
  bool write_ascii(std::string filename) const;

  //! read function data from file filename in ascii Format
  bool read_ascii(std::string filename);

  //! write function data in pgm fromat file
  bool write_pgm(std::string filename) const;

  //! read function data from pgm fromat file
  bool read_pgm(std::string filename); 

  //! return pointer to internal array for use of BLAS routines 
  DofType* leakPointer () { return (DofType *) &dofVec_[0][0];  }
  //! return pointer to internal array for use of BLAS routines 
  const DofType* leakPointer () const { return (DofType *) &dofVec_[0][0];  }

  //! return reference to internal block vector 
  DofStorageType& blockVector () const { return dofVec_; }

private:  
  //! write data to xdr stream 
  bool writeXdrs(XDR * xdrs) const;
  //! read data from xdr stream 
  bool readXdrs(XDR * xdrs);
  
  //! return object pointer of type LocalFunctionImp 
  LocalFunctionImp * newLocalFunctionObject () const;

  //! the name of the function
  std::string name_;

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

}; // end class StaticDiscreteFunction 


//**************************************************************************
//
//  --StaticDiscreteLocalFunction
//
//! Implementation of the local functions 
//
//**************************************************************************
template < class DiscreteFunctionImp > 
class StaticDiscreteLocalFunction 
: public LocalFunctionDefault < typename DiscreteFunctionImp :: DiscreteFunctionSpaceType ,
                                StaticDiscreteLocalFunction < DiscreteFunctionImp > > 
{
  typedef  DiscreteFunctionImp DiscreteFunctionType;
  typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;
  typedef StaticDiscreteLocalFunction < DiscreteFunctionType > ThisType;

  enum { dimrange = DiscreteFunctionSpaceType::DimRange };
  //CompileTimeChecker<dimrange == 1> check; 
  typedef typename DiscreteFunctionSpaceType::Traits::DomainType DomainType;
  typedef typename DiscreteFunctionSpaceType::Traits::RangeType RangeType;
  typedef typename DiscreteFunctionSpaceType::Traits::RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionSpaceType::Traits::JacobianRangeType JacobianRangeType;

  typedef typename DiscreteFunctionType :: MapperType MapperType;

  typedef typename DiscreteFunctionSpaceType::Traits::GridType :: template
    Codim<0> :: Entity EntityType;
  typedef typename DiscreteFunctionType :: DofStorageType DofStorageType;

  typedef typename DofStorageType :: block_type DofBlockType;

  friend class StaticDiscreteFunction <DiscreteFunctionSpaceType,DofStorageType>;
  friend class LocalFunctionWrapper < StaticDiscreteFunction <DiscreteFunctionSpaceType,DofStorageType> >;
public:
  //! Constructor 
  StaticDiscreteLocalFunction ( const DiscreteFunctionSpaceType &f , 
                                const MapperType& mapper, 
                                DofStorageType & dofVec );

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

  //! the corresponding function space which provides the base function set
  const DiscreteFunctionSpaceType &fSpace_;

  const MapperType& mapper_;
  
  //! actual entity 
  mutable const EntityType* en_;

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


