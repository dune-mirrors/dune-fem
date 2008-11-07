#ifndef DUNE_FEM_BLOCKVECTORFUNCTION_HH
#define DUNE_FEM_BLOCKVECTORFUNCTION_HH

//- system includes
#include <iostream>
#include <fstream>

//- Dune inlcudes 
#include <dune/common/exceptions.hh>
#include <dune/fem/space/common/arrays.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/dgspace/dgmapper.hh>

#if HAVE_DUNE_ISTL 
#include <dune/istl/bvector.hh>
#else 
#include <dune/fem/space/common/arrays.hh>
#endif
#include <dune/fem/io/file/xdrio.hh>

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/common/dofiterator.hh>
#include <dune/fem/function/localfunction/standardlocalfunction.hh>

namespace Dune
{

// forward declarations 
template <class DiscreteFunctionSpaceType> class BlockVectorDiscreteFunction;
template <class DofStorageImp,class DofImp> class DofIteratorBlockVectorDiscreteFunction;
template< class Traits > class BlockVectorLocalFunction;
template< class Traits > class BlockVectorLocalFunctionFactory;
template <class BlockVectorImp, class DofImp> class StraightenBlockVector;

template <class DiscreteFunctionSpaceImp>
struct BlockVectorDiscreteFunctionTraits 
{
  enum { localBlockSize = DiscreteFunctionSpaceImp :: localBlockSize };
  typedef typename DiscreteFunctionSpaceImp :: RangeFieldType RangeFieldType;

  typedef FieldVector<RangeFieldType, localBlockSize > DofBlockType;
  typedef const DofBlockType ConstDofBlockType;
  typedef DofBlockType *DofBlockPtrType;
  typedef const DofBlockType *ConstDofBlockPtrType;
  
#if HAVE_DUNE_ISTL
  typedef BlockVector< DofBlockType > DofStorageType;
#else 
  typedef MutableArray < DofBlockType > DofStorageType;
#endif

  typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
  typedef typename DiscreteFunctionSpaceType :: IndexSetType IndexSetType;

  //! needs additional mapper because of block structure 
  typedef typename DiscreteFunctionSpaceType :: BlockMapperType MapperType;

  typedef BlockVectorDiscreteFunction<DiscreteFunctionSpaceType> DiscreteFunctionType;
 
  typedef DofIteratorBlockVectorDiscreteFunction<DofStorageType,RangeFieldType> DofIteratorType;
  typedef ConstDofIteratorDefault<DofIteratorType> ConstDofIteratorType;

  typedef BlockVectorDiscreteFunctionTraits<DiscreteFunctionSpaceImp> ThisType;
  
  typedef StandardLocalFunctionFactory< ThisType > LocalFunctionFactoryType;

  typedef LocalFunctionStack< LocalFunctionFactoryType > LocalFunctionStorageType;
  typedef typename LocalFunctionStorageType :: LocalFunctionType LocalFunctionType;

  typedef StraightenBlockVector<DofStorageType,RangeFieldType> LeakPointerType;
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

template <class BlockVectorImp, class DofImp>
class StraightenBlockVector
{
  typedef BlockVectorImp BlockVectorType;
  typedef DofImp DofType;
  typedef typename BlockVectorType :: block_type BlockType;
  enum { blockSize = BlockType :: dimension };
public:
  StraightenBlockVector(BlockVectorImp& vec) : vec_(vec) {}

  //! return dof
  DofType& operator [] (int i)
  {
    assert((i >=0) && (i < blockSize * vec_.size()));
    const int count = (int) i/blockSize;
    assert( count < vec_.size() );
    const int local = i%blockSize;
    assert( local < blockSize );
    return vec_[count][local];
  }

  //! return dof
  const DofType& operator [] (int i) const
  {
    assert((i >=0) && (i < blockSize * vec_.size()));
    const int count = (int) i/blockSize;
    assert( count < vec_.size() );
    const int local = i%blockSize;
    assert( local < blockSize );
    return vec_[count][local];
  }

  int size () const
  {
    return vec_.size() * blockSize;
  }

  BlockVectorType& blockVector() { return vec_; }
  const BlockVectorType& blockVector() const { return vec_; }
protected:
  BlockVectorType& vec_;
};



//**********************************************************************
//! @ingroup BlockVectorDFunction
//  --BlockVectorDiscreteFunction 
//
//! this is one special implementation of a discrete function using an
//! array for storing the dofs.  
//!
//**********************************************************************
template< class DiscreteFunctionSpaceImp >
class BlockVectorDiscreteFunction 
: public DiscreteFunctionDefault
  < BlockVectorDiscreteFunctionTraits< DiscreteFunctionSpaceImp > > 
{
  typedef BlockVectorDiscreteFunction< DiscreteFunctionSpaceImp > ThisType;
  typedef DiscreteFunctionDefault
    < BlockVectorDiscreteFunctionTraits< DiscreteFunctionSpaceImp > >
    BaseType;

  typedef DiscreteFunctionDefault<BlockVectorDiscreteFunctionTraits <DiscreteFunctionSpaceImp> >
  DiscreteFunctionDefaultType;

public:  
  //! type of discrete functions space 
  typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

  //! traits of this type 
  typedef BlockVectorDiscreteFunctionTraits<DiscreteFunctionSpaceType> Traits;
  
  //! size of local blocks  
  enum { localBlockSize = Traits :: localBlockSize };

  //! needs additional mapper 
  typedef typename Traits :: MapperType MapperType; 

  friend class BlockVectorLocalFunctionFactory< Traits >;

private:
  enum { myId_ = 0};
  
public:
  //! type of underlying array
  typedef typename Traits :: DofStorageType DofStorageType;

  //! my type 
  typedef BlockVectorDiscreteFunction <DiscreteFunctionSpaceType> DiscreteFunctionType;
  
  //! Type of the range field
  typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;

  /** \brief For ISTL-compatibility */
  typedef typename DofStorageType :: block_type block_type; 

  //! Type of the grid
  typedef typename DiscreteFunctionSpaceType::GridType GridType;

  //! type of local function factory 
  typedef typename Traits :: LocalFunctionFactoryType LocalFunctionFactoryType;

  //! LocalFunctionType is the exported lf type 
  typedef typename LocalFunctionFactoryType :: ObjectType
    LocalFunctionImp;

  //! the dof iterator type of this function
  typedef typename Traits :: DofIteratorType DofIteratorType;
  typedef typename Traits :: ConstDofIteratorType ConstDofIteratorType;

  //! The associated discrete function space
  typedef DiscreteFunctionSpaceType FunctionSpaceType;

  //! the type of the unknowns 
  typedef RangeFieldType DofType;

  //! type of block stored in block vector 
  typedef block_type DofBlockType;

  typedef typename Traits :: ConstDofBlockType ConstDofBlockType;
  typedef typename Traits :: DofBlockPtrType DofBlockPtrType;
  typedef typename Traits :: ConstDofBlockPtrType ConstDofBlockPtrType;
  
  //! type of index set 
  typedef typename DiscreteFunctionSpaceType :: IndexSetType IndexSetType; 
  
  //! type of LeakPointer 
  typedef typename Traits :: LeakPointerType LeakPointerType;

public:
  //! \brief Constructor makes Discrete Function  
  BlockVectorDiscreteFunction ( const DiscreteFunctionSpaceType &f );
  
  //! \brief Constructor makes Discrete Function with name 
  BlockVectorDiscreteFunction ( const std :: string &name,
                                const DiscreteFunctionSpaceType &f );
  
  //! \brief Constructor makes Discrete Function  
  BlockVectorDiscreteFunction ( const std :: string &name,
                                const DiscreteFunctionSpaceType &f,
                                const DofStorageType &data );
  
  //! \brief Constructor makes Discrete Function from copy 
  BlockVectorDiscreteFunction ( const ThisType &other ); 

  /**  \brief delete stack of free local functions belonging to this discrete function */
  ~BlockVectorDiscreteFunction ();

  using BaseType :: name;

#if 0
  /** \brief @copydoc DiscreteFunctionInterface::localFunction */ 
  template <class EntityType>
  LocalFunctionType localFunction(const EntityType& en) const;
#endif

  /** \copydoc Dune::DiscreteFunctionInterface::dbegin() */ 
  DofIteratorType dbegin ();
  
  /** \copydoc Dune::DiscreteFunctionInterface::dend() */
  DofIteratorType dend (); 

  /** \copydoc Dune::DiscreteFunctionInterface::dbegin() const */
  ConstDofIteratorType dbegin () const;
  
  /** \copydoc Dune::DiscreteFunctionInterface::dend() const */ 
  ConstDofIteratorType dend () const;

  inline ConstDofBlockPtrType block ( const unsigned int block ) const
  {
    return &(dofVec_[ block ]);
  }
  
  inline DofBlockPtrType block ( const unsigned int block )
  {
    return &(dofVec_[ block ]);
  }

  inline const RangeFieldType &dof ( const unsigned int index ) const
  {
    return leakPtr_[ index ];
  }

  inline RangeFieldType &dof ( const unsigned int index )
  {
    return leakPtr_[ index ];
  }

  /** \copydoc Dune::DiscreteFunctionInterface::size */
  inline int size() const
  {
    return dofVec_.size();
  }

  /** \copydoc Dune::DiscreteFunctionDefault::clear */
  void clear();

  /** \copydoc Dune::DiscreteFunctionDefault::addScaled */
  void addScaled ( const DiscreteFunctionType &g,
                   const RangeFieldType &s ); 
 
  template< class GridIteratorType >
  void DUNE_DEPRECATED addLocal ( GridIteratorType &it,
                                  const DiscreteFunctionType &g );
  
  template< class GridIteratorType >
  void DUNE_DEPRECATED substractLocal ( GridIteratorType &it,
                                        const DiscreteFunctionType &g ); 
  
  template< class GridIteratorType >
  void DUNE_DEPRECATED setLocal ( GridIteratorType &it,
                                  const RangeFieldType &scalar );
  
  /** \copydoc Dune::DiscreteFunctionInterface::print */
  void print( std :: ostream &out ) const;

#if DUNE_FEM_COMPATIBILITY
  /** \copydoc Dune::DiscreteFunctionDefault::write_xdr */
  virtual bool write_xdr( const std::string filename ) const;

  /** \copydoc Dune::DiscreteFunctionDefault::read_xdr */
  virtual bool read_xdr( const std::string filename );

  /** \copydoc Dune::DiscreteFunctionDefault::write_ascii */
  virtual bool write_ascii( const std::string filename ) const;

  /** \copydoc Dune::DiscreteFunctionDefault::read_ascii */
  virtual bool read_ascii( const std::string filename );
#endif

  /** \brief return reference to internal block vector 
      \return reference to blockVector */ 
  DofStorageType& blockVector () const { return dofVec_; }

  /** \brief return reference to leak pointer 
      \return reference to leakPointer */ 
  LeakPointerType& leakPointer() { return leakPtr_; }

  /** \brief return const reference to leak pointer 
      \return constant reference to leakPointer */
  const LeakPointerType& leakPointer() const { return leakPtr_; }

  /** \copydoc Dune::DiscreteFunctionInterface::enableDofCompression() */
  void enableDofCompression();

private:  
  // allocates dof storage 
  DofStorageType& allocateDofStorage();
  
  LocalFunctionFactoryType lfFactory_;

  //! write/read data to/from xdr stream 
  bool processXdrs(XDRStream& xdr) const;
  
  // single mapper for blocks 
  MapperType& mapper_;

  // DofStorage that manages the memory for the dofs of this function
  DofStorageInterface* memObject_;

  //! the dofs stored in an array
  mutable DofStorageType& dofVec_;

  //! leak pointer converting block vector to straight vector 
  LeakPointerType leakPtr_; 

  //! hold one object for addLocal and setLocal and so on 
  LocalFunctionImp localFunc_;
}; // end class BlockVectorDiscreteFunction



//***********************************************************************
//
//  --DofIteratorBlockVectorDiscreteFunction
//
//***********************************************************************
  /** \brief Iterator over an array of dofs 
      \todo Please doc me!
  */
template < class DofStorageImp, class DofImp >
class DofIteratorBlockVectorDiscreteFunction : public
DofIteratorDefault < DofImp , DofIteratorBlockVectorDiscreteFunction < DofStorageImp, DofImp > >
{
public:
  typedef DofImp DofType;
  typedef DofStorageImp DofStorageType;
  typedef DofIteratorBlockVectorDiscreteFunction<DofStorageType,DofType> ThisType;

  typedef typename DofStorageType :: block_type DofBlockType;
  typedef DofTypeWrapper<DofBlockType> DofWrapperType;
  enum { blockSize = DofWrapperType :: blockSize };
  
  //! Default constructor
  DofIteratorBlockVectorDiscreteFunction() :
    dofArray_ (0) ,
    count_(0),
    idx_(0)
  {}

  //! Constructor (const)
  DofIteratorBlockVectorDiscreteFunction ( const DofStorageType & dofArray , int count )
    :  dofArray_ (const_cast<DofStorageType*>(&dofArray)) ,
       count_(count),
       idx_(0) {}
  
  //! Constructor
  DofIteratorBlockVectorDiscreteFunction(DofStorageType& dofArray, int count)
    : dofArray_(&dofArray),
      count_(count),
      idx_(0) {}

  //! Copy Constructor
  DofIteratorBlockVectorDiscreteFunction (const ThisType& other)
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
  DofType& operator *()
  {
    assert((count_ >=0) && (count_ < dofArray_->size()));
    //return DofWrapperType::convert((*dofArray_)[count_],idx_);
    return ((*dofArray_)[count_][idx_]);
  }

  //! return dof read only 
  const DofType& operator * () const
  {
    assert((count_ >=0) && (count_ < dofArray_->size()));
    //return DofWrapperType::convert((*dofArray_)[count_],idx_);
    return ((*dofArray_)[count_][idx_]);
  }

  //! go to next dof
  ThisType& operator++ ()
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

}; // end DofIteratorBlockVectorDiscreteFunction 
template <class DiscreteFunctionSpaceImp>
class ManagedDiscreteFunction<BlockVectorDiscreteFunction<DiscreteFunctionSpaceImp> > :
public BlockVectorDiscreteFunction<DiscreteFunctionSpaceImp> {
  typedef BlockVectorDiscreteFunction<DiscreteFunctionSpaceImp> BaseType;
public:
  typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;
  typedef ManagedDiscreteFunction<BaseType> ThisType;
  //! \brief Constructor makes Discrete Function  
  ManagedDiscreteFunction ( const DiscreteFunctionSpaceType & f ) : BaseType(f) {}
  //! \brief Constructor makes Discrete Function with name 
  ManagedDiscreteFunction ( const std::string name, const DiscreteFunctionSpaceType & f ) : BaseType(name,f) {}
  //! \brief Constructor makes Discrete Function  
  ManagedDiscreteFunction ( const std::string name, const DiscreteFunctionSpaceType & f, const typename BaseType::DofStorageType & data ) : BaseType(name,f,data) {}
  //! \brief Constructor makes Discrete Function from copy 
  ManagedDiscreteFunction (const ThisType & df) : BaseType(df) {}
};

} // end namespace Dune

#include "blockvectorfunction_inline.hh"
#endif
