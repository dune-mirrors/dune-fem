#ifndef DUNE_STATICFUNCTION_HH
#define DUNE_STATICFUNCTION_HH

//- system includes
#include <fstream>

#if ! HAVE_DUNE_ISTL 
#error "Dune-ISTL is needed for this type of discrete function! Re-configure with --with-dune-istl! "
#endif

//- Dune inlcudes 
#include <dune/fem/space/common/arrays.hh>
#include <dune/fem/space/dgspace/dgmapper.hh>
#include <dune/istl/bvector.hh>
#include <dune/fem/io/file/xdrio.hh>

//- local includes 
#include "../common/discretefunction.hh"
#include "../common/localfunction.hh"
#include "../common/dofiterator.hh"
#include "../common/localfunctionwrapper.hh"

namespace Dune{


// forward declarations 
template <class DiscreteFunctionSpaceType> class BlockVectorDiscreteFunction;
template <class DofStorageImp,class DofImp> class DofIteratorBlockVectorDiscreteFunction;
template< class Traits > class BlockVectorLocalFunction;
template< class Traits > class BlockVectorLocalFunctionFactory;


template <class DiscreteFunctionSpaceImp>
struct BlockVectorDiscreteFunctionTraits 
{
  enum { localBlockSize = DiscreteFunctionSpaceImp :: localBlockSize };
  typedef typename DiscreteFunctionSpaceImp :: RangeFieldType RangeFieldType;
  typedef BlockVector< FieldVector<RangeFieldType, localBlockSize > > DofStorageType;
  typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType :: IndexSetType IndexSetType;

  //! needs additional mapper 
  typedef DGMapper<IndexSetType,0,1> MapperType;


  
  typedef BlockVectorDiscreteFunction<DiscreteFunctionSpaceType> DiscreteFunctionType;
  //typedef BlockVectorLocalFunction<DiscreteFunctionType> LocalFunctionImp;
 
  typedef DofIteratorBlockVectorDiscreteFunction<DofStorageType,
            typename DofStorageType::field_type> DofIteratorType;
  typedef ConstDofIteratorDefault<DofIteratorType> ConstDofIteratorType;


  typedef BlockVectorDiscreteFunctionTraits<DiscreteFunctionSpaceImp> ThisType;
  typedef BlockVectorLocalFunctionFactory< ThisType >  LocalFunctionFactoryType; 
  typedef LocalFunctionStack< LocalFunctionFactoryType > LocalFunctionStorageType;
  typedef typename LocalFunctionStorageType :: LocalFunctionType LocalFunctionType;
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

private:
  BlockVectorType& vec_;
};



template< class TraitsImp > 
class BlockVectorLocalFunctionFactory
{
public:
  typedef TraitsImp Traits;
  typedef typename Traits :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename Traits :: MapperType MapperType;
  typedef typename Traits :: DofStorageType DofStorageType;

private:
  typedef BlockVectorLocalFunctionFactory<TraitsImp>  ThisType;

  friend class BlockVectorDiscreteFunction< DiscreteFunctionSpaceType >;

public:
  typedef BlockVectorLocalFunction<Traits> ObjectType;

  typedef BlockVectorDiscreteFunction< DiscreteFunctionSpaceType >
    DiscreteFunctionType;

protected:
  DiscreteFunctionType &discreteFunction_;

protected:
  inline explicit BlockVectorLocalFunctionFactory ( DiscreteFunctionType &df )
  : discreteFunction_( df )
  {
  }

public:
  ObjectType *newObject () const
  {
    return new ObjectType( discreteFunction_.space(),
                           discreteFunction_.mapper_ ,
                           discreteFunction_.dofVec_ );
  }
};


//**********************************************************************
//! @ingroup BlockVectorDFunction
//  --BlockVectorDiscreteFunction 
//
//! this is one special implementation of a discrete function using an
//! array for storing the dofs.  
//!
//**********************************************************************
template<class DiscreteFunctionSpaceImp> 
class BlockVectorDiscreteFunction 
: public DiscreteFunctionDefault <BlockVectorDiscreteFunctionTraits<DiscreteFunctionSpaceImp> > 
{
  typedef DiscreteFunctionDefault<BlockVectorDiscreteFunctionTraits <DiscreteFunctionSpaceImp> >
  DiscreteFunctionDefaultType;

public:  
  //! type of discrete functions space 
  typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

  //! traits of this type 
  typedef BlockVectorDiscreteFunctionTraits<DiscreteFunctionSpaceType> Traits;
  
  //! needs additional mapper 
  typedef typename Traits :: MapperType MapperType; 

  friend class BlockVectorLocalFunctionFactory< Traits > ;
private:

  typedef BlockVectorDiscreteFunction <DiscreteFunctionSpaceImp> ThisType;
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

  //! dof manager 
  typedef DofManager<GridType> DofManagerType;
  typedef DofManagerFactory<DofManagerType> DofManagerFactoryType;


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
  
  //! type of index set 
  typedef typename DiscreteFunctionSpaceType :: IndexSetType IndexSetType; 
  
  //! type of LeakPointer 
  typedef StraightenBlockVector<DofStorageType,DofType> LeakPointerType;

  //! \brief Constructor makes Discrete Function  
  BlockVectorDiscreteFunction ( const DiscreteFunctionSpaceType & f ) ;
  
  //! \brief Constructor makes Discrete Function with name 
  BlockVectorDiscreteFunction ( const std::string name, const DiscreteFunctionSpaceType & f ) ;
  
  //! \brief Constructor makes Discrete Function  
  BlockVectorDiscreteFunction ( const std::string name, const DiscreteFunctionSpaceType & f, const DofStorageType & data ) ;
  
  //! \brief Constructor makes Discrete Function from copy 
  BlockVectorDiscreteFunction (const ThisType & df); 

  /**  \brief delete stack of free local functions belonging to this discrete function */
  ~BlockVectorDiscreteFunction ();

#if 0
  /** \brief @copydoc DiscreteFunctionInterface::localFunction */ 
  template <class EntityType>
  LocalFunctionType localFunction(const EntityType& en) const;
#endif

  /** \brief @copydoc DiscreteFunctionInterface::dbegin */ 
  DofIteratorType dbegin (); 
  
  /** \brief @copydoc DiscreteFunctionInterface::dend */ 
  DofIteratorType dend   (); 

  /** \brief @copydoc DiscreteFunctionInterface::dbegin */ 
  ConstDofIteratorType dbegin () const;
  
  /** \brief @copydoc DiscreteFunctionInterface::dend  */ 
  ConstDofIteratorType dend   () const; 

  /** \brief @copydoc DiscreteFunctionInterface::name  */ 
  const std::string& name() const {return name_;} 

  /** \brief @copydoc DiscreteFunctionInterface::size  */ 
  int size() const { return dofVec_.size(); }

  /** \brief @copydoc DiscreteFunctionDefault::clear */
  void clear();

  /** \brief @copydoc DiscreteFunctionDefault::addScaled */
  void addScaled ( const DiscreteFunctionType & g,
      const RangeFieldType &scalar); 
  
  /** \brief add g to this on local entity
      \param[in] GridIteratorType it 
      \param[in] discrete function that is added 
  */
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
  
  /** \brief @copydoc DiscreteFunctionDefault::print */
  void print(std::ostream& s) const;

  /** \brief @copydoc DiscreteFunctionDefault::write_xdr */
  virtual bool write_xdr( const std::string filename ) const;

  /** \brief @copydoc DiscreteFunctionDefault::read-xdr  */
  virtual bool read_xdr( const std::string filename );

  /** \brief @copydoc DiscreteFunctionDefault::write_ascii  */
  virtual bool write_ascii(const std::string filename) const;

  /** \brief @copydoc DiscreteFunctionDefault::read_ascii */
  virtual bool read_ascii(const std::string filename);

  /** \brief @copydoc DiscreteFunctionDefault::write_pgm  */
  virtual bool write_pgm(const std::string filename) const;

  /** \brief @copydoc DiscreteFunctionDefault::read_pgm  */
  virtual bool read_pgm(const std::string filename); 

  /** \brief return reference to internal block vector 
      \return reference to blockVector */ 
  DofStorageType& blockVector () const { return dofVec_; }

  /** \brief return reference to leak pointer 
      \return reference to leakPointer */ 
  LeakPointerType& leakPointer() { return leakPtr_; }

  /** \brief return const reference to leak pointer 
      \return constant reference to leakPointer */
  const LeakPointerType& leakPointer() const { return leakPtr_; }

private:  
  LocalFunctionFactoryType lfFactory_;

  //! write/read data to/from xdr stream 
  bool processXdrs(XDRStream& xdr) const;
  
  //! the name of the function
  std::string name_;

  // single mapper for blocks 
  MapperType mapper_;

  // dof manager 
  DofManagerType&  dm_;
  
  // MemObject that manages the memory for the dofs of this function
  std::pair<MemObjectInterface*, DofStorageType*> memPair_;

  //! true if memory was allocated
  bool built_;

  //! the dofs stored in an array
  mutable DofStorageType& dofVec_;

  //! hold one object for addLocal and setLocal and so on 
  LocalFunctionImp localFunc_;

  //! leak pointer converting block vector to straight vector 
  LeakPointerType leakPtr_; 
}; // end class BlockVectorDiscreteFunction 


//**************************************************************************
//
//  --BlockVectorLocalFunction
//
//! Implementation of the local functions 
//
//**************************************************************************
template< class TraitsImp >
class BlockVectorLocalFunction
: public LocalFunctionDefault
  < typename TraitsImp :: DiscreteFunctionSpaceType,
    BlockVectorLocalFunction< TraitsImp >   > 
{
public:
  typedef TraitsImp Traits; 
  typedef typename Traits :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename Traits :: MapperType  MapperType;
  typedef typename Traits :: DofStorageType DofStorageType;

private:
  typedef BlockVectorLocalFunction< Traits > ThisType;
  
public:
  typedef typename DiscreteFunctionSpaceType::Traits::GridType GridType;
  typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;

  enum { dimRange = DiscreteFunctionSpaceType::DimRange };
  //CompileTimeChecker<dimrange == 1> check; 
  typedef typename DiscreteFunctionSpaceType::Traits::DomainType DomainType;
  typedef typename DiscreteFunctionSpaceType::Traits::RangeType RangeType;
  typedef typename DiscreteFunctionSpaceType::Traits::RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionSpaceType::Traits::JacobianRangeType JacobianRangeType;

  typedef typename GridType :: template  Codim<0> :: Entity EntityType;

  typedef typename DofStorageType :: block_type DofBlockType;

  friend class BlockVectorDiscreteFunction <DiscreteFunctionSpaceType>;

private:
  typedef typename GridType :: ctype ctype;
  enum { dim = GridType :: dimension };
  typedef FieldMatrix<ctype,dim,dim> JacobianInverseType;

public:
  //! Constructor 
  BlockVectorLocalFunction ( const DiscreteFunctionSpaceType &f , 
                                const MapperType& mapper, 
                                DofStorageType & dofVec );

  //! \brief Destructor 
  ~BlockVectorLocalFunction ();

   /** \brief  @copydoc LocalFunctionDefault::operator []  */
  RangeFieldType & operator [] (int num);

   /** \brief  @copydoc LocalFunctionDefault::operator []  */
  const RangeFieldType & operator [] (int num) const;

  /** \brief  @copydoc LocalFunctionDefault::numDofs  */
  int numDofs () const;

  /** \brief  @copydoc LocalFunctionDefault::evaluate  */
  void evaluate (const DomainType & x, RangeType & ret) const ;
 
   /** \brief  @copydoc LocalFunctionDefault::evaluate  */
  template <class QuadratureType>
  void evaluate (const QuadratureType &quad, const int quadPoint , RangeType & ret) const;

   /** \brief  @copydoc LocalFunctionDefault::jacobian  */
  void jacobian(const DomainType& x, JacobianRangeType& ret) const ;

   /** \brief  @copydoc LocalFunctionDefault::jacobian  */
  template <class QuadratureType>
  void jacobian(const QuadratureType &quad, const int quadPoint , JacobianRangeType & ret) const;

   /** \brief  @copydoc LocalFunctionDefault::baseFunctionSet  */
  const BaseFunctionSetType& baseFunctionSet() const;

   /** \brief  @copydoc LocalFunctionDefault::axpy  */
  template <class QuadratureType>
  inline void axpy(const QuadratureType&, const int qp, const RangeType& factor);

   /** \brief  @copydoc LocalFunctionDefault::axpy  */
  template <class QuadratureType>
  inline void axpy(const QuadratureType&, const int qp, const JacobianRangeType& factor);

   /** \brief  @copydoc LocalFunctionDefault::axpy  */
  template <class QuadratureType>
  inline void axpy(const QuadratureType&, const int qp, const RangeType& factor1, const JacobianRangeType& factor2);

  inline void
  rightMultiply(const JacobianRangeType& factor,
              const JacobianInverseType& jInv,
              JacobianRangeType& result) const;
  
  //! update local function for given Entity  
  void init ( const EntityType &en ) const;

protected:
  //! return reference to entity
  const EntityType& en() const 
  {
    assert( en_ ); 
    return *en_;
  }

  //! the corresponding function space which provides the base function set
  const DiscreteFunctionSpaceType &fSpace_;

  //! special mapper of block vector functions
  const MapperType& mapper_;
  
  //! actual entity 
  mutable const EntityType* en_;

  //! dofVec from all levels of the discrete function 
  DofStorageType & dofVec_;

  //! Array holding pointers to the local dofs 
  mutable MutableArray < RangeFieldType * > values_ ;

  //! needed once 
  mutable RangeType tmp_;

  //! needed once 
  mutable JacobianRangeType tmpGrad_;
  mutable JacobianRangeType factorInv_;

  //! diffVar for evaluate, is empty 
  const DiffVariable<0>::Type diffVar;

  //! number of all dofs 
  mutable int numOfDof_;

  //! do we have the same base function set for all elements
  mutable bool needCheckGeometry_;

  //! corresponding base function set 
  mutable BaseFunctionSetType baseSet_; 
}; // end BlockVectorLocalFunction 


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

} // end namespace Dune

#include "blockvectorfunction_inline.hh"
#endif
