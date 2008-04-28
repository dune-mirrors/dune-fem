#ifndef DUNE_MATRIXSINGLETON_HH
#define DUNE_MATRIXSINGLETON_HH

#include <dune/fem/storage/singletonlist.hh>

namespace Dune {

//! Key for CommManager singleton list 
template <class RowSpaceImp, class ColSpaceImp>
class MatrixObjectSingletonKey
{
  const RowSpaceImp & rowSpace_;
  const ColSpaceImp & colSpace_;
  const std::string paramfile_;
public: 
  //! constructor taking space 
  MatrixObjectSingletonKey(const RowSpaceImp & rowSpace,
                    const ColSpaceImp & colSpace,
                    const std::string& paramfile)
    : rowSpace_(rowSpace)
    , colSpace_(colSpace)
    , paramfile_(paramfile) {}

  //! copy constructor  
  MatrixObjectSingletonKey(const MatrixObjectSingletonKey & org) 
    : rowSpace_(org.rowSpace_) 
    , colSpace_(org.colSpace_) 
    , paramfile_(org.paramfile_)
  {}

  //! returns true if indexSet pointer and numDofs are equal 
  bool operator == (const MatrixObjectSingletonKey & otherKey) const
  {
    // mapper of space is singleton 
    return (  (&(rowSpace_.mapper()) == & (otherKey.rowSpace_.mapper())) && 
              (&(colSpace_.mapper()) == & (otherKey.colSpace_.mapper())) );
  }

  //! return reference to index set 
  const RowSpaceImp & rowSpace() const { return rowSpace_; }
  const ColSpaceImp & colSpace() const { return colSpace_; }
  const std::string& paramfile () const { return paramfile_; }
};

//! Factory class for SingletonList to tell how objects are created and
//! how compared.
template <class KeyImp, class ObjectImp>
class MatrixObjectSingletonFactory
{
  public:
  //! create new matrix object    
  static ObjectImp * createObject( const KeyImp & key )
  {
    return new ObjectImp(key.rowSpace(),key.colSpace(),key.paramfile());
  }

  //! delete matrix object 
  static void deleteObject( ObjectImp * obj )
  {
    delete obj;
  }
};

//! MatrixObject handling an istl matrix 
template <class MatrixObjectImp> 
class MatrixObjectSingleton
{
public: 
  //! type of matrix object 
  typedef MatrixObjectImp MatrixObjectType;
  //! type of space defining row structure
  typedef typename MatrixObjectType :: RowSpaceType RowSpaceType; 
  //! type of space defining column structure
  typedef typename MatrixObjectType :: ColumnSpaceType ColumnSpaceType; 

  //! type of key for singleton list 
  typedef MatrixObjectSingletonKey<RowSpaceType,ColumnSpaceType> KeyType;
  //! factory to create matrix object 
  typedef MatrixObjectSingletonFactory<KeyType, MatrixObjectType> FactoryType;

  //! matrix object provider 
  typedef SingletonList< KeyType , MatrixObjectType , FactoryType > MatrixObjectProviderType;

  //! type of used matrix 
  typedef typename MatrixObjectType :: MatrixType MatrixType; 

  //! type of preconditioner 
  typedef typename MatrixObjectType :: PreconditionMatrixType PreconditionMatrixType; 

  //! type of local matrix 
  typedef typename MatrixObjectType :: LocalMatrixType LocalMatrixType;

  //! constructor recieving reference to matrix object 
  MatrixObjectSingleton(const RowSpaceType & rowSpace,
                        const ColumnSpaceType & colSpace,
                        const std::string& paramfile)
  : key_(rowSpace,colSpace,paramfile)
  , matrixObj_( MatrixObjectProviderType :: getObject( key_ ))
  {}

  //! desctrutor 
  ~MatrixObjectSingleton() 
  { 
    MatrixObjectProviderType :: removeObject( matrixObj_ );
  }

  //! return reference to system matrix
  MatrixType & matrix() const { return matrixObj_.matrix(); }
  //! return true if preconditioner exists 
  bool hasPreconditionMatrix() const { return matrixObj_.hasPreconditionMatrix(); }
  //! return reference to preconditioner 
  const PreconditionMatrixType& preconditionMatrix () const { return matrixObj_.preconditionMatrix(); }

  //! reserve memory for matrix and build structure 
  void reserve(bool verbose = false) { matrixObj_.reserve(verbose); }
  //! set all matrix entitries to zero 
  void clear() { matrixObj_.clear(); }

  //! apply matrix multiplication  
  void multOEM(const double * arg, double * dest) const
  {
    matrixObj_.multOEM(arg,dest);
  }

  //! resort matrix
  void resort() {matrixObj_.resort(); }
  //! print matrix  
  void print(std::ostream & s) const { matrixObj_.print(s); }
  //! create preconditioner, only for some matrices 
  void createPreconditionMatrix() { matrixObj_.createPreconditionMatrix(); }
private:
  // make copy constructor private
  MatrixObjectSingleton(const  MatrixObjectSingleton&);
  // key for obtaining matrix reference
  KeyType key_;
  // reference to matrix object 
  MatrixObjectType& matrixObj_;
};

} // end namespace Dune 
#endif
