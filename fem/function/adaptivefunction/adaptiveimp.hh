#ifndef DUNE_ADAPTIVEFUNCTIONIMP_HH
#define DUNE_ADAPTIVEFUNCTIONIMP_HH

//- System includes
#include <string>
#include <iostream>
#include <fstream>

//- Dune includes
#include <dune/common/exceptions.hh>

#include <dune/fem/io/file/xdrio.hh>
#include <dune/fem/storage/envelope.hh>

namespace Dune {

  //- Forward declarations
  template <class DiscreteFunctionSpaceImp>
  class AdaptiveDiscreteFunction;
  template <class DiscreteFunctionSpaceImp>
  class AdaptiveLocalFunction;
  template <class DiscreteFunctionSpaceImp>
  class AdaptiveDiscreteFunctionTraits;
  template <class DiscreteFunctionSpaceImp>
  class AdaptiveLocalFunctionTraits;
  

  template <class DiscreteFunctionSpaceImp>
  class AdaptiveFunctionImplementation {
  public:
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

    typedef AdaptiveDiscreteFunctionTraits<DiscreteFunctionSpaceImp> Traits;
  private:
    typedef AdaptiveFunctionImplementation<DiscreteFunctionSpaceImp> ThisType;
    
    typedef typename Traits::DofIteratorType DofIteratorType;
    typedef typename Traits::ConstDofIteratorType ConstDofIteratorType;
    
  /*
  public:   
    typedef typename Traits::LocalFunctionImp LocalFunctionImp;
    typedef typename Traits::LocalFunctionType LocalFunctionType;
  */
  private:
    typedef typename Traits::MapperType MapperType;

    typedef typename DiscreteFunctionSpaceImp::Traits SpaceTraits;
    typedef typename SpaceTraits::RangeFieldType RangeFieldType;

    // dof array that is resizeable 
    typedef typename Traits::MutableDofStorageType MutableDofStorageType;
    // static dof array to be used by discrete function 
    typedef typename Traits::DofStorageType DofStorageType;
    typedef typename Traits::GridType GridType;
    typedef DofManager<GridType> DofManagerType;

    typedef typename Traits::DiscreteFunctionType LeafType;

    typedef typename Traits::DofType DofType;
    
  protected:
    template< class Dof, unsigned int Size >
    class DofBlockProxy;

  public:
    enum { blockSize = DiscreteFunctionSpaceType :: localBlockSize };
    typedef DofBlockProxy< DofType, blockSize > DofBlockType;
    typedef DofBlockProxy< const DofType, blockSize > ConstDofBlockType;
    typedef Envelope< DofBlockType > DofBlockPtrType;
    typedef Envelope< ConstDofBlockType > ConstDofBlockPtrType; 

  public:
    const std :: string &name () const;
    int size() const;

    DofIteratorType dbegin();
    DofIteratorType dend();
    ConstDofIteratorType dbegin() const;
    ConstDofIteratorType dend() const;
    
#if 0
    //! return local function object 
    template <class EntityType>
    LocalFunctionType localFunction(const EntityType& en);

    template <class EntityType>
    const LocalFunctionType localFunction(const EntityType& en) const;
#endif
   
    //! Set all elements to zero
    void clear();

    //! daxpy operation
    void addScaled(const ThisType& g,
                   const RangeFieldType& c);

    //! Assignment
    void assignFunction(const ThisType& org);
  
    //! operator += 
    void addFunction(const ThisType& org);
  
    //! operator -= 
    void substractFunction(const ThisType& org);
  
    //! write data of discrete function to file filename|timestep 
    //! with xdr methods 
    virtual bool write_xdr(const std::string filename) const;

    //! write data of discrete function to file filename|timestep 
    //! with xdr methods 
    virtual bool read_xdr(const std::string filename);
    
    //! write function data to file filename|timestep in ascii Format
    virtual bool write_ascii(const std::string filename) const;
    
    //! read function data from file filename|timestep in ascii Format
    virtual bool read_ascii(const std::string filename);

    //! write function data in pgm fromat file
    virtual bool write_pgm(const std::string filename) const;
    
    //! read function data from pgm fromat file
    virtual bool read_pgm(const std::string filename); 
    
    //! return pointer to underlying array 
    DofType       * leakPointer ()       { return dofVec_.leakPointer(); }
    //! return pointer to underlying array 
    const DofType * leakPointer () const { return dofVec_.leakPointer(); }

    inline ConstDofBlockPtrType block ( unsigned int index ) const
    {
      return ConstDofBlockPtrType( leakPointer() + (blockSize * index) );
    }
    
    inline DofBlockPtrType block ( unsigned int index )
    {
      return DofBlockPtrType( leakPointer() + (blockSize * index) );
    } 
    
  protected:
#if 0
    //! return pointer to local function implementation 
    LocalFunctionImp* newObject () const;
#endif

    //! return reference to dof storage 
    DofStorageType& dofStorage() { return dofVec_; }

    //! normal constructor creating discrete function 
    AdaptiveFunctionImplementation(std::string name,
                                   const DiscreteFunctionSpaceType& spc);

    // Constructor getting vector from outside 
    template <class VectorPointerType>
    AdaptiveFunctionImplementation(std::string name,
                                   const DiscreteFunctionSpaceType& spc,
                                   VectorPointerType * vector);
    
    //! create adaptive discrete function with name, space and vector
    AdaptiveFunctionImplementation(std::string name,
                                   const DiscreteFunctionSpaceType& spc,
                                   DofStorageType& dofVec);
    
    //! copy constructor 
    AdaptiveFunctionImplementation(const ThisType& other);
    //! destructor 
    virtual ~AdaptiveFunctionImplementation();

    //! enable dof compressiion for this discrete function 
    void enableDofCompression(); 

  protected:
    virtual const LeafType& interface() const = 0;
    const DiscreteFunctionSpaceType& spc_;
    std::string name_;
    DofManagerType& dm_;
    std::pair<MemObjectInterface*, DofStorageType*> memPair_; 

  protected:
    DofStorageType& dofVec_;
  }; // end class AdaptiveFunctionImplementation


  
  template< class DiscreteFunctionSpace >
  template< class Dof, unsigned int Size >
  class AdaptiveFunctionImplementation< DiscreteFunctionSpace > :: DofBlockProxy
  {
    friend class Envelope< DofBlockProxy >;

  public:
    typedef Dof DofType;

    enum { size = Size };

    typedef std :: size_t size_type;

  protected:
    DofType *const dofBlock_;

  protected:
    inline DofBlockProxy ( DofType *const dofBlock )
    : dofBlock_( dofBlock )
    {}

    inline DofBlockProxy ( const DofBlockProxy &other )
    : dofBlock_( other.dofBlock_ )
    {
    }

  public:
    inline DofBlockProxy &operator= ( const DofBlockProxy &other )
    {
      for( size_type i = 0; i < size; ++i )
        (*this)[ i ] = other[ i ];
      return *this;
    }
    
    inline const DofType &operator[] ( size_type index ) const
    {
      return dofBlock_[ index ];
    }

    inline DofType &operator[] ( size_type index )
    {
      return dofBlock_[ index ];
    }

    inline size_type dim () const
    {
      return size;
    }
  };

} // end namespace Dune

#include "adaptiveimp.cc"

#endif
