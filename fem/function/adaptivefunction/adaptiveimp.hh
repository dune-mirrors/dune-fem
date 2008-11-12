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
#include <dune/fem/space/common/dofmanager.hh>

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
    int size() const;

    DofIteratorType dbegin();
    DofIteratorType dend();
    ConstDofIteratorType dbegin() const;
    ConstDofIteratorType dend() const;
    
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
  
#if DUNE_FEM_COMPATIBILITY
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
#endif

    //! return pointer to underlying array 
    DofType *leakPointer ()
    {
      return dofStorage().leakPointer();
    }

    //! return pointer to underlying array 
    const DofType *leakPointer () const
    {
      return dofStorage().leakPointer();
    }

    inline ConstDofBlockPtrType block ( unsigned int index ) const
    {
      assert( blockSize * (index + 1) <= (unsigned int)size() );
      return ConstDofBlockPtrType( leakPointer() + (blockSize * index) );
    }
    
    inline DofBlockPtrType block ( unsigned int index )
    {
      assert( blockSize * (index + 1) <= (unsigned int)size() );
      return DofBlockPtrType( leakPointer() + (blockSize * index) );
    } 
    
  protected:
    //! return reference to dof storage 
    const DofStorageType &dofStorage () const
    {
      return dofVec_;
    }

    //! return reference to dof storage 
    DofStorageType &dofStorage ()
    {
      return dofVec_;
    }

    //! normal constructor creating discrete function 
    AdaptiveFunctionImplementation ( const std :: string &name,
                                     const DiscreteFunctionSpaceType &spc );

    // Constructor getting vector from outside 
    template <class VectorPointerType>
    AdaptiveFunctionImplementation ( const std :: string &name,
                                     const DiscreteFunctionSpaceType &spc,
                                     VectorPointerType *vector );
    
    //! create adaptive discrete function with name, space and vector
    AdaptiveFunctionImplementation( const std :: string &name,
                                    const DiscreteFunctionSpaceType &spc,
                                    DofStorageType& dofVec );
    
    //! copy constructor 
    AdaptiveFunctionImplementation( const std :: string &name,
                                    const ThisType& other );
    //! destructor 
    virtual ~AdaptiveFunctionImplementation();

    //! enable dof compressiion for this discrete function 
    void enableDofCompression(); 

  protected:

    //! wrapper class to create fake DofStorage from double* 
    template <class VectorPointerType>
    class DofStorageWrapper : public DofStorageInterface
    {
      const std::string name_;
      DofStorageType array_;
    public:
      template <class MapperType>
      DofStorageWrapper(const MapperType& mapper,
                        const std::string& name,
                        const VectorPointerType* v)
        : name_(name),
          array_( mapper.size() , const_cast<VectorPointerType *> (v))
      {}

      //! returns name of this vector 
      const std::string& name () const { return name_; }

      //! return array 
      DofStorageType& getArray() { return array_; }

      //! do nothing here 
      void enableDofCompression () {}

      //! return array's size 
      int size() const { return array_.size(); }
    };

  protected:
    // allocate unmanaged dof storage 
    template <class MapperType, class VectorPointerType>
    DofStorageType& 
    allocateDofStorageWrapper(const MapperType& mapper,
                              const std::string& name, 
                              const VectorPointerType* v)
    {                         
      typedef DofStorageWrapper<VectorPointerType> DofStorageWrapperType;
      DofStorageWrapperType* dsw = new DofStorageWrapperType(mapper,name,v);
      assert( dsw );
      
      // save pointer to object 
      memObject_ = dsw; 
      // return array  
      return dsw->getArray ();
    } 
    
    // allocate managed dof storage 
    DofStorageType& allocateDofStorage(const std::string& name)
    {
      if( memObject_ != 0)
      {
        DUNE_THROW(InvalidStateException,"DofStorage already allocated!");
      } 
      
      // create memory object 
      std::pair< DofStorageInterface* , DofStorageType* > memPair
         = allocateManagedDofStorage( spc_.grid(), spc_.mapper(),
                                      name, (MutableDofStorageType *) 0);
                                      
      // save pointer 
      memObject_ = memPair.first;
      return *(memPair.second);
    }
 
    virtual const LeafType& interface() const = 0;
    const DiscreteFunctionSpaceType& spc_;
    DofStorageInterface* memObject_;

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
