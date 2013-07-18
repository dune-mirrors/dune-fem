#ifndef DUNE_FEM_ADAPTIVEFUNCTIONIMP_HH
#define DUNE_FEM_ADAPTIVEFUNCTIONIMP_HH

//- System includes
#include <string>
#include <iostream>
#include <fstream>

//- Dune includes
#include <dune/common/exceptions.hh>

#include <dune/fem/storage/envelope.hh>
#include <dune/fem/space/common/dofmanager.hh>

namespace Dune 
{

  namespace Fem 
  { 

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
      
    private:
      typedef typename Traits::MapperType MapperType;

      typedef typename DiscreteFunctionSpaceImp::RangeFieldType RangeFieldType;

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

      //! axpy operation
      void axpy(const RangeFieldType& c, const ThisType& g);

      //! Assignment
      void assignFunction(const ThisType& org);
    
      //! operator += 
      void addFunction(const ThisType& org);
    
      //! operator -= 
      void substractFunction(const ThisType& org);
    
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
#ifndef NDEBUG
        const size_t s = size();
        const size_t bs = blockSize ;
        assert( bs * (index + 1) <= (unsigned int)s );
#endif
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

      //! needed for p-adaptation, (inofficial method, do not use) 
      void resize( )
      {
        // cast to ManagedDofStorageInterface, if cast returns NULL we are not allowed to
        // call resize and dofCompress 
        ManagedDofStorageInterface* managedObject = dynamic_cast< ManagedDofStorageInterface* > ( memObject_ );
        if( managedObject ) 
        {
          // resize function 
          managedObject->resize();
        }
      }

    protected:
      //! wrapper class to create fake DofStorage from double* 
      template <class VectorPointerType>
      class DofStorageWrapper : public DofStorageInterface
      {
        const std::string name_;
        DofStorageType array_;
      public:
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
      template <class VectorPointerType>
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
           = allocateManagedDofStorage( spc_.gridPart().grid(), mapper(),
                                        name, (MutableDofStorageType *) 0);
                                        
        // save pointer 
        memObject_ = memPair.first;
        return *(memPair.second);
      }

      MapperType& mapper() { return mapper_ ; }
   
      const DiscreteFunctionSpaceType& spc_;
      MapperType mapper_;
      DofStorageInterface* memObject_;
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

  } // end namespace Fem 

} // end namespace Dune

#include "adaptiveimp.cc"

#endif // #ifndef DUNE_FEM_ADAPTIVEFUNCTIONIMP_HH
