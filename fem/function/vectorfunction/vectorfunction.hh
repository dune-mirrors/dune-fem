#ifndef DUNE_FEM_VECTORFUNCTION_HH
#define DUNE_FEM_VECTORFUNCTION_HH

#include <dune/common/typetraits.hh>

#include <dune/fem/storage/vector.hh>
#include <dune/fem/storage/envelope.hh>
#include <dune/fem/function/common/dofblock.hh>
//#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/localfunction/standardlocalfunction.hh>

namespace Dune
{

#if 0
  namespace VectorDiscreteFunctionHelper
  {
    template< class DofVector, class Dof, unsigned int Size >
    class DofBlockProxy
    {
      friend class Envelope< DofBlockProxy >;

    public:
      typedef DofVector DofVectorType;

      typedef Dof DofType;
      
      enum { size = Size };
      
      typedef unsigned int size_type;

      typedef std :: pair< DofVectorType*, size_type > KeyType;

    protected:
      DofVectorType &dofVector_;
      const size_type first_;

    protected:
      inline DofBlockProxy ( const KeyType key )
      : dofVector_( *(key.first) ),
        first_( size * key.second )
      {}

      inline DofBlockProxy ( const DofBlockProxy &other )
      : dofVector_( other.dofVector_ ),
        first_( other.first_ )
      {}

    public:
      inline DofBlockProxy &operator= ( const DofBlockProxy &other )
      {
        for( size_type i = 0; i < size; ++i )
          (*this)[ i ] = other[ i ];
        return *this;
      }
      
      inline const DofType &operator[] ( size_type index ) const
      {
        return dofVector_[ first_ + index ];
      }
      
      inline DofType &operator[] ( size_type index )
      {
        return dofVector_[ first_ + index ];
      }
     
      inline size_type dim () const
      {
        return size;
      }
    };
  }
#endif



  template< class DiscreteFunctionSpace, class DofVector >
  class VectorDiscreteFunction;



  template< class DiscreteFunctionSpace, class DofVector >
  struct VectorDiscreteFunctionTraits
  {
    typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

    typedef DofVector DofVectorType;

    typedef VectorDiscreteFunctionTraits
      < DiscreteFunctionSpaceType, DofVectorType >
      DiscreteFunctionTraits;
    typedef VectorDiscreteFunction
      < DiscreteFunctionSpaceType, DofVectorType >
      DiscreteFunctionType;

    typedef StandardLocalFunctionFactory< DiscreteFunctionTraits >
      LocalFunctionFactoryType;
    typedef LocalFunctionStack< LocalFunctionFactoryType >
      LocalFunctionStorageType;
    typedef typename LocalFunctionStorageType :: LocalFunctionType
      LocalFunctionType;

    typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;

    typedef typename DiscreteFunctionSpaceType :: DomainFieldType
      DomainFieldType;
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType
      RangeFieldType;

    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType
      JacobianRangeType;

    typedef typename DiscreteFunctionSpaceType :: MapperType MapperType;
    //typedef typename DiscreteFunctionSpaceType :: GridType GridType;

    typedef typename DofVectorType :: FieldType DofType;
    typedef DofVectorType DofStorageType;

    typedef typename DofVectorType :: ConstIteratorType ConstDofIteratorType;
    typedef typename DofVectorType :: IteratorType DofIteratorType;

    enum { blockSize = DiscreteFunctionSpaceType :: localBlockSize };

#if 0
    typedef VectorDiscreteFunctionHelper :: DofBlockProxy
      < DofVectorType, DofType, blockSize >
      DofBlockType;
    typedef VectorDiscreteFunctionHelper :: DofBlockProxy
      < DofVectorType, const DofType, blockSize >
      ConstDofBlockType;
#endif
    typedef DofBlockProxy< DiscreteFunctionType, DofType, blockSize >
      DofBlockType;
    typedef DofBlockProxy
      < const DiscreteFunctionType, const DofType, blockSize >
      ConstDofBlockType;
    typedef Envelope< DofBlockType > DofBlockPtrType;
    typedef Envelope< ConstDofBlockType > ConstDofBlockPtrType;
  };



  template< class DiscreteFunctionSpace, class DofVector >
  class VectorDiscreteFunction
  : public DiscreteFunctionDefault
    < VectorDiscreteFunctionTraits< DiscreteFunctionSpace, DofVector > >
  {
  public:
    //! type of this class's traits
    typedef VectorDiscreteFunctionTraits< DiscreteFunctionSpace, DofVector >
      Traits;

    //! type of the associated discrete function space
    typedef typename Traits :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;
    //! type of the DoF storage array
    typedef typename Traits :: DofVectorType DofVectorType;

  private:
    typedef VectorDiscreteFunction
      < DiscreteFunctionSpaceType, DofVectorType > ThisType;
    typedef DiscreteFunctionDefault< Traits > BaseType;

  public:
    typedef typename Traits :: DiscreteFunctionType DiscreteFunctionType;

    typedef typename Traits :: LocalFunctionType LocalFunctionType;
    typedef typename Traits :: LocalFunctionFactoryType
      LocalFunctionFactoryType;
    typedef typename Traits :: LocalFunctionStorageType
      LocalFunctionStorageType;

    typedef typename Traits :: DomainType DomainType;
    typedef typename Traits :: RangeType RangeType;

    typedef typename Traits :: DomainFieldType DomainFieldType;
    typedef typename Traits :: RangeFieldType RangeFieldType;

    typedef typename Traits :: JacobianRangeType JacobianRangeType;

    typedef typename Traits :: DofType DofType;
    typedef typename Traits :: DofStorageType DofStorageType;

    typedef typename Traits :: ConstDofIteratorType ConstDofIteratorType;
    typedef typename Traits :: DofIteratorType DofIteratorType;

    typedef typename Traits :: DofBlockType DofBlockType;
    typedef typename Traits :: ConstDofBlockType ConstDofBlockType;
    typedef typename Traits :: DofBlockPtrType DofBlockPtrType;
    typedef typename Traits :: ConstDofBlockPtrType ConstDofBlockPtrType;

  private:
    typedef CheckVectorInterface< DofVectorType > CheckDofVectorType;

    typedef CompileTimeChecker
      < Conversion< RangeFieldType, DofType > :: sameType >
      Check_RangeFieldType_and_DofType_are_Equal;

  private:
    const LocalFunctionFactoryType lfFactory_;

    DofVectorType *const dofVector_;
    const bool freeDofVector_;

  public:
    //! Constructor
    inline VectorDiscreteFunction ( const std :: string &name,
                                    const DiscreteFunctionSpaceType &dfSpace,
                                    DofVectorType &dofVector )
    : BaseType( name, dfSpace, lfFactory_ ),
      lfFactory_( *this ),
      dofVector_( &dofVector ),
      freeDofVector_( false )
    {
      assert( dofVector_->size() == (unsigned int)dfSpace.size() );
    }

    inline VectorDiscreteFunction ( const ThisType &other )
    : BaseType( other.name(), other.space(), lfFactory_ ),
      lfFactory_( *this ),
      dofVector_( new DofVectorType( other.dofVector() ) ),
      freeDofVector_( true )
    {}

    inline ~VectorDiscreteFunction ()
    {
      if( freeDofVector_ )
        delete dofVector_;
    }

  private:
    // prohibit assignment
    ThisType &operator= ( const ThisType & );

  public:
    inline ThisType &operator+= ( const ThisType &u )
    {
      dofVector() += u.dofVector();
      return *this;
    }

    inline ThisType &operator-= ( const ThisType &u )
    {
      dofVector() -= u.dofVector();
      return *this;
    }

    inline void addScaled ( const RangeFieldType &s,
                            const ThisType &u )
    {
      dofVector().addScaled( s, u.dofVector() );
    }

    inline void addScaled ( const ThisType &u,
                            const RangeFieldType &s )
    {
      addScaled( s, u );
    }

    inline void assign ( const ThisType &u )
    {
      dofVector().assign( u.dofVector() );
    }

    inline void clear ()
    {
      dofVector().assign( 0 );
    }

    inline ConstDofIteratorType dbegin () const
    {
      return dofVector().begin();
    }

    inline DofIteratorType dbegin ()
    {
      return dofVector().begin();
    }

    inline ConstDofIteratorType dend () const
    {
      return dofVector().end();
    }

    inline DofIteratorType dend ()
    {
      return dofVector().end();
    }

    inline ConstDofBlockPtrType block ( unsigned int index ) const
    {
#if 0
      typename ConstDofBlockType :: KeyType key( &(dofVector()), index );
      return ConstDofBlockPtrType( key );
#endif
      typename ConstDofBlockType :: KeyType key( this, index );
      return ConstDofBlockPtrType( key );
    }

    inline DofBlockPtrType block ( unsigned int index )
    {
#if 0
      typename DofBlockType :: KeyType key( &(dofVector()), index );
      return DofBlockPtrType( key );
#endif
      typename DofBlockType :: KeyType key( this, index );
      return DofBlockPtrType( key );
    }

    inline const RangeFieldType &dof ( unsigned int index ) const
    {
      return dofVector()[ index ];
    }

    inline RangeFieldType &dof ( unsigned int index )
    {
      return dofVector()[ index ];
    }

    inline const RangeFieldType *leakPointer () const
    {
      return dofVector().leakPointer();
    }

    inline RangeFieldType *leakPointer ()
    {
      return dofVector().leakPointer();
    }

    inline RangeFieldType scalarProductDofs ( const ThisType &u ) const
    {
      return dofVector() * u.dofVector();
    }

    inline int size () const
    {
      return dofVector().size();
    }

    const DofVectorType &dofVector () const
    {
      return *dofVector_;
    }
    
    DofVectorType &dofVector ()
    {
      return *dofVector_;
    }
  };



  // Capabilibies
  // ------------

  namespace Capabilities
  {

    template< class DiscreteFunctionSpace, class DofVector >
    struct HasLeakPointer
      < VectorDiscreteFunction< DiscreteFunctionSpace, DofVector > >
    : public HasLeakPointer< DofVector >
    {};
    
  }

};

#include "managedvectorfunction.hh"

#endif
