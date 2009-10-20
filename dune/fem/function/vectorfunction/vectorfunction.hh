#ifndef DUNE_FEM_VECTORFUNCTION_HH
#define DUNE_FEM_VECTORFUNCTION_HH

#include <dune/common/static_assert.hh>
#include <dune/common/typetraits.hh>

#include <dune/fem/storage/vector.hh>
#include <dune/fem/storage/envelope.hh>
#include <dune/fem/function/common/dofblock.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/localfunction/standardlocalfunction.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< class DiscreteFunctionSpace, class DofVector >
  class VectorDiscreteFunction;



  // VectorDiscreteFunctionTraits
  // ----------------------------

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

    typedef typename DofVectorType :: FieldType DofType;
    typedef DofVectorType DofStorageType;

    typedef typename DofVectorType :: ConstIteratorType ConstDofIteratorType;
    typedef typename DofVectorType :: IteratorType DofIteratorType;

    enum { blockSize = DiscreteFunctionSpaceType :: localBlockSize };

    typedef DofBlockProxy< DiscreteFunctionType, DofType, blockSize >
      DofBlockType;
    typedef DofBlockProxy
      < const DiscreteFunctionType, const DofType, blockSize >
      ConstDofBlockType;
    typedef Envelope< DofBlockType > DofBlockPtrType;
    typedef Envelope< ConstDofBlockType > ConstDofBlockPtrType;
  };



  // VectorDiscreteFunction
  // ----------------------

  template< class DiscreteFunctionSpace, class DofVector >
  class VectorDiscreteFunction
  : public DiscreteFunctionDefault< VectorDiscreteFunctionTraits< DiscreteFunctionSpace, DofVector > >
  {
    typedef VectorDiscreteFunction< DiscreteFunctionSpace, DofVector > ThisType;
    typedef DiscreteFunctionDefault< VectorDiscreteFunctionTraits< DiscreteFunctionSpace, DofVector > > BaseType;

    dune_static_assert( SupportsVectorInterface< DofVector >::v, "DofVector must support VectorInterface." );

  public:
    //! type of this class's traits
    typedef VectorDiscreteFunctionTraits< DiscreteFunctionSpace, DofVector > Traits;

    //! type of the associated discrete function space
    typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    //! type of the DoF storage array
    typedef typename Traits::DofVectorType DofVectorType;

    typedef typename Traits::DiscreteFunctionType DiscreteFunctionType;

    typedef typename Traits::LocalFunctionType LocalFunctionType;
    typedef typename Traits::LocalFunctionFactoryType  LocalFunctionFactoryType;
    typedef typename Traits::LocalFunctionStorageType LocalFunctionStorageType;

    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;

    typedef typename Traits::DomainFieldType DomainFieldType;
    typedef typename Traits::RangeFieldType RangeFieldType;

    typedef typename Traits::JacobianRangeType JacobianRangeType;

    typedef typename Traits::DofType DofType;
    typedef typename Traits::DofStorageType DofStorageType;

    typedef typename Traits::ConstDofIteratorType ConstDofIteratorType;
    typedef typename Traits::DofIteratorType DofIteratorType;

    typedef typename Traits::DofBlockType DofBlockType;
    typedef typename Traits::ConstDofBlockType ConstDofBlockType;
    typedef typename Traits::DofBlockPtrType DofBlockPtrType;
    typedef typename Traits::ConstDofBlockPtrType ConstDofBlockPtrType;

  private:
    dune_static_assert( (Conversion< RangeFieldType, DofType >::sameType), "RangeFieldType and DofType must equal." );

  public:
    //! Constructor
    VectorDiscreteFunction ( const std::string &name,
                             const DiscreteFunctionSpaceType &dfSpace,
                             DofVectorType &dofVector )
    : BaseType( name, dfSpace, lfFactory_ ),
      lfFactory_( *this ),
      dofVector_( &dofVector ),
      freeDofVector_( false )
    {
      assert( dofVector_->size() == (unsigned int)dfSpace.size() );
    }

    VectorDiscreteFunction ( const ThisType &other )
    : BaseType( other.name(), other.space(), lfFactory_ ),
      lfFactory_( *this ),
      dofVector_( new DofVectorType( other.dofVector() ) ),
      freeDofVector_( true )
    {}

    ~VectorDiscreteFunction ()
    {
      if( freeDofVector_ )
        delete dofVector_;
    }

  private:
    // prohibit assignment
    ThisType &operator= ( const ThisType & );

  public:
    ThisType &operator+= ( const ThisType &u )
    {
      dofVector() += u.dofVector();
      return *this;
    }

    ThisType &operator-= ( const ThisType &u )
    {
      dofVector() -= u.dofVector();
      return *this;
    }

    void addScaled ( const RangeFieldType &s, const ThisType &u )
    {
      dofVector().addScaled( s, u.dofVector() );
    }

    void addScaled ( const ThisType &u, const RangeFieldType &s )
    {
      addScaled( s, u );
    }

    void assign ( const ThisType &u )
    {
      dofVector().assign( u.dofVector() );
    }

    void clear ()
    {
      dofVector().assign( 0 );
    }

    ConstDofIteratorType dbegin () const
    {
      return dofVector().begin();
    }

    DofIteratorType dbegin ()
    {
      return dofVector().begin();
    }

    ConstDofIteratorType dend () const
    {
      return dofVector().end();
    }

    DofIteratorType dend ()
    {
      return dofVector().end();
    }

    ConstDofBlockPtrType block ( unsigned int index ) const
    {
      typename ConstDofBlockType::KeyType key( this, index );
      return ConstDofBlockPtrType( key );
    }

    DofBlockPtrType block ( unsigned int index )
    {
      typename DofBlockType::KeyType key( this, index );
      return DofBlockPtrType( key );
    }

    const RangeFieldType &dof ( unsigned int index ) const
    {
      return dofVector()[ index ];
    }

    RangeFieldType &dof ( unsigned int index )
    {
      return dofVector()[ index ];
    }

    const RangeFieldType *leakPointer () const
    {
      return dofVector().leakPointer();
    }

    RangeFieldType *leakPointer ()
    {
      return dofVector().leakPointer();
    }

    RangeFieldType scalarProductDofs ( const ThisType &u ) const
    {
      return dofVector() * u.dofVector();
    }

    int size () const
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

  private:
    const LocalFunctionFactoryType lfFactory_;

    DofVectorType *const dofVector_;
    const bool freeDofVector_;
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
