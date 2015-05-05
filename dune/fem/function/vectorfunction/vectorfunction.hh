#ifndef DUNE_FEM_VECTORFUNCTION_HH
#define DUNE_FEM_VECTORFUNCTION_HH

#include <dune/common/typetraits.hh>

#include <dune/fem/common/referencevector.hh>
#include <dune/fem/common/stackallocator.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/common/dofblock.hh>
#include <dune/fem/function/localfunction/mutable.hh>
#include <dune/fem/storage/envelope.hh>
#include <dune/fem/storage/vector.hh>

#include <dune/fem/function/blockvectordiscretefunction/discretefunction.hh>
#include <dune/fem/function/blockvectors/simpleblockvector.hh>

namespace Dune
{

  namespace Fem
  {


#if 0
    // Internal Forward Declarations
    // -----------------------------

    template< class DiscreteFunctionSpace, class DofVector >
    class VectorDiscreteFunction;


    // VectorDiscreteFunctionTraits
    // ----------------------------

    template< class DiscreteFunctionSpace, class DofVector >
    struct DiscreteFunctionTraits< VectorDiscreteFunction< DiscreteFunctionSpace, DofVector > >
    {
      typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

      typedef DofVector DofVectorType;

      typedef VectorDiscreteFunction< DiscreteFunctionSpaceType, DofVectorType >  DiscreteFunctionType;

      typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
      typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;

      typedef typename DiscreteFunctionSpaceType :: DomainFieldType
        DomainFieldType;
      typedef typename DiscreteFunctionSpaceType :: RangeFieldType
        RangeFieldType;

      typedef typename DiscreteFunctionSpaceType :: JacobianRangeType
        JacobianRangeType;

      typedef typename DofVectorType :: value_type  DofType;
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

      typedef ThreadSafeValue< UninitializedObjectStack > LocalDofVectorStackType;
      typedef StackAllocator< DofType, LocalDofVectorStackType* > LocalDofVectorAllocatorType;
      typedef DynamicReferenceVector< DofType, LocalDofVectorAllocatorType > LocalDofVectorType;

      typedef MutableLocalFunction< DiscreteFunctionType > LocalFunctionType;
    };


    // VectorDiscreteFunction
    // ----------------------

    template< class DiscreteFunctionSpace, class DofVector >
    class VectorDiscreteFunction
    : public DiscreteFunctionDefault< VectorDiscreteFunction< DiscreteFunctionSpace, DofVector > >
    {
      typedef VectorDiscreteFunction< DiscreteFunctionSpace, DofVector > ThisType;
      typedef DiscreteFunctionDefault< VectorDiscreteFunction< DiscreteFunctionSpace, DofVector > > BaseType;

      static_assert( SupportsVectorInterface< DofVector >::v, "DofVector must support VectorInterface." );

    public:
      //! type of this class's traits
      typedef DiscreteFunctionTraits< ThisType > Traits;

      //! type of the associated discrete function space
      typedef typename BaseType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

      typedef typename Traits::DiscreteFunctionType DiscreteFunctionType;


      typedef typename BaseType::LocalFunctionType LocalFunctionType;


      typedef typename BaseType::DomainType DomainType;
      typedef typename BaseType::RangeType RangeType;

      typedef typename BaseType::DomainFieldType DomainFieldType;
      typedef typename BaseType::RangeFieldType RangeFieldType;

      typedef typename BaseType::JacobianRangeType JacobianRangeType;

      typedef typename BaseType::DofType DofType;

      //! type of the DoF storage array
      typedef typename Traits::DofVectorType DofVectorType;
      typedef typename Traits::DofStorageType DofStorageType;

      typedef typename BaseType::ConstDofIteratorType ConstDofIteratorType;
      typedef typename BaseType::DofIteratorType DofIteratorType;

      typedef typename BaseType::DofBlockType DofBlockType;
      typedef typename BaseType::ConstDofBlockType ConstDofBlockType;
      typedef typename BaseType::DofBlockPtrType DofBlockPtrType;
      typedef typename BaseType::ConstDofBlockPtrType ConstDofBlockPtrType;

      typedef typename BaseType :: LocalDofVectorAllocatorType LocalDofVectorAllocatorType;

      //! Constructor
      VectorDiscreteFunction ( const std::string &name,
                               const DiscreteFunctionSpaceType &dfSpace,
                               DofVectorType &dofVector )
      : BaseType( name, dfSpace, LocalDofVectorAllocatorType( &ldvStack_ ) ),
        ldvStack_( std::max( sizeof( DofType ), sizeof( DofType* ) ) * space().blockMapper().maxNumDofs() * DiscreteFunctionSpaceType::localBlockSize ),
        dofVector_( &dofVector ),
        freeDofVector_( false )
      {
        // size of dof vector must be size of space in blocks x localBlockSize
        assert( dofVector_->size() == (unsigned int)dfSpace.blockMapper().size() *
            DiscreteFunctionSpaceType :: localBlockSize );
      }

      VectorDiscreteFunction ( const ThisType &other )
      : BaseType( other.name(), other.space(), LocalDofVectorAllocatorType( &ldvStack_ ) ),
        ldvStack_( std::max( sizeof( DofType ), sizeof( DofType* ) ) * space().blockMapper().maxNumDofs() * DiscreteFunctionSpaceType::localBlockSize ),
        dofVector_( new DofVectorType( other.dofVector() ) ),
        freeDofVector_( true )
      {}

      ~VectorDiscreteFunction ()
      {
        if( freeDofVector_ )
          delete dofVector_;
      }

      using BaseType::space;

    private:
      // prohibit assignment
      ThisType &operator= ( const ThisType & );

    public:
      using BaseType :: assign ;
      using BaseType :: operator+=;
      using BaseType :: operator-=;

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

      void axpy ( const RangeFieldType &s, const ThisType &u )
      {
        dofVector().addScaled( s, u.dofVector() );
      }

      void axpy ( const ThisType &u, const RangeFieldType &s )
      {
        axpy( s, u );
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

      const DofType &dof ( unsigned int index ) const
      {
        return dofVector()[ index ];
      }

      DofType &dof ( unsigned int index )
      {
        return dofVector()[ index ];
      }

      const DofType *leakPointer () const
      {
        return dofVector().leakPointer();
      }

      DofType *leakPointer ()
      {
        return dofVector().leakPointer();
      }

      RangeFieldType scalarProductDofs ( const ThisType &other ) const
      {
        return this->scalarProduct_.scalarProductDofs( *this, other );
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
      typename Traits :: LocalDofVectorStackType ldvStack_;

      DofVectorType *const dofVector_;
      const bool freeDofVector_;
    };
#else

    //! @ingroup AdaptiveDFunction
    //! An adaptive discrete function
    //! This class is comparable to DFAdapt, except that it provides a
    //! specialisation for CombinedSpace objects which provides enriched
    //! functionality (access to subfunctions) and runtime optimisations
    template < class DiscreteFunctionSpace, class Vector >
    class VectorDiscreteFunction
    : public DiscreteFunction< DiscreteFunctionSpace,
                               SimpleBlockVector< Vector, DiscreteFunctionSpace::localBlockSize > >
    {
      typedef VectorDiscreteFunction< DiscreteFunctionSpace, Vector > ThisType;
      typedef DiscreteFunction< DiscreteFunctionSpace,
                                SimpleBlockVector< Vector, DiscreteFunctionSpace::localBlockSize > > BaseType;

    public:
      typedef Vector VectorType;
      typedef typename BaseType :: DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;
      typedef typename BaseType :: DofVectorType              DofVectorType;
      typedef typename BaseType :: DofType                    DofType;

      using BaseType::assign;

      VectorDiscreteFunction( const std::string &name,
                              const DiscreteFunctionSpaceType &space,
                              VectorType& vector )
        : BaseType( name, space, dofVector_ ),
          vec_(),
          dofVector_( vector )
      {
      }

      VectorDiscreteFunction( const VectorDiscreteFunction& other )
        : BaseType( "copy of " + other.name(), other.space(), dofVector_ ),
          vec_(),
          dofVector_( allocateDofVector( other.space() ) )
      {
        assign( other );
      }

    protected:
      // allocate managed dof storage
      VectorType& allocateDofVector ( const DiscreteFunctionSpaceType& space )
      {
        vec_.reset( new VectorType( space.size() ) );
        return *vec_;
      }

      // pointer to DofContainer
      std::unique_ptr< VectorType > vec_;
      // pointer to dof vector
      DofVectorType dofVector_;
    };
#endif

    // Capabilibies
    // ------------

    namespace Capabilities
    {

      template< class DiscreteFunctionSpace, class DofVector >
      struct HasLeakPointer
        < Fem :: VectorDiscreteFunction< DiscreteFunctionSpace, DofVector > >
      : public HasLeakPointer< DofVector >
      {};

    }

  } // namespace Fem

} // namespace Dune

#include "managedvectorfunction.hh"

#endif // #ifndef DUNE_FEM_VECTORFUNCTION_HH
