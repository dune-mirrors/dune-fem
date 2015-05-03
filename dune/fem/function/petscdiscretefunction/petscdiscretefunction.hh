// vim: set expandtab ts=2 sw=2 sts=2:
#ifndef DUNE_FEM_PETSCDISCRETEFUNCTION_HH
#define DUNE_FEM_PETSCDISCRETEFUNCTION_HH

#include <string>
#include <algorithm>
#include <map>
#include <vector>
#include <iterator>
#include <utility>


#if HAVE_PETSC

#include <dune/fem/misc/petsc/petsccommon.hh>
#include <dune/fem/misc/petsc/petscdofmappings.hh>
#include <dune/fem/misc/petsc/petscvector.hh>
#include <dune/fem/misc/petsc/petscslavedofprovider.hh>


#include <dune/common/fvector.hh>
#include <dune/common/dynvector.hh>

#include <dune/fem/common/stackallocator.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/common/functor.hh>
#include <dune/fem/function/localfunction/mutable.hh>

#include <dune/common/shared_ptr.hh>

namespace Dune
{

  namespace Fem
  {


    /* ==============================
     * forward declarations
     */
    template< class DiscreteFunctionSpace > class PetscDiscreteFunction;

    // AssignVectorReference
    // ---------------------

    template< class DofProxy, class Allocator >
    struct AssignVectorReference< Dune::DynamicVector< DofProxy, Allocator >  >
    {
      // we need to overload this
      typedef Dune::DynamicVector< DofProxy, Allocator > Vector;

      AssignVectorReference ( Vector &vector )
        : vector_( vector )
      {}

      void operator() ( const std::size_t index, DofProxy value )
      {
        vector_[ index ].assign( value );
      }

    protected:
      Vector &vector_;
    };

    /** \class DiscreteFunctionTraits
     *  \brief Traits class for a DiscreteFunction
     *
     *  \tparam  DiscreteFunctionSpace   space the discrete function lives in
     *  \tparam  DofVector             implementation class of the block vector
     */
    template< typename DiscreteFunctionSpace >
    struct DiscreteFunctionTraits< DiscreteFunction< DiscreteFunctionSpace, PetscVector< DiscreteFunctionSpace > > >
    {
      typedef PetscVector< DiscreteFunctionSpace >  DofVectorType;

      typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
      typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
      typedef DiscreteFunction< DiscreteFunctionSpace, DofVectorType >   DiscreteFunctionType;

      typedef typename DofVectorType::IteratorType DofIteratorType;
      typedef typename DofVectorType::ConstIteratorType ConstDofIteratorType;
      typedef typename DofVectorType::DofBlockType DofBlockType;
      typedef typename DofVectorType::ConstDofBlockType ConstDofBlockType;
      typedef Fem::Envelope< DofBlockType >                          DofBlockPtrType;
      typedef Fem::Envelope< ConstDofBlockType >                     ConstDofBlockPtrType;
      typedef typename DiscreteFunctionSpaceType::BlockMapperType MapperType;
      typedef typename DofVectorType::FieldType DofType;

      typedef typename DofBlockType::DofProxy DofProxyType;

      typedef ThreadSafeValue< UninitializedObjectStack > LocalDofVectorStackType;
      typedef StackAllocator< DofProxyType, LocalDofVectorStackType* > LocalDofVectorAllocatorType;
      typedef Dune::DynamicVector< DofProxyType, LocalDofVectorAllocatorType > LocalDofVectorType;

      typedef MutableLocalFunction< DiscreteFunctionType > LocalFunctionType;
    };


#if 0
    /* ========================================
     * struct PetscDiscreteFunctionTraits
     * =======================================
     */
    template< class DiscreteFunctionSpace >
    struct DiscreteFunctionTraits< PetscDiscreteFunction< DiscreteFunctionSpace > >
    {
      typedef DiscreteFunctionSpace                                     DiscreteFunctionSpaceType;
      typedef PetscDiscreteFunction< DiscreteFunctionSpaceType >        DiscreteFunctionType;

      typedef typename DiscreteFunctionSpaceType::DomainType            DomainType;
      typedef typename DiscreteFunctionSpaceType::RangeType             RangeType;
      typedef typename DiscreteFunctionSpaceType::RangeFieldType        RangeFieldType;
      typedef typename DiscreteFunctionSpaceType::JacobianRangeType     JacobianRangeType;

      typedef typename DiscreteFunctionSpaceType::BlockMapperType       BlockMapperType;
      typedef typename DiscreteFunctionSpaceType::GridPartType          GridPartType;

      typedef PetscVector< DiscreteFunctionSpaceType >                  PetscVectorType;
      typedef typename PetscVectorType::DofBlockType                    DofBlockType;
      typedef typename PetscVectorType::ConstDofBlockType               ConstDofBlockType;
      typedef typename PetscVectorType::DofIteratorType                 DofIteratorType;
      typedef typename PetscVectorType::ConstDofIteratorType            ConstDofIteratorType;
      typedef typename PetscVectorType::DofBlockPtrType                 DofBlockPtrType;
      typedef typename PetscVectorType::ConstDofBlockPtrType            ConstDofBlockPtrType;

      typedef typename DofBlockType::DofProxy DofProxyType;

      typedef RangeFieldType DofType;

      typedef ThreadSafeValue< UninitializedObjectStack > LocalDofVectorStackType;
      typedef StackAllocator< DofProxyType, LocalDofVectorStackType* > LocalDofVectorAllocatorType;
      typedef Dune::DynamicVector< DofProxyType, LocalDofVectorAllocatorType > LocalDofVectorType;

      typedef MutableLocalFunction< DiscreteFunctionType > LocalFunctionType;
    };


    /* ========================================
     * class PetscDiscreteFunction
     * =======================================
     */
    template< class DiscreteFunctionSpace >
    class PetscDiscreteFunction
      : public DiscreteFunctionDefault< PetscDiscreteFunction< DiscreteFunctionSpace > >
    {
      typedef PetscDiscreteFunction< DiscreteFunctionSpace > ThisType;
      typedef DiscreteFunctionDefault< PetscDiscreteFunction< DiscreteFunctionSpace > > BaseType;

    public:

      /*
       * types
       */
      typedef DiscreteFunctionTraits< ThisType > Traits;
      typedef typename Traits::DiscreteFunctionSpaceType                DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType::FunctionSpaceType     FunctionSpaceType;
      typedef typename Traits::DiscreteFunctionSpaceType::EntityType    EntityType ;

      const static size_t localBlockSize = DiscreteFunctionSpaceType::localBlockSize;

      typedef typename Traits :: PetscVectorType                        PetscVectorType;

      typedef typename DiscreteFunctionSpaceType::RangeFieldType        RangeFieldType;
      typedef typename DiscreteFunctionSpaceType::DomainFieldType       DomainFieldType;
      typedef typename Traits::BlockMapperType                          BlockMapperType;
      typedef typename Traits::GridPartType                             GridPartType;

      typedef typename Traits::DomainType                               DomainType;
      typedef typename Traits::RangeType                                RangeType;
      typedef typename Traits::JacobianRangeType                        JacobianRangeType;

      typedef typename PetscVectorType::DofBlockType                    DofBlockType;
      typedef typename PetscVectorType::ConstDofBlockType               ConstDofBlockType;
      typedef typename PetscVectorType::DofIteratorType                 DofIteratorType;
      typedef typename PetscVectorType::ConstDofIteratorType            ConstDofIteratorType;
      typedef typename PetscVectorType::DofBlockPtrType                 DofBlockPtrType;
      typedef typename PetscVectorType::ConstDofBlockPtrType            ConstDofBlockPtrType;

      typedef typename Traits :: DofType                                DofType;
      typedef typename Traits :: DofProxyType                           DofProxyType;

      typedef typename BaseType :: LocalDofVectorAllocatorType LocalDofVectorAllocatorType;

    protected:
      typedef PetscManagedDofStorage< DiscreteFunctionSpace, BlockMapperType > PetscManagedDofStorageType;
    public:
      using BaseType :: space;
      using BaseType :: name;

      /*
       * methods
       */
      PetscDiscreteFunction ( const std::string &name,
                              const DiscreteFunctionSpaceType &dfSpace )
      : BaseType( name, dfSpace, LocalDofVectorAllocatorType( &ldvStack_ ) ),
        ldvStack_( std::max( sizeof( DofType ), sizeof( DofProxyType ) ) * space().blockMapper().maxNumDofs() * DiscreteFunctionSpaceType::localBlockSize ),
        memObject_( space(), space().blockMapper(), name ),
        petscVector_( memObject_.getArray() )
      {}

      PetscDiscreteFunction ( const ThisType &other )
      : BaseType( "copy of " + other.name(), other.space(), LocalDofVectorAllocatorType( &ldvStack_ ) ),
        ldvStack_( other.ldvStack_ ),
        memObject_( space(), space().blockMapper(), name() ),
        petscVector_( memObject_.getArray() )
      {
        // copy data
        assign( other );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::block( unsigned int index ) const */
      ConstDofBlockPtrType block ( unsigned int index ) const
      {
        return petscVector().block( index );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::block( unsigned int index ) const */
      DofBlockPtrType block ( unsigned int index )
      {
        return petscVector().block( index );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::size() const */
      size_t size () const { return petscVector().size(); }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::dbegin() const */
      DofIteratorType dbegin () { return petscVector().dbegin(); }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::dbegin() const */
      ConstDofIteratorType dbegin () const { return petscVector().dbegin(); }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::dend() const */
      DofIteratorType dend() { return petscVector().dend(); }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::dend() const */
      ConstDofIteratorType dend() const { return petscVector().dend(); }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::clear() */
      void clear ()
      {
        petscVector().clear();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::assign( clear() */
      void assign( const ThisType& other )
      {
        petscVector().assign( other.petscVector() );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::communicate() */
      void communicate ()
      {
        petscVector().communicateNow();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::operator+= */
      const ThisType& operator+= ( const ThisType &other )
      {
        petscVector() += other.petscVector();
        return *this;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::operator-= */
      const ThisType& operator-= ( const ThisType &other )
      {
        petscVector() -= other.petscVector();
        return *this;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::operator*= */
      const ThisType& operator*= ( const PetscScalar scalar )
      {
        petscVector() *= scalar;
        return *this;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::operator/= */
      const ThisType& operator/= ( const PetscScalar scalar )
      {
        petscVector() /= scalar;
        return *this;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::axpy */
      void axpy ( const RangeFieldType &scalar, const ThisType &other )
      {
        petscVector().axpy( static_cast< const PetscScalar& >( scalar ), other.petscVector() );
      }

      /** \brief obtain a constand pointer to the underlying PETSc Vec */
      const Vec* petscVec () const { return petscVector().getVector(); }

      /** \brief obtain a pointer to the underlying PETSc Vec */
      Vec* petscVec () { return petscVector().getVector(); }

      void enableDofCompression ()
      {
        memObject_.enableDofCompression();
      }

      void print( std::ostream& out )
      {
        petscVector().printGlobal( true );
      }

    protected:
      PetscDiscreteFunction ();
      ThisType& operator= ( const ThisType &other );

      PetscVectorType& petscVector () { return petscVector_; }
      const PetscVectorType& petscVector () const { return petscVector_; }

      /*
       * data fields
       */
      typename Traits::LocalDofVectorStackType ldvStack_;
      PetscManagedDofStorageType memObject_;
      PetscVectorType&           petscVector_;
    };
#else

    //! @ingroup AdaptiveDFunction
    //! An adaptive discrete function
    //! This class is comparable to DFAdapt, except that it provides a
    //! specialisation for CombinedSpace objects which provides enriched
    //! functionality (access to subfunctions) and runtime optimisations
    template <class DiscreteFunctionSpace>
    class PetscDiscreteFunction
    : public DiscreteFunction< DiscreteFunctionSpace,
                               PetscVector< DiscreteFunctionSpace > >
    {
      typedef ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpace > ThisType;
      typedef DiscreteFunction< DiscreteFunctionSpace,
                                PetscVector< DiscreteFunctionSpace > > BaseType;

    public:
      typedef typename BaseType :: DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;
      typedef typename BaseType :: DofVectorType              DofVectorType;
      typedef typename BaseType :: DofType                    DofType;

      using BaseType::assign;

      PetscDiscreteFunction( const std::string &name,
                             const DiscreteFunctionSpaceType &space )
        : BaseType( name, space, allocateDofStorage( space ) )
      {
      }

      PetscDiscreteFunction( const std::string &name,
                             const DiscreteFunctionSpaceType &space,
                             DofVectorType& dofVector )
        : BaseType( name, space, dofVector )
      {
        // in this case we have no allocated mem object
        memObject_ = 0;
      }

      PetscDiscreteFunction( const PetscDiscreteFunction& other )
        : BaseType( "copy of " + other.name(), other.space(), allocateDofStorage( other.space() ) )
      {
        assign( other );
      }

      ~PetscDiscreteFunction()
      {
        if( memObject_ )
        {
          delete memObject_;
          memObject_ = 0;
        }
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::enableDofCompression()
       */
      void enableDofCompression ()
      {
        if( memObject_ )
          memObject_->enableDofCompression();
      }

    protected:
      typedef typename DiscreteFunctionSpaceType :: BlockMapperType  BlockMapperType;
      typedef PetscManagedDofStorage< DiscreteFunctionSpaceType, BlockMapperType > PetscManagedDofStorageType;

      // allocate managed dof storage
      DofVectorType& allocateDofStorage ( const DiscreteFunctionSpaceType &space )
      {
        std::string name("deprecated");

        memObject_ = new PetscManagedDofStorageType( space, space.blockMapper(), name ) ;

        return memObject_->getArray();
      }

      // pointer to allocated DofVector
      PetscManagedDofStorageType* memObject_;
    };

#endif


  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_PETSC

#endif // #ifndef DUNE_FEM_PETSCDISCRETEFUNCTION_HH
