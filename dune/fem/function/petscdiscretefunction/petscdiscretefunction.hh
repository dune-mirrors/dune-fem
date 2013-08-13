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
#include <dune/fem/misc/petsc/petsclocalfunction.hh>

#include <dune/fem/misc/petsc/petscslavedofprovider.hh>


#include <dune/common/fvector.hh> 

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/localfunction/localfunction.hh>
#include <dune/fem/function/blockvectordiscretefunction/blockvectordiscretefunction.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/common/shared_ptr.hh>

namespace Dune 
{

  namespace Fem 
  {


    /* ==============================
     * forward declarations
     */
    template< class DiscreteFunctionSpace > class PetscDiscreteFunction;


    /* ========================================
     * struct PetscDiscreteFunctionTraits
     * =======================================
     */
    template< class DiscreteFunctionSpace >
    struct PetscDiscreteFunctionTraits
    {
      typedef DiscreteFunctionSpace                                     DiscreteFunctionSpaceType;
      typedef PetscDiscreteFunctionTraits< DiscreteFunctionSpaceType >  Traits;
      typedef typename DiscreteFunctionSpaceType::DomainType            DomainType;
      typedef typename DiscreteFunctionSpaceType::RangeType             RangeType;

      typedef typename DiscreteFunctionSpaceType::JacobianRangeType     JacobianRangeType;
      typedef typename DiscreteFunctionSpaceType::MapperType            MapperType;
      typedef typename DiscreteFunctionSpaceType::BlockMapperType       BlockMapperType;
      typedef typename DiscreteFunctionSpaceType::GridPartType          GridPartType;
      typedef PetscDiscreteFunction< DiscreteFunctionSpaceType >        DiscreteFunctionType;

      typedef PetscLocalFunctionFactory< Traits >                       LocalFunctionFactoryType;
      typedef LocalFunctionStack< LocalFunctionFactoryType >            LocalFunctionStorageType;
      typedef PetscLocalFunction< DiscreteFunctionType >                LocalFunctionType;

      typedef PetscVector< DiscreteFunctionSpaceType  >                 PetscVectorType;
      typedef typename PetscVectorType::DofBlockType                    DofBlockType;
      typedef typename PetscVectorType::ConstDofBlockType               ConstDofBlockType;
      typedef typename PetscVectorType::DofIteratorType                 DofIteratorType;
      typedef typename PetscVectorType::ConstDofIteratorType            ConstDofIteratorType;
      typedef typename PetscVectorType::DofBlockPtrType                 DofBlockPtrType;
      typedef typename PetscVectorType::ConstDofBlockPtrType            ConstDofBlockPtrType;
    };

    /* ========================================
     * class PetscDiscreteFunction
     * =======================================
     */
    template< class DiscreteFunctionSpace >
    class PetscDiscreteFunction
      : public DiscreteFunctionDefault< PetscDiscreteFunctionTraits< DiscreteFunctionSpace > >
    {
      typedef PetscDiscreteFunction< DiscreteFunctionSpace > ThisType;
      typedef DiscreteFunctionDefault< PetscDiscreteFunctionTraits< DiscreteFunctionSpace > > BaseType;
      friend class PetscLocalFunction< ThisType >;

    public:

      /*
       * types
       */
      typedef PetscDiscreteFunctionTraits< DiscreteFunctionSpace > Traits;
      typedef typename Traits::DiscreteFunctionSpaceType                DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType::FunctionSpaceType     FunctionSpaceType;
      typedef typename Traits::DiscreteFunctionSpaceType::EntityType    EntityType ;

      const static size_t localBlockSize = DiscreteFunctionSpaceType::localBlockSize;

      typedef typename Traits :: PetscVectorType                        PetscVectorType;

      typedef typename DiscreteFunctionSpaceType::RangeFieldType        RangeFieldType;
      typedef typename DiscreteFunctionSpaceType::DomainFieldType       DomainFieldType;
      typedef typename Traits::BlockMapperType                          BlockMapperType;
      typedef typename Traits::GridPartType                             GridPartType;

      typedef typename Traits :: LocalFunctionFactoryType               LocalFunctionFactoryType;

      typedef typename Traits::DomainType                               DomainType;
      typedef typename Traits::RangeType                                RangeType;
      typedef typename Traits::JacobianRangeType                        JacobianRangeType;

      typedef typename PetscVectorType::DofBlockType                    DofBlockType;
      typedef typename PetscVectorType::ConstDofBlockType               ConstDofBlockType;
      typedef typename PetscVectorType::DofIteratorType                 DofIteratorType;
      typedef typename PetscVectorType::ConstDofIteratorType            ConstDofIteratorType;
      typedef typename PetscVectorType::DofBlockPtrType                 DofBlockPtrType;
      typedef typename PetscVectorType::ConstDofBlockPtrType            ConstDofBlockPtrType;

      // what is DoFType ?!?!?!? this
      typedef RangeFieldType DofType;

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
      : BaseType( name, dfSpace, lfFactory_ ),
        lfFactory_( *this ),
        memObject_( space(), space().blockMapper(), name ),
        petscVector_( memObject_.getArray() )
      {}

      PetscDiscreteFunction ( const ThisType &other )
      : BaseType( "copy of " + other.name(), other.space(), lfFactory_ ),
        lfFactory_( *this ),
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
      LocalFunctionFactoryType   lfFactory_;
      PetscManagedDofStorageType memObject_;
      PetscVectorType&           petscVector_;
    };


  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_PETSC

#endif // #ifndef DUNE_FEM_PETSCDISCRETEFUNCTION_HH
