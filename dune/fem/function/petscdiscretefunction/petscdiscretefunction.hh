// vim: set expandtab ts=2 sw=2 sts=2:
#ifndef DUNE_FEM_PETSCDISCRETEFUNCTION_HH
#define DUNE_FEM_PETSCDISCRETEFUNCTION_HH

#include <string>
#include <algorithm>
#include <map>
#include <vector>
#include <iterator>
#include <utility>

#if defined HAVE_PETSC

#include <dune/fem/misc/petsc/petsccommon.hh>
#include <dune/fem/misc/petsc/petscdofmappings.hh>
#include <dune/fem/misc/petsc/petscvector.hh>
#include <dune/fem/misc/petsc/petsclocalfunction.hh>

#include <dune/fem/misc/petsc/petscslavedofprovider.hh>


#include <dune/common/fvector.hh> 

// If we don't define this, there is an error when headerchecking this file...
#define WANT_CACHED_COMM_MANAGER 0

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/localfunction/localfunction.hh>
//#include <dune/fem/function/localfunction/standardlocalfunction.hh>
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
    template< typename S > class PetscDiscreteFunction;


    /* ========================================
     * struct PetscDiscreteFunctionTraits
     * =======================================
     */
    template< typename S >
    struct PetscDiscreteFunctionTraits
    {
      typedef PetscDiscreteFunctionTraits< S >                          ThisType;
      typedef S                                                         DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType::DomainType            DomainType;
      typedef typename DiscreteFunctionSpaceType::RangeType             RangeType;
      typedef typename DiscreteFunctionSpaceType::MapperType            MapperType;
      typedef typename DiscreteFunctionSpaceType::BlockMapperType       BlockMapperType;
      typedef typename DiscreteFunctionSpaceType::GridPartType          GridPartType;
      typedef PetscDiscreteFunction< S >                                DiscreteFunctionType;
    };

    /* ========================================
     * class PetscDiscreteFunction
     * =======================================
     */
    template< typename S >
    class PetscDiscreteFunction
    : public HasLocalFunction,
      public IsBlockVectorDiscreteFunction
    {
      typedef PetscDiscreteFunction< S > ThisType;
      friend class PetscLocalFunction< ThisType >;

    public:

      /*
       * types
       */
      typedef PetscDiscreteFunctionTraits< S > Traits;
      typedef typename Traits::DiscreteFunctionSpaceType                DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType::FunctionSpaceType     FunctionSpaceType;

      const static size_t localBlockSize = DiscreteFunctionSpaceType::localBlockSize;

      typedef PetscVector< DiscreteFunctionSpaceType  >                 PetscVectorType;

      typedef typename DiscreteFunctionSpaceType::RangeFieldType        RangeFieldType;
      typedef typename DiscreteFunctionSpaceType::DomainFieldType       DomainFieldType;
      typedef typename Traits::BlockMapperType                          BlockMapperType;
      typedef typename Traits::GridPartType                             GridPartType;
      typedef PetscLocalFunctionFactory< ThisType >                     LocalFunctionFactoryType;
      typedef PetscLocalFunctionStack< ThisType >                       LocalFunctionStorageType;

      typedef PetscLocalFunction< ThisType >                            LocalFunctionType;

      typedef typename Traits::DomainType                               DomainType;
      typedef typename Traits::RangeType                                RangeType;

      typedef typename PetscVectorType::DofBlockType                    DofBlockType;
      typedef typename PetscVectorType::ConstDofBlockType               ConstDofBlockType;
      typedef typename PetscVectorType::DofIteratorType                 DofIteratorType;
      typedef typename PetscVectorType::ConstDofIteratorType            ConstDofIteratorType;
      typedef typename PetscVectorType::DofBlockPtrType                 DofBlockPtrType;
      typedef typename PetscVectorType::ConstDofBlockPtrType            ConstDofBlockPtrType;

      typedef typename DofBlockType::DofProxy DofType;

    public:

      /*
       * methods
       */
      PetscDiscreteFunction ( const std::string &name,
                              const DiscreteFunctionSpaceType &dfSpace )
      : name_( name ),
        dfSpace_( dfSpace ),
        lfFactory_( *this ),
        lfStorage_( lfFactory_ ),
        petscVector_( dfSpace )
      {}

      PetscDiscreteFunction ( const ThisType &other )
      : name_( "copy of " + other.name_ ),
        dfSpace_( other.dfSpace_ ),
        lfFactory_( other.lfFactory_ ),
        lfStorage_( other.lfStorage_ ),
        petscVector_( other.petscVector_ )
      {}

      ~PetscDiscreteFunction () 
      {}

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

        /** \copydoc Dune::Fem::DiscreteFunctionInterface::communicate() */
      void communicate ()
      {
        petscVector().communicateNow();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::scalarProductDofs */
      PetscScalar scalarProductDofs ( const ThisType &other ) const
      {
        return petscVector()*other.petscVector();
      }

      /** \copydoc Dune::Fem::Function::evaluate(const FieldVector<int,diffOrder> &diffVariable,const DomainType &x,RangeType &ret) const */
      template< int diffOrder >
      void evaluate ( const ::Dune::FieldVector< int, diffOrder > &diffVariable,
                      const DomainType &x,
                      RangeType &ret ) const
      {
        typedef typename DiscreteFunctionSpaceType::IteratorType Iterator;
        typedef typename Iterator::Entity Entity;
        typedef typename Entity::Geometry Geometry;
        typedef typename Geometry::LocalCoordinate LocalCoordinateType;

        const int dimLocal = LocalCoordinateType::dimension;
        
        const Iterator end = space().end();
        for( Iterator it = space().begin(); it != end; ++it )
        {
          const Entity &entity = *it;
          const Geometry &geometry = entity.geometry();

          const Dune::ReferenceElement< DomainFieldType, dimLocal > &refElement
            = Dune::ReferenceElements< DomainFieldType, dimLocal >::general( geometry.type() );

          const LocalCoordinateType xlocal = geometry.local( x );
          if( refElement.checkInside( xlocal ) )
          {
            localFunction( entity ).evaluate( diffVariable, xlocal, ret );
            return;
          }
        }
        DUNE_THROW( RangeError, "PetscDiscreteFunction::evaluate: x is not within domain." );
      }

      /** \copydoc Dune::Fem::Function::evaluate(const DomainType &x,RangeType &ret) const */
      void evaluate ( const DomainType &x, RangeType &ret ) const
      {
        Dune::FieldVector< int, 0 > diffVariable;
        evaluate( diffVariable, x, ret );
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

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::space() const */
      const DiscreteFunctionSpaceType& space () const { return dfSpace_; }

      /** \brief obtain a constand pointer to the underlying PETSc Vec */
      const Vec* petscVec () const { return petscVector().getVector(); }

      /** \brief obtain a pointer to the underlying PETSc Vec */
      Vec* petscVec () { return petscVector().getVector(); }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::localFunction(const EntityType &entity) */
      template< class EntityType >
      LocalFunctionType localFunction ( const EntityType &entity )
      {
        return localFunctionStorage().localFunction( entity );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::localFunction(const EntityType &entity) const */
      template< class EntityType >
      LocalFunctionType localFunction ( const EntityType &entity ) const
      {
        return localFunctionStorage().localFunction( entity );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::name() const */
      const std::string& name () const
      {
        return name_;
      }

    private:
      PetscDiscreteFunction ();
      ThisType& operator= ( const ThisType &other );

      /*
       * private methods
       */
      LocalFunctionStorageType& localFunctionStorage () const 
      {
        return lfStorage_;
      }

      PetscVectorType& petscVector () { return petscVector_; }

      const PetscVectorType& petscVector () const { return petscVector_; }

      /*
       * data fields
       */
      std::string name_;
      const DiscreteFunctionSpaceType &dfSpace_;
      LocalFunctionFactoryType lfFactory_;
      mutable LocalFunctionStorageType lfStorage_;
      PetscVectorType petscVector_;
    };


  } // namespace Fem

} // namespace Dune

#endif // #if defined HAVE_PETSC

#endif // #ifndef DUNE_FEM_PETSCDISCRETEFUNCTION_HH
