#ifndef DUNE_FEM_ADAPTIVEFUNCTION_HH
#define DUNE_FEM_ADAPTIVEFUNCTION_HH

//- System includes
#include <string>
#include <vector>

//- Dune includes
#include <dune/common/typetraits.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/storage/subarray.hh>

//- Local includes
#include "adaptiveimp.hh"
#include <dune/fem/function/localfunction/mutable.hh>
#include <dune/fem/common/referencevector.hh>
#include <dune/fem/common/stackallocator.hh>

namespace Dune
{

  namespace Fem
  {

    //- Forward declarations
    template <class DiscreteFunctionSpaceImp>
    class AdaptiveDiscreteFunction;

    template <class DiscreteFunctionImp>
    class SubFunctionStorage;

    //- Class definitions
    //! Traits class for AdaptiveDiscreteFunction
    template< class DiscreteFunctionSpace >
    struct DiscreteFunctionTraits< AdaptiveDiscreteFunction< DiscreteFunctionSpace > >
    {
      typedef AdaptiveDiscreteFunction< DiscreteFunctionSpace > DiscreteFunctionType;
    private:

      typedef DiscreteFunctionTraits< DiscreteFunctionType > Traits;

      typedef AdaptiveFunctionImplementation< DiscreteFunctionSpace > ImplementationType;

    public:
      typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

      typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
      typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;

      typedef typename DiscreteFunctionSpaceType :: DomainFieldType DomainFieldType;
      typedef typename DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType;
      typedef RangeFieldType DofType;

      typedef typename DiscreteFunctionSpaceType :: JacobianRangeType JacobianRangeType;
      typedef typename DiscreteFunctionSpaceType :: BlockMapperType BlockMapperType;
      typedef typename DiscreteFunctionSpaceType :: GridType GridType;

      typedef ThreadSafeValue< UninitializedObjectStack > LocalDofVectorStackType;
      typedef StackAllocator< DofType, LocalDofVectorStackType* > LocalDofVectorAllocatorType;
      typedef DynamicReferenceVector< DofType, LocalDofVectorAllocatorType > LocalDofVectorType;

      typedef MutableLocalFunction< DiscreteFunctionType > LocalFunctionType;

      // type used to determine size etc.
      typedef NonBlockMapper< BlockMapperType,
                              DiscreteFunctionSpaceType :: localBlockSize > MapperType;

      // type of Array seen by functions
      typedef StaticArray<DofType>  DofStorageType;
      // tpye of array created
      typedef MutableArray<DofType> MutableDofStorageType;

      typedef typename DofStorageType :: DofIteratorType DofIteratorType;
      typedef typename DofStorageType :: ConstDofIteratorType ConstDofIteratorType;

      typedef DofManager< GridType > DofManagerType;

      typedef typename ImplementationType :: DofBlockType DofBlockType;
      typedef typename ImplementationType :: ConstDofBlockType ConstDofBlockType;
      typedef typename ImplementationType :: DofBlockPtrType DofBlockPtrType;
      typedef typename ImplementationType :: ConstDofBlockPtrType ConstDofBlockPtrType;
    };



    //! @ingroup AdaptiveDFunction
    //! An adaptive discrete function
    //! This class is comparable to DFAdapt, except that it provides a
    //! specialisation for CombinedSpace objects which provides enriched
    //! functionality (access to subfunctions) and runtime optimisations
    template <class DiscreteFunctionSpaceImp>
    class AdaptiveDiscreteFunction
    : public DiscreteFunctionDefault< AdaptiveDiscreteFunction< DiscreteFunctionSpaceImp > >,
      private AdaptiveFunctionImplementation< DiscreteFunctionSpaceImp >
    {
      typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceImp > ThisType;
      typedef DiscreteFunctionDefault< AdaptiveDiscreteFunction< DiscreteFunctionSpaceImp > > BaseType;

      typedef AdaptiveFunctionImplementation< DiscreteFunctionSpaceImp > Imp;

    public:
      //! Discrete function space this discrete function belongs to
      typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

      //! Traits class with all necessary type definitions
      typedef DiscreteFunctionTraits< AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > > Traits;

      typedef typename BaseType::DiscreteFunctionInterfaceType DiscreteFunctionInterfaceType;

      using BaseType::assign;
      using BaseType::axpy;
      using BaseType::space;
      using BaseType::asImp;
      using BaseType :: operator+=;
      using BaseType :: operator-=;

    public:
      //- Typedefs and enums
      //! Class containing the actual implementation
      typedef Imp ImplementationType;

      //! Local function implementation
      typedef typename BaseType::GridType GridType;

      //! Discrete function type (identical to this type, needed as
      //! Barton-Nackman parameter
      typedef typename BaseType::DiscreteFunctionType DiscreteFunctionType;

      //! Intrinsic type used for the dofs (typically a float type)
      typedef typename BaseType::DofType DofType;
      //! Intrinsic type used for the range field (identical to DofType)
      typedef typename BaseType::RangeFieldType RangeFieldType;
      //! Intrinsic type used for the domain field
      typedef typename BaseType::DomainFieldType DomainFieldType;
      //! Vector type used for the range field
      typedef typename BaseType::RangeType RangeType;
      //! Vector type used for the domain field
      typedef typename BaseType::DomainType DomainType;

      //! Container class type for the dofs (managed by the DofManager)
      typedef typename Traits::MutableDofStorageType MutableDofStorageType;
      //! Container class type for the dofs (managed by the DofManager)
      typedef typename Traits::DofStorageType DofStorageType;
       //! Mapper type (from the space)
      typedef typename Traits::MapperType MapperType;

      //! Iterator over dof container
      typedef typename BaseType::DofIteratorType DofIteratorType;
      //! Read-only iterator over dof container
      typedef typename BaseType::ConstDofIteratorType ConstDofIteratorType;

      typedef typename BaseType :: DofBlockType DofBlockType;
      typedef typename BaseType :: ConstDofBlockType ConstDofBlockType;
      typedef typename BaseType :: DofBlockPtrType DofBlockPtrType;
      typedef typename BaseType :: ConstDofBlockPtrType ConstDofBlockPtrType;

      typedef typename BaseType :: LocalDofVectorAllocatorType LocalDofVectorAllocatorType;

    protected:
      typename Traits :: LocalDofVectorStackType ldvStack_;

    public:
      //- Public methods
      //! Constructor
      AdaptiveDiscreteFunction( const std :: string &name,
                                const DiscreteFunctionSpaceType &spc )
      : BaseType( name, spc, LocalDofVectorAllocatorType( &ldvStack_ ) ),
        Imp( name, spc ),
        ldvStack_( std::max( sizeof( DofType ), sizeof( DofType* ) ) * spc.blockMapper().maxNumDofs() * DiscreteFunctionSpaceType::localBlockSize )
      {
        spc.addFunction( *this );
      }

      //! Constructor
      template <class VectorPointerType>
      AdaptiveDiscreteFunction( const std :: string &name,
                                const DiscreteFunctionSpaceType &spc,
                                VectorPointerType *vector)
      : BaseType( name, spc, LocalDofVectorAllocatorType( &ldvStack_ ) ),
        Imp( name, spc, vector ),
        ldvStack_( std::max( sizeof( DofType ), sizeof( DofType* ) ) * spc.blockMapper().maxNumDofs() * DiscreteFunctionSpaceType::localBlockSize )
      {}

      //! Constructor for SubDiscreteFunctions
      //! This constructor is only called internally
      AdaptiveDiscreteFunction( const std :: string &name,
                                const DiscreteFunctionSpaceType &spc,
                                DofStorageType &dofVec )
      : BaseType( name, spc, LocalDofVectorAllocatorType( &ldvStack_ ) ),
        Imp( name, spc, dofVec ),
        ldvStack_( std::max( sizeof( DofType ), sizeof( DofType* ) ) * spc.blockMapper().maxNumDofs() * DiscreteFunctionSpaceType::localBlockSize )
      {}

      //! Copy constructor
      //! The copy constructor copies the dofs
      AdaptiveDiscreteFunction( const ThisType & other )
      : BaseType( "copy of " + other.name(), other.space(), LocalDofVectorAllocatorType( &ldvStack_ ) ),
        Imp( BaseType :: name(), other ),
        ldvStack_( other.ldvStack_ )
      {
        space().addFunction( *this );
      }

      //! destructor
      ~AdaptiveDiscreteFunction()
      {
        space().removeFunction( *this );
      }

    private:
      ThisType &operator= ( const ThisType &other );

    public:
      /** \copydoc Dune::Fem::DiscreteFunctionInterface::assign(const DiscreteFunctionInterfaceType &g) */
      void assign ( const DiscreteFunctionInterfaceType &g )
      {
        Imp::assignFunction( asImp( g ) );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::operator+=(const DiscreteFunctionInterfaceType &g) */
      ThisType &operator+= ( const DiscreteFunctionInterfaceType &g )
      {
        Imp::addFunction( asImp( g ) );
        return *this;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::operator-=(const DFType &g) */
      ThisType &operator-= ( const DiscreteFunctionInterfaceType &g )
      {
        Imp::substractFunction( asImp( g ) );
        return *this;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::axpy(const RangeFieldType &s,const DiscreteFunctionInterfaceType &g) */
      void axpy( const RangeFieldType &s, const DiscreteFunctionInterfaceType &g )
      {
        Imp::axpy( s, asImp( g ) );
      }

    protected:
      using Imp :: dofVec_;
      using Imp :: dofStorage;

    public:
      inline const DofType &dof ( unsigned int index ) const
      {
        return dofVec_[ index ];
      }

      inline DofType &dof ( unsigned int index )
      {
        return dofVec_[ index ];
      }

      using Imp::resize ;
      using Imp::clear;
      using Imp::size;
      using Imp::dbegin;
      using Imp::dend;

      using Imp::leakPointer;
      using Imp::block;
      using Imp::enableDofCompression;

      friend class SubFunctionStorage < ThisType >;
    }; // end class AdaptiveDiscreteFunction


  } // end namespace Fem

} // end namespace Dune
#endif // #ifndef DUNE_FEM_ADAPTIVEFUNCTION_HH
