#ifndef DUNE_ADAPTIVEFUNCTION_HH
#define DUNE_ADAPTIVEFUNCTION_HH

//- System includes
#include <string>
#include <vector>

//- Dune includes
#include <dune/common/typetraits.hh>
#include <dune/fem/space/common/dofstorage.hh>
#include <dune/fem/space/common/dofmanager.hh>

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/storage/subarray.hh>
#include <dune/fem/function/vectorfunction/vectorfunction.hh>

//- Local includes
#include "adaptiveimp.hh"
#include <dune/fem/function/localfunction/standardlocalfunction.hh>
#include <dune/fem/function/localfunction/genericlocalfunction.hh>

namespace Dune
{
  
  //- Forward declarations
  template <class DiscreteFunctionSpaceImp>
  class AdaptiveDiscreteFunction;

  //- Forward declarations of Combined Space 
  template <class, int , DofStoragePolicy > 
  class CombinedSpace;

  
  //- Class definitions
  //! Traits class for AdaptiveDiscreteFunction and AdaptiveLocalFunction
  template< class DiscreteFunctionSpaceImp >
  struct AdaptiveDiscreteFunctionTraits
  {
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;
 
    typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

    typedef AdaptiveFunctionImplementation< DiscreteFunctionSpaceType >
      ImplementationType;

    typedef AdaptiveDiscreteFunctionTraits< DiscreteFunctionSpaceType > Traits;

    static const bool isGenericSpace = Conversion< DiscreteFunctionSpaceType, GenericDiscreteFunctionSpace >::exists;
    typedef typename SelectType< isGenericSpace, GenericLocalFunctionFactory< Traits >, StandardLocalFunctionFactory< Traits > >::Type
      LocalFunctionFactoryType;

    typedef LocalFunctionStack< LocalFunctionFactoryType > LocalFunctionStorageType;

    typedef typename LocalFunctionStorageType :: LocalFunctionType LocalFunctionType;

    typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;

    typedef typename DiscreteFunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType;
    typedef RangeFieldType DofType;

    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType JacobianRangeType;
    typedef typename DiscreteFunctionSpaceType :: MapperType MapperType;
    typedef typename DiscreteFunctionSpaceType :: GridType GridType;

    // type of Array seen by functions 
    typedef StaticArray<DofType> DofStorageType;
    // tpye of array created 
    typedef MutableArray<DofType> MutableDofStorageType;
     
    typedef typename DofStorageType :: DofIteratorType DofIteratorType;
    typedef typename DofStorageType :: ConstDofIteratorType ConstDofIteratorType;
 
    typedef DofManager< GridType > DofManagerType;

    typedef typename ImplementationType :: DofBlockType DofBlockType;
    typedef typename ImplementationType :: ConstDofBlockType ConstDofBlockType;
    typedef typename ImplementationType :: DofBlockPtrType DofBlockPtrType;
    typedef typename ImplementationType :: ConstDofBlockPtrType ConstDofBlockPtrType;
  }; // end class AdaptiveDiscreteFunctionTraits



  //! @ingroup AdaptiveDFunction
  //! An adaptive discrete function
  //! This class is comparable to DFAdapt, except that it provides a 
  //! specialisation for CombinedSpace objects which provides enriched 
  //! functionality (access to subfunctions) and runtime optimisations
  template <class DiscreteFunctionSpaceImp>
  class AdaptiveDiscreteFunction
  : public DiscreteFunctionDefault< AdaptiveDiscreteFunctionTraits< DiscreteFunctionSpaceImp > >,
    private AdaptiveFunctionImplementation< DiscreteFunctionSpaceImp >
  {
    typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceImp > ThisType;
    typedef DiscreteFunctionDefault< AdaptiveDiscreteFunctionTraits< DiscreteFunctionSpaceImp > >
      BaseType;

    typedef AdaptiveFunctionImplementation< DiscreteFunctionSpaceImp > Imp;

  public:
    //! Discrete function space this discrete function belongs to
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

    //! Traits class with all necessary type definitions
    typedef AdaptiveDiscreteFunctionTraits< DiscreteFunctionSpaceType > Traits;

    typedef typename BaseType::DiscreteFunctionInterfaceType DiscreteFunctionInterfaceType;

    using BaseType::assign;

  public:
    //- Typedefs and enums
    //! Class containing the actual implementation
    typedef Imp ImplementationType;

    //! Local function implementation 
    typedef typename Traits::GridType GridType;
    
    //! Discrete function type (identical to this type, needed as 
    //! Barton-Nackman parameter
    typedef typename Traits::DiscreteFunctionType DiscreteFunctionType;
    
    //! Intrinsic type used for the dofs (typically a float type)
    typedef typename Traits::DofType DofType;
    //! Intrinsic type used for the range field (identical to DofType)
    typedef typename Traits::RangeFieldType RangeFieldType;
    //! Intrinsic type used for the domain field
    typedef typename Traits::DomainFieldType DomainFieldType;
    //! Vector type used for the range field
    typedef typename Traits::RangeType RangeType;
    //! Vector type used for the domain field
    typedef typename Traits::DomainType DomainType;
    //! Mapper type (from the space)
    typedef typename Traits::MapperType MapperType;
    
    //! Container class type for the dofs (managed by the DofManager)
    typedef typename Traits :: MutableDofStorageType MutableDofStorageType;
    //! Container class type for the dofs (managed by the DofManager)
    typedef typename Traits::DofStorageType DofStorageType;
 
    //! Iterator over dof container
    typedef typename Traits::DofIteratorType DofIteratorType;
    //! Read-only iterator over dof container
    typedef typename Traits::ConstDofIteratorType ConstDofIteratorType;
    
    typedef typename Traits :: DofBlockType DofBlockType;
    typedef typename Traits :: ConstDofBlockType ConstDofBlockType;
    typedef typename Traits :: DofBlockPtrType DofBlockPtrType;
    typedef typename Traits :: ConstDofBlockPtrType ConstDofBlockPtrType;

    typedef Mapping<DomainFieldType, RangeFieldType,
                    DomainType, RangeType> MappingType;

    typedef typename Traits :: LocalFunctionFactoryType LocalFunctionFactoryType;

  protected:
    const LocalFunctionFactoryType lfFactory_;

  public:
    //- Public methods
    //! Constructor
    AdaptiveDiscreteFunction( const std :: string &name,
                              const DiscreteFunctionSpaceType &spc )
    : BaseType( name, spc, lfFactory_ ),
      Imp( name, spc ),
      lfFactory_( *this )
    {}

    //! Constructor
    template <class VectorPointerType>
    AdaptiveDiscreteFunction( const std :: string &name,
                              const DiscreteFunctionSpaceType &spc,
                              VectorPointerType *vector)
    : BaseType( name, spc, lfFactory_ ),
      Imp( name, spc, vector ),
      lfFactory_( *this )
    {}

    //! Constructor for SubDiscreteFunctions
    //! This constructor is only called internally
    AdaptiveDiscreteFunction( const std :: string &name,
                              const DiscreteFunctionSpaceType &spc,
                              DofStorageType &dofVec )
    : BaseType( name, spc, lfFactory_ ),
      Imp( name, spc, dofVec ),
      lfFactory_( *this )
    {}

    //! Copy constructor
    //! The copy constructor copies the dofs
    AdaptiveDiscreteFunction( const ThisType & other )
    : BaseType( "copy of " + other.name(), other.space(), lfFactory_ ),
      Imp( BaseType :: name(), other ),
      lfFactory_( *this )
    {}

  private:
    ThisType &operator= ( const ThisType &other );
    
  public:
    /** \copydoc Dune::DiscreteFunctionInterface::assign(const DiscreteFunctionInterfaceType &g) */
    void assign ( const DiscreteFunctionInterfaceType &g )
    {
      Imp::assignFunction( asImp( g ) );
    }

    /** \copydoc Dune::DiscreteFunctionInterface::operator+=(const DiscreteFunctionInterfaceType &g) */ 
    ThisType &operator+= ( const DiscreteFunctionInterfaceType &g )
    {
      Imp::addFunction( asImp( g ) );
      return *this;
    }

    /** \copydoc Dune::DiscreteFunctionInterface::operator-=(const DFType &g) */ 
    ThisType &operator-= ( const DiscreteFunctionInterfaceType &g )
    {
      Imp::substractFunction( asImp( g ) );
      return *this;
    }

    /** \copydoc Dune::DiscreteFunctionDefault::addScaled(const DiscreteFunctionInterfaceType &g,RangeFieldType &s) */
    void addScaled( const DiscreteFunctionInterfaceType &g, const RangeFieldType &s )
    {
      Imp::addScaled( asImp( g ), s );
    }

  protected:
    using Imp :: dofVec_;

  public:
    inline const RangeFieldType &dof ( unsigned int index ) const
    {
      return dofVec_[ index ];
    }

    inline RangeFieldType &dof ( unsigned int index )
    {
      return dofVec_[ index ];
    }

    using Imp::clear;
    using Imp::size;
    using Imp::dbegin;
    using Imp::dend;

#if DUNE_FEM_COMPATIBILITY
    using Imp::write_xdr;
    using Imp::read_xdr;
    using Imp::write_ascii;
    using Imp::read_ascii;
#endif

    using Imp::leakPointer;
    using Imp::block;
    using Imp::enableDofCompression;

#if 0
    //- Forbidden members
  private:
    const ThisType &interface () const
    {
      return *this;
    }
#endif
  }; // end class AdaptiveDiscreteFunction
 


  //- Specialisations
  //! Specialised version of AdaptiveDiscreteFunction for CombinedSpace

  template< class ContainedFunctionSpaceImp, int N, DofStoragePolicy p >
  class AdaptiveDiscreteFunction< CombinedSpace<ContainedFunctionSpaceImp, N, p > >
  : public DiscreteFunctionDefault< AdaptiveDiscreteFunctionTraits< CombinedSpace< ContainedFunctionSpaceImp, N, p > > >,
    private AdaptiveFunctionImplementation< CombinedSpace< ContainedFunctionSpaceImp, N, p > >
  {
    typedef AdaptiveDiscreteFunction< CombinedSpace< ContainedFunctionSpaceImp, N, p > > ThisType;
    typedef DiscreteFunctionDefault< AdaptiveDiscreteFunctionTraits< CombinedSpace< ContainedFunctionSpaceImp, N, p > > >
      BaseType;

    typedef AdaptiveFunctionImplementation< CombinedSpace< ContainedFunctionSpaceImp, N, p > > Imp;

  public:
    typedef CombinedSpace< ContainedFunctionSpaceImp, N, p > DiscreteFunctionSpaceType;

    using BaseType::assign;

    //- Typedefs and enums
    typedef AdaptiveDiscreteFunctionTraits< DiscreteFunctionSpaceType > Traits;

    typedef typename BaseType::DiscreteFunctionInterfaceType DiscreteFunctionInterfaceType;

    typedef Imp ImplementationType;

    typedef typename Traits::LocalFunctionType LocalFunctionType;
    typedef typename Traits::DiscreteFunctionType DiscreteFunctionType;
    
    typedef typename Traits::DofType DofType;
    typedef typename Traits::RangeFieldType RangeFieldType;
    typedef typename Traits::DomainFieldType DomainFieldType;
    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::MapperType MapperType;
    
    typedef typename Traits::DofStorageType DofStorageType;
     
    typedef typename Traits::DofIteratorType DofIteratorType;
    typedef typename Traits::ConstDofIteratorType ConstDofIteratorType;

    typedef typename Traits :: DofBlockType DofBlockType;
    typedef typename Traits :: ConstDofBlockType ConstDofBlockType;
    typedef typename Traits :: DofBlockPtrType DofBlockPtrType;
    typedef typename Traits :: ConstDofBlockPtrType ConstDofBlockPtrType;

    //! type of local functions factory 
    typedef typename Traits :: LocalFunctionFactoryType  LocalFunctionFactoryType;
    //- Additional typedefs
    typedef typename DiscreteFunctionSpaceType::ContainedDiscreteFunctionSpaceType SubSpaceType;
    typedef typename DiscreteFunctionSpaceType::SubMapperType SubMapperType;
    typedef SubVector<DofStorageType,SubMapperType> SubDofVectorType;
    typedef VectorDiscreteFunction<SubSpaceType,SubDofVectorType> SubDiscreteFunctionType;
   
    typedef Mapping<DomainFieldType, RangeFieldType,
                    DomainType, RangeType> MappingType;

  public:
    //- Public methods
    //! Constructor
    AdaptiveDiscreteFunction ( const std :: string &name,
                               const DiscreteFunctionSpaceType &spc)
    : BaseType( name, spc, lfFactory_ ),
      Imp( name, spc ),
      lfFactory_( *this )
    {
      initializeSubFunctions();
    }

    //! Constructor
    template <class VectorPointerType>
    AdaptiveDiscreteFunction ( const std :: string &name,
                               const DiscreteFunctionSpaceType &spc,
                               VectorPointerType *vector )
    : BaseType( name, spc, lfFactory_ ),
      Imp( name, spc , vector ),
      lfFactory_( *this ) 
    {
      initializeSubFunctions();
    }
    
    //! Constructor
    AdaptiveDiscreteFunction ( const std :: string &name,
                               const DiscreteFunctionSpaceType &spc,
                               DofStorageType &dofVec ) 
    : BaseType( name, spc, lfFactory_ ),
      Imp( name, spc , dofVec ),
      lfFactory_( *this )
    {
      initializeSubFunctions();
    }

    //! Copy constructor
    AdaptiveDiscreteFunction ( const ThisType &other )
    : BaseType( "copy of " + other.name(), other.space(), lfFactory_ ),
      Imp( BaseType :: name(), other ),
      lfFactory_( *this )
    {
      initializeSubFunctions();
    }
    
    //! destructor 
    ~AdaptiveDiscreteFunction();

  protected:  
    // initialize sub functions 
    void initializeSubFunctions() ;

  private:
    // local function factory 
    const LocalFunctionFactoryType lfFactory_;

    ThisType &operator= ( const ThisType &other );
 
  public:
    /** \copydoc Dune::DiscreteFunctionDefault::assign(const DiscreteFunctionType &g) */
    void assign ( const DiscreteFunctionInterfaceType &g )
    {
      Imp::assignFunction( asImp( g ) );
    }

    /** \copydoc Dune::DiscreteFunctionInterface::operator+=(const DiscreteFunctionType &g) */
    ThisType &operator+= ( const DiscreteFunctionInterfaceType &g )
    {
      Imp::addFunction( asImp( g ) );
      return *this;
    }

    /** \copydoc Dune::DiscreteFunctionInterface::operator-=(const DFType &g) */
    ThisType &operator-= ( const DiscreteFunctionInterfaceType &g )
    {
      Imp::substractFunction( asImp( g ) );
      return *this;
    }

    /** \copydoc Dune::DiscreteFunctionDefault::addScaled(const DiscreteFunctionInterfaceType &g,RangeFieldType &s) */
    void addScaled( const DiscreteFunctionInterfaceType &g, const RangeFieldType &s )
    {
      Imp::addScaled( asImp( g ), s );
    }

  protected:
    using Imp :: dofVec_;

  public:
    inline const RangeFieldType &dof ( unsigned int index ) const
    {
      return dofVec_[ index ];
    }

    inline RangeFieldType &dof ( unsigned int index )
    {
      return dofVec_[ index ];
    }

    using Imp::clear;
    using Imp::size;
    using Imp::dbegin;
    using Imp::dend;

#if DUNE_FEM_COMPATIBILITY
    using Imp::write_xdr;
    using Imp::read_xdr;
    using Imp::write_ascii;
    using Imp::read_ascii;
#endif
    using Imp::leakPointer;
    using Imp::block;
    using Imp::enableDofCompression;
    
    //- Additional methods
    SubDiscreteFunctionType& subFunction(int component) {
      return *subDiscFunc_[component];
    }
    const SubDiscreteFunctionType& 
      subFunction(int component) const {
      return *subDiscFunc_[component];
    }

    int numComponents() const { return N; }

#if 0
  private:
    const ThisType &interface () const
    {
      return *this;
    }
#endif

    SubMapperType* subDofMapper_[N];
    SubDofVectorType* subDofVector_[N]; 
    SubDiscreteFunctionType* subDiscFunc_[N];
  }; // end class AdaptiveDiscreteFunction (specialised for CombinedSpace)
  
template <class DiscreteFunctionSpaceImp>
class ManagedDiscreteFunction<AdaptiveDiscreteFunction<DiscreteFunctionSpaceImp> > :
public AdaptiveDiscreteFunction<DiscreteFunctionSpaceImp> {
  typedef AdaptiveDiscreteFunction<DiscreteFunctionSpaceImp> BaseType;
public:
  typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;
  typedef ManagedDiscreteFunction<BaseType> ThisType;
  //! Constructor
  ManagedDiscreteFunction(std::string name,
                           const DiscreteFunctionSpaceType& spc)
  : BaseType( name, spc) {}
  //! Constructor
  template <class VectorPointerType>
  ManagedDiscreteFunction(std::string name,
                           const DiscreteFunctionSpaceType& spc,
                           VectorPointerType * vector)
  : BaseType( name,spc,vector) {}
  //! Constructor
  ManagedDiscreteFunction(std::string name,
                           const DiscreteFunctionSpaceType& spc,
                           typename BaseType::DofStorageType& dofVec) 
  : BaseType( name,spc,dofVec) {}
  //! Copy constructor
  ManagedDiscreteFunction(const ThisType& other)
  : BaseType( other ) {} 
};
} // end namespace Dune

#include "adaptivefunction.cc"

#endif // #ifndef DUNE_ADAPTIVEFUNCTION_HH
