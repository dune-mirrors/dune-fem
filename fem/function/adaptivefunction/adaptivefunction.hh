#ifndef DUNE_ADAPTIVEFUNCTION_HH
#define DUNE_ADAPTIVEFUNCTION_HH

//- System includes
#include <string>
#include <vector>

//- Dune includes
#include <dune/fem/space/common/dofstorage.hh>
#include <dune/fem/space/common/dofmanager.hh>

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/storage/subarray.hh>
#include <dune/fem/function/vectorfunction/vectorfunction.hh>

//- Local includes
#include "adaptiveimp.hh"
#include <dune/fem/function/localfunction/standardlocalfunction.hh>

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

    typedef StandardLocalFunctionFactory
      < AdaptiveDiscreteFunctionTraits< DiscreteFunctionSpaceType > >
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
  }; // end class AdaptiveDiscreteFunctionTraits



  //! @ingroup AdaptiveDFunction
  //! An adaptive discrete function
  //! This class is comparable to DFAdapt, except that it provides a 
  //! specialisation for CombinedSpace objects which provides enriched 
  //! functionality (access to subfunctions) and runtime optimisations
  template <class DiscreteFunctionSpaceImp>
  class AdaptiveDiscreteFunction : 
    public DiscreteFunctionDefault<
    AdaptiveDiscreteFunctionTraits<DiscreteFunctionSpaceImp > >,
    private AdaptiveFunctionImplementation<DiscreteFunctionSpaceImp > 
  {
  public:
    //! Discrete function space this discrete function belongs to
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

    //! Traits class with all necessary type definitions
    typedef AdaptiveDiscreteFunctionTraits< DiscreteFunctionSpaceType > Traits;

  private:
    typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > MyType;
    typedef AdaptiveFunctionImplementation< DiscreteFunctionSpaceType > Imp;
    typedef DiscreteFunctionDefault< Traits > BaseType;
    typedef AdaptiveDiscreteFunction<DiscreteFunctionSpaceImp> ThisType;

  public:
    using BaseType :: assign;

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

    typedef Mapping<DomainFieldType, RangeFieldType,
                    DomainType, RangeType> MappingType;

    typedef typename Traits :: LocalFunctionFactoryType LocalFunctionFactoryType;

  protected:
    const LocalFunctionFactoryType lfFactory_;

  public:
    //- Public methods
    //! Constructor
    AdaptiveDiscreteFunction( std :: string name,
                              const DiscreteFunctionSpaceType &spc )
    : BaseType( spc, lfFactory_ ),
      Imp( name, spc ),
      lfFactory_( *this )
    {
    }

    //! Constructor
    template <class VectorPointerType>
    AdaptiveDiscreteFunction( std :: string name,
                              const DiscreteFunctionSpaceType &spc,
                              VectorPointerType *vector)
    : BaseType( spc, lfFactory_ ),
      Imp( name, spc, vector ),
      lfFactory_( *this )
    {
    }

    //! Constructor for SubDiscreteFunctions
    //! This constructor is only called internally
    AdaptiveDiscreteFunction( std :: string name,
                              const DiscreteFunctionSpaceType &spc,
                              DofStorageType &dofVec )
    : BaseType( spc, lfFactory_ ),
      Imp(name, spc, dofVec),
      lfFactory_( *this )
    {
    }

    //! Copy constructor
    //! The copy constructor copies the dofs
    AdaptiveDiscreteFunction(const MyType &other)
    : BaseType( other.space(), lfFactory_ ),
      Imp( other ),
      lfFactory_( *this )
    {
    }

  private:
    ThisType &operator= ( const ThisType &other );
    
  public:
    /** \copydoc Dune::DiscreteFunctionInterface::assign(const DiscreteFunctionType &g) */
    inline void assign( const DiscreteFunctionType &g )
    {
      Imp :: assignFunction( g );
    }

    /** \copydoc Dune::DiscreteFunctionDefault::operator+=
     */ 
    inline ThisType &operator += ( const ThisType &g )
    {
      Imp :: addFunction( g );
      return *this;
    }

    /** \copydoc Dune::DiscreteFunctionDefault::operator-=
     */ 
    inline BaseType &operator-= ( const ThisType &g )
    {
      Imp::substractFunction(g);
      return *this;
    }

    /** \copydoc Dune::DiscreteFunctionDefault::addScaled
     */
    inline void addScaled( const ThisType &g,
                           const RangeFieldType &s )
    {
      Imp :: addScaled( g, s );
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
    using Imp::name;
    using Imp::size;
    using Imp::dbegin;
    using Imp::dend;
    using Imp::write_xdr;
    using Imp::read_xdr;
    using Imp::write_ascii;
    using Imp::read_ascii;
    using Imp::write_pgm;
    using Imp::read_pgm;

    using Imp::leakPointer;

#if 0
  protected:  
    using Imp::newObject;
#endif
    //- Forbidden members
  private:
    const MyType& interface() const { return *this; }
  }; // end class AdaptiveDiscreteFunction



  //- Specialisations
  //! Specialised version of AdaptiveDiscreteFunction for CombinedSpace

  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  class AdaptiveDiscreteFunction<
    CombinedSpace<ContainedFunctionSpaceImp, N, p> > : 
    public DiscreteFunctionDefault<AdaptiveDiscreteFunctionTraits<CombinedSpace<ContainedFunctionSpaceImp, N, p> > >,
    private AdaptiveFunctionImplementation<CombinedSpace<ContainedFunctionSpaceImp, N, p> >
  {
  private:
    typedef CombinedSpace<
      ContainedFunctionSpaceImp, N, p> DiscreteFunctionSpaceImp;
    typedef AdaptiveDiscreteFunction<DiscreteFunctionSpaceImp> MyType;
    typedef AdaptiveFunctionImplementation<DiscreteFunctionSpaceImp> Imp;
    typedef AdaptiveDiscreteFunctionTraits<DiscreteFunctionSpaceImp> MyTraits;
    typedef DiscreteFunctionDefault<MyTraits> BaseType;
    typedef AdaptiveDiscreteFunction<CombinedSpace<ContainedFunctionSpaceImp,N,p> > ThisType;

    using BaseType :: assign;

  public:
    //- Typedefs and enums
    typedef MyTraits Traits;
    typedef Imp ImplementationType;
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;


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
    AdaptiveDiscreteFunction(std::string name,
                             const DiscreteFunctionSpaceType& spc)
    : BaseType( spc, lfFactory_ ),
      Imp( name, spc ),
      lfFactory_( *this )
    {
      const SubSpaceType& subSpace = spc.containedSpace();
      for (int i=0;i<N;i++) {
        subDofMapper_[i] = new SubMapperType(this->spc_,i);
        subDofVector_[i] = new SubDofVectorType(this->dofStorage(), *subDofMapper_[i]);
        subDiscFunc_[i]  = new SubDiscreteFunctionType(
                               std::string("Subfunction of ")+this->name(),
                               subSpace,*(subDofVector_[i]));
      }
    }

    //! Constructor
    template <class VectorPointerType>
    AdaptiveDiscreteFunction(std::string name,
                             const DiscreteFunctionSpaceType& spc,
                             VectorPointerType * vector)
    : BaseType( spc, lfFactory_ ),
      Imp( name, spc , vector ),
      lfFactory_( *this ) {
      const SubSpaceType& subSpace = spc.containedSpace();
      for (int i=0;i<N;i++) {
        subDofMapper_[i] = new SubMapperType(this->spc_,i);
        subDofVector_[i] = new SubDofVectorType(this->dofStorage(), *subDofMapper_[i]);
        subDiscFunc_[i]  = new SubDiscreteFunctionType(
                               std::string("Subfunction of ")+this->name(),
                               subSpace,*(subDofVector_[i]));
      }
    }
    
    //! Constructor
    AdaptiveDiscreteFunction(std::string name,
                             const DiscreteFunctionSpaceType& spc,
                             DofStorageType& dofVec) 
    : BaseType( spc, lfFactory_ ),
      Imp( name, spc , dofVec ),
      lfFactory_( *this )
    {
      const SubSpaceType& subSpace = spc.containedSpace();
      for (int i=0;i<N;i++) {
        subDofMapper_[i] = new SubMapperType(this->spc_,i);
        subDofVector_[i] = new SubDofVectorType(this->dofStorage(), *subDofMapper_[i]);
        subDiscFunc_[i]  = new SubDiscreteFunctionType(
                               std::string("Subfunction of ")+this->name(),
                               subSpace,*(subDofVector_[i]));
      }
    }

    //! Copy constructor
    AdaptiveDiscreteFunction(const MyType& other)
    : BaseType( other.space(), lfFactory_ ),
      Imp(other),
      lfFactory_( *this )
    {
      const SubSpaceType& subSpace = other.space().containedSpace();
      for (int i=0;i<N;i++) {
        subDofMapper_[i] = new SubMapperType(this->spc_,i);
        subDofVector_[i] = new SubDofVectorType(this->dofStorage(), *subDofMapper_[i]);
        subDiscFunc_[i]  = new SubDiscreteFunctionType(
                               std::string("Subfunction of ")+this->name(),
                               subSpace,*(subDofVector_[i]));
      }
    }
    
    ~AdaptiveDiscreteFunction();
    

  private:
    // local function factory 
    const LocalFunctionFactoryType lfFactory_;

    ThisType &operator= ( const ThisType &other );
 
  public:
    /** \copydoc Dune::DiscreteFunctionDefault::assign(const DiscreteFunctionType &g) */
    inline void assign( const DiscreteFunctionType &g )
    {
      Imp :: assignFunction( g );
    }

    /** \copydoc Dune::DiscreteFunctionDefault::operator+= */
    inline ThisType &operator += ( const ThisType &g )
    {
      Imp::addFunction(g);
      return *this;
    }

    /** \copydoc Dune::DiscreteFunctionDefault::operator-= */
    inline ThisType &operator-= ( const ThisType &g )
    {
      Imp::substractFunction(g);
      return *this;
    }

    /** \copydoc Dune::DiscreteFunctionDefault::addScaled
     */
    inline void addScaled( const ThisType &g,
                           const RangeFieldType &s )
    {
      Imp :: addScaled( g, s );
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
    using Imp::name;
    using Imp::size;
    using Imp::dbegin;
    using Imp::dend;
#if 0
    /** \brief  @copydoc DiscreteFunctionDefault::localFunction */
    template <class EntityType> 
    LocalFunctionType localFunction (const EntityType &en) { return LocalFunctionType(en,*this); }
    /** \brief  @copydoc DiscreteFunctionDefault::localFunction */
    template <class EntityType> 
    const LocalFunctionType localFunction (const EntityType &en) const { return LocalFunctionType(en,*this); }
    // using Imp::localFunction;
#endif
    using Imp::write_xdr;
    using Imp::read_xdr;
    using Imp::write_ascii;
    using Imp::read_ascii;
    using Imp::write_pgm;
    using Imp::read_pgm;

    using Imp::leakPointer;
    
    //- Additional methods
    SubDiscreteFunctionType& subFunction(int component);

    int numComponents() const { return N; }

#if 0
  public:
    friend class DiscreteFunctionInterface<MyTraits>;

  protected:
    using Imp::newObject;
#endif
    
  private:
    const MyType& interface() const { return *this; }
    SubMapperType* subDofMapper_[N];
    SubDofVectorType* subDofVector_[N]; 
    SubDiscreteFunctionType* subDiscFunc_[N];
  }; // end class AdaptiveDiscreteFunction (specialised for CombinedSpace)
  
} // end namespace Dune

#include "adaptivefunction.cc"

#endif
