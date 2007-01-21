#ifndef DUNE_ADAPTIVEFUNCTION_HH
#define DUNE_ADAPTIVEFUNCTION_HH

//- System includes
#include <string>
#include <vector>

//- Dune includes
#include <dune/common/array.hh>
#include "../common/discretefunction.hh"
#include "../common/localfunction.hh"
#include "../../space/combinedspace/combinedspace.hh"
#include "../../space/combinedspace/subspace.hh"
#include "../../space/common/dofmanager.hh"

//- Local includes
#include "adaptiveimp.hh"

namespace Dune {
  
  //- Forward declarations
  template <class DiscreteFunctionSpaceImp>
  class AdaptiveDiscreteFunction;
  template <class DiscreteFunctionSpaceImp>
  class AdaptiveLocalFunction;
  
  //- Class definitions
  //! Traits class for AdaptiveDiscreteFunction and AdaptiveLocalFunction
  template <class DiscreteFunctionSpaceImp>
  struct AdaptiveDiscreteFunctionTraits {
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;
 
    typedef AdaptiveDiscreteFunction<
      DiscreteFunctionSpaceImp> DiscreteFunctionType;
    // the local functions implementation 
    typedef AdaptiveLocalFunction<
      DiscreteFunctionSpaceImp> LocalFunctionImp;
    
    // local function type 
    typedef LocalFunctionWrapper<
      DiscreteFunctionType> LocalFunctionType;

    typedef typename DiscreteFunctionSpaceType::Traits::RangeFieldType DofType;
    typedef typename DiscreteFunctionSpaceType::Traits::RangeFieldType RangeFieldType;
    typedef typename DiscreteFunctionSpaceType::Traits::RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType::Traits::DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType::Traits::JacobianRangeType JacobianRangeType;
    typedef typename DiscreteFunctionSpaceType::Traits::MapperType MapperType;
    typedef typename DiscreteFunctionSpaceType::Traits::GridType GridType;

    typedef DofArray<DofType> DofStorageType;
     
    typedef typename DofStorageType::DofIteratorType DofIteratorType;
    typedef typename DofStorageType::ConstDofIteratorType ConstDofIteratorType;
 
    typedef DofManager<GridType> DofManagerType;
   
  }; // end class AdaptiveDiscreteFunctionTraits


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
    //- friends
    
  private:
    typedef AdaptiveDiscreteFunction<
      DiscreteFunctionSpaceImp > MyType;
    typedef AdaptiveFunctionImplementation<
      DiscreteFunctionSpaceImp > Imp;
    typedef AdaptiveDiscreteFunctionTraits<
      DiscreteFunctionSpaceImp > MyTraits;
    typedef DiscreteFunctionDefault<MyTraits> BaseType;

  public:
    //- Typedefs and enums
    //! Traits class with all necessary type definitions
    typedef MyTraits Traits;
    //! Class containing the actual implementation
    typedef Imp ImplementationType;
    //! Discrete function space this discrete function belongs to
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

    //! Local function implementation 
    typedef typename Traits::LocalFunctionImp LocalFunctionImp;
    typedef typename Traits::GridType GridType;
    
    //! Local function type
    typedef typename Traits::LocalFunctionType LocalFunctionType;
    
    //! Discrete function type (identical to this type, needed as 
    //! Barton-Nackman parameter
    typedef typename Traits::DiscreteFunctionType DiscreteFunctionType;
    
    //! Intrinsic type used for the dofs (typically a float type)
    typedef typename Traits::DofType DofType;
    //! Intrinsic type used for the range field (identical to DofType)
    typedef typename Traits::RangeFieldType RangeFieldType;
    //! Vector type used for the range field
    typedef typename Traits::RangeType RangeType;
    //! Vector type used for the domain field
    typedef typename Traits::DomainType DomainType;
    //! Mapper type (from the space)
    typedef typename Traits::MapperType MapperType;
    
    //! Container class type for the dofs (managed by the DofManager)
    typedef typename Traits::DofStorageType DofStorageType;
 
    //! Iterator over dof container
    typedef typename Traits::DofIteratorType DofIteratorType;
    //! Read-only iterator over dof container
    typedef typename Traits::ConstDofIteratorType ConstDofIteratorType;
  public:
    //- Public methods
    //! Constructor
    AdaptiveDiscreteFunction(std::string name,
                             const DiscreteFunctionSpaceType& spc) :
      BaseType(spc),
      Imp(name, spc) 
    {}

    //! Constructor
    template <class VectorPointerType>
    AdaptiveDiscreteFunction(std::string name,
                             const DiscreteFunctionSpaceType& spc,
                             VectorPointerType * vector) :
      BaseType(spc),
      Imp(name, spc, vector) 
    {}

    //! Constructor for SubDiscreteFunctions
    //! This constructor is only called internally
    AdaptiveDiscreteFunction(std::string name,
                             const DiscreteFunctionSpaceType& spc,
                             DofStorageType& dofVec) :
      BaseType(spc),
      Imp(name, spc, dofVec)
    {}

    //! Copy constructor
    //! The copy constructor copies the dofs
    AdaptiveDiscreteFunction(const MyType& other) :
      BaseType(other.space()),
      Imp(other)
    {}

    using Imp::name;
    using Imp::size;
    using Imp::dbegin;
    using Imp::dend;
    using Imp::newLocalFunctionObject;
    //! return empty local function 
    LocalFunctionType newLocalFunction () DUNE_DEPRECATED { return LocalFunctionType(*this); }
    using Imp::localFunction;
    using Imp::write_xdr;
    using Imp::read_xdr;
    using Imp::write_ascii;
    using Imp::read_ascii;
    using Imp::write_pgm;
    using Imp::read_pgm;

    using Imp::leakPointer;
  private:
    //- Forbidden members

    const MyType& interface() const { return *this; }
  }; // end class AdaptiveDiscreteFunction

  // Note: could use Traits class for Barton-Nackman instead
  //! Local function belonging to AdaptiveDiscreteFunction
  template <class DiscreteFunctionSpaceImp>
  class AdaptiveLocalFunction : 
    public LocalFunctionDefault<
    DiscreteFunctionSpaceImp,
    AdaptiveLocalFunction<DiscreteFunctionSpaceImp > >
  {
  public:
    friend class AdaptiveFunctionImplementation<DiscreteFunctionSpaceImp>;
    friend class LocalFunctionWrapper<
      AdaptiveDiscreteFunction< DiscreteFunctionSpaceImp > > ;

  private:
    typedef AdaptiveLocalFunction<DiscreteFunctionSpaceImp> ThisType;
  public:
    //- Public typedefs and enums
    //! The discrete function space this local function belongs to
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;
    //! The discrete function this local function belongs to
    typedef AdaptiveDiscreteFunction<
      DiscreteFunctionSpaceImp > DiscreteFunctionType;
    //! Traits class with type definitions for AdaptiveDiscreteFunction and
    //! AdaptiveLocalFunction
    typedef AdaptiveDiscreteFunctionTraits<
      DiscreteFunctionSpaceType > Traits;
    //! Traits class of DiscreteFunctionSpaceType
    typedef typename DiscreteFunctionSpaceType::Traits SpaceTraits;

    //! Function space type
    typedef typename SpaceTraits::FunctionSpaceType FunctionSpaceType;
    //! The base function set of DiscreteFunctionSpaceType
    
    //! Intrinsic data type for range field
    typedef typename Traits::RangeFieldType RangeFieldType;
    //! Vector type for the domain field
    typedef typename Traits::DomainType DomainType;
    //! Vector type for the range field
    typedef typename Traits::RangeType RangeType;
    //! Tensor type for the jacobian
    typedef typename Traits::JacobianRangeType JacobianRangeType;
    //! Intrinsic data type for the degrees of freedom (dof)
    typedef RangeFieldType DofType;

    //! Container class type for the dofs
    typedef typename Traits::DofStorageType DofStorageType;

    typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;
    //! Dimension of the range field
    enum { dimRange = DiscreteFunctionSpaceType::DimRange };

    //! type of codim 0 entity
    typedef typename DiscreteFunctionSpaceType :: GridType :: template
      Codim<0> :: Entity EntityType;

  public:
    //- Public methods
    //- Constructors and destructors
    
    //! Constructor
    AdaptiveLocalFunction(const DiscreteFunctionSpaceType& spc,
                          DofStorageType& dofVec);
    
    //! Copy constructor
    AdaptiveLocalFunction(const ThisType& other);

    //! Destructor
    ~AdaptiveLocalFunction();

    //- Operators
    //! Random access operator
    inline
    DofType& operator[] (const int num);

     //! Cosnt random access operator
    inline
    const DofType& operator[] (const int num) const;

    //- Methods

    //! Number of dofs on this element
    inline
    int numDofs() const;

    //! Evaluation of the discrete function
    inline
    void evaluate(const DomainType& x, 
                  RangeType & ret) const;

    //! Evaluation of the discrete function
    template <class QuadratureType>
    inline
    void evaluate(const QuadratureType& quad,
                  const int quadPoint,
                  RangeType& ret) const;

    //! Jacobian of the discrete function
    inline
    void jacobian(const DomainType& x, 
            		  JacobianRangeType& ret) const; 

    //! Jacobian of the discrete function
    template <class QuadratureType>
    inline
    void jacobian(const QuadratureType& quad,
                  const int quadPoint,
                  JacobianRangeType& ret) const;

    #if OLDFEM
    //! Evaluation of the discrete function
    inline
    void evaluate(EntityType& en, 
                  const DomainType& x, 
                  RangeType & ret) const {evaluate(en.gemoetry().local(x),ret);}
    //! Evaluation of the discrete function
    template <class EntityType>
    inline
    void evaluateLocal(EntityType& en, 
                       const DomainType& x, 
                       RangeType & ret) const {evaluate(x,ret);}
    //! Evaluation of the discrete function
    template <class QuadratureType>
    inline
    void evaluate(EntityType& en, 
                  QuadratureType& quad,
                  int quadPoint,
                  RangeType& ret) const {evaluate(quad,quadPoint,ret);}
    //! Jacobian of the discrete function
    inline
    void jacobianLocal(EntityType& en, 
                       const DomainType& x, 
                       JacobianRangeType& ret) const; 

    //! Jacobian of the discrete function
    inline
    void jacobian(EntityType& en, 
		  const DomainType& x, 
		  JacobianRangeType& ret) const; 
    //! Jacobian of the discrete function
    template <class QuadratureType>
    inline
    void jacobian(EntityType& en,
                  QuadratureType& quad,
                  int quadPoint,
                  JacobianRangeType& ret) const;

    #endif
    //! get the base function set
    const BaseFunctionSetType& baseFunctionSet() const;
    const BaseFunctionSetType& getBaseFunctionSet() const DUNE_DEPRECATED {return baseFunctionSet();}
  private:
    //- Forbidden methods
    //! assignment operator
    ThisType& operator=(const ThisType& other);

    //! init local function 
    inline
    void init(const EntityType& en);

  private:
    // return reference to actual entity
    const EntityType& en() const;

    //- Data members
    const DiscreteFunctionSpaceType& spc_;
    DofStorageType& dofVec_;
    
    // vector holding pointer to local dofs 
    mutable Array<RangeFieldType*> values_;
 
    mutable RangeType tmp_;
    mutable JacobianRangeType tmpGrad_;

    mutable bool init_;
    const bool multipleGeometryTypes_;

    mutable const BaseFunctionSetType* baseSet_;

    mutable const EntityType* en_;

    mutable GeometryType geoType_;
  }; // end class AdaptiveLocalFunction

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

  public:
    //- Typedefs and enums
    typedef MyTraits Traits;
    typedef Imp ImplementationType;
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

    typedef typename Traits::LocalFunctionImp LocalFunctionImp;
    typedef typename Traits::LocalFunctionType LocalFunctionType;
    typedef typename Traits::DiscreteFunctionType DiscreteFunctionType;
    
    typedef typename Traits::DofType DofType;
    typedef typename Traits::RangeFieldType RangeFieldType;
    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::MapperType MapperType;
    
    typedef typename Traits::DofStorageType DofStorageType;
     
    typedef typename Traits::DofIteratorType DofIteratorType;
    typedef typename Traits::ConstDofIteratorType ConstDofIteratorType;

    //- Additional typedefs
    typedef SubSpace<DiscreteFunctionSpaceType> SubSpaceType;
    typedef AdaptiveDiscreteFunction<SubSpaceType> SubDiscreteFunctionType;
   
  public:
    //- Public methods
    //! Constructor
    AdaptiveDiscreteFunction(std::string name,
                             const DiscreteFunctionSpaceType& spc) :
      BaseType(spc),
      Imp(name, spc),
      subSpaces_(0)
    {}

    //! Constructor
    template <class VectorPointerType>
    AdaptiveDiscreteFunction(std::string name,
                             const DiscreteFunctionSpaceType& spc,
                             VectorPointerType * vector) :
      BaseType(spc),
      Imp(name, spc, vector)
    {}
    
    //! Constructor
    AdaptiveDiscreteFunction(std::string name,
                             const DiscreteFunctionSpaceType& spc,
                             DofStorageType& dofVec) :
      BaseType(spc),
      Imp(name, spc, dofVec),
      subSpaces_(0)
    {}

    //! Copy constructor
    AdaptiveDiscreteFunction(const MyType& other) :
      BaseType(other.space()),
      Imp(other),
      subSpaces_(0) // Don't copy them
    {}
    
    ~AdaptiveDiscreteFunction();

    using Imp::name;
    using Imp::size;
    using Imp::dbegin;
    using Imp::dend;
    using Imp::newLocalFunctionObject;
    //! return empty local function 
    LocalFunctionType newLocalFunction () DUNE_DEPRECATED { return LocalFunctionType(*this); }
    //! return local function for given entity
    template <class EntityType> 
    LocalFunctionType localFunction (const EntityType &en) { return LocalFunctionType(en,*this); }
    using Imp::localFunction;
    using Imp::write_xdr;
    using Imp::read_xdr;
    using Imp::write_ascii;
    using Imp::read_ascii;
    using Imp::write_pgm;
    using Imp::read_pgm;

    using Imp::leakPointer;
    
    //- Additional methods
    SubDiscreteFunctionType subFunction(int component);

    int numComponents() const { return N; }

  private:
    std::vector<SubSpaceType*> subSpaces_;

    const MyType& interface() const { return *this; }
  }; // end class AdaptiveDiscreteFunction (specialised for CombinedSpace)

  //- class AdaptiveLocalFunction (specialised)
  //! Specialised version of AdaptiveLocalFunction for CombinedSpace
  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  class AdaptiveLocalFunction<
    CombinedSpace<ContainedFunctionSpaceImp, N, p> >
      : public LocalFunctionDefault<
    CombinedSpace<ContainedFunctionSpaceImp, N, p>,
    AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> > >
  {
  public:
    //- Friends
    friend class AdaptiveFunctionImplementation<
      CombinedSpace<ContainedFunctionSpaceImp, N, p> >;
    friend class LocalFunctionWrapper<
      AdaptiveDiscreteFunction<
      CombinedSpace<ContainedFunctionSpaceImp, N, p> > >;
  public:
    //- Public typedefs and enums
    typedef CombinedSpace<
      ContainedFunctionSpaceImp, N, p> DiscreteFunctionSpaceType;
    typedef AdaptiveLocalFunction<
      DiscreteFunctionSpaceType > ThisType;
    typedef AdaptiveDiscreteFunctionTraits<
      DiscreteFunctionSpaceType > Traits;
    typedef typename DiscreteFunctionSpaceType::Traits SpaceTraits;

    enum { dimRange = DiscreteFunctionSpaceType::DimRange };
    
    typedef typename SpaceTraits::ContainedRangeType ContainedRangeType;
    typedef typename SpaceTraits::ContainedJacobianRangeType 
    ContainedJacobianRangeType;
    typedef typename SpaceTraits::BaseFunctionSetType BaseFunctionSetType;
       
    typedef typename Traits::RangeFieldType RangeFieldType;
    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::JacobianRangeType JacobianRangeType;
    typedef typename Traits::DofType DofType;
    
    typedef typename Traits::DiscreteFunctionType DiscreteFunctionType;
    typedef typename Traits::DofStorageType DofStorageType;
  public:
    //- Public methods
    //- Constructors and destructors
    
    //! Constructor
    AdaptiveLocalFunction(const DiscreteFunctionSpaceType& spc,
                          DofStorageType& dofVec);

    //! Copy constructor
    AdaptiveLocalFunction(const ThisType& other);
    
    //! Destructor
    ~AdaptiveLocalFunction();

    //- Operators
    //! Random access operator
    inline
    DofType& operator[] (int num);

     //! Cosnt random access operator
     inline
   const DofType& operator[] (int num) const;

    //! Number of degrees of freedom
    inline
    int numDofs() const;

    //! Evaluation
    inline
    void evaluate(const DomainType& x, 
                  RangeType & ret) const;
    #if OLDFEM
    //! Evaluation
    template <class EntityType>
    inline
    void evaluate(EntityType& en, 
                  const DomainType& x, 
                  RangeType & ret) const {evaluate(en.geometry().local(x),ret);}
    //! Evaluation
    template <class EntityType>
    inline
    void evaluateLocal(EntityType& en, 
                       const DomainType& x, 
                       RangeType & ret) const {evaluate(en,x,ret);}
    //! Evaluation
    template <class EntityType, class QuadratureType>
    inline
    void evaluate(EntityType& en,
                  QuadratureType& quad,
                  int quadPoint, 
                  RangeType & ret) const {evaluate(quad,quadPoint,ret);}
    #endif

    //! Evaluation
    template <class QuadratureType>
    inline
    void evaluate(QuadratureType& quad,
                  int quadPoint, 
                  RangeType & ret) const;
    #if OLDFEM
    //! Evaluation
    template <class EntityType>
    inline
    void jacobianLocal(EntityType& en, 
                       const DomainType& x, 
                       JacobianRangeType& ret) const;
    #endif
    //! Evaluation
    template <class EntityType>
    inline
    void jacobian(EntityType& en, 
		  const DomainType& x, 
		  JacobianRangeType& ret) const;
    //! Evaluation
    template <class EntityType, class QuadratureType>
    inline
    void jacobian(EntityType& en,
                  QuadratureType& quad,
                  int quadPoint,
                  JacobianRangeType& ret) const;

    //- Additional methods for specialisation
    //! Assign a vector of dofs
    inline
    void assign(int dofNum, const RangeType& dofs);
    
    //! Number of contained scalar base functions
    inline
    int numDifferentBaseFunctions() const;

  private:
    //- Private methods
    template <class EntityType>
    inline
    void init(const EntityType& en);

  private:
    //- Typedefs
    typedef typename FieldVector<DofType*, N>::size_type SizeType;

  private:
    //- Member data
    const DiscreteFunctionSpaceType& spc_;
    DofStorageType& dofVec_;
    
    mutable Array<FieldVector<DofType*, N> > values_;
 
    mutable ContainedRangeType cTmp_;
    mutable ContainedJacobianRangeType cTmpGradRef_;
    mutable ContainedJacobianRangeType cTmpGradReal_;
    mutable RangeType tmp_;

    mutable bool init_;
    mutable const BaseFunctionSetType* baseSet_;

    mutable GeometryType geoType_;
  }; // end class AdaptiveLocalFunction (specialised for CombinedSpace)
  
} // end namespace Dune

#include "adaptivefunction.cc"

#endif
