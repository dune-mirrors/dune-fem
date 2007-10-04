#ifndef DUNE_ADAPTIVEFUNCTION_HH
#define DUNE_ADAPTIVEFUNCTION_HH

//- System includes
#include <string>
#include <vector>

//- Dune includes
#include <dune/fem/space/common/dofstorage.hh>
#include <dune/fem/space/common/dofmanager.hh>

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/common/localfunction.hh>

//- Local includes
#include "adaptiveimp.hh"

namespace Dune
{
  
  //- Forward declarations
  template <class DiscreteFunctionSpaceImp>
  class AdaptiveDiscreteFunction;
  template <class DiscreteFunctionSpaceImp>
  class AdaptiveLocalFunction;

  //- Forward declarations of Combined Space 
  template <class, int , DofStoragePolicy > 
  class CombinedSpace;
  template <class CombinedSpaceImp>
  class SubSpace;



  template< class DiscreteFunctionSpaceImp >
  class AdaptiveLocalFunctionFactory
  {
  public:
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

  private:
    typedef AdaptiveLocalFunctionFactory< DiscreteFunctionSpaceType > ThisType;

    friend class AdaptiveDiscreteFunction< DiscreteFunctionSpaceType >;

  public:
    typedef AdaptiveLocalFunction< DiscreteFunctionSpaceType > ObjectType;

    typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType >
      DiscreteFunctionType;

  protected:
    DiscreteFunctionType &discreteFunction_;

  protected:
    inline explicit AdaptiveLocalFunctionFactory ( DiscreteFunctionType &df )
    : discreteFunction_( df )
    {
    }

  public:
    ObjectType *newObject () const
    {
      return new ObjectType( discreteFunction_.space(),
                             discreteFunction_.dofVec_ );
    }
  };


  
  //- Class definitions
  //! Traits class for AdaptiveDiscreteFunction and AdaptiveLocalFunction
  template< class DiscreteFunctionSpaceImp >
  struct AdaptiveDiscreteFunctionTraits
  {
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;
 
    typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

    typedef AdaptiveLocalFunctionFactory< DiscreteFunctionSpaceType >
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

    friend class AdaptiveLocalFunctionFactory< DiscreteFunctionSpaceType >;
  
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
    /** \copydoc Dune::DiscreteFunctionDefault::assign
     */
    inline void assign( const ThisType &g )
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

    //! GridType
    typedef typename Traits::GridType GridType;
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

  private:  
    typedef typename GridType :: ctype ctype;
    enum { dim = GridType :: dimension };
    typedef FieldMatrix<ctype,dim,dim> JacobianInverseType;
  public:
    //- Public methods
    //- Constructors and destructors
    
    //! Constructor
    inline
    AdaptiveLocalFunction(const DiscreteFunctionSpaceType& spc,
                          DofStorageType& dofVec);
    
    //! Copy constructor
    inline
    AdaptiveLocalFunction(const ThisType& other);

    //! Destructor
    ~AdaptiveLocalFunction();

    //- Operators
    /** \copydoc Dune::LocalFunctionDefault::operator[](const int num)
     */
    inline DofType &operator[] ( const int num );

    /** \copydoc Dune::LocalFunctionDefault::operator[](const int num) const
     */
    inline const DofType &operator[] ( const int num ) const;

    //- Methods

    /** \copydoc LocalFunctionInterface::numDofs
     */
    inline int numDofs() const;

    /** \copydoc Dune::LocalFunctionInterface::evaluate(const DomainType &x,RangeType &ret) const
     */
    inline void evaluate( const DomainType &x,
                          RangeType &ret ) const;

    /** \copydoc Dune::LocalFunctionInterface::evaluate(const QuadratureType &quadrature,const int quadPoint,RangeType &ret) const
     */
    template< class QuadratureType >
    inline void evaluate ( const QuadratureType &quadrature,
                           const int quadPoint,
                           RangeType &ret ) const;

    /** \copydoc Dune::LocalFunctionInterface::jacobian(const DomainType &x,JacobianRangeType &ret) const
     */
    inline void jacobian( const DomainType &x,
                          JacobianRangeType &ret ) const; 

    /** \copydoc Dune::LocalFunctionInterface::jacobian( const QuadratureType &quadrature,const int quadPoint,JacobianRangeType &ret) const
     */
    template< class QuadratureType >
    inline void jacobian( const QuadratureType &quadrature,
                          const int quadPoint,
                          JacobianRangeType &ret ) const;
    
    /** \brief evaluate jacobian of the discrete function on point x
        \param[in] entity en
        \param[in] domain type x
        \param[out] jacobian range type ret
    */
    inline
    void jacobian(EntityType& en, 
      const DomainType& x, 
      JacobianRangeType& ret) const; 
    
    /** \brief evaluate jacobian of the discrete function on quadrature point quadPoint
        \param[in] entity en
        \param[in] quadrature type quad
        \param[in] constant quadrature point quadPoint 
        \param[out] jacobian range type
    */
    template <class QuadratureType>
    inline
    void jacobian(EntityType& en,
                  QuadratureType& quad,
                  int quadPoint,
                  JacobianRangeType& ret) const;

    /** \copydoc Dune::LocalFunctionDefault::baseFunctionSet */
    const BaseFunctionSetType& baseFunctionSet() const;

    /** \brief @copydoc LocalFunctionDefault::axpy */
    template <class QuadratureType>
    inline void axpy(const QuadratureType&, const int qp, const RangeType& factor);

    /** \brief @copydoc LocalFunctionDefault::axpy */
    template <class QuadratureType>
    inline void axpy(const QuadratureType&, const int qp, const JacobianRangeType& factor);

    /** \brief @copydoc LocalFunctionDefault::axpy */
    template <class QuadratureType>
    inline void axpy(const QuadratureType&, const int qp, const RangeType& factor1, const JacobianRangeType& factor2);

  private:
    //- Forbidden methods
    //! assignment operator
    ThisType& operator=(const ThisType& other);

  public:
    //! init local function 
    inline void init ( const EntityType &entity );

  private:
    inline
    void rightMultiply(const JacobianRangeType& factor, 
           const JacobianInverseType& jInv, 
                       JacobianRangeType& result) const;
    // return reference to actual entity
    const EntityType& en() const;

    //- Data members
    const DiscreteFunctionSpaceType& spc_;
  protected:
    DofStorageType& dofVec_;
    
    // vector holding pointer to local dofs 
    mutable MutableArray<RangeFieldType*> values_;
    
      //! number of local dofs 
    int numDofs_;
 
    mutable RangeType tmp_;
    mutable JacobianRangeType tmpGrad_;
    mutable JacobianRangeType factorInv_;

    mutable bool needCheckGeometry_;

    // base function set 
    mutable BaseFunctionSetType baseSet_;

    // actual entity
    mutable const EntityType* en_;
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
    typedef SubSpace<DiscreteFunctionSpaceType> SubSpaceType;
    typedef AdaptiveDiscreteFunction<SubSpaceType> SubDiscreteFunctionType;
   
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
    {}

    //! Constructor
    template <class VectorPointerType>
    AdaptiveDiscreteFunction(std::string name,
                             const DiscreteFunctionSpaceType& spc,
                             VectorPointerType * vector)
    : BaseType( spc, lfFactory_ ),
      Imp( name, spc , vector ),
      lfFactory_( *this )
    {}
    
    //! Constructor
    AdaptiveDiscreteFunction(std::string name,
                             const DiscreteFunctionSpaceType& spc,
                             DofStorageType& dofVec) 
    : BaseType( spc, lfFactory_ ),
      Imp( name, spc , dofVec ),
      lfFactory_( *this )
    {}

    //! Copy constructor
    AdaptiveDiscreteFunction(const MyType& other)
    : BaseType( other.space(), lfFactory_ ),
      Imp(other),
      lfFactory_( *this )
    {}
    
    ~AdaptiveDiscreteFunction();
    

  private:
    friend class AdaptiveLocalFunctionFactory< DiscreteFunctionSpaceType >;

    // local function factory 
    const LocalFunctionFactoryType lfFactory_;

    ThisType &operator= ( const ThisType &other );
 
  public:
    /** \brief @copydoc DiscreteFunctionDefault::assign */
    void assign(const ThisType& g)
    {
      Imp::assignFunction(g);
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
    SubDiscreteFunctionType subFunction(int component);

    int numComponents() const { return N; }

#if 0
  public:
    friend class DiscreteFunctionInterface<MyTraits>;

  protected:
    using Imp::newObject;
#endif
    
  private:
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

  public:
    //- Public typedefs and enums
    typedef CombinedSpace<
      ContainedFunctionSpaceImp, N, p> DiscreteFunctionSpaceType;
    typedef AdaptiveLocalFunction<
      DiscreteFunctionSpaceType > ThisType;
    typedef AdaptiveDiscreteFunctionTraits<
      DiscreteFunctionSpaceType > Traits;
    typedef typename DiscreteFunctionSpaceType::Traits SpaceTraits;

    typedef typename Traits::GridType GridType;
    
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

    typedef typename DiscreteFunctionSpaceType::GridType::
    template Codim<0>:: Entity EntityType;

  private:
    //- Typedefs
    typedef typename FieldVector<DofType*, N>::size_type SizeType;

    // store pointer to dofs belonging to local function 
    template <class Dummy, DofStoragePolicy policy>
    struct MapLocalDofs
    {
      template <class LocalDofVector>
      static void map(const DiscreteFunctionSpaceType& spc,
                      const EntityType& en, 
                      const int i, 
                      DofStorageType& dofVec,
                      LocalDofVector& values)
      {
        for (SizeType j = 0; j < N; ++j )
        {
          values[i][j] = &(dofVec[spc.mapToGlobal(en, i*N+j)]);
        } // end for j
      }
    };
    
    // specialisation for point based storage 
    template <class Dummy>
    struct MapLocalDofs<Dummy,PointBased>
    {
      template <class LocalDofVector>
      static void map(const DiscreteFunctionSpaceType& spc,
                      const EntityType& en, 
                      const int i, 
                      DofStorageType& dofVec,
                      LocalDofVector& values)
      {
        const int idx = spc.mapToGlobal(en, i*N );
        for (SizeType j = 0; j < N; ++j )
        {
          // only works for point base mappings 
          assert( (idx+ (int) j) == spc.mapToGlobal(en, i*N+j) );
          values[i][j] = &(dofVec[ idx+j ]);
        } // end for j
      }
    };
    
  private:  
    typedef typename GridType :: ctype ctype;
    enum { dim = GridType :: dimension };
    typedef FieldMatrix<ctype,dim,dim> JacobianInverseType;
  public:
    //- Public methods
    //- Constructors and destructors
    
    //! Constructor
    inline
    AdaptiveLocalFunction(const DiscreteFunctionSpaceType& spc,
                          DofStorageType& dofVec);

    //! Copy constructor
    inline
    AdaptiveLocalFunction(const ThisType& other);
    
    //! Destructor
    ~AdaptiveLocalFunction();

    //- Operators
    /** \copydoc Dune::LocalFunctionInterface::operator[](const int num) const */
    inline const RangeFieldType &operator[]( const int num ) const;

    /** \copydoc Dune::LocalFunctionInterface::operator[](const int num) */
    inline RangeFieldType &operator[]( const int num );

    /** \copydoc Dune::LocalFunctionInterface::numDofs */
    inline int numDofs() const;

    /** \copydoc Dune::LocalFunctionInterface::evaluate(const DomainType &x,RangeType &ret) const */
    inline void evaluate ( const DomainType &x,
                           RangeType &ret ) const;

    /** \copydoc Dune::LocalFunctionInterface::evaluate(const QuadratureType &quadrature,const int quadPoint,RangeType &ret) const */
    template< class QuadratureType >
    inline void evaluate ( const QuadratureType &quadrature,
                           const int quadPoint, 
                           RangeType &ret ) const;
    
    /** \brief evaluate jacobian of the discrete function on point x
        \param[in] entity en 
        \param[in] domain type x
        \param[out] jacobian range type ret
    */ 
    inline
    void jacobian(EntityType& en, 
      const DomainType& x, 
      JacobianRangeType& ret) const;
    
    
    /** \brief evaluate jacobian of the discrete function on quadrature point quadPoint
        \param[in] entity en
        \param[in] quadrature type quad
        \param[in] constant quadrature point quadPoint 
        \param[out] jacobian range type
    */
    template <class QuadratureType>
    inline
    void jacobian(EntityType& en,
                  QuadratureType& quad,
                  int quadPoint,
                  JacobianRangeType& ret) const;

    
    /** \copydoc Dune::LocalFunctionInterface::jacobian(const DomainType &x,JacobianRangeType &ret) const */
    inline void jacobian ( const DomainType &x,
                           JacobianRangeType &ret ) const;
        
    /** \copydoc Dune::LocalFunctionInterface::jacobian(const QuadratureType &quadrature,const int quadPoint,JacobianRangeType &ret) const */
    template< class QuadratureType >
    inline void jacobian ( const QuadratureType &quadrature,
                           const int quadPoint,
                           JacobianRangeType &ret ) const;

    //- Additional methods for specialisation
    /** \todo Please doc me!
     *
     *  \note Neither LocalFunctionInterface nor LocalFunctionDefault have this
     *        method.
     */
    inline
    void assign(int dofNum, const RangeType& dofs);
    
    /** \brief Number of contained scalar base functions */
    inline
    int numDifferentBaseFunctions() const;

    /** \copydoc Dune::LocalFunctionInterface::baseFunctionSet */
    inline const BaseFunctionSetType &baseFunctionSet() const;
   
     /** \brief @copydoc LocalFunctionDefault::axpy */
    template <class QuadratureType>
    inline void axpy(const QuadratureType&, const int qp, const RangeType& factor);


    /** \brief @copydoc LocalFunctionDefault::axpy */
    template <class QuadratureType>
    inline void axpy(const QuadratureType&, const int qp, const JacobianRangeType& factor);

   
     /** \brief @copydoc LocalFunctionDefault::axpy */
    template <class QuadratureType>
    inline void axpy(const QuadratureType&, const int qp, const RangeType& factor1, const JacobianRangeType& factor2);

  public:
    inline void init ( const EntityType &entity );

  private:
    //- Private methods
    // return reference to entity 
    inline
    const EntityType& en() const;

    // apply factor.rightmultiply(jInv) and strore in result 
    inline
    void rightMultiply(const JacobianRangeType& factor, const JacobianInverseType& jInv, 
                       JacobianRangeType& result) const;
  private:
    //- Member data
    const DiscreteFunctionSpaceType& spc_;

  protected:
    DofStorageType& dofVec_;
    
  private:
    // local dofs 
    typedef DofType* DofPointerType; 
    mutable MutableArray< DofPointerType[N] > values_;

    // number of local dofs 
    int numDofs_;
 
    // temporary variables 
    mutable ContainedRangeType cTmp_;
    mutable ContainedJacobianRangeType cTmpGradRef_;
    mutable ContainedJacobianRangeType cTmpGradReal_;
    mutable JacobianRangeType factorInv_;
    mutable RangeType tmp_;

    mutable bool needCheckGeometry_;
    mutable BaseFunctionSetType baseSet_;

    mutable const EntityType* en_;
  }; // end class AdaptiveLocalFunction (specialised for CombinedSpace)
  
} // end namespace Dune

#include "adaptivefunction.cc"

#endif
