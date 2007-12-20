#ifndef DUNE_COMBINEDFUNCTION_HH
#define DUNE_COMBINEDFUNCTION_HH

//- System includes
#include <string>
#include <vector>

//- Dune includes
#include <dune/fem/space/common/dofstorage.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/combinedspace/combinedspace.hh>

#include <dune/fem/function/common/discretefunction.hh>

//- Local includes
#include <dune/fem/function/localfunction/standardlocalfunction.hh>

namespace Dune
{
  
  //- Forward declarations
  template <class ContainedDiscreteFunctionImp,int N>
  class CombinedDiscreteFunction;

  
  //- Class definitions
  //! Traits class for AdaptiveDiscreteFunction and 
  //! AdaptiveLocalFunction
  template< class ContainedDiscreteFunctionImp,int N>
  struct CombinedDiscreteFunctionTraits
  {
    typedef ContainedDiscreteFunctionImp 
    ContainedDiscreteFunctionType;
    typedef typename ContainedDiscreteFunctionType::
    DiscreteFunctionSpaceType ContainedDiscreteFunctionSpaceType;
    typedef CombinedSpace<ContainedDiscreteFunctionSpaceType,N,VariableBased>
    DiscreteFunctionSpaceType;
    typedef CombinedDiscreteFunction<ContainedDiscreteFunctionType,N> DiscreteFunctionType;
    typedef StandardLocalFunctionFactory
    <CombinedDiscreteFunctionTraits<ContainedDiscreteFunctionType,N> >
    LocalFunctionFactoryType;
    typedef LocalFunctionStack<LocalFunctionFactoryType> 
    LocalFunctionStorageType;
    typedef typename LocalFunctionStorageType :: 
    LocalFunctionType LocalFunctionType;
    typedef typename DiscreteFunctionSpaceType :: 
    DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType :: 
    RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType :: 
    DomainFieldType DomainFieldType;
    typedef typename DiscreteFunctionSpaceType :: 
    RangeFieldType RangeFieldType;
    typedef RangeFieldType DofType;
    typedef typename DiscreteFunctionSpaceType :: 
    JacobianRangeType JacobianRangeType;
    typedef typename DiscreteFunctionSpaceType :: 
    MapperType MapperType;
    typedef typename DiscreteFunctionSpaceType :: 
    GridType GridType;
    
    typedef typename ContainedDiscreteFunctionType::DofIteratorType DofIteratorType;
    typedef typename ContainedDiscreteFunctionType::ConstDofIteratorType ConstDofIteratorType;
  }; // end class CombinedDiscreteFunctionTraits

  //! @ingroup CombinedDFunction
  //! A class for combining N discrete function of the same
  //! type to a vector valued function
  template <class ContainedDiscreteFunctionImp,int N >
  class CombinedDiscreteFunction : 
      public DiscreteFunctionDefault<
    CombinedDiscreteFunctionTraits<ContainedDiscreteFunctionImp,N> >
  {
  public:
    //! Discrete function this discrete function belongs to
    typedef ContainedDiscreteFunctionImp 
    ContainedDiscreteFunctionType;
    typedef ContainedDiscreteFunctionType 
    SubDiscreteFunctionType;
    //! Traits class with all necessary type definitions
    typedef CombinedDiscreteFunctionTraits<ContainedDiscreteFunctionType,N> Traits;

  private:
    typedef CombinedDiscreteFunction< ContainedDiscreteFunctionType,N > ThisType;
    typedef DiscreteFunctionDefault< Traits > BaseType;

  public:
    using BaseType :: assign; // needs DofIterator!

    //! Grid implementation 
    typedef typename Traits::GridType GridType;
    
    //! Discrete function type (identical to this type, 
    //! needed as Barton-Nackman parameter
    typedef typename Traits::DiscreteFunctionType 
    DiscreteFunctionType;
    //! the combined discrete function type
    typedef typename Traits::DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;
    //! Contained discrete function space 
    typedef typename Traits::ContainedDiscreteFunctionSpaceType
    ContainedDiscreteFunctionSpaceType;
    typedef ContainedDiscreteFunctionSpaceType 
    SubDiscreteFunctionSpaceType; 
    //! Intrinsic type used for dofs (typically a float type)
    typedef typename Traits::DofType DofType;
    //! Intrinsic type used for range field (like DofType)
    typedef typename Traits::RangeFieldType RangeFieldType;
    //! Intrinsic type used for the domain field
    typedef typename Traits::DomainFieldType DomainFieldType;
    //! Vector type used for the range field
    typedef typename Traits::RangeType RangeType;
    //! Vector type used for the domain field
    typedef typename Traits::DomainType DomainType;
    //! Mapper type (from the space)
    typedef typename Traits::MapperType MapperType;
     
    //! Iterator over dof container
    typedef typename Traits::DofIteratorType DofIteratorType;
    //! Read-only iterator over dof container
    typedef typename Traits::ConstDofIteratorType 
    ConstDofIteratorType;

    typedef Mapping<DomainFieldType, RangeFieldType,
                    DomainType, RangeType> MappingType;
    typedef typename Traits :: LocalFunctionFactoryType 
    LocalFunctionFactoryType;

  public:
    //- Public methods
    //! Constructor 
    //! WARNING: here we have to use a const cast for the
    //! function space!
    CombinedDiscreteFunction(ContainedDiscreteFunctionType& func) 
      : BaseType( spc_, lfFactory_ ),
	spc_(const_cast<ContainedDiscreteFunctionSpaceType&>(func.space()).gridPart()),
	lfFactory_( *this )
    {
      for (int i=0;i<N;i++) {
	func_[i] = new ContainedDiscreteFunctionType(func);
      }
    }
    //! Copy constructor
    //! The copy constructor copies the dofs
    CombinedDiscreteFunction(const ThisType &other)
      : BaseType( other.space(), lfFactory_ ),
	spc_(other.spc_),
	lfFactory_( *this )
    {
      for (int i=0;i<N;i++) {
	func_[i] = new ContainedDiscreteFunctionType(other.subFunction(i));
      }
    }
    //! Destructor
    ~CombinedDiscreteFunction() {
      for (int i=0;i<N;i++) 
	delete func_[i];
    }

  private:
    ThisType &operator= ( const ThisType &other );
    
  public:
    /** \copydoc Dune::DiscreteFunctionInterface::name */
    const std :: string &name () const
    {
      return func_[0]->name() + std::string("<N>");
    }
    /** \copydoc Dune::DiscreteFunctionInterface::clear */
    inline void clear() {
      for (int i=0;i<N;i++)
	      func_[i]->clear();
    }
    /** \copydoc Dune::DiscreteFunctionInterface::assign(const DiscreteFunctionType &g) */
    inline void assign( const DiscreteFunctionType &g )
    {
      for (int i=0;i<N;i++)
	      func_[i]->assign(g.subFunction(i));
    }
    /** \copydoc Dune::DiscreteFunctionDefault::size */ 
    inline int size() const
    {
      return func_[0]->size()*N;
    }
    /** \copydoc Dune::DiscreteFunctionDefault::operator+=
     */ 
    inline ThisType &operator += ( const ThisType &g )
    {
      for (int i=0;i<N;i++)
	      *func_[i] += g.subFunction(i);
      return *this;
    }
    /** \copydoc Dune::DiscreteFunctionDefault::operator-=
     */ 
    inline BaseType &operator-= ( const ThisType &g )
    {
      for (int i=0;i<N;i++)
	*func_[i] -= g.subFunction(i);
     return *this;
    }
    /** \copydoc Dune::DiscreteFunctionDefault::operator*=(const RangeFieldType &scalar) */    
    DiscreteFunctionType& operator *= (const RangeFieldType &scalar) {
      for (int i=0;i<N;i++)
	*func_[i] *= scalar;
      return *this;
    }
    /** \copydoc Dune::DiscreteFunctionDefault::operator*=(const RangeFieldType &scalar) */    
    DiscreteFunctionType& operator /= (const RangeFieldType &scalar) {
      for (int i=0;i<N;i++)
	*func_[i] /= scalar;
      return *this;
    }
    /** \copydoc Dune::DiscreteFunctionDefault::addScaled
     */
    inline void addScaled( const ThisType &g,
                           const RangeFieldType &s )
    {
      for (int i=0;i<N;i++)
	func_[i]->addScaled(g.subFunction(i),s);
    }
    /** \copydoc Dune::DiscreteFunctionDefault::scalarProductDofs(const DiscreteFunctionType &g) */
    RangeFieldType scalarProductDofs ( const DiscreteFunctionType &g ) const
    {
      double ret=func_[0]->scalarProductDofs(g.subFunction(0));
      for (int i=1;i<N;i++)
	ret += func_[i]->scalarProductDofs(g.subFunction(i));
    }

    /** \copydoc Dune::DiscreteFunctionDefault::read(InStreamInterface< StreamTraits>& in) */
    template< class StreamTraits >
    inline void read ( InStreamInterface< StreamTraits >& in)
    {
      for (int i=0;i<N;i++)
	func_[i]->read(in);
    }
    /** \copydoc Dune::DiscreteFunctionDefault::write(OutStreamInterface< StreamTraits>& out) */
    template< class StreamTraits >
    inline void write ( OutStreamInterface< StreamTraits>& out) const
    {
      for (int i=0;i<N;i++)
	func_[i]->write(out);
    }
    /** \copydoc Dune::DiscreteFunctionDefault::print(std :: ostream &out) */
    inline void print( std :: ostream &out ) const {
      for (int i=0;i<N;i++)
	func_[i]->print(out);
    }
    /** \copydoc Dune::DiscreteFunctionDefault::dofValid */
    inline bool dofsValid () const {
      bool ret = func_[0]->dofsValid();
      for (int i=1;i<N;i++)
	ret |= func_[i]->dofsValid();
      return ret;
    }

    /** \copydoc Dune::DiscreteFunctionDefault::dofs(unsigned int index) const */
    inline const RangeFieldType &dof(unsigned int index) const
    {
      /*
      int variable = index % N;
      int point = index / N;
      */
      int variable = index / func_[0]->size();
      int point    = index % func_[0]->size();
      return func_[variable]->dof(point);
    }
    /** \copydoc Dune::DiscreteFunctionDefault::dofs(unsigned int index) */
    inline RangeFieldType &dof ( unsigned int index )
    {
      /*
      int variable = index % N;
      int point = index / N;
      */
      int variable = index / func_[0]->size();
      int point    = index % func_[0]->size();
      return func_[variable]->dof(point);
    }
    
    /** \copydoc Dune::DiscreteFunctionDefault::dbegin */
    inline ConstDofIteratorType dbegin () const
    {
      return func_[0]->dbegin();
    }
    /** \copydoc Dune::DiscreteFunctionDefault::dend */
    inline ConstDofIteratorType dend () const
    {
      return func_[0]->dend();
    }
    /** \copydoc Dune::DiscreteFunctionDefault::dbegin */
    inline DofIteratorType dbegin ()
    {
      return func_[0]->dbegin();
    }
    /** \copydoc Dune::DiscreteFunctionDefault::dend */
    inline DofIteratorType dend ()
    {
      return func_[0]->dend();
    }
    
    inline ContainedDiscreteFunctionType& subFunction(int i) {
      return *(func_[i]);
    }
    inline ContainedDiscreteFunctionSpaceType& subSpace() {
      return spc_.containedSpace();
    }
    //- Forbidden members
  private:
    typedef ThisType MyType;
    const MyType& interface() const { return *this; }
    
    DiscreteFunctionSpaceType spc_;
    const LocalFunctionFactoryType lfFactory_;
    ContainedDiscreteFunctionType* func_[N];
  }; // end class AdaptiveDiscreteFunction
  
} // end namespace Dune

#endif
