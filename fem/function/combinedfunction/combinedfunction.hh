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
  template< class ContainedDiscreteFunction, int N >
  class CombinedDiscreteFunction;
  
  template< class ContainedDiscreteFunction, int N >
  class CombinedDiscreteFunctionDofIterator;

  
  //- Class definitions
  //! Traits class for AdaptiveDiscreteFunction and 
  //! AdaptiveLocalFunction
  template< class ContainedDiscreteFunction, int N >
  struct CombinedDiscreteFunctionTraits
  {
    typedef ContainedDiscreteFunction ContainedDiscreteFunctionType;
    typedef CombinedDiscreteFunction< ContainedDiscreteFunctionType, N >
      DiscreteFunctionType;

    typedef typename ContainedDiscreteFunctionType :: DiscreteFunctionSpaceType
      ContainedDiscreteFunctionSpaceType;
    typedef CombinedSpace< ContainedDiscreteFunctionSpaceType, N, VariableBased >
      DiscreteFunctionSpaceType;

    typedef StandardLocalFunctionFactory
      < CombinedDiscreteFunctionTraits< ContainedDiscreteFunctionType, N > >
      LocalFunctionFactoryType;
    typedef LocalFunctionStack< LocalFunctionFactoryType >
      LocalFunctionStorageType;
    typedef typename LocalFunctionStorageType :: LocalFunctionType
      LocalFunctionType;

    typedef typename DiscreteFunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType;
    typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
    
    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType JacobianRangeType;

    typedef RangeFieldType DofType;
    typedef typename DiscreteFunctionSpaceType :: MapperType MapperType;
    typedef typename DiscreteFunctionSpaceType :: GridType GridType;
    
    typedef CombinedDiscreteFunctionDofIterator< ContainedDiscreteFunctionType, N >
      DofIteratorType;
    typedef ConstDofIteratorDefault< DofIteratorType > ConstDofIteratorType;

    typedef typename ContainedDiscreteFunctionType :: DofBlockType DofBlockType;
    typedef typename ContainedDiscreteFunctionType :: ConstDofBlockType
      ConstDofBlockType;
    typedef typename ContainedDiscreteFunctionType :: DofBlockPtrType
      DofBlockPtrType;
    typedef typename ContainedDiscreteFunctionType :: ConstDofBlockPtrType
      ConstDofBlockPtrType;
  }; // end class CombinedDiscreteFunctionTraits



  //! @ingroup CombinedDFunction
  //! A class for combining N discrete function of the same
  //! type to a vector valued function
  template <class ContainedDiscreteFunctionImp,int N >
  class CombinedDiscreteFunction
  : public DiscreteFunctionDefault
    < CombinedDiscreteFunctionTraits< ContainedDiscreteFunctionImp, N > >
  {
  public:
    //! Discrete function this discrete function belongs to
    typedef ContainedDiscreteFunctionImp 
    ContainedDiscreteFunctionType;
    typedef ContainedDiscreteFunctionType 
    SubDiscreteFunctionType;
    //! Traits class with all necessary type definitions
    typedef CombinedDiscreteFunctionTraits<ContainedDiscreteFunctionType,N> Traits;

    //! type of the interface of this discrete function
    typedef DiscreteFunctionInterface< Traits > DiscreteFunctionInterfaceType;

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

    typedef typename Traits :: DofBlockPtrType DofBlockPtrType;
    typedef typename Traits :: ConstDofBlockPtrType ConstDofBlockPtrType;

  public:
    //- Public methods
    //! Constructor 
    //! WARNING: here we have to use a const cast for the
    //! function space!
    CombinedDiscreteFunction(ContainedDiscreteFunctionType& func) 
      : BaseType( "", spc_, lfFactory_ ),
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
      : BaseType( "", other.space(), lfFactory_ ),
	spc_(other.spc_),
	lfFactory_( *this )
    {
      for (int i=0;i<N;i++) {
	func_[i] = new 
	  ContainedDiscreteFunctionType(other.subFunction(i));
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
    /** \copydoc Dune::DiscreteFunctionInterface::size() const */ 
    inline int size() const
    {
      return func_[0]->size()*N;
    }

    /** \copydoc Dune::DiscreteFunctionInterface::operator+=(const DiscreteFunctionType &g) */
    inline ThisType &operator+= ( const ThisType &g )
    {
      for (int i=0;i<N;i++)
	      *func_[i] += g.subFunction(i);
      return *this;
    }

    /** \copydoc Dune::DiscreteFunctionInterface::operator-=
     */ 
    using BaseType::operator-=;
    inline BaseType &operator-= ( const ThisType &g )
    {
      // std::cout << "     special operator -= in combineddf" 
      //	<< std::endl;
      for (int i=0;i<N;i++)
	*func_[i] -= g.subFunction(i);
     return *this;
    }
    /** \copydoc Dune::DiscreteFunctionInterface::operator*=(const RangeFieldType &scalar) */    
    DiscreteFunctionType& operator *= (const RangeFieldType &scalar) {
      for (int i=0;i<N;i++)
	*func_[i] *= scalar;
      return *this;
    }
    /** \copydoc Dune::DiscreteFunctionInterface::operator*=(const RangeFieldType &scalar) */    
    DiscreteFunctionType& operator /= (const RangeFieldType &scalar) {
      for (int i=0;i<N;i++)
	*func_[i] /= scalar;
      return *this;
    }
    /** \copydoc Dune::DiscreteFunctionInterface::addScaled
     */
    inline void addScaled( const ThisType &g,
                           const RangeFieldType &s )
    {
      for (int i=0;i<N;i++)
	func_[i]->addScaled(g.subFunction(i),s);
    }
    /** \copydoc Dune::DiscreteFunctionInterface::scalarProductDofs(const DiscreteFunctionInterfaceType &other) const */
    RangeFieldType scalarProductDofs ( const DiscreteFunctionInterfaceType &other ) const
    {
      RangeFieldType ret( 0 );
      for( int i = 0; i < N; ++i )
	ret += func_[ i ]->scalarProductDofs( other.asImp().subFunction( i ) );
      return ret;
    }

    /** \copydoc Dune::DiscreteFunctionInterface::read */
    template< class StreamTraits >
    inline void read ( InStreamInterface< StreamTraits >& in)
    {
      for (int i=0;i<N;i++)
	func_[i]->read(in);
    }
    /** \copydoc Dune::DiscreteFunctionInterface::write */
    template< class StreamTraits >
    inline void write ( OutStreamInterface< StreamTraits >& out) const
    {
      for (int i=0;i<N;i++)
	func_[i]->write(out);
    }
    /** \copydoc Dune::DiscreteFunctionInterface::print(std::ostream &out) const */
    inline void print( std :: ostream &out ) const {
      for (int i=0;i<N;i++)
	func_[i]->print(out);
    }
    /** \copydoc Dune::DiscreteFunctionInterface::dofsValid() const */
    inline bool dofsValid () const {
      bool ret = func_[0]->dofsValid();
      for (int i=1;i<N;i++)
	ret |= func_[i]->dofsValid();
      return ret;
    }

    inline ConstDofBlockPtrType block ( unsigned int index ) const
    {
      // This is wrong with the current implementation of CombinedSpace
      const int containedSize = func_[ 0 ]->space().blockMapper().size();
      const int component = index / containedSize;
      const int containedIndex = index % containedSize;
      return func_[ component ]->block( containedIndex );
    }
    
    inline DofBlockPtrType block ( unsigned int index )
    {
      // This is wrong with the current implementation of CombinedSpace
      const int containedSize = func_[ 0 ]->space().blockMapper().size();
      const int component = index / containedSize;
      const int containedIndex = index % containedSize;
      return func_[ component ]->block( containedIndex );
    }

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
    
    /** \copydoc Dune::DiscreteFunctionInterface::dbegin() const */
    inline ConstDofIteratorType dbegin () const
    {
      return ConstDofIteratorType(DofIteratorType(*this));
    }
    /** \copydoc Dune::DiscreteFunctionInterface::dend() const */
    inline ConstDofIteratorType dend () const
    {
      return ConstDofIteratorType(DofIteratorType(false,*this));
    }
    /** \copydoc Dune::DiscreteFunctionInterface::dbegin() */
    inline DofIteratorType dbegin ()
    {
      return DofIteratorType(*this);
    }
    /** \copydoc Dune::DiscreteFunctionInterface::dend() */
    inline DofIteratorType dend ()
    {
      return DofIteratorType(false,*this);
    }
    
    inline ContainedDiscreteFunctionType& subFunction(int i) {
      return *(func_[i]);
    }
    inline const ContainedDiscreteFunctionType& 
    subFunction(int i) const {
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
    friend class CombinedDiscreteFunctionDofIterator<ContainedDiscreteFunctionType,N>;
  }; // end class AdaptiveDiscreteFunction
  
  /** \brief Iterator over an array of dofs 
      \todo Please doc me!
  */
  template <class ContainedDiscreteFunctionImp,int N>
  class CombinedDiscreteFunctionDofIterator
    : public DofIteratorDefault < 
    typename ContainedDiscreteFunctionImp::DofType , 
    CombinedDiscreteFunctionDofIterator<ContainedDiscreteFunctionImp,N> >
  {
  public:
    typedef CombinedDiscreteFunctionDofIterator<ContainedDiscreteFunctionImp,N> ThisType;
    typedef CombinedDiscreteFunctionTraits<ContainedDiscreteFunctionImp,N> Traits;
    typedef typename Traits::DiscreteFunctionType DiscreteFunctionType;
    typedef typename Traits::ContainedDiscreteFunctionType ContainedDiscreteFunctionType;
    typedef typename ContainedDiscreteFunctionType::DofIteratorType ContainedDofIteratorType;
    typedef typename ContainedDiscreteFunctionType::ConstDofIteratorType ContainedConstDofIteratorType;
    typedef typename Traits::DofType DofType;
    
    //! End constructor
    CombinedDiscreteFunctionDofIterator
    (bool end,const DiscreteFunctionType& df) :
      df_(const_cast<DiscreteFunctionType&>(df)),
      comp_(N-1),
      iter_(df.func_[N-1]->dend()),
      endIter_(df.func_[N-1]->dend())
    {}
    //! Constructor (const)
    CombinedDiscreteFunctionDofIterator
    (const DiscreteFunctionType& df) :
      df_(const_cast<DiscreteFunctionType&>(df)),
      comp_(0),
      iter_(df.func_[0]->dbegin()),
      endIter_(df.func_[0]->dend())
    {}
    //! End constructor
    CombinedDiscreteFunctionDofIterator
    (bool end,DiscreteFunctionType& df) :
      df_(df),
      comp_(N-1),
      iter_(df.func_[N-1]->dend()),
      endIter_(df.func_[N-1]->dend())
    {}
    //! Constructor
    CombinedDiscreteFunctionDofIterator
    (DiscreteFunctionType& df) :
      df_(df),
      comp_(0),
      iter_(df.func_[0]->dbegin()),
      endIter_(df.func_[0]->dend())
    {}
    //! Copy Constructor
    CombinedDiscreteFunctionDofIterator(const ThisType& other):
      df_(other.df_),
      comp_(other.comp_),
      iter_(other.iter_),
      endIter_(other.endIter_)
    {}

    //! Assignment operator
    ThisType& operator=(const ThisType& other) {
      df_ = other.df_;
      comp_ = other.comp_;
      iter_ = other.iter_;
      endIter_ = other.endIter_;
      return *this;
    }
    //! return dof
    DofType& operator *() {
      return *iter_;
    }
    //! return dof read only 
    const DofType& operator * () const {
      return *iter_;
    }
    //! go to next dof
    ThisType& operator++ () {
      ++iter_;
      if (iter_==endIter_ && comp_<N-1) {
	++comp_;
	iter_ = df_.func_[comp_]->dbegin();
	endIter_ = df_.func_[comp_]->dend();
      } 
      return *this;
    }
  
    //! compare
    bool operator == (const ThisType & I ) const
    {
      return (comp_ == I.comp_) && (iter_ == I.iter_);
    }
    //! compare 
    bool operator != (const ThisType & I ) const
    {
      return !((*this) == I);
    }
    
private:
    DiscreteFunctionType& df_;
    //! index 
    mutable int comp_;
    mutable ContainedDofIteratorType iter_,endIter_;
}; // end DofIteratorCombinedDiscreteFunction 


} // end namespace Dune

#endif
