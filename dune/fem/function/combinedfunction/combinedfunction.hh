#ifndef DUNE_FEM_FUNCTION_COMBINEDFUNCTION_COMBINEDFUNCTION_HH
#define DUNE_FEM_FUNCTION_COMBINEDFUNCTION_COMBINEDFUNCTION_HH

#include <algorithm>
#include <iostream>
#include <string>

#include <dune/fem/common/stackallocator.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/localfunction/mutable.hh>
#include <dune/fem/space/combinedspace/combinedspace.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/common/dofstorage.hh>
#include <dune/fem/storage/referencevector.hh>

namespace Dune
{

  namespace Fem
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
    struct DiscreteFunctionTraits< CombinedDiscreteFunction< ContainedDiscreteFunction, N > >
    {
      typedef ContainedDiscreteFunction ContainedDiscreteFunctionType;
      typedef CombinedDiscreteFunction< ContainedDiscreteFunctionType, N >
        DiscreteFunctionType;

      typedef typename ContainedDiscreteFunctionType :: DiscreteFunctionSpaceType
        ContainedDiscreteFunctionSpaceType;
      typedef CombinedSpace< ContainedDiscreteFunctionSpaceType, N, VariableBased >
        DiscreteFunctionSpaceType;

      typedef typename DiscreteFunctionSpaceType :: DomainFieldType DomainFieldType;
      typedef typename DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType;
      typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
      typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;

      typedef typename DiscreteFunctionSpaceType :: JacobianRangeType JacobianRangeType;

      typedef RangeFieldType DofType;
      typedef typename DiscreteFunctionSpaceType :: MapperType MapperType;
      typedef typename DiscreteFunctionSpaceType :: GridType GridType;
      typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;

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

      typedef ThreadSafeValue< UninitializedObjectStack > LocalDofVectorStackType;
      typedef StackAllocator< DofType, LocalDofVectorStackType* > LocalDofVectorAllocatorType;
      typedef DynamicReferenceVector< DofType, LocalDofVectorAllocatorType > LocalDofVectorType;

      typedef MutableLocalFunction< DiscreteFunctionType > LocalFunctionType;
    };


    //! @ingroup CombinedDFunction
    //! A class for combining N discrete function of the same
    //! type to a vector valued function
    template <class ContainedDiscreteFunctionImp,int N >
    class CombinedDiscreteFunction
    : public DiscreteFunctionDefault< CombinedDiscreteFunction< ContainedDiscreteFunctionImp, N > >
    {
      typedef CombinedDiscreteFunction< ContainedDiscreteFunctionImp, N > ThisType;
      typedef DiscreteFunctionDefault< CombinedDiscreteFunction< ContainedDiscreteFunctionImp, N > >
        BaseType;

    public:
      //! Discrete function this discrete function belongs to
      typedef ContainedDiscreteFunctionImp ContainedDiscreteFunctionType;
      typedef ContainedDiscreteFunctionType SubDiscreteFunctionType;

      //! Traits class with all necessary type definitions
      typedef DiscreteFunctionTraits< ThisType > Traits;

      //! Grid implementation
      typedef typename BaseType::GridType GridType;

      //! GridPart implementation
      typedef typename BaseType::GridPartType GridPartType;

      //! Discrete function type (identical to this type,
      //! needed as Barton-Nackman parameter
      typedef typename BaseType::DiscreteFunctionType
      DiscreteFunctionType;
      //! the combined discrete function type
      typedef typename BaseType::DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;
      //! Contained discrete function space
      typedef typename Traits::ContainedDiscreteFunctionSpaceType
      ContainedDiscreteFunctionSpaceType;
      typedef ContainedDiscreteFunctionSpaceType
      SubDiscreteFunctionSpaceType;
      //! Intrinsic type used for dofs (typically a float type)
      typedef typename BaseType::DofType DofType;
      //! Intrinsic type used for range field (like DofType)
      typedef typename BaseType::RangeFieldType RangeFieldType;
      //! Intrinsic type used for the domain field
      typedef typename BaseType::DomainFieldType DomainFieldType;
      //! Vector type used for the range field
      typedef typename BaseType::RangeType RangeType;
      //! Vector type used for the domain field
      typedef typename BaseType::DomainType DomainType;
      //! Mapper type (from the space)
      typedef typename Traits::MapperType MapperType;

      //! Iterator over dof container
      typedef typename BaseType::DofIteratorType DofIteratorType;
      //! Read-only iterator over dof container
      typedef typename BaseType::ConstDofIteratorType
      ConstDofIteratorType;

      typedef typename BaseType :: DofBlockPtrType DofBlockPtrType;
      typedef typename BaseType :: ConstDofBlockPtrType ConstDofBlockPtrType;

      typedef typename BaseType :: LocalDofVectorAllocatorType LocalDofVectorAllocatorType;

      using BaseType :: assign; // needs DofIterator!
      using BaseType :: axpy;
      using BaseType :: space;

      //! Constructor
      CombinedDiscreteFunction( const ContainedDiscreteFunctionType& func )
        : BaseType( "combined_"+func.name(), createSpace( func.space().gridPart() ), LocalDofVectorAllocatorType( &ldvStack_ ) ),
          ldvStack_( std::max( sizeof( DofType ), sizeof( DofType* ) ) * space().blockMapper().maxNumDofs() * DiscreteFunctionSpaceType::localBlockSize )
      {
        for (int i=0; i<N; ++i)
          func_[i] = new ContainedDiscreteFunctionType(func);
      }

      CombinedDiscreteFunction(const std::string &name,
                               const ContainedDiscreteFunctionSpaceType& spc)
        : BaseType( "combined_"+name, createSpace( spc.gridPart() ), LocalDofVectorAllocatorType( &ldvStack_ ) ),
          ldvStack_( std::max( sizeof( DofType ), sizeof( DofType* ) ) * space().blockMapper().maxNumDofs() * DiscreteFunctionSpaceType::localBlockSize )
      {
        for (int i=0; i<N; ++i)
          func_[i] = new ContainedDiscreteFunctionType(name,spc);
      }

      CombinedDiscreteFunction(const std::string &name,
                               const DiscreteFunctionSpaceType& spc )
        : BaseType( "combined_"+name, createSpace( spc.gridPart() ), LocalDofVectorAllocatorType( &ldvStack_) ),
          ldvStack_( std::max( sizeof( DofType ), sizeof( DofType* ) ) * space().blockMapper().maxNumDofs() * DiscreteFunctionSpaceType::localBlockSize )
      {
        for (int i=0; i<N; ++i)
        {
          func_[i] = new ContainedDiscreteFunctionType(name,space().containedSpace());
        }
      }

      //! Copy constructor
      //! The copy constructor copies the dofs
      CombinedDiscreteFunction(const ThisType &other)
        : BaseType( other.name()+"_copy", createSpace( other.space().gridPart() ), LocalDofVectorAllocatorType( &ldvStack_ ) ),
          ldvStack_( other.ldvStack_ )
      {
        for (int i=0; i<N; ++i)
          func_[i] = new ContainedDiscreteFunctionType(other.subFunction(i));
      }

      //! Destructor
      ~CombinedDiscreteFunction()
      {
        for (int i=0; i<N; ++i)
          delete func_[i];

        delete &space();
      }

      CombinedDiscreteFunction() = delete;
      ThisType& operator= ( const ThisType& ) = delete;
      ThisType& operator= ( ThisType&& ) = delete;

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::clear */
      void clear()
      {
        for (int i=0; i<N; ++i)
          func_[i]->clear();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::assign(const DiscreteFunctionInterfaceType &g) */
      void assign( const ThisType &g )
      {
        for( int i=0; i<N; ++i)
          func_[i]->assign( g.subFunction( i ) );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::size() const */
      int size() const
      {
        return func_[0]->size()*N;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::operator+=(const DiscreteFunctionInterfaceType &g) */
      ThisType &operator+= ( const ThisType &g )
      {
        for (int i=0; i<N; ++i)
          *func_[ i ] +=  g.subFunction( i );
        return *this;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::operator-=
       */
      using BaseType::operator-=;
      ThisType &operator-= ( const ThisType &g )
      {
        for( int i = 0; i < N; ++i )
          *func_[ i ] -=  g.subFunction( i );
        return *this;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::operator*=(const RangeFieldType &scalar) */
      DiscreteFunctionType& operator *= (const RangeFieldType &scalar)
      {
        for (int i=0; i<N; ++i)
          *func_[i] *= scalar;
        return *this;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::operator*=(const RangeFieldType &scalar) */
      DiscreteFunctionType& operator /= (const RangeFieldType &scalar)
      {
        for (int i=0; i<N; ++i)
          *func_[i] /= scalar;
        return *this;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::addScaled
       */
      void addScaled( const ThisType &g, const RangeFieldType &s )
      {
        axpy( g, s );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::axpy
       */
      void axpy( const RangeFieldType &s, const ThisType &g )
      {
        for (int i=0; i<N; ++i)
          func_[i]->axpy( s, g.subFunction( i ) );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::scalarProductDofs(const DiscreteFunctionInterfaceType &other) const */
      RangeFieldType scalarProductDofs ( const ThisType &other ) const
      {
        RangeFieldType ret( 0 );
        for( int i = 0; i < N; ++i )
          ret += func_[ i ]->scalarProductDofs( other.subFunction( i ) );
        return ret;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::read */
      template< class StreamTraits >
      void read ( InStreamInterface< StreamTraits >& in)
      {
        for (int i=0; i<N; ++i)
          func_[i]->read(in);
      }
      /** \copydoc Dune::Fem::DiscreteFunctionInterface::write */
      template< class StreamTraits >
      void write ( OutStreamInterface< StreamTraits >& out) const
      {
        for (int i=0; i<N; ++i)
          func_[i]->write(out);
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::print(std::ostream &out) const */
      void print( std :: ostream &out ) const
      {
        for (int i=0; i<N; ++i)
          func_[i]->print(out);
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::dofsValid() const */
      bool dofsValid () const
      {
        bool ret = func_[0]->dofsValid();
        for (int i=1;i<N;i++)
          ret |= func_[i]->dofsValid();
        return ret;
      }

      ConstDofBlockPtrType block ( unsigned int index ) const
      {
        // This is wrong with the current implementation of CombinedSpace
        const int containedSize = func_[ 0 ]->space().blockMapper().size();
        const int component = index / containedSize;
        const int containedIndex = index % containedSize;
        const ContainedDiscreteFunctionType& func = *(func_[ component ]);
        return func.block( containedIndex );
      }

      DofBlockPtrType block ( unsigned int index )
      {
        // This is wrong with the current implementation of CombinedSpace
        const int containedSize = func_[ 0 ]->space().blockMapper().size();
        const int component = index / containedSize;
        const int containedIndex = index % containedSize;
        return func_[ component ]->block( containedIndex );
      }

      const RangeFieldType &dof(unsigned int index) const
      {
        int variable = index / func_[0]->size();
        int point    = index % func_[0]->size();
        return func_[variable]->dof(point);
      }

      RangeFieldType &dof ( unsigned int index )
      {
        int variable = index / func_[0]->size();
        int point    = index % func_[0]->size();
        return func_[variable]->dof(point);
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::dbegin() const */
      ConstDofIteratorType dbegin () const
      {
        return ConstDofIteratorType(DofIteratorType(*this));
      }
      /** \copydoc Dune::Fem::DiscreteFunctionInterface::dend() const */
      ConstDofIteratorType dend () const
      {
        return ConstDofIteratorType(DofIteratorType(false,*this));
      }
      /** \copydoc Dune::Fem::DiscreteFunctionInterface::dbegin() */
      DofIteratorType dbegin ()
      {
        return DofIteratorType(*this);
      }
      /** \copydoc Dune::Fem::DiscreteFunctionInterface::dend() */
      DofIteratorType dend ()
      {
        return DofIteratorType(false,*this);
      }

      ContainedDiscreteFunctionType& subFunction( const int i )
      {
        return *(func_[i]);
      }

      const ContainedDiscreteFunctionType& subFunction( const int i ) const
      {
        return *(func_[i]);
      }

      ContainedDiscreteFunctionSpaceType& subSpace()
      {
        return space().containedSpace();
      }

    private:
      const ThisType& interface() const
      {
        return *this;
      }

      DiscreteFunctionSpaceType& createSpace( GridPartType& gp )
      {
        // we need to delete the space in the destructor
        return *(new DiscreteFunctionSpaceType( gp ));
      }

      typename Traits :: LocalDofVectorStackType ldvStack_;
      ContainedDiscreteFunctionType* func_[N];
      friend class CombinedDiscreteFunctionDofIterator<ContainedDiscreteFunctionType,N>;
    };

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
      typedef DiscreteFunctionTraits<CombinedDiscreteFunction< ContainedDiscreteFunctionImp,N> > Traits;
      typedef typename Traits::DiscreteFunctionType DiscreteFunctionType;
      typedef typename Traits::ContainedDiscreteFunctionType ContainedDiscreteFunctionType;
      typedef typename ContainedDiscreteFunctionType::DofIteratorType ContainedDofIteratorType;
      typedef typename ContainedDiscreteFunctionType::ConstDofIteratorType ContainedConstDofIteratorType;
      typedef typename Traits::DofType DofType;

      //! End constructor
      CombinedDiscreteFunctionDofIterator(bool end,const DiscreteFunctionType& df) :
        df_(const_cast<DiscreteFunctionType&>(df)),
        comp_(N-1),
        iter_(df.func_[N-1]->dend()),
        endIter_(df.func_[N-1]->dend())
      {}
      //! Constructor (const)
      CombinedDiscreteFunctionDofIterator(const DiscreteFunctionType& df) :
        df_(const_cast<DiscreteFunctionType&>(df)),
        comp_(0),
        iter_(df.func_[0]->dbegin()),
        endIter_(df.func_[0]->dend())
      {}
      //! End constructor
      CombinedDiscreteFunctionDofIterator(bool end,DiscreteFunctionType& df) :
        df_(df),
        comp_(N-1),
        iter_(df.func_[N-1]->dend()),
        endIter_(df.func_[N-1]->dend())
      {}
      //! Constructor
      CombinedDiscreteFunctionDofIterator(DiscreteFunctionType& df) :
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
      ThisType& operator=(const ThisType& other)
      {
        df_ = other.df_;
        comp_ = other.comp_;
        iter_ = other.iter_;
        endIter_ = other.endIter_;
        return *this;
      }

      //! return dof
      DofType& operator *()
      {
        return *iter_;
      }

      //! return dof read only
      const DofType& operator * () const
      {
        return *iter_;
      }

      //! go to next dof
      ThisType& operator++ ()
      {
        ++iter_;
        if (iter_==endIter_ && comp_<N-1)
        {
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

      void reset()
      {
        comp_    = 0;
        iter_    = df_.func_[ 0 ]->dbegin();
        endIter_ = df_.func_[ 0 ]->dend();
      }

  private:
      DiscreteFunctionType& df_;
      //! index
      mutable int comp_;
      mutable ContainedDofIteratorType iter_,endIter_;

    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_COMBINEDFUNCTION_COMBINEDFUNCTION_HH
