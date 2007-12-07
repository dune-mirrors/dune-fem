#ifndef DUNE_FEM_VECTORFUNCTION_VECTORFUNCTION_HH
#define DUNE_FEM_VECTORFUNCTION_VECTORFUNCTION_HH

#include <dune/common/typetraits.hh>

#include <dune/fem/storage/vector.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/localfunction/standardlocalfunction.hh>

namespace Dune
{

  template< class DiscreteFunctionSpaceImp, class DofVectorImp >
  class VectorDiscreteFunction;



  template< class DiscreteFunctionSpaceImp, class DofVectorImp >
  struct VectorDiscreteFunctionTraits
  {
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

    typedef DofVectorImp DofVectorType;

    typedef VectorDiscreteFunction< DiscreteFunctionSpaceType, DofVectorType >
      DiscreteFunctionType;

    typedef StandardLocalFunctionFactory
      < VectorDiscreteFunctionTraits< DiscreteFunctionSpaceType,
	                              DofVectorType > >
      LocalFunctionFactoryType;

    typedef LocalFunctionStack< LocalFunctionFactoryType >
      LocalFunctionStorageType;

    typedef typename LocalFunctionStorageType :: LocalFunctionType
      LocalFunctionType;

    typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;

    typedef typename DiscreteFunctionSpaceType :: DomainFieldType
      DomainFieldType;
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType
      RangeFieldType;

    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType
      JacobianRangeType;

    typedef typename DiscreteFunctionSpaceType :: MapperType MapperType;
    typedef typename DiscreteFunctionSpaceType :: GridType GridType;

    typedef typename DofVectorType :: FieldType DofType;
    typedef DofVectorType DofStorageType;

    typedef typename DofVectorType :: ConstIteratorType ConstDofIteratorType;
    typedef typename DofVectorType :: IteratorType DofIteratorType;

    typedef DofManager< GridType > DofManagerType;
  };



  template< class DiscreteFunctionSpaceImp, class DofVectorImp >
  class VectorDiscreteFunction
  : public DiscreteFunctionDefault
    < VectorDiscreteFunctionTraits< DiscreteFunctionSpaceImp, DofVectorImp > >
  {
  public:
    //! type of this class's traits
    typedef VectorDiscreteFunctionTraits
      < DiscreteFunctionSpaceImp, DofVectorImp >
      Traits;

    //! type of the associated discrete function space
    typedef typename Traits :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;

    //! type of the DoF storage array
    typedef typename Traits :: DofVectorType DofVectorType;

  private:
    typedef VectorDiscreteFunction< DiscreteFunctionSpaceType, DofVectorType >
      ThisType;
    typedef DiscreteFunctionDefault< Traits > BaseType;

  public:
    typedef typename Traits :: DiscreteFunctionType DiscreteFunctionType;

    typedef typename Traits :: LocalFunctionType LocalFunctionType;
    typedef typename Traits :: LocalFunctionFactoryType
      LocalFunctionFactoryType;
    typedef typename Traits :: LocalFunctionStorageType
      LocalFunctionStorageType;

    typedef typename Traits :: DomainType DomainType;
    typedef typename Traits :: RangeType RangeType;

    typedef typename Traits :: DomainFieldType DomainFieldType;
    typedef typename Traits :: RangeFieldType RangeFieldType;

    typedef typename Traits :: JacobianRangeType JacobianRangeType;

    typedef typename Traits :: DofType DofType;
    typedef typename Traits :: DofStorageType DofStorageType;

    typedef typename Traits :: ConstDofIteratorType ConstDofIteratorType;
    typedef typename Traits :: DofIteratorType DofIteratorType;

  private:
    typedef CheckVectorInterface< DofVectorType >
      CheckDofVectorType;

    typedef CompileTimeChecker
      < Conversion< RangeFieldType, DofType > :: sameType >
      Check_RangeFieldType_and_DofType_are_Equal;

  private:
    const LocalFunctionFactoryType lfFactory_;

    const std :: string name_;
    DofVectorType *const dofVector_;
    const bool freeDofVector_;

  public:
    //! Constructor
    inline VectorDiscreteFunction ( const std :: string name,
                                    const DiscreteFunctionSpaceType &dfSpace,
                                    DofVectorType &dofVector )
    : BaseType( dfSpace, lfFactory_ ),
      lfFactory_( *this ),
      name_( name ),
      dofVector_( &dofVector ),
      freeDofVector_( false )
    {
      assert( dofVector_->size() == (unsigned int)dfSpace.size() );
    }

    inline VectorDiscreteFunction ( const ThisType &other )
    : BaseType( other.space(), lfFactory_ ),
      lfFactory_( *this ),
      name_( other.name() ),
      dofVector_( new DofVectorType( other.dofVector() ) ),
      freeDofVector_( true )
    {
    }

    inline ~VectorDiscreteFunction ()
    {
      if( freeDofVector_ )
        delete dofVector_;
    }

  private:
    // prohibit assignment
    ThisType &operator= ( const ThisType & );

  public:
    inline ThisType &operator+= ( const ThisType &u )
    {
      dofVector() += u.dofVector();
      return *this;
    }

    inline ThisType &operator-= ( const ThisType &u )
    {
      dofVector() -= u.dofVector();
      return *this;
    }

    inline void addScaled ( const RangeFieldType &s,
                            const ThisType &u )
    {
      dofVector().addScaled( s, u.dofVector() );
    }

    inline void addScaled ( const ThisType &u,
                            const RangeFieldType &s )
    {
      addScaled( s, u );
    }

    inline void assign ( const ThisType &u )
    {
      dofVector().assign( u.dofVector() );
    }

    inline void clear ()
    {
      dofVector().assign( 0 );
    }

    inline ConstDofIteratorType dbegin () const
    {
      return dofVector().begin();
    }

    inline DofIteratorType dbegin ()
    {
      return dofVector().begin();
    }

    inline ConstDofIteratorType dend () const
    {
      return dofVector().end();
    }

    inline DofIteratorType dend ()
    {
      return dofVector().end();
    }

    inline const RangeFieldType &dof ( unsigned int index ) const
    {
      return dofVector()[ index ];
    }

    inline RangeFieldType &dof ( unsigned int index )
    {
      return dofVector()[ index ];
    }

    inline const std :: string &name () const
    {
      return name_;
    }

    inline RangeFieldType scalarProductDofs ( const ThisType &u ) const
    {
      return dofVector() * u.dofVector();
    }

    inline int size () const
    {
      return dofVector().size();
    }

    const DofVectorType &dofVector () const
    {
      return *dofVector_;
    }
    
    DofVectorType &dofVector ()
    {
      return *dofVector_;
    }
  };

};

#endif
